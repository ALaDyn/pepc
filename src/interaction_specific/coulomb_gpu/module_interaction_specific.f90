! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_interaction_specific
     use module_pepc_types
     use module_interaction_specific_types
     use module_accelerator
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public t_particle_thread
      real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8 !< TODO: currently unused
      real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8

      integer, public :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      logical, public :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
      real*8, public  :: theta2       = 0.36  !< square of multipole opening angle
      real*8, public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
      real*8, public  :: kelbg_invsqrttemp = 0.0 !< inverse square root of temperature for kelbg potential

      !> debug stuff for interaction partners (currently only used by pepc-f frontend) 
      !> see module_treediags::write_interaction_partners_to_vtk() and module_interaction_specific::calc_force_per_interaction(force_law==6) 
      integer(kind_node), allocatable,public :: interaction_nodelist(:,:)
      integer(kind_node), allocatable,public :: no_interaction_partners(:)
      real*8, allocatable,public :: interaction_vbox(:,:,:)
      namelist /calc_force_coulomb/ force_law, mac_select, include_far_field_if_periodic, theta2, eps2, kelbg_invsqrttemp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! currently, all public functions in module_interaction_specific are obligatory
      public multipole_from_particle
      public shift_multipoles_up
      public results_add
      interface calc_force_per_interaction
         module procedure calc_force_per_interaction_global, calc_force_per_interaction_thread ! these stay private and are for thread-(non)-local types
      end interface calc_force_per_interaction
      public calc_force_per_interaction
      public compute_iact_list
      public compute_iact_list_leaf
      public calc_force_per_particle
      interface mac
         module procedure mac_global, mac_thread ! these stay private and are for thread-(non)-local types
      end interface mac
      public mac
      public particleresults_clear
      public calc_force_read_parameters
      public calc_force_write_parameters
      public calc_force_finalize
      public calc_force_prepare
      public calc_force_after_grow
      public get_number_of_interactions_per_particle

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Computes multipole properties of a single particle
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        implicit none
        real*8, intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole

        multipole = t_tree_node_interaction_data(particle_pos, &
                                     particle%q,   &
                                 abs(particle%q),  &
                                     [0., 0., 0.], &
                                     [0., 0., 0.], &
                                      0., 0., 0., 0. )
      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: nchild, j

        real*8 :: shift(1:3)

        nchild = size(children)

        parent%charge     = SUM( children(1:nchild)%charge )
        parent%abs_charge = SUM( children(1:nchild)%abs_charge )

        ! centre of charge
        parent%coc        = [0., 0., 0.]

        if (parent%abs_charge .ne. 0.) then
          ! use center-of-charge because we may divide by abs_charge
          do j=1,nchild
            parent%coc(1:3) = parent%coc(1:3) + ( children(j)%coc(1:3) * children(j)%abs_charge )
          end do

          parent%coc(1:3) = parent%coc(1:3) / parent%abs_charge
        else
          ! use geometric center
          do j=1,nchild
            parent%coc(1:3) = parent%coc(1:3) +   children(j)%coc(1:3)
          end do

         parent%coc(1:3) = parent%coc(1:3) / nchild
        endif

        ! multipole properties
        parent%dip    = [0., 0., 0.]
        parent%quad   = [0., 0., 0.]
        parent%xyquad = 0.
        parent%yzquad = 0.
        parent%zxquad = 0.

        do j=1,nchild
          shift(1:3) = parent%coc(1:3) - children(j)%coc

          ! dipole moment
          parent%dip = parent%dip + children(j)%dip - children(j)%charge*shift(1:3)

          ! quadrupole moment
          parent%quad(1:3) = parent%quad(1:3) + children(j)%quad(1:3) - 2*children(j)%dip(1:3)*shift(1:3) + children(j)%charge*shift(1:3)**2

          parent%xyquad = parent%xyquad + children(j)%xyquad - children(j)%dip(1)*shift(2) - children(j)%dip(2)*shift(1) + children(j)%charge*shift(1)*shift(2)
          parent%yzquad = parent%yzquad + children(j)%yzquad - children(j)%dip(2)*shift(3) - children(j)%dip(3)*shift(2) + children(j)%charge*shift(2)*shift(3)
          parent%zxquad = parent%zxquad + children(j)%zxquad - children(j)%dip(3)*shift(1) - children(j)%dip(1)*shift(3) + children(j)%charge*shift(3)*shift(1)
        end do

        parent%bmax = maxval(sqrt((parent%coc(1)-children(1:nchild)%coc(1))**2+(parent%coc(2)-children(1:nchild)%coc(2))**2+(parent%coc(3)-children(1:nchild)%coc(3))**2) + children(1:nchild)%bmax)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> adds res2 to res1
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in) :: res2

        res1%e    = res1%e    + res2%e
        res1%pot  = res1%pot  + res2%pot
      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> reads interaction-specific parameters from file
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section calc_force_coulomb")
        read(filehandle, NML=calc_force_coulomb)

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> writes interaction-specific parameters to file
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_write_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=calc_force_coulomb)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> computes derived parameters for calc force module
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_prepare()
        use treevars, only : me, MPI_COMM_lpepc
        use module_fmm_framework, only : fmm_framework_init
        use pthreads_stuff, only: pthreads_createthread, pthreads_sched_yield, THREAD_TYPE_ACCELERATOR
        use module_atomic_ops, only: atomic_load_int, atomic_store_int, atomic_allocate_int
        !use module_tree_communicator, only: tree_comm_thread_counter

        implicit none

        call fmm_framework_init(me, MPI_COMM_lpepc)

        ! start ACCELERATOR thread (GPUs, MIC, BG/Q Quad, ...)
        if (.not. associated(acc%thread_status)) then 
           call atomic_allocate_int(acc%thread_status)
           call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_WAITING)
        end if

        ! do that only once...
        if (      atomic_load_int(acc%thread_status) /= ACC_THREAD_STATUS_STARTED &
            .and. atomic_load_int(acc%thread_status) /= ACC_THREAD_STATUS_STARTING ) then
           call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTING)
           !tree_comm_thread_counter = tree_comm_thread_counter + 1
           ERROR_ON_FAIL(pthreads_createthread(acc%acc_thread, c_funloc(acc_loop), c_null_ptr, thread_type = THREAD_TYPE_ACCELERATOR, counter = 99))

           ! we have to wait here until the communicator has really started to find out its processor id
           do while (atomic_load_int(acc%thread_status) /= ACC_THREAD_STATUS_STARTED)
              ERROR_ON_FAIL(pthreads_sched_yield())
           end do
        end if

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> initializes static variables of calc force module that depend 
      !> on particle data and might be reused on subsequent traversals
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_after_grow(particles)
        use module_pepc_types
        use module_fmm_framework, only : fmm_framework_timestep
        use module_mirror_boxes, only : do_periodic
        implicit none
        type(t_particle), dimension(:), intent(in) :: particles

        ! calculate spherical multipole expansion of central box
        ! this cannot be done in calc_force_per_particle() since there, possibly
        ! other particles are used than we need for the multipoles
        ! e.g. in the case of a second traverse for test/grid particles
        if ((do_periodic) .and. (include_far_field_if_periodic)) then
          call fmm_framework_timestep(particles)
        end if

      end subroutine      


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> subroutine must return the estimated number of iteractions per particle
      !> for the current mac and/or parameters and the supplied total number of particles
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
        implicit none
        integer*8, intent(in) :: npart_total !< total number of particles
        integer*8, intent(out) :: nintmax !< maximum number of interactions per particle

        real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

        ! Estimate of interaction list length - Hernquist expression
        ! applies for BH-MAC
        invnintmax = max(theta2 / (35._8*log(1._8*npart_total)) , 1._8/npart_total)
        nintmax    = int(1._8/invnintmax)

      end subroutine



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> finalizes the calc force module at end of simulation
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_finalize()
        use module_atomic_ops, only: atomic_load_int, atomic_store_int, atomic_allocate_int
        implicit none

        ! stop GPU thread, or at least signal it to stop
        call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STOPPED)

        return

      end subroutine calc_force_finalize


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> generic Multipole Acceptance Criterion
      !>
      !> we have two version here for calls with a thread-local particle type and a
      !> global particle type
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function mac_(node, dist2, boxlength2)
        implicit none

        logical :: mac_
        type(t_tree_node_interaction_data), intent(in) :: node
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2

        select case (mac_select)
            case (0)
              ! Barnes-Hut-MAC
              mac_ = (theta2 * dist2 > boxlength2)
            case (1)
               ! Bmax-MAC
              mac_ = (theta2 * dist2 > min(node%bmax**2, 3.0 * boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
            case default
              ! N^2 code
              mac_ = .false.
        end select

      end function
      function mac_global(particle, node, dist2, boxlength2)
        implicit none

        logical :: mac_global
        type(t_tree_node_interaction_data), intent(inout) :: node
        type(t_particle), intent(in) :: particle
        real*8, intent(inout) :: dist2
        real*8, intent(inout) :: boxlength2

        mac_global = mac_(node, dist2, boxlength2)

      end function
      function mac_thread(particle, node, dist2, boxlength2)
        implicit none

        logical :: mac_thread
        type(t_tree_node_interaction_data), intent(inout) :: node
        type(t_particle_thread), intent(in) :: particle
        real*8, intent(inout) :: dist2
        real*8, intent(inout) :: boxlength2

        mac_thread = mac_(node, dist2, boxlength2)

      end function

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> clears result in t_particle datatype - usually, this function does not need to be touched
      !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
      !> function cannot reside in module_interaction_specific that may not include module_pepc_types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine particleresults_clear(particles)
        use module_pepc_types
        implicit none
        type(t_particle), intent(inout) :: particles(:)

        particles(:)%results = EMPTY_PARTICLE_RESULTS

      end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !> It comes in a thread-local and a 'global' version to be called for
        !> single particles (global) and particle lists (thread-local)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction_global(particle, node, node_idx, delta, dist2, vbox, node_is_leaf)
          use module_pepc_types
          use treevars
          use module_coulomb_kernels
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: node
          integer(kind_node), intent(in) :: node_idx
          type(t_particle), intent(inout) :: particle
          logical, intent(in) :: node_is_leaf
          real*8, intent(in) :: vbox(3), delta(3), dist2


          real*8 :: exyz(3), phic

          select case (force_law)
            case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list

                if (node_is_leaf) then
                  call calc_force_coulomb_2D_direct(node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
                else
                  call calc_force_coulomb_2D(       node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
                end if
                exyz(3) = 0.

            case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list

                if (node_is_leaf) then
                    call calc_force_coulomb_3D_direct(node, delta, dist2 + eps2, exyz, phic)
                else
                    call calc_force_coulomb_3D(       node, delta, dist2 + eps2, exyz, phic)
                end if

            case (4)  ! LJ potential for quiet start
                call calc_force_LJ(node, delta, dist2, eps2, exyz, phic)

            case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction
                      !  and Kelbg for particle-particle interaction

                if (node_is_leaf) then
                    ! It's a leaf, do direct summation with kelbg
                    call calc_force_kelbg_3D_direct(particle, node, delta, dist2, kelbg_invsqrttemp, exyz, phic)
                else
                    ! It's a twig, do ME with coulomb
                    call calc_force_coulomb_3D(node, delta, dist2, exyz, phic)
                end if
            case (6)  !  used to save interaction partners
                no_interaction_partners(particle%label)=no_interaction_partners(particle%label)+1
                interaction_nodelist(particle%label,no_interaction_partners(particle%label))=node_idx
                interaction_vbox(particle%label,no_interaction_partners(particle%label),1:3)=vbox(1:3)
            case default
              exyz = 0.
              phic = 0.
          end select

          particle%results%e         = particle%results%e    + exyz
          particle%results%pot       = particle%results%pot  + phic
          particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

        end subroutine calc_force_per_interaction_global

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for thread-local list.
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction_thread(particle, node, node_idx, delta, dist2, vbox, node_is_leaf)
          implicit none

          type(t_tree_node_interaction_data), intent(inout) :: node
          integer(kind_node), intent(in) :: node_idx
          type(t_particle_thread), intent(inout) :: particle
          logical, intent(inout) :: node_is_leaf
          real*8, intent(in) :: vbox(3), delta(3), dist2

          !! this will be an array of interaction partners, to be associated with the thread-local particle data
          !!type(t_iact_partner), pointer :: iact_partner(:)

#ifdef LEAFY
          if(node_is_leaf) then
             select case (particle%queued_l)

             case (-1) ! allocate interaction list and start filling it

                allocate(particle%partner_l(MAX_IACT_PARTNERS))
                particle%queued_l = 0
                call add_to_iact_list_leaf(particle, node, delta)

             case (MAX_IACT_PARTNERS) ! compute list and add current partner to empty list

                ! compute list, which leaves us we NO list
                call compute_iact_list_leaf(particle)
                ! create new list
                allocate(particle%partner_l(MAX_IACT_PARTNERS))
                particle%queued_l = 0
                ! add to new list
                call add_to_iact_list_leaf(particle, node, delta)

             case default ! add to list

                call add_to_iact_list_leaf(particle, node, delta)

             end select
          else
#endif
             select case (particle%queued)

             case (-1) ! allocate interaction list and start filling it

                allocate(particle%partner(MAX_IACT_PARTNERS))
                particle%queued = 0
                call add_to_iact_list(particle, node, delta, node_is_leaf)

             case (MAX_IACT_PARTNERS) ! compute list and add current partner to empty list

                ! compute list, which leaves us we NO list
                call compute_iact_list(particle)
                ! create new list
                allocate(particle%partner(MAX_IACT_PARTNERS))
                particle%queued = 0
                ! add to new list
                call add_to_iact_list(particle, node, delta, node_is_leaf)

             case default ! add to list

                call add_to_iact_list(particle, node, delta, node_is_leaf)

             end select
#ifdef LEAFY
          endif
#endif

          return

          contains

             subroutine add_to_iact_list(particle, node, delta, node_is_leaf)
                implicit none

                type(t_particle_thread), intent(inout) :: particle
                type(t_tree_node_interaction_data), intent(in) :: node
                real*8, intent(in) :: delta(3)
                logical, intent(in) :: node_is_leaf

                integer :: queued

                queued = particle%queued + 1

                particle%partner(queued)%delta = delta
                particle%partner(queued)%node  = node
                particle%partner(queued)%leaf  = node_is_leaf
                particle%queued = queued

             end subroutine add_to_iact_list

             subroutine add_to_iact_list_leaf(particle, node, delta)
                implicit none

                type(t_particle_thread), intent(inout) :: particle
                type(t_tree_node_interaction_data), intent(in) :: node
                real*8, intent(in) :: delta(3)

                integer :: queued

                queued = particle%queued_l + 1

                particle%partner_l(queued)%delta = delta
                particle%partner_l(queued)%charge  = node%charge
                particle%queued_l = queued

             end subroutine add_to_iact_list_leaf

        end subroutine calc_force_per_interaction_thread

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> compute the list of interactions
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine compute_iact_list(particle)
           use module_coulomb_kernels
           implicit none

           type(t_particle_thread), intent(inout) :: particle

           integer :: idx
           real*8 :: exyz(3), phic
           real*8 :: delta(3), dist2
           type(t_tree_node_interaction_data) :: node
#ifndef LEAFY
           logical :: node_is_leaf
#endif

!           write(*,*) ' LLLLLLLL computing a list of length', particle%queued

           select case (force_law)
           case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list

              do idx = 1, particle%queued

#ifndef LEAFY
                 node_is_leaf = particle%partner(idx)%leaf
#endif
                 node  = particle%partner(idx)%node
                 delta = particle%partner(idx)%delta
                 dist2 = dot_product(delta, delta)

#ifndef LEAFY
                 if (node_is_leaf) then
                    call calc_force_coulomb_2D_direct(node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
                 else
#endif
                    call calc_force_coulomb_2D(       node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
#ifndef LEAFY
                 end if
#endif
                 exyz(3) = 0.

                 particle%results%e         = particle%results%e    + exyz
                 particle%results%pot       = particle%results%pot  + phic
                 particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

              end do

           case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list

              call dispatch_list(particle, eps2, WORKLOAD_PENALTY_INTERACTION)

           case (4)  ! LJ potential for quiet start

              do idx = 1, particle%queued

                 node  = particle%partner(idx)%node
                 delta = particle%partner(idx)%delta
                 dist2 = dot_product(delta, delta)

                 call calc_force_LJ(node, delta, dist2, eps2, exyz, phic)

                 particle%results%e         = particle%results%e    + exyz
                 particle%results%pot       = particle%results%pot  + phic
                 particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

              end do

           case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction

              do idx = 1, particle%queued

#ifndef LEAFY
                 node_is_leaf = particle%partner(idx)%leaf
#endif
                 node  = particle%partner(idx)%node
                 delta = particle%partner(idx)%delta
                 dist2 = dot_product(delta, delta)

                 !  and Kelbg for particle-particle interaction
#ifndef LEAFY
                 if (node_is_leaf) then
                    ! It's a leaf, do direct summation with kelbg
                    call calc_force_kelbg_3D_direct(particle, node, delta, dist2, kelbg_invsqrttemp, exyz, phic)
                 else
#endif
                    ! It's a twig, do ME with coulomb
                    call calc_force_coulomb_3D(node, delta, dist2, exyz, phic)
#ifndef LEAFY
                 end if
#endif

                 particle%results%e         = particle%results%e    + exyz
                 particle%results%pot       = particle%results%pot  + phic
                 particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

              end do

           case default

              exyz = 0.
              phic = 0.

              particle%results%e         = particle%results%e    + exyz
              particle%results%pot       = particle%results%pot  + phic
              particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

           end select

           ! free space and refresh list
           particle%queued = -1
           nullify(particle%partner)

           return

        end subroutine compute_iact_list
        subroutine compute_iact_list_leaf(particle)
           use module_coulomb_kernels
           implicit none

           type(t_particle_thread), intent(inout) :: particle

           integer :: idx
           real*8 :: exyz(3), phic
           real*8 :: delta(3), dist2

           interface 
              subroutine kernel_leaf(particle, eps2, WORKLOAD_PENALTY_INTERACTION)
              use module_pepc_types
              use module_interaction_specific_types
              implicit none
              type(t_particle_thread), intent(inout) :: particle
              real*8, intent(in) :: eps2, WORKLOAD_PENALTY_INTERACTION
              end subroutine
           end interface

!           write(*,*) ' LLLL computing a list of length', particle%queued_l

           select case (force_law)
           case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list

              do idx = 1, particle%queued_l

                 delta = particle%partner_l(idx)%delta
                 dist2 = dot_product(delta, delta)

!!!                 call calc_force_coulomb_2D_leaf(particle%partner_l(idx)%charge, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
                 exyz(3) = 0.

                 particle%results%e         = particle%results%e    + exyz
                 particle%results%pot       = particle%results%pot  + phic
                 particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

              end do

           case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list

              call kernel_leaf(particle, eps2, WORKLOAD_PENALTY_INTERACTION)

           case (4)  ! LJ potential for quiet start

              do idx = 1, particle%queued_l

!!!                 node  = particle%partner(idx)%node
                 delta = particle%partner_l(idx)%delta
                 dist2 = dot_product(delta, delta)

!!!                 call calc_force_LJ(node, delta, dist2, eps2, exyz, phic)

                 particle%results%e         = particle%results%e    + exyz
                 particle%results%pot       = particle%results%pot  + phic
                 particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

              end do

           case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction

              do idx = 1, particle%queued_l

!!!                 node  = particle%partner(idx)%node
                 delta = particle%partner_l(idx)%delta
                 dist2 = dot_product(delta, delta)

!!!                 !  and Kelbg for particle-particle interaction
!!!                 if (node_is_leaf) then
!!!                    ! It's a leaf, do direct summation with kelbg
!!!                    call calc_force_kelbg_3D_direct(particle, node, delta, dist2, kelbg_invsqrttemp, exyz, phic)
!!!                 else
!!!                    ! It's a twig, do ME with coulomb
!!!                    call calc_force_coulomb_3D(node, delta, dist2, exyz, phic)
!!!                 end if

                 particle%results%e         = particle%results%e    + exyz
                 particle%results%pot       = particle%results%pot  + phic
                 particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

              end do

           case default

              exyz = 0.
              phic = 0.

              particle%results%e         = particle%results%e    + exyz
              particle%results%pot       = particle%results%pot  + phic
              particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

           end select

           ! free space and refresh list
           particle%queued_l = -1
           nullify(particle%partner_l)

           return

        end subroutine compute_iact_list_leaf


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles)
          use treevars, only: num_threads
          use module_debug, only : pepc_status
          use module_pepc_types
          use treevars, only : me
          use module_fmm_framework
          use module_mirror_boxes
          implicit none

          type(t_particle), intent(inout) :: particles(:)
          real*8 :: e_lattice(3), phi_lattice
          integer(kind_particle) :: p

          call pepc_status('CALC FORCE PER PARTICLE')

          potfarfield  = 0.
          potnearfield = 0.

          if ((do_periodic) .and. (include_far_field_if_periodic)) then

             if ((me==0) .and. (force_law .ne. 3)) write(*,*) "Warning: far-field lattice contribution is currently only supported for force_law==3"
          !$ call omp_set_num_threads(num_threads)
          !$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(particles) SCHEDULE(RUNTIME) REDUCTION(+:potfarfield,potnearfield)
             do p=1,size(particles)
                call fmm_sum_lattice_force(particles(p)%x, e_lattice, phi_lattice)

                potfarfield  = potfarfield  + phi_lattice               * particles(p)%data%q
                potnearfield = potnearfield + particles(p)%results%pot  * particles(p)%data%q

                particles(p)%results%e     = particles(p)%results%e     + e_lattice
                particles(p)%results%pot   = particles(p)%results%pot   +  phi_lattice
             end do
          !$OMP  END PARALLEL DO
          !$ call omp_set_num_threads(1)

          end if

          call pepc_status('CALC FORCE PER PARTICLE DONE')

        end subroutine calc_force_per_particle

end module module_interaction_specific
