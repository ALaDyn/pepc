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

#include "multipole.h"

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
     use module_multipole
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8 !< TODO: currently unused
      real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8

      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      logical, public :: include_far_field_if_periodic = .false. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
      real*8, public  :: theta2       = 0.36  !< square of multipole opening angle
      real*8, public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential

      ! these are required by module_treediags, for now
      integer*8, allocatable, public :: interaction_keylist(:,:)
      integer, allocatable, public :: no_interaction_partners(:)
      real*8, allocatable, public :: interaction_vbox(:,:,:)

      namelist /calc_force_log2d/ mac_select, include_far_field_if_periodic, theta2, eps2


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! currently, all public functions in module_interaction_specific are obligatory
      public multipole_from_particle
      public shift_multipoles_up
      public results_add
      public calc_force_per_interaction
      public calc_force_per_particle
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
      !> Q = q and omega_k = 0, k = 1, ..., as the moments are centered
      !> on the particle
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        implicit none
        real*8, intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole
        complex(kind = 8), parameter :: omega0(pMultipole) = 0
        
        multipole = t_tree_node_interaction_data(particle_pos, &
                                     particle%q,   &
                                 abs(particle%q),  &
                                     omega0, 0. )
      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Accumulates multipole properties of child nodes to parent node.
      !> The monopole is just summed:
      !>
      !> Q_parent = Sum_children Q_child
      !>
      !> Higher moments are shifted by A (M2M) and summed:
      !>
      !> omega_parent,k = Sum_children ( -Q_child O_k(z0) / k
      !>   + Sum_l=1,k A_k,l(z0) omega_l )
      !>
      !> where z0 points from parent -> child.
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: nchild, j, k, l

        real*8, dimension(3) :: shift
        complex(kind = 8), parameter :: ic = (0, 1)
        complex(kind = 8) :: z0

        nchild = size(children)

        ! monopole moment Q
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

        ! higher order multipole moments omega
        parent%omega = 0

        do j=1,nchild
          shift = children(j)%coc - parent%coc
          z0 = shift(1) + ic * shift(2)

          do k = 1, pMultipole
            parent%omega(k) = parent%omega(k) - children(j)%charge * OMultipole(k, z0) / k

            do l = 1, k
              parent%omega(k) = parent%omega(k) + children(j)%omega(l) * ATranslate(k, l, z0)
            end do
          end do
        end do

        parent%bmax = maxval(sqrt((parent%coc(1)-children(1:nchild)%coc(1))**2+(parent%coc(2)-children(1:nchild)%coc(2))**2) + children(1:nchild)%bmax)

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

        call pepc_status("READ PARAMETERS, section calc_force_log2d")
        read(filehandle, NML=calc_force_log2d)

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

        write(filehandle, NML=calc_force_log2d)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> computes derived parameters for calc force module
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_prepare()
        use treevars, only : me, MPI_COMM_lpepc
        use module_fmm_periodicity, only : fmm_periodicity_init
        implicit none

        call fmm_periodicity_init(me, MPI_COMM_lpepc)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> initializes static variables of calc force module that depend 
      !> on particle data and might be reused on subsequent traversals
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_after_grow(particles, nparticles)
        use module_pepc_types
        use module_fmm_periodicity, only : fmm_periodicity_timestep
        use module_mirror_boxes, only : do_periodic
        implicit none
        type(t_particle), dimension(:), intent(in) :: particles
        integer, intent(in) :: nparticles

        ! calculate spherical multipole expansion of central box
        ! this cannot be done in calc_force_per_particle() since there, possibly
        ! other particles are used than we need for the multipoles
        ! e.g. in the case of a second traverse for test/grid particles
        if ((do_periodic) .and. (include_far_field_if_periodic)) then
          call fmm_periodicity_timestep(particles, nparticles)
        end if

      end subroutine      


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> subroutine must return the estimated number of iteractions per particle
      !> for the current mac and/or parameters and the supplied total number of particles
      !>
      !> TODO: is this correct for 2D?
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
        implicit none
        integer, intent(in) :: npart_total !< total number of particles
        integer, intent(out) :: nintmax !< maximum number of interactions per particle

        real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

        ! Estimate of interaction list length - Hernquist expression
        ! applies for BH-MAC
        invnintmax = max(theta2 / (35.*log(1.*npart_total)) , 1._8/npart_total)
        nintmax    = int(1._8/invnintmax)

      end subroutine



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> finalizes the calc force module at end of simulation
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_finalize()
        implicit none
        ! nothing to do here
      end subroutine calc_force_finalize


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> generic Multipole Acceptance Criterion
      !>
      !> TODO: is this correct for 2D?
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function mac(particle, node, dist2, boxlength2)

        use treevars, only : tree_nodes
        implicit none

        logical :: mac
        integer, intent(in) :: node
        type(t_particle), intent(in) :: particle
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2

        select case (mac_select)
            case (0)
              ! Barnes-Hut-MAC
              mac = (theta2 * dist2 > boxlength2)
            case (1)
               ! Bmax-MAC
              mac = (theta2 * dist2 > min(tree_nodes(node)%bmax**2,3.0*boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
            case default
              ! N^2 code
              mac = .false.
        end select

      end function

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> clears result in t_particle datatype - usually, this function does not need to be touched
      !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
      !> function cannot reside in module_interaction_specific that may not include module_pepc_types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine particleresults_clear(particles, nparticles)
        use module_pepc_types
        implicit none
        type(t_particle), intent(inout) :: particles(nparticles)
        integer, intent(in) :: nparticles

        particles(1:nparticles)%results = EMPTY_PARTICLE_RESULTS

      end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction(particle, node, key, delta, dist2, vbox, node_is_leaf)
          use module_pepc_types
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: node
          integer*8, intent(in) :: key
          type(t_particle), intent(inout) :: particle
          logical, intent(in) :: node_is_leaf
          real*8, intent(in) :: vbox(3), delta(3), dist2


          real*8 :: exy(2), phic

          !  compute 2D-Coulomb fields and potential of particle p with node

          if (node_is_leaf) then
            ! It's a leaf, do direct summation
            call calc_force_log_2D_direct(node, delta(1:2), exy(1), exy(2), phic)
          else
            ! It's a twig, do ME
            call calc_force_log_2D(node, delta(1:2), exy(1), exy(2), phic)
          end if

          particle%results%e         = particle%results%e    + exy
          particle%results%pot       = particle%results%pot  + phic
          particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

        end subroutine calc_force_per_interaction

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles)
          use module_debug, only : pepc_status
          use module_pepc_types
          use treevars, only : me
          use module_fmm_periodicity
          use module_mirror_boxes
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)
          real*8 :: e_lattice(2), phi_lattice
          integer :: p

          call pepc_status('CALC FORCE PER PARTICLE')

          potfarfield  = 0.
          potnearfield = 0.

          if ((do_periodic) .and. (include_far_field_if_periodic)) then

             do p=1,nparticles
                call fmm_periodicity_sum_lattice_force(particles(p)%x, e_lattice, phi_lattice)

                potfarfield  = potfarfield  + phi_lattice               * particles(p)%data%q
                potnearfield = potnearfield + particles(p)%results%pot  * particles(p)%data%q

                particles(p)%results%e     = particles(p)%results%e     + e_lattice
                particles(p)%results%pot   = particles(p)%results%pot   + phi_lattice
             end do

          end if

          call pepc_status('CALC FORCE PER PARTICLE DONE')

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Computes the particle multipole interaction potential and force
        !> field:
        !>
        !> Phi = - Re { Q log(z) + Sum_k=1^p omega_k M_k(z) }
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_log_2D(t, d, ex, ey, phi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t
          real*8, intent(in) :: d(2) !< separation vector precomputed in walk_single_particle
          real*8, intent(out) ::  ex, ey, phi

          complex(kind = 8), parameter :: ic = (0, 1)
          complex(kind = 8) :: z, cphi, cf, rz, mtaylor_

          integer :: k

          z = d(1) + ic * d(2)
          rz = 1. / z

          cphi = -t%charge * log(z)
          cf   = t%charge * rz

          !
          ! This uses the fact that
          !
          !   MTaylor(0, z) = 1
          !   MTaylor(k, z) = MTaylor(k - 1, z) / z
          !
          ! and
          !
          !   MTaylorPrime(k, z) = -k MTaylor(k, z) / z 
          !
          ! to save time. A naive implementation would be:
          !
          !   do k = 1, pMultipole
          !     cphi = cphi - t%omega(k) * MTaylor(k, z)
          !     cf   = cf + t%omega(k) * MTaylorPrime(k, z)
          !   end do
          !
          mtaylor_ = rz

          do k = 1, pMultipole
            cphi     = cphi - t%omega(k) * mtaylor_
            mtaylor_ = mtaylor_ * rz
            cf       = cf - k * t%omega(k) * mtaylor_
          end do

          phi = real(cphi, kind = 8)
          ex  = real(cf, kind = 8)
          ey  = -real(aimag(cf), kind = 8)

        end subroutine calc_force_log_2D
        

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Computes the direct particle particle interaction potential and
      !> force field:
      !>
      !> Phi = -q log(||r - r0||)
      !> E = - grad Phi = q / ||r - r0||^2 (x - x0, y - y0)
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_log_2D_direct(t, d, ex, ey, phi)
          use module_pepc_types
          use treevars
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t
          real*8, intent(in) :: d(2) !< separation vector precomputed in walk_single_particle
          real*8, intent(out) ::  ex, ey, phi

          real*8 :: dx, dy, d2, rd2, charge

          dx = d(1)
          dy = d(2)

          d2  = dot_product(d, d) + eps2
          !d2  = dist2
          rd2 = 1./d2

          charge = t%charge

          phi = - 0.5 * charge * log(d2)

          ex = charge * dx * rd2
          ey = charge * dy * rd2

      end subroutine calc_force_log_2D_direct


  end module module_interaction_specific