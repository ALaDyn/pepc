! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>


module module_interaction_specific
     use module_pepc_kinds, only: kind_physics
     use module_pepc_types
     use module_interaction_specific_types
     use module_multipole
     implicit none
     save
     private

      integer, public :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      logical, public :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
      real(kind_physics), public  :: theta2       = 0.36  !< square of multipole opening angle
      real(kind_physics), public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
      real(kind_physics), public  :: kelbg_invsqrttemp = 0.0 !< inverse square root of temperature for kelbg potential

      namelist /calc_force/ force_law, mac_select, include_far_field_if_periodic, theta2, eps2, kelbg_invsqrttemp

      ! currently, all public functions in module_interaction_specific are obligatory
      public multipole_from_particle
      public shift_multipoles_up
      public results_add
      public calc_force_per_interaction_with_self
      public calc_force_per_interaction_with_leaf
      public calc_force_per_interaction_with_twig
      public calc_force_per_particle
      public mac
      public particleresults_clear
      public calc_force_read_parameters
      public calc_force_write_parameters
      public calc_force_finalize
      public calc_force_prepare
      public calc_force_after_grow
      public get_number_of_interactions_per_particle

      contains

      !>
      !> Computes multipole properties of a single particle
      !>
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        use module_globals , only: vtilde
        use module_shortcut, only: zero,one  
        implicit none
        real(kind_physics), intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole

        real(kind_physics) vx, vy, vz, q,g

        
        q  = particle%q
        g  = particle%g
        
        ! particle velocity
        vx = q*particle%v(1)
        vy = q*particle%v(2)
        vz = q*particle%v(3)

        multipole = t_tree_node_interaction_data(particle_pos, &            !position
                                     q,   &                                 !charge
                                     abs(q),  &                             !absolute charge/mass
                                     (/zero, zero, zero/), &                      !dipole
                                     (/zero, zero, zero/), &                      !diagonal quadrupole
                                       zero, zero, zero, &                        !mixed quadrupole
                                     (/ vx/g, vy/g, vz/g /), &                    !current density monopole - qj*vxj - qj*vxj - qj*vxj
                                     (/zero, zero, zero/), &                      !current density dipole - qj*vxj*x - qj*vxj*y - qj*vxj*z
                                     (/zero, zero, zero/), &                      !current density dipole - qj*vyj*x - qj*vyj*y - qj*vyj*z
                                     (/zero, zero, zero/), &                      !current density dipole - qj*vzj*x - qj*vzj*y - qj*vzj*z
                                     (/zero, zero, zero/), &                      !current density diagonal quadrupole - qj*vxj*x**2 - qj*vxj*y**2 - qj*vxj*z**2
                                     (/zero, zero, zero/), &                      !current density diagonal quadrupole - qj*vyj*x**2 - qj*vyj*y**2 - qj*vyj*z**2
                                     (/zero, zero, zero/), &                      !current density diagonal quadrupole - qj*vzj*x**2 - qj*vzj*y**2 - qj*vzj*z**2
                                     (/zero, zero, zero/), &                      !current density mixed quadrupole
                                     (/zero, zero, zero/), &                      !current density mixed quadrupole
                                     (/zero, zero, zero/), zero )                   !current density mixed quadrupole /  bmax
 

      end subroutine


      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      subroutine shift_multipoles_up(parent, children)
        use module_shortcut, only: zero,half  
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: nchild, j

        real(kind_physics) :: shift(1:3)
!        real(kind_physics), parameter :: half = 0.5_kind_physics

        nchild = size(children)

        parent%charge     = SUM( children(1:nchild)%charge )
        parent%abs_charge = SUM( children(1:nchild)%abs_charge )

        do j=1,3
            parent%monoj(j)                 = SUM( children(1:nchild)%monoj(j) )
        enddo


        ! centre of charge
        parent%coc(1:3)     = zero

        if (parent%abs_charge .ne. 0.) then
          ! use center-of-charge because we may divide by abs_charge
          do j=1,nchild
            parent%coc(1:3) = parent%coc(1:3) + ( children(j)%coc(1:3) * children(j)%abs_charge )
          end do

          parent%coc(1:3) = parent%coc(1:3) / parent%abs_charge
        else
          ! use geometric center
          do j=1,nchild
            parent%coc(1:3) = parent%coc(1:3) + children(j)%coc(1:3)
          end do

         parent%coc(1:3) = parent%coc(1:3) / nchild
        endif

        
        ! multipole properties
        parent%dip(1:3)    = zero
        parent%quad(1:3)   = zero
        parent%xyquad      = zero
        parent%yzquad      = zero
        parent%zxquad      = zero
        parent%dipjx(1:3)  = zero
        parent%dipjy(1:3)  = zero
        parent%dipjz(1:3)  = zero
        parent%quadjx(1:3) = zero
        parent%quadjy(1:3) = zero
        parent%quadjz(1:3) = zero
        parent%quadjxy(1:3)= zero
        parent%quadjyz(1:3)= zero
        parent%quadjzx(1:3)= zero
        

        do j=1,nchild

          shift(1:3) = parent%coc(1:3) - children(j)%coc

          ! dipole moment
          parent%dip = parent%dip + children(j)%dip + children(j)%charge * shift(1:3)

          ! quadrupole moment
          parent%quad(1:3) = parent%quad(1:3) + children(j)%quad(1:3) + children(j)%dip(1:3) * shift(1:3) + half * children(j)%charge * shift(1:3)**2

          parent%xyquad = parent%xyquad + children(j)%xyquad + children(j)%dip(1) * shift(2) + children(j)%dip(2) * shift(1) + children(j)%charge * shift(1) * shift(2)
          parent%yzquad = parent%yzquad + children(j)%yzquad + children(j)%dip(2) * shift(3) + children(j)%dip(3) * shift(2) + children(j)%charge * shift(2) * shift(3)
          parent%zxquad = parent%zxquad + children(j)%zxquad + children(j)%dip(3) * shift(1) + children(j)%dip(1) * shift(3) + children(j)%charge * shift(3) * shift(1)


          ! current dipole moments
          parent%dipjx(1:3) = parent%dipjx(1:3) + children(j)%dipjx(1:3) + children(j)%monoj(1:3)* shift(1) !!!! CHECK FROM HERE
          parent%dipjy(1:3) = parent%dipjy(1:3) + children(j)%dipjy(1:3) + children(j)%monoj(1:3)* shift(2)
          parent%dipjz(1:3) = parent%dipjz(1:3) + children(j)%dipjz(1:3) + children(j)%monoj(1:3)* shift(3)

          !Current quadrupole moments

          parent%quadjx(1:3) = parent%quadjx(1:3) + children(j)%quadjx(1:3) + children(j)%dipjx(1:3) * shift(1) + half * children(j)%monoj(1:3) * shift(1)**2
          parent%quadjy(1:3) = parent%quadjy(1:3) + children(j)%quadjy(1:3) + children(j)%dipjy(1:3) * shift(2) + half * children(j)%monoj(1:3) * shift(2)**2
          parent%quadjz(1:3) = parent%quadjz(1:3) + children(j)%quadjz(1:3) + children(j)%dipjz(1:3) * shift(3) + half * children(j)%monoj(1:3) * shift(3)**2

          parent%quadjxy(1:3) = parent%quadjxy(1:3) + children(j)%quadjxy(1:3) + children(j)%dipjx(1:3) * shift(2) + children(j)%dipjy(1:3) * shift(1) + children(j)%monoj(1:3) * shift(1) * shift(2)
          parent%quadjyz(1:3) = parent%quadjyz(1:3) + children(j)%quadjyz(1:3) + children(j)%dipjy(1:3) * shift(3) + children(j)%dipjz(1:3) * shift(2) + children(j)%monoj(1:3) * shift(2) * shift(3)
          parent%quadjzx(1:3) = parent%quadjzx(1:3) + children(j)%quadjzx(1:3) + children(j)%dipjz(1:3) * shift(1) + children(j)%dipjx(1:3) * shift(3) + children(j)%monoj(1:3) * shift(3) * shift(1)

        end do

        parent%bmax = maxval(sqrt((parent%coc(1)-children(1:nchild)%coc(1))**2+(parent%coc(2)-children(1:nchild)%coc(2))**2+(parent%coc(3)-children(1:nchild)%coc(3))**2) + children(1:nchild)%bmax)
      end subroutine


      !>
      !> adds res2 to res1
      !>
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in)    :: res2

        res1%E    = res1%E    + res2%E
        res1%pot  = res1%pot  + res2%pot
!        res1%rho  = res1%rho  + res2%rho
        res1%A    = res1%A    + res2%A
        res1%dxA  = res1%dxA  + res2%dxA
        res1%dyA  = res1%dyA  + res2%dyA
        res1%B    = res1%B    + res2%B
        res1%J    = res1%J    + res2%J
        res1%Jirr = res1%Jirr + res2%Jirr
        
      end subroutine


      !>
      !> reads interaction-specific parameters from file
      !>
      subroutine calc_force_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section calc_force")
        read(filehandle, NML=calc_force)
      end subroutine


      !>
      !> writes interaction-specific parameters to file
      !>
      subroutine calc_force_write_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=calc_force)
      end subroutine


      !>
      !> computes derived parameters for calc force module
      !>
      subroutine calc_force_prepare()
        use treevars, only : me, MPI_COMM_lpepc
        use module_fmm_framework, only : fmm_framework_prepare
        use module_mirror_boxes, only : do_periodic
        implicit none

        if (do_periodic .and. include_far_field_if_periodic) then
          call fmm_framework_prepare(me, MPI_COMM_lpepc)
        end if
      end subroutine


      !>
      !> initializes static variables of calc force module that depend
      !> on particle data and might be reused on subsequent traversals
      !>
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
        if (do_periodic .and. include_far_field_if_periodic) then
          call fmm_framework_timestep(particles)
        end if
      end subroutine


      !>
      !> subroutine must return the estimated number of iteractions per particle
      !> for the current mac and/or parameters and the supplied total number of particles
      !>
      subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
        use module_pepc_types
        implicit none
        integer(kind_particle), intent(in) :: npart_total !< total number of particles
        integer(kind_node), intent(out) :: nintmax !< maximum number of interactions per particle

        real(kind_physics) :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

        ! Estimate of interaction list length - Hernquist expression
        ! applies for BH-MAC
        invnintmax = max(theta2 / (35._8*log(1._8*npart_total)) , 1._8/npart_total)
        nintmax    = int(1._8/invnintmax)
      end subroutine


      !>
      !> finalizes the calc force module at end of simulation
      !>
      subroutine calc_force_finalize()
        implicit none
        ! nothing to do here
      end subroutine calc_force_finalize


      !>
      !> generic Multipole Acceptance Criterion
      !>
      function mac(node, dist2, boxlength2)
        implicit none

        logical :: mac
        type(t_tree_node_interaction_data), intent(in) :: node
        real(kind_physics), intent(in) :: dist2
        real(kind_physics), intent(in) :: boxlength2

        select case (mac_select)
            case (0)
              ! Barnes-Hut-MAC
              mac = (theta2 * dist2 > boxlength2)
            case (1)
               ! Bmax-MAC
              mac = (theta2 * dist2 > min(node%bmax**2, 3.0 * boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
            case default
              ! N^2 code
              mac = .false.
        end select
      end function


      !>
      !> clears result in t_particle datatype - usually, this function does not need to be touched
      !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
      !> function cannot reside in module_interaction_specific that may not include module_pepc_types
      !>
      subroutine particleresults_clear(particles)
        use module_pepc_types
        implicit none
        type(t_particle), intent(inout) :: particles(:)

        particles(:)%results = EMPTY_PARTICLE_RESULTS
      end subroutine


      !>
      !> Force calculation wrapper.
      !> This function is thought for pre- and postprocessing of
      !> calculated fields, and for being able to call several
      !> (different) force calculation routines
      !>
      subroutine calc_force_per_interaction_with_self(particle, node, node_idx, delta, dist2, vbox)
        use module_pepc_types
        use treevars
        use module_coulomb_kernels
        implicit none

        type(t_tree_node_interaction_data), intent(in) :: node
        integer(kind_node), intent(in) :: node_idx
        type(t_particle), intent(inout) :: particle
        real(kind_physics), intent(in) :: vbox(3), delta(3), dist2
      end subroutine


      !>
      !> Force calculation wrapper.
      !> This function is thought for pre- and postprocessing of
      !> calculated fields, and for being able to call several
      !> (different) force calculation routines
      !>
      subroutine calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
        use module_pepc_types
        use treevars
        use module_coulomb_kernels
        use module_shortcut, only:zero
        implicit none

        type(t_tree_node_interaction_data), intent(in) :: node
        integer(kind_node), intent(in) :: node_idx
        type(t_particle), intent(inout) :: particle
        real(kind_physics), intent(in) :: vbox(3), delta(3), dist2

        real(kind_physics) :: exyz(3),elxyz(3), Bxyz(3), Jxyz(3), Jirrxyz(3), phic, Axyz(3), gradxA(1:3), gradyA(1:3)

              !  compute 2D3V-Coulomb/Darwin fields and potential
              ! It's a leaf, do direct summation with calc_force_darwin_3D_direct
              call calc_force_darwin_2D3V_direct(node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, eps2, elxyz(1:2), &
                        Axyz(1:3), Jxyz(1:3),  Jirrxyz(1:3) , Bxyz(1:3), phic, gradxA(1:3), gradyA(1:3) )
                        
              particle%results%A         = particle%results%A         + Axyz           
              particle%results%B         = particle%results%B         + Bxyz           
              particle%results%J         = particle%results%J         + Jxyz           
              particle%results%Jirr(1:2) = particle%results%Jirr(1:2) + Jirrxyz(1:2)   
              particle%results%E(1:2)    = particle%results%E(1:2)    + elxyz(1:2)   
              particle%results%dxA(1:3)  = particle%results%dxA(1:3)  + gradxA(1:3)
              particle%results%dyA(1:3)  = particle%results%dyA(1:3)  + gradyA(1:3)
              particle%results%pot       = particle%results%pot       + phic     
!              particle%results%rho       = particle%results%rho       + rhoc 

              particle%results%Jirr(3)   = zero
              particle%results%E(3)      = zero


      end subroutine

      !>
      !> Force calculation wrapper.
      !> This function is thought for pre- and postprocessing of
      !> calculated fields, and for being able to call several
      !> (different) force calculation routines
      !>
      subroutine calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
        use module_pepc_types
        use treevars
        use module_coulomb_kernels
        use module_shortcut, only:zero
        implicit none

        type(t_tree_node_interaction_data), intent(in) :: node
        integer(kind_node), intent(in) :: node_idx
        type(t_particle), intent(inout) :: particle
        real(kind_physics), intent(in) :: vbox(3), delta(3), dist2

        real(kind_physics) :: exyz(3),elxyz(3), Bxyz(3), Jxyz(3), Jirrxyz(3), phic, Axyz(3), gradxA(1:3), gradyA(1:3) 

            !  compute 2D3V-Coulomb/Darwin fields and potential
              call calc_force_darwin_2D3V(node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, eps2, elxyz(1:2), &
                            Axyz(1:3), Jxyz(1:3), Jirrxyz(1:3), Bxyz(1:3), phic, gradxA(1:3), gradyA(1:3) )
              particle%results%A         = particle%results%A         + Axyz            
              particle%results%B         = particle%results%B         + Bxyz           
              particle%results%J         = particle%results%J         + Jxyz             
              particle%results%Jirr(1:2) = particle%results%Jirr(1:2) + Jirrxyz(1:2)     
              particle%results%E(1:2)    = particle%results%E(1:2)    + elxyz(1:2)   
              particle%results%dxA(1:3)  = particle%results%dxA(1:3)  + gradxA(1:3)
              particle%results%dyA(1:3)  = particle%results%dyA(1:3)  + gradyA(1:3)       
              particle%results%pot       = particle%results%pot       + phic    

              particle%results%Jirr(3)   = zero
              particle%results%E(3)      = zero

              
      end subroutine


        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        subroutine calc_force_per_particle(particles)
          use treevars, only: num_threads
          use module_debug, only : pepc_status
          use module_pepc_types
          use treevars, only : me
          use module_fmm_framework
          use module_mirror_boxes
          implicit none

          type(t_particle), intent(inout) :: particles(:)
          real(kind_physics) :: e_lattice(3), phi_lattice
          integer(kind_particle) :: p

          call pepc_status('CALC FORCE PER PARTICLE')

          potfarfield  = 0.
          potnearfield = 0.

          if (do_periodic .and. include_far_field_if_periodic) then
!             if ((me==0) .and. (force_law .ne. 3)) write(*,*) "Warning: far-field lattice contribution is currently only supported for force_law==3"
!             !$ call omp_set_num_threads(num_threads)
!             !$OMP  PARALLEL DO DEFAULT(PRIVATE) SHARED(particles) SCHEDULE(RUNTIME) REDUCTION(+:potfarfield,potnearfield)
             do p=1,size(particles)
                call fmm_sum_lattice_force(particles(p)%x, e_lattice, phi_lattice)

                potfarfield  = potfarfield  + phi_lattice               * particles(p)%data%q
                potnearfield = potnearfield + particles(p)%results%pot  * particles(p)%data%q

                particles(p)%results%e     = particles(p)%results%e     + e_lattice
                particles(p)%results%pot   = particles(p)%results%pot   + phi_lattice
             end do
!             !$OMP  END PARALLEL DO
!             !$ call omp_set_num_threads(1)
          end if

          call pepc_status('CALC FORCE PER PARTICLE DONE')
        end subroutine calc_force_per_particle
end module module_interaction_specific
