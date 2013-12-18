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

!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
module module_interaction_specific
     use module_pepc_types
     use module_interaction_specific_types
     implicit none
     save
     private

      integer, public :: force_law    = 5      !< 5=NN-list "interaction"
      integer, public :: mac_select   = 1      !< selector for multipole acceptance criterion, 1 = NN-MAC

      namelist /calc_force_nearestneighbour/ force_law, mac_select

      ! currently, all public functions in module_interaction_specific are obligatory
      public multipole_from_particle
      public shift_multipoles_up
      public results_add
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
      public pack_particle_list
      public unpack_particle_list

      contains

      !>
      !> Computes multipole properties of a single particle
      !>
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        use module_spacefilling, only: veccoord_to_key_lastlevel
        implicit none
        real*8, intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole

        ! use velocity (v) at same time step as coordinate, not v_minus_half
        multipole = t_tree_node_interaction_data(particle_pos, particle%particle_id)
      end subroutine


      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: i

        parent%coc         = [0._8, 0._8, 0._8]
        parent%particle_id = 0_8

        do i=1,size(children)
          parent%coc = parent%coc + children(i)%coc
        end do

        parent%coc = parent%coc / size(children)

      end subroutine


      !>
      !> adds res2 to res1
      !>
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in) :: res2

        ! TODONN: ist das wirklich alles?
        res1 = res2
      end subroutine


      !>
      !> reads interaction-specific parameters from file
      !>
      subroutine calc_force_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section calc_force_nearestneighbour")
        read(filehandle, NML=calc_force_nearestneighbour)
      end subroutine


      !>
      !> writes interaction-specific parameters to file
      !>
      subroutine calc_force_write_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=calc_force_nearestneighbour)
      end subroutine


      !>
      !> computes derived parameters for calc force module
      !>
      subroutine calc_force_prepare()
        implicit none
        ! nothing to do here
      end subroutine


      !>
      !> initializes static variables of calc force module that depend
      !> on particle data and might be reused on subsequent traversals
      !>
      subroutine calc_force_after_grow(particles)
        use module_pepc_types
        implicit none
        type(t_particle), dimension(:), intent(in) :: particles

        ! nothing to be done here for now
      end subroutine


      !>
      !> subroutine must return the estimated number of iteractions per particle
      !> for the current mac and/or parameters and the supplied total number of particles
      !>
      subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
        implicit none
        integer(kind_particle), intent(in) :: npart_total !< total number of particles
        integer(kind_node), intent(out) :: nintmax !< maximum number of interactions per particle

        real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

        real*8, parameter :: theta2 = 0.3**2 !TODO: this function must be adapted to the actual MAC

        ! Estimate of interaction list length - Hernquist expression
        ! applies for BH-MAC
        invnintmax = max(theta2 / (35.*log(1.*npart_total)) , 1._8/npart_total)
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
      function mac(particle, node, dist2, boxlength2)
        implicit none

        logical :: mac
        type(t_tree_node_interaction_data), intent(in) :: node
        type(t_particle), intent(in) :: particle
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2

        select case (mac_select)
            case (0)
              ! Barnes-Hut-MAC
              ! mac = (theta2 * dist2 > boxlength2)
            case (1)
              ! NN-MAC: we may "interact" with the node if it is further away than maxdist2 --> this leads to the node *not* being put onto the NN-list (strange, i know)
              ! first line: original formulation, last line: after transition to formulation with only one square root
              ! mac = sqrt(dist2) - sqrt(3.* boxlength2)  >  sqrt(results%maxdist2)                                ! + sqrt(3.*boxlength2)
              !     = sqrt(dist2)                         >  sqrt(results%maxdist2) + sqrt(3.*boxlength2)          ! ^2
              !     =      dist2                          > (sqrt(results%maxdist2) + sqrt(3.*boxlength2))**2
              !     =      dist2                          > results%maxdist2 + 2.*sqrt( 3.*results%maxdist2*boxlength2) + 3.*boxlength2
                mac =      dist2                          > particle%results%maxdist2 +    sqrt(12.*particle%results%maxdist2*boxlength2) + 3.*boxlength2
              ! TODO NN: this estimation should be evaluated without (!!) any square roots for performance reasons (which does not seem to be trivial)
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
         use module_pepc_types, only: t_particle
         implicit none
         type(t_particle), intent(inout) :: particles(:)

         integer(kind_particle) :: i, j

         do i=1,size(particles, kind=kind(i))
            particles(i)%results%maxdist2           = huge(0._8)
            do j=1, size(particles(i)%results%neighbour_nodes)
              particles(i)%results%neighbour_nodes(j)%p => null()
            end do
            particles(i)%results%maxidx             = 1
         end do
       end subroutine particleresults_clear


        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        subroutine calc_force_per_interaction_with_leaf(delta, dist2, particle_pack, node_data)
          use module_pepc_types
          use module_debug
          use treevars
          implicit none

          real*8, intent(in) :: delta(:,:)
          real*8, intent(in) :: dist2(:)
          type(t_particle_pack), intent(inout) :: particle_pack
          type(t_tree_node_interaction_data), target, intent(in) :: node_data

          select case (force_law)
            case (5)
                call update_nn_list(particle_pack, node_data, delta, dist2)
            case default
               DEBUG_ERROR(*, "value of force_law is not allowed in calc_force_per_interaction:", force_law)
          end select
        end subroutine


        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        subroutine calc_force_per_interaction_with_twig(delta, dist2, particle_pack, node_data)
          use module_pepc_types
          use module_debug
          use treevars
          implicit none

          real*8, intent(in) :: delta(:,:)
          real*8, intent(in) :: dist2(:)
          type(t_particle_pack), intent(inout) :: particle_pack
          type(t_tree_node_interaction_data), target, intent(in) :: node_data

          select case (force_law)
            case (5)
                call update_nn_list(particle_pack, node_data, delta, dist2)
            case default
               DEBUG_ERROR(*, "value of force_law is not allowed in calc_force_per_interaction:", force_law)
          end select
        end subroutine


        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        subroutine calc_force_per_particle(particles)
          use module_interaction_specific_types
          implicit none

          type(t_particle), intent(inout) :: particles(:)

          ! currently nothing to do here
        end subroutine calc_force_per_particle


        subroutine update_nn_list(particle_pack, node_data, delta, dist2)
          use module_pepc_types
          use treevars
          implicit none
          include 'mpif.h'

          real*8, intent(in) :: delta(:,:)
          real*8, intent(in) :: dist2(:)
          type(t_particle_pack), intent(inout) :: particle_pack
          type(t_tree_node_interaction_data), target, intent(in) :: node_data

          integer :: tmp(1), p

          do p=1,size(dist2)
            if (dist2(p) < particle_pack%maxdist2(p)) then
              ! add node to NN_list
              particle_pack%neighbour_nodes(p, particle_pack%maxidx(p))%p => node_data
              particle_pack%dist2(p, particle_pack%maxidx(p))           = dist2(p)
              particle_pack%dist_vector(p, :,particle_pack%maxidx(p))   = delta(p,:)
              tmp = maxloc(particle_pack%dist2(p, 1:num_neighbour_particles)) ! this is really ugly, but maxloc returns a 1-by-1 vector instead of the expected scalar
              particle_pack%maxidx(p)   = tmp(1)
              particle_pack%maxdist2(p) = particle_pack%dist2(p, particle_pack%maxidx(p))
            else
              ! node is further away than farest particle in nn-list --> can be ignored
            endif
          end do
        end subroutine update_nn_list


        subroutine pack_particle_list(particles, packed)
          use module_pepc_types, only: t_particle, kind_particle
          implicit none

          type(t_particle), intent(in) :: particles(:)
          type(t_particle_pack), intent(inout) :: packed

          integer(kind_particle) :: ip, np

          np = size(particles, kind = kind_particle)

          allocate(packed%maxdist2(np), &
                   packed%maxidx(np), &
                   packed%neighbour_nodes(np, max_neighbour_particles), &
                   packed%dist2(np, max_neighbour_particles), &
                   packed%dist_vector(np, 3, max_neighbour_particles))

          do ip = 1, np
            packed%maxdist2(ip) = particles(ip)%results%maxdist2
            packed%maxidx(ip)   = particles(ip)%results%maxidx
            packed%neighbour_nodes(ip,:) = particles(ip)%results%neighbour_nodes
            packed%dist2(ip,:)  = particles(ip)%results%dist2
            packed%dist_vector(ip,:,:)   = particles(ip)%results%dist_vector(:,:)
          end do
        end subroutine pack_particle_list


        subroutine unpack_particle_list(packed, particles)
          use module_pepc_types, only: t_particle, kind_particle
          use module_debug
          implicit none

          type(t_particle_pack), intent(inout) :: packed
          type(t_particle), intent(inout) :: particles(:)

          integer(kind_particle) :: ip, np

          np = size(particles, kind = kind_particle)
          DEBUG_ASSERT(np == size(packed%maxdist2, kind = kind_particle))

          do ip = 1, np
            particles(ip)%results%maxdist2 = packed%maxdist2(ip)
            particles(ip)%results%maxidx = packed%maxidx(ip)
            particles(ip)%results%neighbour_nodes = packed%neighbour_nodes(ip,:)
            particles(ip)%results%dist2 = packed%dist2(ip,:)
            particles(ip)%results%dist_vector(:,:) = packed%dist_vector(ip,:,:)
          end do

          deallocate(packed%maxdist2, &
                     packed%maxidx, &
                     packed%neighbour_nodes, &
                     packed%dist2, &
                     packed%dist_vector)
        end subroutine unpack_particle_list
end module module_interaction_specific
