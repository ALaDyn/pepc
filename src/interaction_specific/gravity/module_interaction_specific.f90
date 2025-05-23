! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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
   use module_pepc_kinds
   use module_pepc_types
   use module_interaction_specific_types
   implicit none
   save
   private

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer, public :: force_law = 3                          !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
   integer, public :: mac_select = 0                         !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
   logical, public :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
   real*8, public  :: theta2 = 0.36                          !< square of multipole opening angle
   real*8, public  :: eps2 = 0.0                             !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
   real*8, public  :: kelbg_invsqrttemp = 0.0                !< inverse square root of temperature for kelbg potential

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

      multipole = t_tree_node_interaction_data(particle_pos, particle%q, 0.)
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

      parent%charge = SUM(children(1:nchild)%charge)

      ! centre of charge
      parent%coc = [0., 0., 0.]

      ! use center-of-charge because we may divide by charge
      do j = 1, nchild
         parent%coc(1:3) = parent%coc(1:3) + (children(j)%coc(1:3) * children(j)%charge)
      end do

      parent%coc(1:3) = parent%coc(1:3) / parent%charge
      parent%bmax = maxval(sqrt((parent%coc(1) - children(1:nchild)%coc(1))**2 + (parent%coc(2) - children(1:nchild)%coc(2))**2 + (parent%coc(3) - children(1:nchild)%coc(3))**2) + children(1:nchild)%bmax)

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

      res1%e   = res1%e   + res2%e !&
      res1%pot = res1%pot + res2%pot !&
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
      read (filehandle, NML=calc_force_coulomb)

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

      write (filehandle, NML=calc_force_coulomb)

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> computes derived parameters for calc force module
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_force_prepare()
      implicit none
      !nothing to do here
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> initializes static variables of calc force module that depend
   !> on particle data and might be reused on subsequent traversals
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_force_after_grow(particles)
      implicit none
      type(t_particle), dimension(:), intent(in) :: particles

      ! nothing to be done here for now

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> subroutine must return the estimated number of iteractions per particle
   !> for the current mac and/or parameters and the supplied total number of particles
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
      use module_pepc_types
      implicit none
      integer(kind_particle), intent(in) :: npart_total !< total number of particles
      integer(kind_node), intent(out) :: nintmax !< maximum number of interactions per particle

      real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

      ! Estimate of interaction list length - Hernquist expression
      ! applies for BH-MAC
      invnintmax = max(theta2 / (35._8 * log(1._8 * npart_total)), 1._8 / npart_total)
      nintmax = int(1._8 / invnintmax)

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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mac(node, dist2, boxlength2)
      implicit none

      logical :: mac
      type(t_tree_node_interaction_data), intent(in) :: node
      real*8, intent(in) :: dist2
      real*8, intent(in) :: boxlength2

      select case (mac_select)
      case (0)
         ! Barnes-Hut-MAC
         mac = (theta2 * dist2 .gt. boxlength2)
      case (1)
         ! Bmax-MAC
         mac = (theta2 * dist2 .gt. min(node%bmax**2, 3.0 * boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_force_per_interaction_with_self(particle, node, key, delta, dist2, vbox)
      use module_pepc_types
      use treevars
      use module_coulomb_kernels
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: node
      integer(kind_key), intent(in) :: key
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: vbox(3), delta(3), dist2

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> Force calculation wrapper.
   !> This function is thought for pre- and postprocessing of
   !> calculated fields, and for being able to call several
   !> (different) force calculation routines
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_force_per_interaction_with_leaf(particle, node, key, delta, dist2, vbox)
      use module_pepc_types
      use treevars
      use module_coulomb_kernels
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: node
      integer(kind_key), intent(in) :: key
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: vbox(3), delta(3), dist2

      real*8 :: exyz(3), phic

      select case (force_law)
      case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
         call calc_force_coulomb_2D_direct(node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
         exyz(3) = 0.
      case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
         call calc_force_coulomb_3D_direct(node, delta, dist2 + eps2, exyz, phic)
      case (4)  ! LJ potential for quiet start
         call calc_force_LJ(node, delta, dist2, eps2, exyz, phic)
      case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction
         !  and Kelbg for particle-particle interaction
         call calc_force_kelbg_3D_direct(particle, node, delta, dist2, kelbg_invsqrttemp, exyz, phic)
      case default
         exyz = 0.
         phic = 0.
      end select

      particle%results%e   = particle%results%e   + exyz !&
      particle%results%pot = particle%results%pot + phic !&
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> Force calculation wrapper.
   !> This function is thought for pre- and postprocessing of
   !> calculated fields, and for being able to call several
   !> (different) force calculation routines
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_force_per_interaction_with_twig(particle, node, key, delta, dist2, vbox)
      use module_pepc_types
      use treevars
      use module_coulomb_kernels
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: node
      integer(kind_key), intent(in) :: key
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: vbox(3), delta(3), dist2

      real*8 :: exyz(3), phic

      select case (force_law)
      case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
         call calc_force_coulomb_2D(node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
         exyz(3) = 0.
      case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
         call calc_force_coulomb_3D(node, delta, dist2 + eps2, exyz, phic)
      case (4)  ! LJ potential for quiet start
         call calc_force_LJ(node, delta, dist2, eps2, exyz, phic)
      case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction
         !  and Kelbg for particle-particle interaction
         call calc_force_coulomb_3D(node, delta, dist2, exyz, phic)
      case default
         exyz = 0.
         phic = 0.
      end select

      particle%results%e   = particle%results%e   + exyz !&
      particle%results%pot = particle%results%pot + phic !&
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> Force calculation wrapper for contributions that only have
   !> to be added once per particle
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_force_per_particle(particles)
      implicit none
      type(t_particle), intent(in) :: particles(:)
      ! nothing to do here
   end subroutine calc_force_per_particle

end module module_interaction_specific
