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

!>
!> Facilitates different code/interfaces for the interaction backend.
!>
submodule (module_interaction_specific) force !&
   implicit none

contains

   !>
   !> adds res2 to res1
   !>
   subroutine results_add(res1, res2)
      implicit none
      type(t_particle_results), intent(inout) :: res1
      type(t_particle_results), intent(in) :: res2

      res1%u = res1%u + res2%u
      res1%af = res1%af + res2%af
      res1%psi = res1%psi + res2%psi
      res1%div = res1%div + res2%div
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
      use mpi
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: node
      integer(kind_node), intent(in) :: node_idx
      type(t_particle), intent(inout) :: particle
      real(kind_physics), intent(in) :: vbox(3), delta(3), dist2

      integer :: ierr
      real(kind_physics) :: u(3), af(3), psi(3), div

      u = 0.
      af = 0.
      psi = 0.
      div = 0.

      select case (force_law)
      case (21)  !  use 2nd order Gaussian kernel, transposed scheme
         call calc_2nd_gaussian_transposed_direct(particle, node, delta, dist2, u, af, div)
      case (22)  !  use 2nd order algebraic kernel, transposed scheme
         call calc_2nd_algebraic_transposed_direct(particle, node, delta, dist2, u, af, psi, div)
      case (61)  ! use 6th order algebraic kernel, classical scheme
         call calc_6th_gaussian_transposed_direct(particle, node, delta, dist2, u, af, div) !TODO: 6xth order direct summation
      case (62)  ! use 6th order algebraic kernel, transposed scheme
         call calc_6th_algebraic_transposed_direct(particle, node, delta, dist2, u, af, div) !TODO: 6xth order direct summation
      case default
         write (*, *) 'What force law is this?', force_law
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end select

      particle%results%u(1:3) = particle%results%u(1:3) - u(1:3)
      particle%results%af(1:3) = particle%results%af(1:3) + af(1:3)
      particle%results%psi(1:3) = particle%results%psi(1:3) + psi(1:3)
      particle%results%div = particle%results%div + div
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
      use mpi
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: node
      integer(kind_node), intent(in) :: node_idx
      type(t_particle), intent(inout) :: particle
      real(kind_physics), intent(in) :: vbox(3), delta(3), dist2

      integer :: ierr
      real(kind_physics) :: u(3), af(3), psi(3), div

      u = 0.
      af = 0.
      psi = 0.
      div = 0.

      select case (force_law)
      case (21)  !  use 2nd order Gaussian kernel, transposed scheme
         ! TODO: ME 2nd order classical scheme
         write (*, *) 'ME not implemented for Gaussian kernels, aborting ...'
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      case (22)  !  use 2nd order algebraic kernel, transposed scheme
         call calc_2nd_algebraic_transposed(particle, node, delta, dist2, u, af, psi)
      case (61)  ! use 6th order algebraic kernel, classical scheme
         ! TODO: ME 6th order classical scheme
         write (*, *) 'ME not implemented for Gaussian kernels, aborting ...'
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      case (62)  ! use 6th order algebraic kernel, transposed scheme
         call calc_6th_algebraic_transposed(particle, node, delta, dist2, u, af)
      case default
         write (*, *) 'What force law is this?', force_law
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end select

      !&<
      particle%results%u(1:3)   = particle%results%u(1:3)   -   u(1:3)
      particle%results%af(1:3)  = particle%results%af(1:3)  +  af(1:3)
      particle%results%psi(1:3) = particle%results%psi(1:3) + psi(1:3)
      particle%results%div      = particle%results%div      + div
      !&>
   end subroutine

   !>
   !> Calculates 3D 2nd order condensed algebraic kernel interaction
   !> of particle p with tree node t, results are returned in u and af
   !>
   subroutine calc_2nd_algebraic_transposed(particle, t, d, dist2, u, af, psi)
      use module_pepc_types
      use module_interaction_specific_types
      implicit none

      type(t_particle), intent(in) :: particle
      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  u(1:3), af(1:3), psi(1:3)

      integer :: i1, i2, i3 !< helper variables for the tensor structures

      real(kind_physics) :: dx, dy, dz !< temp variables for distance
      real(kind_physics) :: Gc15, Gc25, Gc35, Gc45, Gc55, MPa1, DPa1, DPa2, QPa1, QPa2 !< prefactors for the multipole expansion
      real(kind_physics) :: Gc(0:4)
      real(kind_physics), dimension(3) :: vort !< temp variables for vorticity (or better: alpha)

      ! tensors allow nice and short code and better comparison with (my!) theory
      real(kind_physics), dimension(3) :: m0, CP0 !< data structures for the monopole moments
      real(kind_physics), dimension(3, 3) :: m1, CP1 !< data structures for the dipole moments
      real(kind_physics), dimension(3, 3, 3) :: m2, CP2 !< data structures for the quadrupole moments

      dx = d(1)
      dy = d(2)
      dz = d(3)

      vort = [particle%data%alpha(1), particle%data%alpha(2), particle%data%alpha(3)]  ! need particle`s vorticity for cross-product here

      m0 = [t%chargex, t%chargey, t%chargez]          ! monopole moment tensor
      CP0 = cross_prod(m0, vort)                      ! cross-product for 1st expansion term

      m1(:, 1) = [t%xdip1, t%xdip2, t%xdip3]          ! dipole moment tensor
      m1(:, 2) = [t%ydip1, t%ydip2, t%ydip3]
      m1(:, 3) = [t%zdip1, t%zdip2, t%zdip3]
      CP1(:, 1) = cross_prod(m1(:, 1), vort)          ! cross-product for 2nd expansion term
      CP1(:, 2) = cross_prod(m1(:, 2), vort)
      CP1(:, 3) = cross_prod(m1(:, 3), vort)

      m2(:, 1, 1) = [t%xxquad1, t%xxquad2, t%xxquad3] ! quadrupole moment tensor
      m2(:, 2, 1) = [t%xyquad1, t%xyquad2, t%xyquad3]
      m2(:, 3, 1) = [t%xzquad1, t%xzquad2, t%xzquad3]
      m2(:, 1, 2) = [t%xyquad1, t%xyquad2, t%xyquad3]
      m2(:, 2, 2) = [t%yyquad1, t%yyquad2, t%yyquad3]
      m2(:, 3, 2) = [t%yzquad1, t%yzquad2, t%yzquad3]
      m2(:, 1, 3) = [t%xzquad1, t%xzquad2, t%xzquad3]
      m2(:, 2, 3) = [t%yzquad1, t%yzquad2, t%yzquad3]
      m2(:, 3, 3) = [t%zzquad1, t%zzquad2, t%zzquad3]
      CP2(:, 1, 1) = cross_prod(m2(:, 1, 1), vort)     ! cross-product for 3rd expansion term
      CP2(:, 2, 1) = cross_prod(m2(:, 2, 1), vort)
      CP2(:, 3, 1) = cross_prod(m2(:, 3, 1), vort)
      CP2(:, 1, 2) = cross_prod(m2(:, 1, 2), vort)
      CP2(:, 2, 2) = cross_prod(m2(:, 2, 2), vort)
      CP2(:, 3, 2) = cross_prod(m2(:, 3, 2), vort)
      CP2(:, 1, 3) = cross_prod(m2(:, 1, 3), vort)
      CP2(:, 2, 3) = cross_prod(m2(:, 2, 3), vort)
      CP2(:, 3, 3) = cross_prod(m2(:, 3, 3), vort)

      ! precompute kernel function evaluations of various order
      Gc = G_core([dist2, dist2, dist2, dist2, dist2], & !&
                  [sig2,  sig2, sig2, sig2, sig2], & !&
                  [1.5D0, 2.5D0, 3.5D0, 4.5D0, 5.5d0])   !&
      Gc15 = Gc(0)
      Gc25 = Gc(1)
      Gc35 = Gc(2)
      Gc45 = Gc(3)
      Gc55 = Gc(4)

      MPa1 = 3.0D00 * Gc35 * dot_product(d, CP0)   ! monopole prefactor for af
      DPa1 = 15.0D00 * Gc45 * sum((/(sum((/(CP1(i2, i1) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/))  ! dipole prefators for af
      DPa2 = (CP1(1, 1) + CP1(2, 2) + CP1(3, 3))
      QPa1 = 52.5D00 * Gc55 * sum((/(sum((/(sum((/(CP2(i3, i2, i1) * d(i3), i3=1, 3)/)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/)) ! quadrupole prefactors for af
      QPa2 = dot_product(d, CP2(1, 1, :) + CP2(2, 2, :) + CP2(3, 3, :) + CP2(1, :, 1) + CP2(2, :, 2) + CP2(3, :, 3) + CP2(:, 1, 1) + CP2(:, 2, 2) + CP2(:, 3, 3))

      u(1) = Gc25 * (dy * m0(3) - dz * m0(2)) &                                                                                ! MONOPOLE
             + 3.0D00 * Gc35 * sum((/((m1(3, i1) * dy - m1(2, i1) * dz) * d(i1), i1=1, 3)/)) - Gc25 * (m1(3, 2) - m1(2, 3)) &  ! DIPOLE
             - 1.5D00 * Gc35 * (sum((/(m2(3, i1, i1) * dy - m2(2, i1, i1) * dz, i1=1, 3)/)) + 2.0 * sum((/((m2(3, i1, 2) - m2(2, i1, 3)) * d(i1), i1=1, 3)/))) &
             + 7.5D00 * Gc45 * sum((/(sum((/((dy * m2(3, i2, i1) - dz * m2(2, i2, i1)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/))! QUADRUPOLE

      u(2) = Gc25 * (dz * m0(1) - dx * m0(3)) &                                                                                ! MONOPOLE
             + 3.0D00 * Gc35 * sum((/((m1(1, i1) * dz - m1(3, i1) * dx) * d(i1), i1=1, 3)/)) - Gc25 * (m1(1, 3) - m1(3, 1)) &  ! DIPOLE
             - 1.5D00 * Gc35 * (sum((/(m2(1, i1, i1) * dz - m2(3, i1, i1) * dx, i1=1, 3)/)) + 2.0 * sum((/((m2(1, i1, 3) - m2(3, i1, 1)) * d(i1), i1=1, 3)/))) &
             + 7.5D00 * Gc45 * sum((/(sum((/((dz * m2(1, i2, i1) - dx * m2(3, i2, i1)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/))! QUADRUPOLE

      u(3) = Gc25 * (dx * m0(2) - dy * m0(1)) &                                                                                ! MONOPOLE
             + 3.0D00 * Gc35 * sum((/((m1(2, i1) * dx - m1(1, i1) * dy) * d(i1), i1=1, 3)/)) - Gc25 * (m1(2, 1) - m1(1, 2)) &  ! DIPOLE
             - 1.5D00 * Gc35 * (sum((/(m2(2, i1, i1) * dx - m2(1, i1, i1) * dy, i1=1, 3)/)) + 2.0 * sum((/((m2(2, i1, 1) - m2(1, i1, 2)) * d(i1), i1=1, 3)/))) &
             + 7.5D00 * Gc45 * sum((/(sum((/((dx * m2(2, i2, i1) - dy * m2(1, i2, i1)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/))! QUADRUPOLE

      af(1) = Mpa1 * dx - Gc25 * CP0(1) &                                                                                      ! MONOPOLE
              + DPa1 * dx - 3.0D00 * Gc35 * (dot_product(CP1(:, 1) + CP1(1, :), d) + dx * DPa2) &                              ! DIPOLE
              + QPa1 * dx &
              - 7.5D00 * Gc45 * (sum((/(sum((/((CP2(i2, i1, 1) + CP2(i2, 1, i1) + CP2(1, i2, i1)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/)) + dx * QPa2) &
              + 1.5D00 * Gc35 * (CP2(1, 1, 1) + CP2(2, 2, 1) + CP2(3, 3, 1) + CP2(1, 1, 1) + CP2(2, 1, 2) + CP2(3, 1, 3) + CP2(1, 1, 1) + CP2(1, 2, 2) + CP2(1, 3, 3)) ! QUADRUPOLE

      af(2) = Mpa1 * dy - Gc25 * CP0(2) &                                                                                      ! MONOPOLE
              + DPa1 * dy - 3.0D00 * Gc35 * (dot_product(CP1(:, 2) + CP1(2, :), d) + dy * DPa2) &                              ! DIPOLE
              + QPa1 * dy &
              - 7.5D00 * Gc45 * (sum((/(sum((/((CP2(i2, i1, 2) + CP2(i2, 2, i1) + CP2(2, i2, i1)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/)) + dy * QPa2) &
              + 1.5D00 * Gc35 * (CP2(1, 1, 2) + CP2(2, 2, 2) + CP2(3, 3, 2) + CP2(1, 2, 1) + CP2(2, 2, 2) + CP2(3, 2, 3) + CP2(2, 1, 1) + CP2(2, 2, 2) + CP2(2, 3, 3)) ! QUADRUPOLE

      af(3) = Mpa1 * dz - Gc25 * CP0(3) &                                                                                      ! MONOPOLE
              + DPa1 * dz - 3.0D00 * Gc35 * (dot_product(CP1(:, 3) + CP1(3, :), d) + dz * DPa2) &                              ! DIPOLE
              + QPa1 * dz &
              - 7.5D00 * Gc45 * (sum((/(sum((/((CP2(i2, i1, 3) + CP2(i2, 3, i1) + CP2(3, i2, i1)) * d(i2), i2=1, 3)/)) * d(i1), i1=1, 3)/)) + dz * QPa2) &
              + 1.5D00 * Gc35 * (CP2(1, 1, 3) + CP2(2, 2, 3) + CP2(3, 3, 3) + CP2(1, 3, 1) + CP2(2, 3, 2) + CP2(3, 3, 3) + CP2(3, 1, 1) + CP2(3, 2, 2) + CP2(3, 3, 3)) ! QUADRUPOLE

      psi(1) = Gc15 * m0(1) & ! MONOPOLE
               + Gc25 * m1(1, 1) * dx + Gc25 * m1(1, 2) * dy + Gc25 * m1(1, 3) * dz & ! DIPOLE
               + 1.5D0 * Gc35 * (sum((/(sum((/(m2(1, i2, i1) * d(i1), i1=1, 3)/)) * d(i2), i2=1, 3)/))) &
               - 0.5d0 * Gc25 * (sum((/(m2(1, i1, i1), i1=1, 3)/)))                                         ! QUADRUPOLE

      psi(2) = Gc15 * m0(2) & ! MONOPOLE
               + Gc25 * m1(2, 1) * dx + Gc25 * m1(2, 2) * dy + Gc25 * m1(2, 3) * dz & ! DIPOLE
               + 1.5D0 * Gc35 * (sum((/(sum((/(m2(2, i2, i1) * d(i1), i1=1, 3)/)) * d(i2), i2=1, 3)/))) &
               - 0.5d0 * Gc25 * (sum((/(m2(2, i1, i1), i1=1, 3)/)))                                         ! QUADRUPOLE

      psi(3) = Gc15 * m0(3) & ! MONOPOLE
               + Gc25 * m1(3, 1) * dx + Gc25 * m1(3, 2) * dy + Gc25 * m1(3, 3) * dz & ! DIPOLE
               + 1.5D0 * Gc35 * (sum((/(sum((/(m2(3, i2, i1) * d(i1), i1=1, 3)/)) * d(i2), i2=1, 3)/))) &
               - 0.5d0 * Gc25 * (sum((/(m2(3, i1, i1), i1=1, 3)/)))                                         ! QUADRUPOLE

   end subroutine calc_2nd_algebraic_transposed

   !>
   !> Calculates 3D 2nd order algebraic kernel interaction, transposed scheme
   !> of particle p with tree *particle*, results are returned in u, psi and af
   !>
   subroutine calc_2nd_algebraic_transposed_direct(particle, t, d, dist2, u, af, psi, div)
      use module_pepc_types
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t
      type(t_particle), intent(inout) :: particle
      real(kind_physics), intent(in) :: d(3), dist2
      real(kind_physics), intent(out) :: u(1:3), af(1:3), psi(1:3), div

      real(kind_physics), dimension(3) :: m0, CP0 !< data structures for the monopole moments
      real(kind_physics) :: Gc15, Gc25, MPa1, nom, nom45, nom35, nom25, nom15
      real(kind_physics), dimension(3) :: vort !< temp variables for vorticity (or better: alpha)

      m0 = [t%chargex, t%chargey, t%chargez]    ! monopole moment tensor
      vort = [particle%data%alpha(1), particle%data%alpha(2), particle%data%alpha(3)]  ! need particle`s vorticity for cross-product here
      CP0 = cross_prod(m0, vort)                ! cross-product for 1st expansion term

      nom = dist2 + sig2
      nom45 = nom**(-4.5)
      nom35 = nom45 * nom
      nom25 = nom35 * nom
      nom15 = nom25 * nom

      Gc25 = (dist2 + 2.5 * sig2) * nom25
      Gc15 = (dist2 + 1.5 * sig2) * nom15

      MPa1 = 3.0 * (dist2 + 3.5 * sig2) * nom35 * dot_product(d, CP0)

      u = Gc25 * cross_prod(d, m0) ! monopole

      af = Mpa1 * d - Gc25 * CP0   ! monopole

      ! Stream function = sum_p G( d ) alpha_p(t)
      psi = Gc15 * m0

      !The divergence of omega = div ( sum_p Zeta( d ) alpha_p(t) )
      !                        = sum_p alpha_p(t) \dot grad_Zeta( d )
      div = -52.5 * nom45 * sig2**2 * dot_product(d, m0)
   end subroutine calc_2nd_algebraic_transposed_direct

end submodule force
