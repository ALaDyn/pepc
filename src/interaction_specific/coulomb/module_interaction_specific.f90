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

!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
module module_interaction_specific
     use module_pepc_kinds, only: kind_physics
     use module_pepc_types
     use module_interaction_specific_types
     implicit none
     save
     private

      integer, public :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      logical, public :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
      real(kind_physics), public  :: theta2       = 0.36  !< square of multipole opening angle
      real(kind_physics), public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
      real(kind_physics), public  :: kelbg_invsqrttemp = 0.0 !< inverse square root of temperature for kelbg potential

      namelist /calc_force_coulomb/ force_law, mac_select, include_far_field_if_periodic, theta2, eps2, kelbg_invsqrttemp

      ! currently, all public functions in module_interaction_specific are obligatory
      public multipole_from_particle
      public shift_multipoles_up
      public results_add
      public calc_force_per_interaction_with_self
      public calc_force_per_interaction_with_leaf
      public calc_force_per_interaction_with_twig
      public calc_force_per_particle
      public mac
      public dual_mac
      public particleresults_clear
      public calc_force_read_parameters
      public calc_force_write_parameters
      public calc_force_finalize
      public calc_force_prepare
      public calc_force_after_grow
      public get_number_of_interactions_per_particle
      public multipole_to_local
      public shift_coefficients_down
      public evaluate_at_particle

      contains

      !>
      !> Computes multipole properties of a single particle
      !>
      subroutine multipole_from_particle(particle_pos, particle, multipole_center, multipole)
        implicit none
        real(kind_physics), intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        real(kind_physics), intent(out) :: multipole_center(3)
        type(t_multipole_moments), intent(out) :: multipole

        multipole_center = particle_pos
        multipole = t_multipole_moments(particle%q,   &
                                 abs(particle%q),  &
                                     (/0., 0., 0./), &
                                     (/0., 0., 0./), &
                                       0., 0., 0., 0. )
      end subroutine


      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      subroutine shift_multipoles_up(parent_center, parent, children_centers, children)
        implicit none
        real(kind_physics), intent(out) :: parent_center(3)
        type(t_multipole_moments), intent(out) :: parent
        real(kind_physics), intent(in) :: children_centers(:, :)
        type(t_multipole_moments), intent(in) :: children(:)

        integer :: nchild, j

        real(kind_physics) :: shift(1:3)
        real(kind_physics), parameter :: half = 0.5_kind_physics

        nchild = size(children)

        parent%charge     = SUM( children(1:nchild)%charge )
        parent%abs_charge = SUM( children(1:nchild)%abs_charge )

        ! centre of charge
        parent_center = [0., 0., 0.]

        if (parent%abs_charge .ne. 0.) then
          ! use center-of-charge because we may divide by abs_charge
          do j=1,nchild
            parent_center(1:3) = parent_center(1:3) + ( children_centers(1:3, j) * children(j)%abs_charge )
          end do

          parent_center(1:3) = parent_center(1:3) / parent%abs_charge
        else
          ! use geometric center
          do j=1,nchild
            parent_center(1:3) = parent_center(1:3) + children_centers(1:3, j)
          end do

          parent_center(1:3) = parent_center(1:3) / nchild
        endif

        ! multipole properties
        parent%dip    = [0., 0., 0.]
        parent%quad   = [0., 0., 0.]
        parent%xyquad = 0.
        parent%yzquad = 0.
        parent%zxquad = 0.

        do j=1,nchild
          shift(1:3) = parent_center(1:3) - children_centers(1:3, j)

          ! dipole moment
          parent%dip = parent%dip + children(j)%dip + children(j)%charge * shift(1:3)

          ! quadrupole moment
          parent%quad(1:3) = parent%quad(1:3) + children(j)%quad(1:3) + children(j)%dip(1:3) * shift(1:3) + half * children(j)%charge * shift(1:3)**2

          parent%xyquad = parent%xyquad + children(j)%xyquad + children(j)%dip(1) * shift(2) + children(j)%dip(2) * shift(1) + children(j)%charge * shift(1) * shift(2)
          parent%yzquad = parent%yzquad + children(j)%yzquad + children(j)%dip(2) * shift(3) + children(j)%dip(3) * shift(2) + children(j)%charge * shift(2) * shift(3)
          parent%zxquad = parent%zxquad + children(j)%zxquad + children(j)%dip(3) * shift(1) + children(j)%dip(1) * shift(3) + children(j)%charge * shift(3) * shift(1)
        end do

        parent%bmax = maxval(sqrt((parent_center(1)-children_centers(1, 1:nchild))**2+(parent_center(2)-children_centers(2, 1:nchild))**2+(parent_center(3)-children_centers(3, 1:nchild))**2) + children(1:nchild)%bmax)
      end subroutine

      !>
      !> Uses the M2L operator to convert-shift a set of multipole moments `m` along the separation vector `d`
      !> into a set of local coefficients `t`.
      !>
      subroutine multipole_to_local(d, d2, m, t)
        implicit none
        real(kind_physics), intent(in) :: d(3), d2
        type(t_multipole_moments), intent(in) :: m
        type(t_local_coefficients), intent(inout) :: t

        real(kind_physics) :: invr, invr2, invr3
        real(kind_physics) :: invr5t3, invr5t3dx, invr5t3dy, invr5t3dz
        real(kind_physics) :: invr7t15, invr7t15dx2, invr7t15dy2, invr7t15dz2
        real(kind_physics) :: dx, dy, dz
        real(kind_physics) :: kx, ky, kz
        real(kind_physics) :: kxx, kyy, kzz
        real(kind_physics) :: kxy, kyz, kzx
        real(kind_physics) :: kxxx, kyyy, kzzz
        real(kind_physics) :: kxxy, kxxz
        real(kind_physics) :: kyyx, kyyz
        real(kind_physics) :: kzzx, kzzy
        real(kind_physics) :: kxyz


        real(kind_physics), parameter :: one = 1._kind_physics
        real(kind_physics), parameter :: three = 3._kind_physics
        real(kind_physics), parameter :: five = 5._kind_physics

        dx = d(1)
        dy = d(2)
        dz = d(3)

        invr = one / sqrt(d2 + eps2)
        invr2 = invr * invr
        invr3 = invr * invr2

        kx = -invr3 * dx
        ky = -invr3 * dy
        kz = -invr3 * dz

        invr5t3 = invr2 * invr3 * three

        invr5t3dx = invr5t3 * dx
        invr5t3dy = invr5t3 * dy
        invr5t3dz = invr5t3 * dz

        kxx = dx * invr5t3dx - invr3
        kyy = dy * invr5t3dy - invr3
        kzz = dz * invr5t3dz - invr3

        kxy = invr5t3dx * dy
        kyz = invr5t3dy * dz
        kzx = invr5t3dz * dx

        invr7t15 = invr5t3 * invr2 * five
        invr7t15dx2 = invr7t15 * dx * dx
        invr7t15dy2 = invr7t15 * dy * dy
        invr7t15dz2 = invr7t15 * dz * dz

        kxxx = -invr7t15dx2 * dx + three * invr5t3dx
        kyyy = -invr7t15dy2 * dy + three * invr5t3dy
        kzzz = -invr7t15dz2 * dz + three * invr5t3dz

        kxxy = -invr7t15dx2 * dy + invr5t3dy
        kxxz = -invr7t15dx2 * dz + invr5t3dz

        kyyx = -invr7t15dy2 * dx + invr5t3dx
        kyyz = -invr7t15dy2 * dz + invr5t3dz

        kzzx = -invr7t15dz2 * dx + invr5t3dx
        kzzy = -invr7t15dz2 * dy + invr5t3dy

        kxyz = -invr7t15 * dx * dy * dz

        t%f = t%f + m%charge * invr &
          + m%dip(1) * kx + m%dip(2) * ky + m%dip(3) * kz &
          + m%quad(1) * kxx + m%quad(2) * kyy + m%quad(3) * kzz &
          + m%xyquad * kxy + m%yzquad * kyz + m%zxquad * kzx

        t%fx = t%fx + m%charge * kx &
          + m%dip(1) * kxx + m%dip(2) * kxy + m%dip(3) * kzx &
          + m%quad(1) * kxxx + m%quad(2) * kyyx + m%quad(3) * kzzx &
          + m%xyquad * kxxy + m%yzquad * kxyz + m%zxquad * kxxz
        t%fy = t%fy + m%charge * ky &
          + m%dip(1) * kxy + m%dip(2) * kyy + m%dip(3) * kyz &
          + m%quad(1) * kxxy + m%quad(2) * kyyy + m%quad(3) * kzzy &
          + m%xyquad * kyyx + m%yzquad * kyyz + m%zxquad * kxyz
        t%fz = t%fz + m%charge * kz &
          + m%dip(1) * kzx + m%dip(2) * kyz + m%dip(3) * kzz &
          + m%quad(1) * kxxz + m%quad(2) * kyyz + m%quad(3) * kzzz &
          + m%xyquad * kxyz + m%yzquad * kzzy + m%zxquad * kzzx

        t%fxx = t%fxx + m%charge * kxx &
          + m%dip(1) * kxxx + m%dip(2) * kxxy + m%dip(3) * kxxz
        t%fyy = t%fyy + m%charge * kyy &
          + m%dip(1) * kyyx + m%dip(2) * kyyy + m%dip(3) * kyyz
        t%fzz = t%fzz + m%charge * kzz &
          + m%dip(1) * kzzx + m%dip(2) * kzzy + m%dip(3) * kzzz
        t%fxy = t%fxy + m%charge * kxy &
          + m%dip(1) * kxxy + m%dip(2) * kyyx + m%dip(3) * kxyz
        t%fyz = t%fyz + m%charge * kyz &
          + m%dip(1) * kxyz + m%dip(2) * kyyz + m%dip(3) * kzzy
        t%fzx = t%fzx + m%charge * kzx &
          + m%dip(1) * kxxz + m%dip(2) * kxyz + m%dip(3) * kzzx

        t%fxxx = t%fxxx + m%charge * kxxx
        t%fyyy = t%fyyy + m%charge * kyyy
        t%fzzz = t%fzzz + m%charge * kzzz

        t%fxxy = t%fxxy + m%charge * kxxy
        t%fxxz = t%fxxz + m%charge * kxxz

        t%fyyx = t%fyyx + m%charge * kyyx
        t%fyyz = t%fyyz + m%charge * kyyz

        t%fzzx = t%fzzx + m%charge * kzzx
        t%fzzy = t%fzzy + m%charge * kzzy

        t%fxyz = t%fxyz + m%charge * kxyz

      end subroutine

      !>
      !> Uses the L2L operator to translate the coefficients `p` along the separation vector `d` into new coefficients `c`.
      !>
      subroutine shift_coefficients_down(d, p, c)
        implicit none
        real(kind_physics), intent(in) :: d(3)
        type(t_local_coefficients), intent(in) :: p
        type(t_local_coefficients), intent(inout) :: c

        real(kind_physics) :: dx, dy, dz, dx2, dy2, dz2, dxy, dyz, dzx

        real(kind_physics), parameter :: one = 1.0_kind_physics
        real(kind_physics), parameter :: two = 2.0_kind_physics
        real(kind_physics), parameter :: six = 6.0_kind_physics
        real(kind_physics), parameter :: half = one / two
        real(kind_physics), parameter :: sixth = one / six

        dx = d(1)
        dy = d(2)
        dz = d(3)

        dx2 = dx * dx
        dy2 = dy * dy
        dz2 = dz * dz
        dxy = dx * dy
        dyz = dy * dz
        dzx = dx * dz

        c%f = c%f + p%f &
          + dx * p%fx + dy * p%fy + dz * p%fz &
          + half * (dx2 * p%fxx + dy2 * p%fyy + dz2 * p%fzz) &
          + dxy * p%fxy + dyz * p%fyz + dzx * p%fzx &
          + sixth * (dx2 * dx * p%fxxx + dy2 * dy * p%fyyy + dz2 * dz * p%fzzz) &
          + half * ( dx2 * dy * p%fxxy + dx2 * dz * p%fxxz &
            + dy2 * dx * p%fyyx + dy2 * dz * p%fyyz &
            + dz2 * dx * p%fzzx + dz2 * dy * p%fzzy) &
          + dxy * dz * p%fxyz

        c%fx = c%fx + p%fx &
          + dx * p%fxx + dy * p%fxy + dz * p%fzx &
          + half * (dx2 * p%fxxx + dy2 * p%fyyx + dz2 * p%fzzx) &
          + dxy * p%fxxy + dzx * p%fxxz + dyz * p%fxyz
        c%fy = c%fy + p%fy &
          + dx * p%fxy + dy * p%fyy + dz * p%fyz &
          + half * (dx2 * p%fxxy + dy2 * p%fyyy + dz2 * p%fzzy) &
          + dxy * p%fyyx + dzx * p%fxyz + dyz * p%fyyz
        c%fz = c%fz + p%fz &
          + dx * p%fzx + dy * p%fyz + dz * p%fzz &
          + half * (dx2 * p%fxxz + dy2 * p%fyyz + dz2 * p%fzzz) &
          + dxy * p%fxyz + dzx * p%fzzx + dyz * p%fzzy

        c%fxx = c%fxx + p%fxx + dx * p%fxxx + dy * p%fxxy + dz * p%fxxz
        c%fyy = c%fyy + p%fyy + dx * p%fyyx + dy * p%fyyy + dz * p%fyyz
        c%fzz = c%fzz + p%fzz + dx * p%fzzx + dy * p%fzzy + dz * p%fzzz
        c%fxy = c%fxy + p%fxy + dx * p%fxxy + dy * p%fyyx + dz * p%fxyz
        c%fyz = c%fyz + p%fyz + dx * p%fxyz + dy * p%fyyz + dz * p%fzzy
        c%fzx = c%fzx + p%fzx + dx * p%fxxz + dy * p%fxyz + dz * p%fzzx

        c%fxxx = c%fxxx + p%fxxx
        c%fyyy = c%fyyy + p%fyyy
        c%fzzz = c%fzzz + p%fzzz

        c%fxxy = c%fxxy + p%fxxy
        c%fxxz = c%fxxz + p%fxxz

        c%fyyx = c%fyyx + p%fyyx
        c%fyyz = c%fyyz + p%fyyz

        c%fzzx = c%fzzx + p%fzzx
        c%fzzy = c%fzzy + p%fzzy

        c%fxyz = c%fxyz + p%fxyz
      end subroutine


      !>
      !> Evaluates local coefficients `c` into results `r` at their center of expansion.
      !>
      subroutine evaluate_at_particle(c, r)
        implicit none
        type(t_local_coefficients), intent(in) :: c
        type(t_particle_results), intent(inout) :: r

        r%pot = r%pot + c%f
        r%e(1) = r%e(1) - c%fx
        r%e(2) = r%e(2) - c%fy
        r%e(3) = r%e(3) - c%fz

      end subroutine


      !>
      !> adds res2 to res1
      !>
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in) :: res2

        res1%e    = res1%e    + res2%e
        res1%pot  = res1%pot  + res2%pot
      end subroutine


      !>
      !> reads interaction-specific parameters from file
      !>
      subroutine calc_force_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section calc_force_coulomb")
        read(filehandle, NML=calc_force_coulomb)
      end subroutine


      !>
      !> writes interaction-specific parameters to file
      !>
      subroutine calc_force_write_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=calc_force_coulomb)
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
        type(t_multipole_moments), intent(in) :: node
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
      !> Multipole Acceptance Criterion for dual tree traversal
      !>
      function dual_mac(d, d2, r1, m1, r2, m2)
        implicit none

        logical :: dual_mac
        real(kind_physics), intent(in) :: d(3), d2
        real(kind_physics), intent(in) :: r1, r2
        type(t_multipole_moments), intent(in) :: m1, m2

        real(kind_physics) :: sr, sr2

        select case (mac_select)
          case (0)
            ! Barnes-Hut-MAC
            sr = r1 + r2
            sr2 = sr * sr
            dual_mac = (theta2 * d2 > sr2)
          case (1)
            ! Bmax-MAC
            sr = m1%bmax + m2%bmax
            sr2 = sr * sr
            dual_mac = (theta2 * d2 > sr2)
          case default
            ! N^2 code
            dual_mac = .false.
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

        type(t_multipole_moments), intent(in) :: node
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
        implicit none

        type(t_multipole_moments), intent(in) :: node
        integer(kind_node), intent(in) :: node_idx
        type(t_particle), intent(inout) :: particle
        real(kind_physics), intent(in) :: vbox(3), delta(3), dist2

        real(kind_physics) :: exyz(3), phic

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
              ! It's a leaf, do direct summation with kelbg
              call calc_force_kelbg_3D_direct(particle, node, delta, dist2, kelbg_invsqrttemp, exyz, phic)
          case default
            exyz = 0.
            phic = 0.
        end select

        particle%results%e         = particle%results%e    + exyz
        particle%results%pot       = particle%results%pot  + phic
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
        implicit none

        type(t_multipole_moments), intent(in) :: node
        integer(kind_node), intent(in) :: node_idx
        type(t_particle), intent(inout) :: particle
        real(kind_physics), intent(in) :: vbox(3), delta(3), dist2

        real(kind_physics) :: exyz(3), phic

        select case (force_law)
          case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
              call calc_force_coulomb_2D(       node, delta(1:2), dot_product(delta(1:2), delta(1:2)) + eps2, exyz(1:2), phic)
              exyz(3) = 0.
          case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
              call calc_force_coulomb_3D(       node, delta, dist2 + eps2, exyz, phic)
          case (4)  ! LJ potential for quiet start
              call calc_force_LJ(node, delta, dist2, eps2, exyz, phic)
          case (5)  !  compute 3D-Coulomb fields and potential for particle-cluster interaction
                    !  and Kelbg for particle-particle interaction

              ! It's a twig, do ME with coulomb
              call calc_force_coulomb_3D(node, delta, dist2, exyz, phic)
          case default
            exyz = 0.
            phic = 0.
        end select

        particle%results%e         = particle%results%e    + exyz
        particle%results%pot       = particle%results%pot  + phic
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
