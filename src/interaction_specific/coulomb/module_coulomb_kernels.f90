! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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
!> Encapsulates the low-level kernels for Coulomb- and similar interactions
!>
module module_coulomb_kernels
  use module_pepc_kinds
  use module_pepc_types
  implicit none
  save
  private

    ! shortcut notations
    real(kind_physics), parameter :: zero    =  0._kind_physics
    real(kind_physics), parameter :: one     =  1._kind_physics
    real(kind_physics), parameter :: two     =  2._kind_physics
    real(kind_physics), parameter :: three   =  3._kind_physics
    real(kind_physics), parameter :: four    =  4._kind_physics
    real(kind_physics), parameter :: five    =  5._kind_physics
    real(kind_physics), parameter :: six     =  6._kind_physics
    real(kind_physics), parameter :: eight   =  8._kind_physics
    real(kind_physics), parameter :: nine    =  9._kind_physics
    real(kind_physics), parameter :: half    =  0.5_kind_physics

    interface calc_force_coulomb_3D_direct
      module procedure calc_force_coulomb_3D_direct_s, calc_force_coulomb_3D_direct_v
    end interface
    interface calc_force_coulomb_3D
      module procedure calc_force_coulomb_3D_s, calc_force_coulomb_3D_v
    end interface
    public calc_force_coulomb_3D
    public calc_force_coulomb_3D_direct
    public calc_force_coulomb_2D
    public calc_force_coulomb_2D_direct
    public calc_force_LJ
    public calc_force_kelbg_3D_direct

  contains

    !>
    !> Calculates 3D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_s(t, d, dist2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi

      real(kind_physics) :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6, m2rd5, m5rd7, &
        pre1, pre2x, pre2y, pre2z, preQ1, preQ2, preQ3

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r  = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd = one/r
      rd2 = rd *rd
      rd3 = rd *rd2
      rd5 = rd3*rd2
      rd7 = rd5*rd2

      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz
      dx3 = dx*dx2
      dy3 = dy*dy2
      dz3 = dz*dz2

      fd1 = three*dx2*rd5 - rd3
      fd2 = three*dy2*rd5 - rd3
      fd3 = three*dz2*rd5 - rd3
      fd4 = three*dx*dy*rd5
      fd5 = three*dy*dz*rd5
      fd6 = three*dx*dz*rd5

      m2rd5 = two*rd5
      m5rd7 = five*rd7
      pre1 = m5rd7*dx*dy*dz
      pre2x = m5rd7*dx2 - rd5
      pre2y = m5rd7*dy2 - rd5
      pre2z = m5rd7*dz2 - rd5
      preQ1 = pre2x*t%quad(1)
      preQ2 = pre2y*t%quad(2)
      preQ3 = pre2z*t%quad(3)

      phi = t%charge*rd                                             &  !  monopole term
            - (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3         &  !  dipole
            + fd1*t%quad(1) + fd2*t%quad(2) + fd3*t%quad(3)         &  !  quadrupole
            + fd4*t%xyquad  + fd5*t%yzquad  + fd6*t%zxquad

      exyz(1) = t%charge*dx*rd3                                     &  ! monopole term
                - (fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3))      &  ! dipole term
                + three * (                                         &  ! quadrupole term
                   dx * (                                           &
                       ( pre2x - m2rd5 )*t%quad(1)                  &
                     + preQ2                                        &
                     + preQ3                                        &
                   )                                                &
                   + dy*pre2x*t%xyquad                              &
                   + dz*pre2x*t%zxquad                              &
                   + pre1*t%yzquad                                  &
                  )

      exyz(2) = t%charge*dy*rd3                                     &
                - (fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3))      &
                + three * (                                         &
                   dy * (                                           &
                       ( pre2y - m2rd5 )*t%quad(2)                  &
                     + preQ1                                        &
                     + preQ3                                        &
                   )                                                &
                   + dx*pre2y*t%xyquad                              &
                   + dz*pre2y*t%yzquad                              &
                   + pre1*t%zxquad                                  &
                  )

      exyz(3) = t%charge*dz*rd3                                     &
                - (fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1))      &
                + three * (                                         &
                   dz * (                                           &
                     + ( pre2z - m2rd5 )*t%quad(3)                  &
                     + preQ2                                        &
                     + preQ1                                        &
                   )                                                &
                   + dx*pre2z*t%zxquad                              &
                   + dy*pre2z*t%yzquad                              &
                   + pre1*t%xyquad                                  &
                  )
    end subroutine calc_force_coulomb_3D_s
    !>
    !> Calculates 3D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_v(np, t, d, dist2, exyz, phi, mask)
      implicit none

      integer(kind_particle), intent(in) :: np !< number of particles to process
      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3, np), dist2(np) !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3, np), phi(np)
      integer, intent(in) :: mask(np) !< vector mask to skip entries which should not be computed

      integer(kind_particle) :: i
      real(kind_physics) :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6, m2rd5, m5rd7, &
        pre1, pre2x, pre2y, pre2z, preQ1, preQ2, preQ3

      do i = 1, np
        dx = d(1, i)
        dy = d(2, i)
        dz = d(3, i)

        r  = sqrt(dist2(i)) ! eps2 is added in calling routine to have plummer intead of coulomb here
        rd = one/r
        rd2 = rd *rd
        rd3 = rd *rd2
        rd5 = rd3*rd2
        rd7 = rd5*rd2

        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        dx3 = dx*dx2
        dy3 = dy*dy2
        dz3 = dz*dz2

        fd1 = three*dx2*rd5 - rd3
        fd2 = three*dy2*rd5 - rd3
        fd3 = three*dz2*rd5 - rd3
        fd4 = three*dx*dy*rd5
        fd5 = three*dy*dz*rd5
        fd6 = three*dx*dz*rd5

        m2rd5 = two*rd5
        m5rd7 = five*rd7
        pre1 = m5rd7*dx*dy*dz
        pre2x = m5rd7*dx2 - rd5
        pre2y = m5rd7*dy2 - rd5
        pre2z = m5rd7*dz2 - rd5
        preQ1 = pre2x*t%quad(1)
        preQ2 = pre2y*t%quad(2)
        preQ3 = pre2z*t%quad(3)

        phi(i) = t%charge*rd                                          &  !  monopole term
              - (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3         &  !  dipole
              + fd1*t%quad(1) + fd2*t%quad(2) + fd3*t%quad(3)         &  !  quadrupole
              + fd4*t%xyquad  + fd5*t%yzquad  + fd6*t%zxquad

        exyz(1, i) = t%charge*dx*rd3                                  &  ! monopole term
                  - (fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3))      &  ! dipole term
                  + three * (                                         &  ! quadrupole term
                     dx * (                                           &
                         ( pre2x - m2rd5 )*t%quad(1)                  &
                       + preQ2                                        &
                       + preQ3                                        &
                     )                                                &
                     + dy*pre2x*t%xyquad                              &
                     + dz*pre2x*t%zxquad                              &
                     + pre1*t%yzquad                                  &
                    )

        exyz(2, i) = t%charge*dy*rd3                                  &
                  - (fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3))      &
                  + three * (                                         &
                     dy * (                                           &
                         ( pre2y - m2rd5 )*t%quad(2)                  &
                       + preQ1                                        &
                       + preQ3                                        &
                     )                                                &
                     + dx*pre2y*t%xyquad                              &
                     + dz*pre2y*t%yzquad                              &
                     + pre1*t%zxquad                                  &
                    )

        exyz(3, i) = t%charge*dz*rd3                                  &
                  - (fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1))      &
                  + three * (                                         &
                     dz * (                                           &
                       + ( pre2z - m2rd5 )*t%quad(3)                  &
                       + preQ2                                        &
                       + preQ1                                        &
                     )                                                &
                     + dx*pre2z*t%zxquad                              &
                     + dy*pre2z*t%yzquad                              &
                     + pre1*t%xyquad                                  &
                    )

        ! add mask to correct for cutoff or self-interaction
        phi(i) = mask(i) * phi(i)
        exyz(:, i) = mask(i) * exyz(:, i)
      end do
    end subroutine calc_force_coulomb_3D_v


    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D(t, d, d2, exy, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(2), d2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(1:2),phi

      real(kind_physics) :: dx,dy,rd2,rd4,rd6,dx2,dy2,dx3,dy3

      dx = d(1)
      dy = d(2)

      rd2 = one/d2
      rd4 = rd2*rd2
      rd6 = rd4*rd2

      dx2 = dx *dx
      dy2 = dy *dy
      dx3 = dx2*dx
      dy3 = dy2*dy

      phi = - half * t%charge * log(d2)            & !  monopole term
           - (dx * t%dip(1) + dy * t%dip(2)) * rd2 & !  dipole
           + t%quad(1) * (two * dx2 * rd4 - rd2)   & ! quadrupole
           + t%quad(2) * (two * dy2 * rd4 - rd2)   &
           + two * t%xyquad * dx * dy * rd4

      exy(1) = t%charge * dx * rd2                                 & ! monopole
           - t%dip(1) * (two * dx2 * rd4 - rd2)                    & ! dipole
           - t%dip(2) * two * dx * dy * rd4                        &
           + t%quad(1) * (eight * dx3 * rd6 - six * dx * rd4)      & ! quadrupole
           + t%quad(2) * (eight * dx * dy2 * rd6 - two * dx * rd4) &
           + t%xyquad * (eight * dx2 * dy * rd6 - two * dy * rd4)

      exy(2) = t%charge * dy * rd2                                 & ! monopole
           - t%dip(2) * (two * dy2 * rd4 - rd2)                    & ! dipole
           - t%dip(1) * two * dx * dy * rd4                        &
           + t%quad(2) * (eight * dy3 * rd6 - six * dy * rd4)      & ! quadrupole
           + t%quad(1) * (eight * dy * dx2 * rd6 - two * dy * rd4) &
           + t%xyquad * (eight * dy2 * dx * rd6 - two * dx * rd4)
    end subroutine calc_force_coulomb_2D


    !>
    !> CALC_FORCE_LJ
    !>
    !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
    !> shifted by the lattice vector vbox
    !> results are returned exyz, phi
    !>
    subroutine calc_force_LJ(t, d, r2, aii2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), r2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi
      real(kind_physics), intent(in) :: aii2
      real(kind_physics) :: flj, epsc2, aii2_r2,aii4_r4, r, fljrd

      ! epsc should be > a_ii to get evenly spaced ions
      epsc2 = 0.8_kind_physics * aii2

      ! Force is repulsive up to and just beyond aii
      if (r2 > epsc2) then
          aii2_r2 = aii2/r2
      else
          aii2_r2 = aii2/epsc2
      endif

      aii4_r4 = aii2_r2*aii2_r2

      flj = two*(aii4_r4*aii4_r4) - aii4_r4

      ! potential
      phi = zero

      !  forces
      r     = sqrt(r2)
      fljrd = flj/r

      exyz = d*fljrd
    end subroutine calc_force_LJ


    !>
    !> Calculates 3D Coulomb interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_direct_s(t, d, dist2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi

      real(kind_physics) :: rd,r,rd3charge

      r         = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd        = one/r
      rd3charge = t%charge*rd*rd*rd

      phi  = t%charge*rd
      exyz = rd3charge*d
    end subroutine calc_force_coulomb_3D_direct_s
    !>
    !> Calculates 3D Coulomb interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_direct_v(np, t, d, dist2, exyz, phi, mask)
      implicit none

      integer(kind_particle), intent(in) :: np !< number of particles to process
      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3, np), dist2(np) !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3, np), phi(np)
      integer, intent(in) :: mask(np) !< vector mask to skip entries which should not be computed

      integer(kind_particle) :: i
      real(kind_physics) :: rd,r,rd3charge
      real(kind_physics) :: charge

      charge = t%charge

      do i = 1, np
         r         = sqrt(dist2(i)) ! eps2 is added in calling routine to have plummer intead of coulomb here
         rd        = one/r
         rd3charge = charge*rd*rd*rd

         phi(i)       = mask(i)*charge*rd
         exyz(:,i)    = mask(i)*rd3charge*d(:,i)
      end do
    end subroutine calc_force_coulomb_3D_direct_v


    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D_direct(t, d, d2, exy, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(2), d2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(2), phi

      real(kind_physics) :: rd2charge

      phi       = - half*t%charge*log(d2)
      rd2charge = t%charge/d2
      exy       = rd2charge*d
    end subroutine calc_force_coulomb_2D_direct


    !>
    !> Calculates 3D Kelbg interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_kelbg_3D_direct(particle, t, d, dist2, kelbg_invsqrttemp, exyz, phi)
      implicit none

      type(t_particle), intent(inout) :: particle
      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi
      real(kind_physics), intent(in) :: kelbg_invsqrttemp
      real(kind_physics) :: rd,r,rd3
      real(kind_physics), parameter :: sqrtpi = sqrt(acos(-1.0_8))
      real(kind_physics) :: ome, rol, lambda, q, fprefac

      q = t%charge

      ! TODO: lambda must be adjusted depending on mass and temperature of interacting partners - currently it is fixed for electron-proton interactions
      if (particle%data%q * q < 0.) then
        ! e-i or i-e interaction
        lambda = 1.00027227_8 * kelbg_invsqrttemp
      else
        if ( q > 0. ) then
          ! i-i interaction
          lambda = 0.03300355_8 * kelbg_invsqrttemp
        else
          ! e-e interaction
          lambda = 1.41421356_8 * kelbg_invsqrttemp
        endif
      endif

      r   = sqrt(dist2)
      rd  = one / r
      rd3 = rd*rd*rd
      rol = r  / lambda        !< "r over lambda"
      ome = 1  - exp(-rol*rol) !< "one minus exp(stuff)"

      ! potential
      phi = q * rd  * (ome + sqrtpi*rol*(1-erf(rol)))
      !  forces
      fprefac = q * rd3 * ome
      exyz    = fprefac * d
    end subroutine calc_force_kelbg_3D_direct
end module
