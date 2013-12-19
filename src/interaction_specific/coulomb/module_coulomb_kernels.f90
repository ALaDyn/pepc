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
!> Encapsulates the low-level kernels for Coulomb- and similar interactions
!>
module module_coulomb_kernels
  use module_pepc_types
  use module_interaction_specific_types
  #ifndef NO_SPATIAL_INTERACTION_CUTOFF
  use module_mirror_boxes, only: spatial_interaction_cutoff
  #endif
  implicit none
  save
  private

    integer, parameter :: kfp       =  8 ! numeric precision (kind value)
    ! shortcut notations
    real(kfp), parameter :: zero    =  0._kfp
    real(kfp), parameter :: one     =  1._kfp
    real(kfp), parameter :: two     =  2._kfp
    real(kfp), parameter :: three   =  3._kfp
    real(kfp), parameter :: four    =  4._kfp
    real(kfp), parameter :: five    =  5._kfp
    real(kfp), parameter :: eight   =  8._kfp
    real(kfp), parameter :: nine    =  9._kfp
    real(kfp), parameter :: half    =  0.5_kfp

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
    subroutine calc_force_coulomb_3D(np, delta, dist2, ex, ey, ez, pot, t, eps2)
      implicit none

      integer(kind_particle), intent(in) :: np
      real(kfp), intent(in) :: delta(np,3)
      real(kfp), intent(in) :: dist2(np)
      real(kfp), intent(inout) :: ex(np), ey(np), ez(np), pot(np)
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

      real(kfp) :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
      integer(kind_particle) :: ip

      do ip = 1, np
        dx = delta(ip, 1)
        dy = delta(ip, 2)
        dz = delta(ip, 3)

        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (dx >= spatial_interaction_cutoff(1) .or. dy >= spatial_interaction_cutoff(2) .or. dz >= spatial_interaction_cutoff(3)) &
          cycle
        #endif

        r  = sqrt(dist2(ip) + eps2)
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

        pot(ip) = pot(ip) &
              + t%charge*rd                                           &  !  monopole term
              + (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3         &  !  dipole
              + half*(fd1*t%quad(1) + fd2*t%quad(2) + fd3*t%quad(3))  &  !  quadrupole
              +       fd4*t%xyquad  + fd5*t%yzquad  + fd6*t%zxquad

        ex(ip) = ex(ip) &
                  + t%charge*dx*rd3                                  &  ! monopole term
                  + fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3)       &  ! dipole term
                  + three * (                                        &  ! quadrupole term
                    half * (                                         &
                        ( five*dx3   *rd7 - three*dx*rd5 )*t%quad(1) &
                      + ( five*dx*dy2*rd7 -       dx*rd5 )*t%quad(2) &
                      + ( five*dx*dz2*rd7 -       dx*rd5 )*t%quad(3) &
                    )                                                &
                    + ( five*dy*dx2  *rd7 - dy*rd5 )*t%xyquad        &
                    + ( five*dz*dx2  *rd7 - dz*rd5 )*t%zxquad        &
                    + ( five*dx*dy*dz*rd7          )*t%yzquad        &
                    )

        ey(ip) = ey(ip) &
                  + t%charge*dy*rd3                                  &
                  + fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3)       &
                  + three * (                                        &
                    half * (                                         &
                        ( five*dy3*rd7    - three*dy*rd5 )*t%quad(2) &
                      + ( five*dy*dx2*rd7 -       dy*rd5 )*t%quad(1) &
                      + ( five*dy*dz2*rd7 -       dy*rd5 )*t%quad(3) &
                    )                                                &
                    + ( five*dx*dy2  *rd7 - dx*rd5 )*t%xyquad        &
                    + ( five*dz*dy2  *rd7 - dz*rd5 )*t%yzquad        &
                    + ( five*dx*dy*dz*rd7          )*t%zxquad        &
                    )

        ez(ip) = ez(ip) &
                  + t%charge*dz*rd3                                  &
                  + fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1)       &
                  + three * (                                        &
                    half * (                                         &
                      + ( five*dz3   *rd7 - three*dz*rd5 )*t%quad(3) &
                      + ( five*dz*dy2*rd7 -       dz*rd5 )*t%quad(2) &
                      + ( five*dz*dx2*rd7 -       dz*rd5 )*t%quad(1) &
                    )                                                &
                    + ( five*dx*dz2  *rd7 - dx*rd5 )*t%zxquad        &
                    + ( five*dy*dz2  *rd7 - dy*rd5 )*t%yzquad        &
                    + ( five*dx*dy*dz*rd7          )*t%xyquad        &
                    )
      end do
    end subroutine calc_force_coulomb_3D


    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D(np, delta, ex, ey, pot, t, eps2)
      implicit none

      integer(kind_particle), intent(in) :: np
      real(kfp), intent(in) :: delta(np,3)
      real(kfp), intent(inout) :: ex(np), ey(np), pot(np)
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

      real(kfp) :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3
      integer(kind_particle) :: ip

      do ip = 1, np
        dx = delta(ip,1)
        dy = delta(ip,2)

        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (dx >= spatial_interaction_cutoff(1) .or. dy >= spatial_interaction_cutoff(2)) cycle
        #endif

        dx2 = dx *dx
        dy2 = dy *dy
        dx3 = dx2*dx
        dy3 = dy2*dy

        d2 = dx2 + dy2 + eps2

        rd2 = one/d2
        rd4 = rd2*rd2
        rd6 = rd4*rd2

        pot(ip) = pot(ip) &
            - half*t%charge*log(d2)               & !  monopole term
            + (dx*t%dip(1) + dy*t%dip(2))*rd2     & !  dipole
            + half*(t%quad(1)*(dx2*rd4 - rd2)     & ! quadrupole
            +       t%quad(2)*(dy2*rd4 - rd2))    &
            + t%xyquad*dx*dy*rd4

        ex(ip) = ex(ip) &
            + t%charge*dx*rd2                              & ! monopole
            + t%dip(1)*(two*dx2  *rd4 - rd2)               & ! dipole
            + t%dip(2)* two*dx*dy*rd4                      &
            + t%quad(1)*(four *dx3    *rd6 - three*dx*rd4) & ! quadrupole
            + t%quad(2)*(four *dx *dy2*rd6 -       dx*rd4) &
            + t%xyquad *(eight*dx2*dy *rd6 -   two*dy*rd4)

        ey(ip) = ey(ip) &
            + t%charge*dy*rd2                              & ! monopole
            + t%dip(2)*(two*dy2  *rd4 - rd2)               & ! dipole
            + t%dip(1)* two*dx*dy*rd4                      &
            + t%quad(2)*(four *dy3    *rd6 - three*dy*rd4) & ! quadrupole
            + t%quad(1)*(four *dy *dx2*rd6 -       dy*rd4) &
            + t%xyquad *(eight*dy2*dx *rd6 -   two*dx*rd4)
      end do
    end subroutine calc_force_coulomb_2D


    !>
    !> CALC_FORCE_LJ
    !>
    !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
    !> shifted by the lattice vector vbox
    !> results are returned exyz, phi
    !>
    subroutine calc_force_LJ(np, delta, dist2, ex, ey, ez, t, aii2)
      implicit none

      integer(kind_particle), intent(in) :: np
      real(kfp), intent(in) :: delta(np,3)
      real(kfp), intent(in) :: dist2(np)
      real(kfp), intent(inout) :: ex(np), ey(np), ez(np)
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: aii2

      real(kfp) :: flj, epsc2, aii2_r2,aii4_r4, r, fljrd
      integer(kind_particle) :: ip

      ! epsc should be > a_ii to get evenly spaced ions
      epsc2 = 0.8_kfp * aii2

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if ( &
          delta(ip,1) >= spatial_interaction_cutoff(1) .or. &
          delta(ip,2) >= spatial_interaction_cutoff(2) .or. &
          delta(ip,3) >= spatial_interaction_cutoff(3) &
        ) cycle
        #endif

        ! Force is repulsive up to and just beyond aii
        aii2_r2 = aii2 / max(dist2(ip), epsc2)
        aii4_r4 = aii2_r2*aii2_r2
        flj = two*(aii4_r4*aii4_r4) - aii4_r4

        !  forces
        r     = sqrt(dist2(ip))
        fljrd = flj/r

        ex(ip) = ex(ip) + delta(ip,1) * fljrd
        ey(ip) = ey(ip) + delta(ip,2) * fljrd
        ez(ip) = ez(ip) + delta(ip,3) * fljrd
      end do
    end subroutine calc_force_LJ


    !>
    !> Calculates 3D Coulomb interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_direct(np, delta, dist2, ex, ey, ez, pot, t, eps2)
      implicit none

      integer(kind_particle), intent(in) :: np
      real(kfp), intent(in) :: delta(np,3)
      real(kfp), intent(in) :: dist2(np)
      real(kfp), intent(inout) :: ex(np), ey(np), ez(np), pot(np)
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

      real(kfp) :: rd,r,rd3charge
      integer(kind_particle) :: ip
      
      real(kfp) :: mask(np)
      
      mask = 1._kfp
      
      where (dist2 <= 0.0_kfp)
        mask = 0.0_kfp
      end where

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if ( &
          delta(ip,1) >= spatial_interaction_cutoff(1) .or. &
          delta(ip,2) >= spatial_interaction_cutoff(2) .or. &
          delta(ip,3) >= spatial_interaction_cutoff(3) &
        ) cycle
        #endif

        r         = sqrt(dist2(ip) + eps2)
        rd        = mask(ip) / r
        rd3charge = t%charge*rd*rd*rd

        pot(ip) = pot(ip) + t%charge*rd
        ex(ip)  = ex(ip) + rd3charge*delta(ip,1)
        ey(ip)  = ey(ip) + rd3charge*delta(ip,2)
        ez(ip)  = ez(ip) + rd3charge*delta(ip,3)
      end do
    end subroutine calc_force_coulomb_3D_direct


    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D_direct(np, delta, ex, ey, pot, t, eps2)
      implicit none

      integer(kind_particle), intent(in) :: np
      real(kfp), intent(in) :: delta(np,3)
      real(kfp), intent(inout) :: ex(np), ey(np), pot(np)
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

      real(kfp) :: rd2charge
      integer(kind_particle) :: ip
      
      real(kfp) :: mask(np), d2(np)
      
      d2 = delta(:,1)*delta(:,1) + delta(:,2)*delta(:,2)

      mask = 1._kfp
      
      where (d2 <= 0.0_kfp)
        mask = 0.0_kfp
      end where
      
      d2 = d2 + eps2

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (delta(ip,1) >= spatial_interaction_cutoff(1) .or. delta(ip,2) >= spatial_interaction_cutoff(2)) cycle
        #endif

        rd2charge = t%charge/d2(ip) * mask(ip)

        pot(ip) = pot(ip) - half*t%charge*log(d2(ip)) * mask(ip)
        ex(ip) = ex(ip) + rd2charge * delta(ip,1)
        ey(ip) = ey(ip) + rd2charge * delta(ip,2)
      end do
    end subroutine calc_force_coulomb_2D_direct


    !>
    !> Calculates 3D Kelbg interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_kelbg_3D_direct(delta, dist2, particle_pack, t, kelbg_invsqrttemp)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: kelbg_invsqrttemp

      real(kfp) :: rd,r,rd3
      real(kfp), parameter :: sqrtpi = sqrt(acos(-1.0_8))
      real(kfp) :: ome, rol, lambda, q, fprefac
      integer(kind_particle) :: ip, np

      np = size(dist2, kind = kind_particle)

      q = t%charge

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if ( &
          delta(ip,1) >= spatial_interaction_cutoff(1) .or. &
          delta(ip,2) >= spatial_interaction_cutoff(2) .or. &
          delta(ip,3) >= spatial_interaction_cutoff(3) &
        ) cycle
        #endif

        if (dist2(ip) > 0.0_kfp) then
          ! TODO: lambda must be adjusted depending on mass and temperature of interacting partners - currently it is fixed for electron-proton interactions
          if (particle_pack%q(ip) * q < 0.) then
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

          r   = sqrt(dist2(ip))
          rd  = one / r
          rd3 = rd*rd*rd
          rol = r  / lambda        !< "r over lambda"
          ome = 1  - exp(-rol*rol) !< "one minus exp(stuff)"

          ! potential
          particle_pack%pot(ip) = particle_pack%pot(ip) + q * rd  * (ome + sqrtpi*rol*(1-erf(rol)))
          !  forces
          fprefac = q * rd3 * ome
          particle_pack%ex(ip) = particle_pack%ex(ip) + fprefac * delta(ip,1)
          particle_pack%ey(ip) = particle_pack%ey(ip) + fprefac * delta(ip,2)
          particle_pack%ez(ip) = particle_pack%ez(ip) + fprefac * delta(ip,3)
        end if
      end do
    end subroutine calc_force_kelbg_3D_direct
end module
