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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates the low-level kernels for Coulomb- and similar interactions
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_coulomb_kernels
  use module_pepc_types
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
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates 3D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_force_coulomb_3D(t, d, dist2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kfp), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kfp), intent(out) ::  exyz(3), phi

      real(kfp) :: dx,dy,dz,r,rd,rd2,rd3

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r  = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd = one/r
      rd2 = rd *rd
      rd3 = rd *rd2

      phi = -t%charge*rd

      exyz(1) = -t%charge*dx*rd3
      exyz(2) = -t%charge*dy*rd3
      exyz(3) = -t%charge*dz*rd3

    end subroutine calc_force_coulomb_3D

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is: 
    !>   Phi = -2q log R 
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_force_coulomb_2D(t, d, d2, exy, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kfp), intent(in) :: d(2), d2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kfp), intent(out) ::  exy(1:2),phi

     ! real(kfp) :: dx,dy,rd2,rd4,rd6,dx2,dy2,dx3,dy3

     ! Should abort here, because we always use 3D

    end subroutine calc_force_coulomb_2D
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> CALC_FORCE_LJ
    !>
    !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
    !> shifted by the lattice vector vbox
    !> results are returned exyz, phi
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_force_LJ(t, d, r2, aii2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kfp), intent(in) :: d(3), r2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kfp), intent(out) ::  exyz(3), phi
      real(kfp), intent(in) :: aii2
      real(kfp) :: flj, epsc2, aii2_r2,aii4_r4, r, fljrd

      ! epsc should be > a_ii to get evenly spaced ions
      epsc2 = 0.8_kfp * aii2

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates 3D Coulomb interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_force_coulomb_3D_direct(t, d, dist2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kfp), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kfp), intent(out) ::  exyz(3), phi

      real(kfp) :: rd,r,rd3charge

      r         = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd        = one/r
      rd3charge = t%charge*rd*rd*rd

      !phi  = t%charge*rd
      !exyz = rd3charge*d

      phi  = -t%charge*rd
      exyz = -rd3charge*d

    end subroutine calc_force_coulomb_3D_direct

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_force_coulomb_2D_direct(t, d, d2, exy, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kfp), intent(in) :: d(2), d2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kfp), intent(out) ::  exy(2), phi

      real(kfp) :: rd2charge

      phi       = - half*t%charge*log(d2)
      rd2charge = t%charge/d2
      exy       = rd2charge*d

    end subroutine calc_force_coulomb_2D_direct

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates 3D Kelbg interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_force_kelbg_3D_direct(particle, t, d, dist2, kelbg_invsqrttemp, exyz, phi)
      implicit none

      type(t_particle), intent(inout) :: particle
      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kfp), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kfp), intent(out) ::  exyz(3), phi
      real(kfp), intent(in) :: kelbg_invsqrttemp
      real(kfp) :: rd,r,rd3
      real(kfp), parameter :: sqrtpi = sqrt(acos(-1.0_8))
      real(kfp) :: ome, rol, lambda, q, fprefac

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
