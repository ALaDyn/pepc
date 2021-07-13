! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2019 Juelich Supercomputing Centre,
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
!>  Encapsulates functions for setting up particle velocities with different models
!>
module module_velocity_setup
   use module_pepc_kinds
   use physvars
   use module_pepc_types
   implicit none
   private

   public maxwell1
   public maxwell2
   public maxwell3
   public scramble_v
   public cold_start
   public rano

contains

   !>
   !>   COLD_START
   !>   initialises velocity array slice to 0.
   !>   @param u array of velocities to be initialized
   !>   @param nmax total number of entries in u
   !>   @param i1 minimal index in u to be used
   !>   @param n maximum index in u to be used
   !>
   subroutine cold_start(ux, uy, uz)
      implicit none
      real*8 :: ux(:), uy(:), uz(:)

      ux(:) = 0.
      uy(:) = 0.
      uz(:) = 0.
   end subroutine cold_start

   !>
   !>   MAXWELL1
   !>   initialises 1D Maxwellian velocity distribution
   !>   @param u array of velocities to be initialized
   !>   @param nmax total number of entries in u
   !>   @param i1 minimal index in u to be used
   !>   @param n maximum index in u to be used
   !>   @param vt desired average velocity of particles
   !>
   subroutine maxwell1(u, nmax, i1, n, vt)
      implicit none
      integer(kind_particle), intent(in) ::  nmax, i1, n
      real*8, intent(in) :: vt
      real*8 :: u(nmax)
      real, parameter :: pi = 3.141592654

      integer(kind_particle) :: ip1, ip2, i, cntr, nv
      real :: f0, df, v, dv, vmax, finf, deltv
      real*8 :: vip

      if (n .eq. 0) return
      nv = 30 * n
      vmax = 4.0
      dv = vmax / nv
      f0 = 0.5
      finf = sqrt(pi / 2.0)
      cntr = 1
      ip1 = n / 2 + i1 - 1
      ip2 = n / 2 + i1

      do i = 1, nv
         v = vmax - (i - 0.5) * dv
         deltv = vmax / 10.
         df = exp(max(-30.0, -0.5 * v**2)) * dv / finf * n / 2.
         f0 = f0 + df       ! integrate dist. fn.
         if (f0 .ge. cntr) then
            vip = vt * (v - dv * (f0 - cntr) / df)
            u(ip1) = vip
            u(ip2) = -vip
            cntr = cntr + 1
            ip1 = ip1 - 1
            ip2 = ip2 + 1
         endif
      end do

      !  odd one out
      if (mod(n, 2_kind_particle) .eq. 1) then
         u(n) = 0.
      endif
   end subroutine maxwell1

   !>
   !>   MAXWELL2
   !>   initialises 2D Maxwellian velocity distribution by direct inversion
   !>   assumes thermal distribution exp{ -(vx**2+vy**2)/2v_t**2 }
   !>   @param ux,uy array of velocities to be initialized
   !>   @param n maximum index in ux to be used
   !>   @param vt desired thermal velocity of particles
   !>
   subroutine maxwell2(u1, u2, n, vt)
      implicit none
      integer(kind_particle), intent(in) :: n
      real, intent(in) :: vt
      real*8 :: u1(n), u2(n)
      real, parameter :: pi = 3.141592654

      integer(kind_particle) :: i
      real*8 :: theta, u0
      integer :: dum1 = -319

      if (n .eq. 0) return

      do i = 1, n
         u0 = vt * sqrt(-2.*log((i - 0.5) / n))
         theta = 2 * pi * rano(dum1)
         u1(i) = u0 * cos(theta)
         u2(i) = u0 * sin(theta)
      end do
   end subroutine maxwell2

   !>
   !>   MAXWELL3
   !>   initialises 3D Maxwellian velocity distribution
   !>   @param ux,uy,uz arrays of velocities to be initialized
   !>   @param nmax total number of entries in u
   !>   @param i1 minimal index in u to be used
   !>   @param n maximum index in u to be used
   !>   @param vt desired average velocity of particles
   !>
   subroutine maxwell3(ux, uy, uz, nmax, i1, n, vt)
      implicit none
      integer(kind_particle), intent(in) ::  nmax, i1, n
      real*8, intent(in) :: vt
      real*8 :: ux(nmax), uy(nmax), uz(nmax), vc

      vc = vt / sqrt(3.) ! homogeneous: vx=vy=vz = sqrt(|v|^2/3)

      call maxwell1(ux, nmax, i1, n, vc)
      call maxwell1(uy, nmax, i1, n, vc)
      call maxwell1(uz, nmax, i1, n, vc)
      call scramble_v(ux, uy, uz, nmax, i1, n)   ! remove x,y,z correlations
   end subroutine maxwell3

   !>
   !>   SCRAMBLE_V
   !>
   subroutine scramble_v(ux, uy, uz, nmax, i1, n)
      implicit none

      integer :: dum1, dum2, dum3
      real*8 :: uxt, uyt, uzt
      integer(kind_particle) :: i, j, k, kk, p, i1, n, n1, nmax
      real*8 :: ux(nmax), uy(nmax), uz(nmax)

      dum1 = -71 - 10 * my_rank
      dum2 = -113301 - 10 * my_rank
      dum3 = -8651 - 10 * my_rank
      !  exclude odd one out
      if (mod(n, 2_kind_particle) .ne. 0) then
         n1 = n - 1
      else
         n1 = n
      endif

      !  scramble indices to remove correlation between ux,uy,uz
      do i = 1, n1
         p = i + i1 - 1
         j = int(n1 * rano(dum1) + i1)
         k = int(n1 * rano(dum2) + i1)
         kk = int(n1 * rano(dum3) + i1)
         uxt = ux(p)
         uyt = uy(p)
         uzt = uz(p)
         ux(p) = ux(kk)
         uy(p) = uy(j)
         uz(p) = uz(k)
         ux(kk) = uxt
         uy(j) = uyt
         uz(k) = uzt
      end do
   end subroutine scramble_v

   ! Random number scrambler
   ! =======================
   !
   !  called with:
   !               x=rano(iseed)
   !
   !  returns real number in the interval (0.0 - 1.0)
   !  set iseed = -ve integer before using call
   !
   !  Example of calling routine:
   !
   !      subroutine coords
   !      include 'common.h'
   !
   !      iseed1 = -11
   !      iseed2 = -7
   !
   !
   !      do i=1,n
   ! x(i)=xlen*rano(iseed1)
   !       y(i)=ylen*rano(iseed2)
   !      end do
   !
   !
   !      end

   ! - routine taken from Numerical Recipies, p195
   !
   real * 8 function rano(idum)
      implicit none
      integer :: idum
      real*8, save :: dseed, dum
      real*8, save :: v(97), y
      integer, save :: iff, icall, i, j
      data iff, icall/0, 0/
      if (idum .lt. 0 .or. iff .eq. 0) then
         !  first call generates new sequence
         iff = 1
         dseed = abs(idum) * 1.0
         idum = 1
         do j = 1, 97
            dum = genran(dseed)
         end do
         do j = 1, 97
            v(j) = genran(dseed)
         end do
         y = genran(dseed)
      endif

      !  next index - make sure we don`t overstep array bounds if
      !  generator returns a 0.0 or 1.0

      j = max(mod(1 + int(97.*y), 98), 1)
      if (j .gt. 97 .or. j .lt. 1) then
         write (6, *) 'Call: ', icall
         write (6, *) 'idum = ', idum, 'j = ', j, ' y= ', y
         write (6, *) 'Random No. generator not initialised properly'
         write (6, *) 'dummy =', dum, ' dseed=', dseed
         write (6, 100) (i, v(i), i=1, 97)
100      format(i4, f10.6)
         stop
      endif
      !  get next variate and generate new one to fill gap

      y = v(j)
      rano = y
      v(j) = genran(dseed)
      icall = icall + 1

      return
   end function rano

   real * 8 function genran(dseed)
      real*8 ::  dseed
      real*8 ::  d2p31m, d2p31
      ! d2p31m=(2**31) - 1
      ! d2p31 =(2**31)(or an adjusted value)
      data d2p31m/2147483647.0/
      data d2p31/2147483648.0/

      dseed = mod(16807.0 * dseed, d2p31m)
      genran = dseed / d2p31
   end function genran
end module module_velocity_setup
