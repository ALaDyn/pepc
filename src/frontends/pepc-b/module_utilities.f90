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

!
!      General utility module
!   *   I/O routines
!   *   Random number gen.
!

module module_utilities

   interface blank
      module procedure blankn, blank6
   end interface

   interface cputime
      module procedure cput
   end interface

   public rano
   public phase

contains

   subroutine blankn(ichan)
      integer :: ichan
      write (ichan, '(/)')
   end subroutine blankn

   subroutine blank6
      write (6, '(/)')
   end subroutine blank6

   subroutine cput(sec)
      real :: sec
      !    integer :: ic1, ir1, im1
      !    CALL SYSTEM_CLOCK(COUNT=IC1, COUNT_RATE=IR1, COUNT_MAX=IM1)
      !    sec = 1.*ic1/ir1
      call cpu_time(sec)

      !  timef returns elapsed wall-clock time in ms
      !      sec=1.e-3*timef()
      !  IBM function rtc returns elapsed wall-clock time in micro-s - compile with -lxlf90
      !    sec = rtc()

   end subroutine cput

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
   !        x(i)=xlen*rano(iseed1)
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
      end if

      !  next index - make sure we do not overstep array bounds if
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
      end if
      !  get next variate and generate new one to fill gap

      y = v(j)
      rano = y
      v(j) = genran(dseed)
      icall = icall + 1

      return
   end function rano

   real * 8 function genran(dseed)
      !                                  specifications for arguments
      real*8 ::  dseed
      !                                  specifications for local variables
      !    real ::  d2p31m,d2p31
      real*8 ::  d2p31m, d2p31
      !                                  d2p31m=(2**31) - 1
      !                                  d2p31 =(2**31)(or an adjusted value)
      data d2p31m/2147483647.0/
      data d2p31/2147483648.0/
      !                                  first executable statement
      dseed = mod(16807.0 * dseed, d2p31m)
      genran = dseed / d2p31
      return
   end function genran

   !  ========================================
   !
   !      function   PHASE
   !
   !  Returns -pi -> pi (4-quadrant) phase of complex pair (x,y)
   !
   !  ========================================

   real function phase(x, y)
      implicit none
      real, parameter :: pi = 3.1415926536
      real, intent(in) :: x, y
      integer :: sx, sy, itx

      sx = sign(1., x)
      sy = sign(1., y)
      itx = (1 - sx) / 2  ! 0 or 1

      ! special cases first
      if (x .eq. 0 .and. y .eq. 0) then
         phase = 0.
      else if (x .eq. 0.) then
         phase = sy * pi / 2.
      else if (y .eq. 0) then
         phase = itx * pi
      else
         phase = atan(y / x) + itx * sy * pi
      end if

   end function phase

end module module_utilities
