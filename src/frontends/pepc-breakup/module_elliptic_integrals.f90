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
!> elliptic integrals module
!>

module module_elliptic_integrals
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use module_helper
   implicit none

contains
   function elliptic_fk(k) result(rf_out)

      !*****************************************************************************80
      !
  !! ELLIPTIC_FK evaluates the complete elliptic integral of first kind, F(K).
      !
      !  Discussion:
      !
      !    The value is computed using Carlson elliptic integrals:
      !
      !      F(k) = RF ( 0, 1-k^2, 1 ).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    29 May 2018
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) K, the argument.
      !
      !    Output, real ( kind = 8 ) ELLIPTIC_FK, the function value.
      !
      implicit none

      real(kind_physics) :: errtol
      integer :: ierr
      real(kind_physics) :: k
      real(kind_physics) :: rf_out
      real(kind_physics) :: x
      real(kind_physics) :: y
      real(kind_physics) :: z

      x = 0.0D+00
      y = (1.0D+00 - k) * (1.0D+00 + k)
      z = 1.0D+00
      errtol = 1.0D-03

      rf_out = rf(x, y, z, errtol, ierr)
   end

   function elliptic_ek(k, rf) result(value)

      !*****************************************************************************80
      !
  !! ELLIPTIC_EK evaluates the complete elliptic integral of second kind,  E(K).
      !
      !  Discussion:
      !
      !    The value is computed using Carlson elliptic integrals:
      !
      !      E(k) = RF ( 0, 1-k^2, 1 ) - 1/3 k^2 RD ( 0, 1-k^2, 1 ).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    30 May 2018
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) K, the argument.
      !
      !    Output, real ( kind = 8 ) ELLIPTIC_EK, the function value.
      !
      implicit none

      real(kind_physics) :: errtol
      integer :: ierr
      real(kind_physics) :: k
      real(kind_physics), intent(in) :: rf
      real(kind_physics) :: rd_out
      ! real(kind_physics) :: rf_out
      real(kind_physics) :: value
      real(kind_physics) :: x
      real(kind_physics) :: y
      real(kind_physics) :: z

      x = 0.0D+00
      y = (1.0D+00 - k) * (1.0D+00 + k)
      z = 1.0D+00
      errtol = 1.0D-03

      rd_out = rd(x, y, z, errtol, ierr)
      ! rf_out = rf(x, y, z, errtol, ierr)
      value = rf - k * k * rd_out / 3.0D+00
   end

   function rf(x, y, z, errtol, ierr)

      !*****************************************************************************80
      !
  !! RF computes an incomplete elliptic integral of the first kind, RF(X,Y,Z).
      !
      !  Discussion:
      !
      !    This function computes the incomplete elliptic integral of the first kind.
      !
      !    RF(X,Y,Z) = Integral ( 0 <= T < oo )
      !
      !                    -1/2     -1/2     -1/2
      !          (1/2)(T+X)    (T+Y)    (T+Z)    DT,
      !
      !    where X, Y, and Z are nonnegative and at most one of them is zero.
      !
      !    If X or Y or Z is zero, the integral is complete.
      !
      !    The duplication theorem is iterated until the variables are
      !    nearly equal, and the function is then expanded in Taylor
      !    series to fifth order.
      !
      !    Check by addition theorem:
      !
      !      RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W) = RF(0,Z,W),
      !      where X, Y, Z, W are positive and X * Y = Z * W.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    30 May 2018
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
      !    This FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Bille Carlson,
      !    Computing Elliptic Integrals by Duplication,
      !    Numerische Mathematik,
      !    Volume 33, 1979, pages 1-16.
      !
      !    Bille Carlson, Elaine Notis,
      !    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
      !    ACM Transactions on Mathematical Software,
      !    Volume 7, Number 3, pages 398-403, September 1981.
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) X, Y, Z, the arguments in the integral.
      !
      !    Input, real ( kind = 8 ) ERRTOL, the error tolerance.
      !    Relative error due to truncation is less than
      !      ERRTOL ^ 6 / (4 * (1 - ERRTOL)).
      !    Sample choices:
      !      ERRTOL   Relative truncation error less than
      !      1.D-3    3.D-19
      !      3.D-3    2.D-16
      !      1.D-2    3.D-13
      !      3.D-2    2.D-10
      !      1.D-1    3.D-7
      !
      !    Output, integer ( kind = 4 ) IERR, the error flag.
      !    0, no error occurred.
      !    1, abnormal termination.
      !
      implicit none

      real(kind_physics), parameter :: c1 = 0.0416666666666666666
      real(kind_physics), parameter :: c2 = 0.0681818181818181818
      real(kind_physics), parameter :: c3 = 0.0714285714285714285
      real(kind_physics) :: e2
      real(kind_physics) :: e3
      real(kind_physics) :: epslon
      real(kind_physics) :: errtol
      integer :: ierr, rf_iter, loop_count, l
      real(kind_physics) :: lamda
      real(kind_physics) :: lolim
      real(kind_physics) :: mu
      real(kind_physics) :: rf
      real(kind_physics) :: s
      real(kind_physics) :: uplim
      real(kind_physics) :: x
      real(kind_physics) :: xn
      real(kind_physics) :: xndev
      real(kind_physics) :: xnroot
      real(kind_physics) :: y
      real(kind_physics) :: yn
      real(kind_physics) :: yndev
      real(kind_physics) :: ynroot
      real(kind_physics) :: z
      real(kind_physics) :: zn
      real(kind_physics) :: zndev
      real(kind_physics) :: znroot
      !
      !  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
      !  LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
      !  UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
      !
      save lolim
      save uplim

      data lolim/3.D-78/
      data uplim/1.D+75/

      ! if ( &
      !   x < 0.0D+00 .or. &
      !   y < 0.0D+00 .or. &
      !   z < 0.0D+00 .or. &
      !   x + y < lolim .or. &
      !   x + z < lolim .or. &
      !   y + z < lolim .or. &
      !   uplim <= x .or. &
      !   uplim <= y .or. &
      !   uplim <= z ) then
      !   write ( *, '(a)' ) ''
      !   write ( *, '(a)' ) 'RF - Error!'
      !   write ( *, '(a)' ) '  Invalid input arguments.'
      !   write ( *, '(a,d23.16)' ) '  X = ', x
      !   write ( *, '(a,d23.16)' ) '  Y = ', y
      !   write ( *, '(a,d23.16)' ) '  Z = ', z
      !   write ( *, '(a)' ) ''
      !   ierr = 1
      !   rf = 0.0D+00
      !   return
      ! end if

      if (y .lt. 0.0 .or. y .lt. lolim .or. uplim .le. y) then
         write (*, '(a)') 'RF - Error!'
         write (*, '(a)') '  Invalid input arguments.'
         write (*, '(a,d23.16)') '  Y = ', y
         ierr = 1
         rf = 0.0D+00
         return
      end if

      ierr = 0
      xn = x
      yn = y
      zn = z
      rf_iter = 0
      loop_count = 8

      ! do

      !   mu = ( xn + yn + zn ) / 3.0d0
      !   xndev = 2.0d0 - ( mu + xn ) / mu
      !   yndev = 2.0d0 - ( mu + yn ) / mu
      !   zndev = 2.0d0 - ( mu + zn ) / mu
      !   epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )

      !   if ( epslon < errtol ) then
      !     ! c1 = 1.0d0 / 24.0d0
      !     ! c2 = 3.0d0 / 44.0d0
      !     ! c3 = 1.0d0 / 14.0d0
      !     e2 = xndev * yndev - zndev * zndev
      !     e3 = xndev * yndev * zndev
      !     s = 1.0d0 + ( c1 * e2 - 0.1d0 - c2 * e3 ) * e2 + c3 * e3
      !     rf = s / sqrt ( mu )
      !     print *, "rf iteration: ", rf_iter
      !     return
      !   end if

      !   xnroot = sqrt ( xn )
      !   ynroot = sqrt ( yn )
      !   znroot = sqrt ( zn )
      !   lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
      !   xn = ( xn + lamda ) * 0.25d0
      !   yn = ( yn + lamda ) * 0.25d0
      !   zn = ( zn + lamda ) * 0.25d0
      !   rf_iter = rf_iter + 1
      ! end do

      do l = 1, loop_count
         xnroot = sqrt(xn)
         ynroot = sqrt(yn)
         znroot = sqrt(zn)
         lamda = xnroot * (ynroot + znroot) + ynroot * znroot
         xn = (xn + lamda) * 0.25d0
         yn = (yn + lamda) * 0.25d0
         zn = (zn + lamda) * 0.25d0
      end do

      mu = (xn + yn + zn) / 3.0d0
      xndev = 2.0d0 - (mu + xn) / mu
      yndev = 2.0d0 - (mu + yn) / mu
      zndev = 2.0d0 - (mu + zn) / mu
      epslon = max(abs(xndev), abs(yndev), abs(zndev))

      if (epslon .gt. errtol) print *, "Unsatisfactory rf epslon."
      e2 = xndev * yndev - zndev * zndev
      e3 = xndev * yndev * zndev
      s = 1.0d0 + (c1 * e2 - 0.1d0 - c2 * e3) * e2 + c3 * e3
      rf = s / sqrt(mu)
   end

   function rd(x, y, z, errtol, ierr)

      !*****************************************************************************80
      !
  !! RD computes an incomplete elliptic integral of the second kind, RD(X,Y,Z).
      !
      !  Discussion:
      !
      !    This function computes an incomplete elliptic integral of the second kind.
      !
      !    RD(X,Y,Z) = Integral ( 0 <= T < oo )
      !
      !                    -1/2     -1/2     -3/2
      !          (3/2)(T+X)    (T+Y)    (T+Z)    DT,
      !
      !    where X and Y are nonnegative, X + Y is positive, and Z is positive.
      !
      !    If X or Y is zero, the integral is complete.
      !
      !    The duplication theorem is iterated until the variables are
      !    nearly equal, and the function is then expanded in Taylor
      !    series to fifth order.
      !
      !    Check:
      !
      !      RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y) = 3 / sqrt ( X * Y * Z ),
      !      where X, Y, and Z are positive.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    30 May 2018
      !
      !  Author:
      !
      !    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
      !    This FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Bille Carlson,
      !    Computing Elliptic Integrals by Duplication,
      !    Numerische Mathematik,
      !    Volume 33, 1979, pages 1-16.
      !
      !    Bille Carlson, Elaine Notis,
      !    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
      !    ACM Transactions on Mathematical Software,
      !    Volume 7, Number 3, pages 398-403, September 1981.
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) X, Y, Z, the arguments in the integral.
      !
      !    Input, real ( kind = 8 ) ERRTOL, the error tolerance.
      !    The relative error due to truncation is less than
      !      3 * ERRTOL ^ 6 / (1-ERRTOL) ^ 3/2.
      !    Sample choices:
      !      ERRTOL   Relative truncation error less than
      !      1.D-3    4.D-18
      !      3.D-3    3.D-15
      !      1.D-2    4.D-12
      !      3.D-2    3.D-9
      !      1.D-1    4.D-6
      !
      !    Output, integer ( kind = 4 ) IERR, the error flag.
      !    0, no error occurred.
      !    1, abnormal termination.
      !
      implicit none

      real(kind_physics), parameter :: c1 = 0.21428571428571428571
      real(kind_physics), parameter :: c2 = 0.16666666666666666666
      real(kind_physics), parameter :: c3 = 0.40909090909090909090
      real(kind_physics), parameter :: c4 = 0.11538461538461538461
      real(kind_physics) :: ea
      real(kind_physics) :: eb
      real(kind_physics) :: ec
      real(kind_physics) :: ed
      real(kind_physics) :: ef
      real(kind_physics) :: epslon
      real(kind_physics) :: errtol
      integer :: ierr, rd_iter, loop_count, l
      real(kind_physics) :: lamda
      real(kind_physics) :: lolim
      real(kind_physics) :: mu
      real(kind_physics) :: power4
      real(kind_physics) :: rd
      real(kind_physics) :: sigma
      real(kind_physics) :: s1
      real(kind_physics) :: s2
      real(kind_physics) :: uplim
      real(kind_physics) :: x
      real(kind_physics) :: xn
      real(kind_physics) :: xndev
      real(kind_physics) :: xnroot
      real(kind_physics) :: y
      real(kind_physics) :: yn
      real(kind_physics) :: yndev
      real(kind_physics) :: ynroot
      real(kind_physics) :: z
      real(kind_physics) :: zn
      real(kind_physics) :: zndev
      real(kind_physics) :: znroot
      !
      !  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
      !  LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ^ (2/3).
      !  UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
      !  MINIMUM) ^ (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
      !  IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
      !  NEVER BE CHOSEN SMALLER THAN 1.D-5.
      !
      save lolim
      save uplim

      data lolim/6.D-51/
      data uplim/1.D+48/

      ! if ( &
      !   x < 0.0D+00 .or. &
      !   y < 0.0D+00 .or. &
      !   x + y < lolim .or. &
      !   z < lolim .or. &
      !   uplim < x .or. &
      !   uplim < y .or. &
      !   uplim < z ) then
      !   write ( *, '(a)' ) ''
      !   write ( *, '(a)' ) 'RD - Error!'
      !   write ( *, '(a)' ) '  Invalid input arguments.'
      !   write ( *, '(a,d23.16)' ) '  X = ', x
      !   write ( *, '(a,d23.16)' ) '  Y = ', y
      !   write ( *, '(a,d23.16)' ) '  Z = ', z
      !   write ( *, '(a)' ) ''
      !   ierr = 1
      !   rd = 0.0D+00
      !   return
      ! end if

      if (y .lt. 0.0 .or. x + y .lt. lolim .or. uplim .lt. y) then
         write (*, '(a)') 'RD - Error!'
         write (*, '(a)') '  Invalid input arguments.'
         write (*, '(a,d23.16)') '  Y = ', y
         ierr = 1
         rd = 0.0D+00
         return
      end if

      ierr = 0
      xn = x
      yn = y
      zn = z
      sigma = 0.0d0
      power4 = 1.0d0
      rd_iter = 0
      loop_count = 8

      ! do
      !
      !   mu = ( xn + yn + 3.0d0 * zn ) * 0.2d0
      !   xndev = ( mu - xn ) / mu
      !   yndev = ( mu - yn ) / mu
      !   zndev = ( mu - zn ) / mu
      !   epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )

      !   if ( epslon < errtol ) then
      !     ! c1 = 3.0d0 / 14.0d0
      !     ! c2 = 1.0d0 / 6.0d0
      !     ! c3 = 9.0d0 / 22.0d0
      !     ! c4 = 3.0d0 / 26.0d0
      !     ea = xndev * yndev
      !     eb = zndev * zndev
      !     ec = ea - eb
      !     ed = ea - 6.0d0 * eb
      !     ef = ed + ec + ec
      !     s1 = ed * ( - c1 + 0.25d0 * c3 * ed - 1.5d0 * c4 * zndev * ef )
      !     s2 = zndev * ( c2 * ef + zndev * ( - c3 * ec + zndev * c4 * ea ) )
      !     rd = 3.0d0 * sigma + power4 * ( 1.0d0 + s1 + s2 ) / ( mu * sqrt ( mu ) )
      !     print *, "rd iteration: ", rd_iter
      !     return
      !   end if

      !   xnroot = sqrt ( xn )
      !   ynroot = sqrt ( yn )
      !   znroot = sqrt ( zn )
      !   lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
      !   sigma = sigma + power4 / ( znroot * ( zn + lamda ) )
      !   power4 = power4 * 0.25d0
      !   xn = ( xn + lamda ) * 0.25d0
      !   yn = ( yn + lamda ) * 0.25d0
      !   zn = ( zn + lamda ) * 0.25d0
      !   rd_iter = rd_iter + 1
      ! end do

      do l = 1, loop_count
         xnroot = sqrt(xn)
         ynroot = sqrt(yn)
         znroot = sqrt(zn)
         lamda = xnroot * (ynroot + znroot) + ynroot * znroot
         sigma = sigma + power4 / (znroot * (zn + lamda))
         power4 = power4 * 0.25d0
         xn = (xn + lamda) * 0.25d0
         yn = (yn + lamda) * 0.25d0
         zn = (zn + lamda) * 0.25d0
      end do

      mu = (xn + yn + 3.0d0 * zn) * 0.2d0
      xndev = (mu - xn) / mu
      yndev = (mu - yn) / mu
      zndev = (mu - zn) / mu
      epslon = max(abs(xndev), abs(yndev), abs(zndev))

      if (epslon .gt. errtol) print *, "Unsatisfactory rd epslon."
      ea = xndev * yndev
      eb = zndev * zndev
      ec = ea - eb
      ed = ea - 6.0d0 * eb
      ef = ed + ec + ec
      s1 = ed * (-c1 + 0.25d0 * c3 * ed - 1.5d0 * c4 * zndev * ef)
      s2 = zndev * (c2 * ef + zndev * (-c3 * ec + zndev * c4 * ea))
      rd = 3.0d0 * sigma + power4 * (1.0d0 + s1 + s2) / (mu * sqrt(mu))

   end
end module
