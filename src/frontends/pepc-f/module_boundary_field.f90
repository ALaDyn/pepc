! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
!                         Forschungszentrum Juelich GmbH,
!                         Germany
!
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

!!!!!!!!!!!!!!!!!!!!
!! boundary field module module
!! This module uses an analytical solution for the field
!! and the potential of a homogeneously charged rectangular
!! surface in the xy-plane (-x0 < x < x0; -y0 < y < y0) and
!! adds as an external field
!!
!! A force constant (e.g. 1/(4*pi*eps0) for SI system) has to be multiplied to the results
!!!!!!!!!!!!!!!!!!!!

MODULE module_boundary_field
    use module_pepc_types
    use variables

    implicit none

    contains

!======================================================================================

    real*8 FUNCTION I1(x,y,z,x0,y0)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0
        real(KIND=8) :: t1,t2,t3,t4

        IF (z == 0) THEN
            I1 = 0._8
        ELSE
            t1 = atan( ((y-y0)*(x-x0)) / ( z * sqrt( (x-x0)**2 + (y-y0)**2 + z**2 ) ) )
            t2 = atan( ((y-y0)*(x+x0)) / ( z * sqrt( (x+x0)**2 + (y-y0)**2 + z**2 ) ) )
            t3 = atan( ((y+y0)*(x-x0)) / ( z * sqrt( (x-x0)**2 + (y+y0)**2 + z**2 ) ) )
            t4 = atan( ((y+y0)*(x+x0)) / ( z * sqrt( (x+x0)**2 + (y+y0)**2 + z**2 ) ) )
            I1 = t1 - t2 - t3 + t4
        END IF

    END FUNCTION I1

!======================================================================================

    real*8 FUNCTION I2(x,y,z,x0,y0)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0
        real(KIND=8) :: t1,t2,t3,t4
        real(KIND=8) :: xp,xm,yp,ym
        real(KIND=8) :: eps = 1e-20

        xp = x+x0
        xm = x-x0
        yp = y+y0
        ym = y-y0
        t1 = 0._8
        t2 = 0._8
        t3 = 0._8
        t4 = 0._8

        IF ((xm**2 + z**2) > eps) THEN
            t1 = log( sqrt( xm**2 + ym**2 + z**2 ) - ym )
            t2 = log( sqrt( xm**2 + yp**2 + z**2 ) - yp )
        END IF
        IF ((xp**2 + z**2) > eps) THEN
            t3 = log( sqrt( xp**2 + ym**2 + z**2 ) - ym )
            t4 = log( sqrt( xp**2 + yp**2 + z**2 ) - yp )
        END IF
        I2 = t1 - t2 - t3 + t4

    END FUNCTION I2

 !======================================================================================

    real*8 FUNCTION I3(x,y,z,x0,y0)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0
        real(KIND=8) :: t1,t2,t3,t4
        real(KIND=8) :: xp,xm,yp,ym
        real(KIND=8) :: eps = 1e-20

        xp = x+x0
        xm = x-x0
        yp = y+y0
        ym = y-y0
        t1 = 0._8
        t2 = 0._8
        t3 = 0._8
        t4 = 0._8

        IF ((ym**2 + z**2) > eps) THEN
            t1 = log( sqrt( xm**2 + ym**2 + z**2 ) - xm )
            t3 = log( sqrt( xp**2 + ym**2 + z**2 ) - xp )
        END IF
        IF ((yp**2 + z**2) > eps) THEN
            t2 = log( sqrt( xm**2 + yp**2 + z**2 ) - xm )
            t4 = log( sqrt( xp**2 + yp**2 + z**2 ) - xp )
        END IF
        I3 = t1 - t2 - t3 + t4

    END FUNCTION I3

 !======================================================================================

    SUBROUTINE E_bnd(x,y,z,x0,y0,Q,A,E)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0,Q,A
        real(KIND=8), intent(out):: E(3)

        IF (Q == 0) THEN
            E = 0
            RETURN
        END IF

        E(1) = I2(x,y,z,x0,y0)
        E(2) = I3(x,y,z,x0,y0)
        E(3) = I1(x,y,z,x0,y0)
        E = E*Q/A

    END SUBROUTINE E_bnd

 !======================================================================================

    real*8 FUNCTION Phi_bnd(x,y,z,x0,y0,Q,A)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0,Q,A
        real(KIND=8) :: xp,xm,yp,ym
        real(KIND=8) :: xpym,xmyp,xmym,xpyp
        real(KIND=8) :: xpymx,xmypx,xmymx,xpypx,xpymy,xmypy,xmymy,xpypy
        real(KIND=8) :: Ao,Ao1,Ao2,Au,Au1,Au2,Bo,Bo1,Bo2,Bu,Bu1,Bu2
        real(KIND=8) :: eps = 1e-20

        IF (Q == 0) THEN
            Phi_bnd = 0
            RETURN
        END IF

        xp = x+x0
        xm = x-x0
        yp = y+y0
        ym = y-y0

        xpym = sqrt(xp**2 + ym**2 + z**2)
        xmyp = sqrt(xm**2 + yp**2 + z**2)
        xmym = sqrt(xm**2 + ym**2 + z**2)
        xpyp = sqrt(xp**2 + yp**2 + z**2)

        xpymx = xpym - xp
        xmypx = xmyp - xm
        xmymx = xmym - xm
        xpypx = xpyp - xp
        xpymy = xpym - ym
        xmypy = xmyp - yp
        xmymy = xmym - ym
        xpypy = xpyp - yp

        Ao1 = 0
        Ao2 = 0
        Au1 = 0
        Au2 = 0
        Bo1 = 0
        Bo2 = 0
        Bu1 = 0
        Bu2 = 0

        IF (ym**2 > eps) THEN
            Ao1 = -ym*log(xmymx)
            Au1 = -ym*log(xpymx)
        END IF
        IF (yp**2 > eps) THEN
            Bo1 = -yp*log(xmypx)
            Bu1 = -yp*log(xpypx)
        END IF
        IF (xm**2 > eps) THEN
            Ao2 = -xm*log(xmymy)
            Bo2 = -xm*log(xmypy)
        END IF
        IF (xp**2 > eps) THEN
            Au2 = -xp*log(xpymy)
            Bu2 = -xp*log(xpypy)
        END IF

        Ao = Ao1 + Ao2
        Au = Au1 + Au2
        Bo = Bo1 + Bo2
        Bu = Bu1 + Bu2

        IF (z**2 > eps) THEN
            Ao = Ao - z*atan(xm*ym/xmym/z)
            Au = Au - z*atan(xp*ym/xpym/z)
            Bo = Bo - z*atan(xm*yp/xmyp/z)
            Bu = Bu - z*atan(xp*yp/xpyp/z)
        END IF

        Phi_bnd = Q/A*(Ao-Au-Bo+Bu)

    END FUNCTION Phi_bnd

!======================================================================================
END MODULE module_boundary_field
