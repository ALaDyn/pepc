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

    real(KIND=8), private :: xp,xm,yp,ym
    real(KIND=8), private :: xpym,xmyp,xmym,xpyp
    real(KIND=8), private :: xpymx,xmypx,xmymx,xpypx,xpymy,xmypy,xmymy,xpypy
    real(KIND=8), private :: lnxpymx,lnxmypx,lnxmymx,lnxpypx,lnxpymy,lnxmypy,lnxmymy,lnxpypy
    real(KIND=8), private :: atanxpym,atanxmyp,atanxmym,atanxpyp
    real(KIND=8), private :: I1,I2,I3
    real(KIND=8), private :: eps = 1e-20

    private :: calculate_common_terms

    contains

!======================================================================================

    SUBROUTINE calculate_common_terms(x,y,z,x0,y0)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0
        real(KIND=8) :: t1,t2,t3,t4

        xp = x+x0
        xm = x-x0
        yp = y+y0
        ym = y-y0
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

        lnxpymx = 0
        lnxmypx = 0
        lnxmymx = 0
        lnxpypx = 0
        lnxpymy = 0
        lnxmypy = 0
        lnxmymy = 0
        lnxpypy = 0

        IF (xpymx**2 > eps) lnxpymx = log(xpymx)
        IF (xmypx**2 > eps) lnxmypx = log(xmypx)
        IF (xmymx**2 > eps) lnxmymx = log(xmymx)
        IF (xpypx**2 > eps) lnxpypx = log(xpypx)
        IF (xpymy**2 > eps) lnxpymy = log(xpymy)
        IF (xmypy**2 > eps) lnxmypy = log(xmypy)
        IF (xmymy**2 > eps) lnxmymy = log(xmymy)
        IF (xpypy**2 > eps) lnxpypy = log(xpypy)

        IF (z**2 > eps) THEN
            atanxmym = atan( ((ym)*(xm)) / ( z * xmym ) )
            atanxpym = atan( ((ym)*(xp)) / ( z * xpym ) )
            atanxmyp = atan( ((yp)*(xm)) / ( z * xmyp ) )
            atanxpyp = atan( ((yp)*(xp)) / ( z * xpyp ) )
        ELSE
            atanxmym = 0
            atanxpym = 0
            atanxmyp = 0
            atanxpyp = 0
        END IF

        t1 = atanxmym
        t2 = atanxpym
        t3 = atanxmyp
        t4 = atanxpyp
        I1 = t1 - t2 - t3 + t4

        t1 = lnxmymy
        t2 = lnxmypy
        t3 = lnxpymy
        t4 = lnxpypy
        I2 = t1 - t2 - t3 + t4

        t1 = lnxmymx
        t3 = lnxpymx
        t2 = lnxmypx
        t4 = lnxpypx
        I3 = t1 - t2 - t3 + t4

    END SUBROUTINE calculate_common_terms


 !======================================================================================

    SUBROUTINE E_bnd(x,y,z,x0,y0,Q,A,E,opt_calc_common_terms)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0,Q,A
        real(KIND=8), intent(out):: E(3)
        logical, intent(in),optional :: opt_calc_common_terms
        logical :: calc_common_terms

        IF (Q == 0) THEN
            E = 0
            RETURN
        END IF

        IF (present(opt_calc_common_terms)) THEN
            calc_common_terms = opt_calc_common_terms
        ELSE
            calc_common_terms = .true.
        END IF
        IF (calc_common_terms) call calculate_common_terms(x,y,z,x0,y0)
        E(1) = I2*Q/A
        E(2) = I3*Q/A
        E(3) = I1*Q/A

    END SUBROUTINE E_bnd

 !======================================================================================

    real*8 FUNCTION Phi_bnd(x,y,z,x0,y0,Q,A,opt_calc_common_terms)
        implicit none
        real(KIND=8), intent(in) :: x,y,z,x0,y0,Q,A
        logical, intent(in), optional :: opt_calc_common_terms
        logical :: calc_common_terms
        real(KIND=8) :: Ao,Ao1,Ao2,Au,Au1,Au2,Bo,Bo1,Bo2,Bu,Bu1,Bu2

        IF (Q == 0) THEN
            Phi_bnd = 0
            RETURN
        END IF

        IF (present(opt_calc_common_terms)) THEN
            calc_common_terms = opt_calc_common_terms
        ELSE
            calc_common_terms = .true.
        END IF
        IF (calc_common_terms) call calculate_common_terms(x,y,z,x0,y0)

        Ao1 = -ym*lnxmymx
        Au1 = -ym*lnxpymx
        Bo1 = -yp*lnxmypx
        Bu1 = -yp*lnxpypx
        Ao2 = -xm*lnxmymy
        Bo2 = -xm*lnxmypy
        Au2 = -xp*lnxpymy
        Bu2 = -xp*lnxpypy

        Ao = Ao1 + Ao2 - z*atanxmym
        Au = Au1 + Au2 - z*atanxpym
        Bo = Bo1 + Bo2 - z*atanxmyp
        Bu = Bu1 + Bu2 - z*atanxpyp

        Phi_bnd = Q/A*(Ao-Au-Bo+Bu)

    END FUNCTION Phi_bnd

!======================================================================================
END MODULE module_boundary_field
