! ==============================================
!
!                ADD_RAMP
!
!   $Revision: 1.0
!
!  Add ramp region to front of profile
! 
! ==============================================

subroutine add_ramp(d_layer)

    use physvars
    use treevars
    use utils

    implicit none

    real, intent(in) :: d_layer
    integer              :: i, j, ierr
    integer              :: idum, iseed1, iseed2, iseed3, i1, n1,p, k, nramp, nx, ny
    real                 :: qtot, qramp, s, xi, xedge, rind

    
! Transform particles in layer at front of target
! to form exponential density ramp

    qtot = rho0*d_layer ! total line density
    s = lolam*(1-rho_min)  ! required layer thickness for stretching
    qramp = rho0*s
    nramp = npart*s/d_layer
    xedge = plasma_center(1)-d_layer/2.  ! leading edge (slab & disc)

    do i=1,npp
       if (x(i).le.xedge+s) then
          rind = (x(i)-xedge)/s  ! fractional index ~ i/nramp
          xi = -lolam*log(1. - rind*(1.-rho_min))
          x(i) = xedge+s - xi   ! transform coordinate out into vacuum
       endif
    end do



end subroutine add_ramp


