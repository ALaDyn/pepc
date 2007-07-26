! ==============================================
!
!               STRETCH_SPHERE 
!
!   $Revision: 1.0
!
!  Convert uniform sphere to tapered density profile
! 
! ==============================================

subroutine stretch_sphere(r0)

    use physvars
    use treevars
    use utils

    implicit none

    real, intent(in) :: r0 
    integer              :: i, j, ierr
    integer              :: idum, iseed1, iseed2, iseed3, i1, n1,p, k, nramp, nx, ny
    real                 :: qtot, qramp, xi, xedge, rcut 
    real*8 :: rt, xt, yt, zt, s
    
! Transform particles in layer at front of target
! to form exponential density ramp


    do i=1,npp
	xt = x(i)-plasma_centre(1)
	yt = y(i)-plasma_centre(2)
	zt = z(i)-plasma_centre(3)
	xi = pelabel(i)/1.001/npart ! impose cutoff
	rt = sqrt(xt**2+yt**2+zt**2)
! inverted  coord for 1/(1+s^3)^2 distribution
	s = (1./(1.-xi) - 1.)**(1./3.)
	x(i) = plasma_centre(1) + xt/rt*r0*s  ! keep direction vector; scale by sphere radius 
	y(i) = plasma_centre(2) + yt/rt*r0*s
	z(i) = plasma_centre(3) + zt/rt*r0*s
    end do



end subroutine stretch_sphere 


