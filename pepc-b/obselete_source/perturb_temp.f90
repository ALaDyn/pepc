 
!  ===============================================================
!
!                       PERTURB_TEMP 
!
!   $Revision: 332 $
!
!   Add temperature perturbation to isotropic thermal distrib. 
!
!  ===============================================================


subroutine perturb_temp 

  use physvars
  use treevars
  use utils

  integer :: i
  real :: Te0, deltaT, k_therm, s
  real :: lambdaT, xpert, xpe


! Scale velocities by space-dep. temperature perturbation 

Te0=Te_keV
xpert=xl*0.8
lambdaT = xpert/kpert
k_therm = 2*pi/lambdaT   ! Perturbation wavenumber - leave 10% buffer at either end 
deltaT0 = tpert*Te0  ! 50% temperature variation 
if (me==0) then 
  write(*,*) 'PEPC-B | Doing electron temperature perturbation'
  write(*,*) 'PEPC-B | Wavelength: ',lambdaT
endif

do i=1,npp
  if (q(i)<0 .and. x(i) >= 0.1*xl .and. x(i) < 0.9*xl) then
      xpe = x(i)-xl/10. 
     deltaT = deltaT0 * sin(k_therm*xpe)
     s = sqrt((deltaT+Te0)/Te0)
     ux(i) = s*ux(i)
     uy(i) = s*uy(i)
     uz(i) = s*uz(i)
  else
! buffer region at Te0
  endif
enddo
 

end subroutine perturb_temp
