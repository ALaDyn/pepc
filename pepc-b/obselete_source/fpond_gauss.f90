
! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute relativistic fpond for standing wave field
!
! ==================================================================

subroutine fpond(t,tpulse,sigma0,vosc,omega,xd,yd,zd,epon_x,epon_y,epon_z,phipon)


  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration
  real, intent(in) :: vosc ! quiver strength
  real, intent(in) :: sigma0 ! pulse width (1/e)
  real, intent(in) :: omega ! laser frequency
  real, intent(in) :: xd,yd,zd ! position to evaluate force; xd is away from target

  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real :: xf, yf, zf, Tpon, Rpon, Xpon, gamma, sigma, atten

  !  linear rise
  Tpon = 2*vosc**2*min(1.,t/tpulse) * (sin(omega*t))**2

  sigma = sigma0*sqrt(1.+abs(xd)**2/4./sigma0**2) ! take Rayleigh length 2*sigma0

  Rpon = exp((-yd**2-zd**2)/sigma**2)

  !  Standing wave in vacuum; evanescent in overdense plasma
  if (xd.ge.0) then
     xf = sin(2*omega*xd)
     Xpon = cos(omega*xd)**2
  else
     xf = -2*exp(2*xd)
     Xpon = exp(2*xd)
  endif

  yf = -yd/sigma**2     ! suppress radial fpond for scale-model
  zf = -zd/sigma**2
  atten = sigma0**2/sigma**2
  phipon = Tpon*Xpon*Rpon*atten  ! include attenuation factor

  gamma = sqrt(1.+abs(phipon)/2.)

  Epon_x = Tpon*Rpon/4/gamma*xf*atten
  Epon_y = phipon/2/gamma*yf
  Epon_z = phipon/2/gamma*zf
! Epon_y = 0.
! Epon_z = 0.

end subroutine fpond
