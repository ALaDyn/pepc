
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

  real :: xf, yf, zf, Tpon, Rpon, Xpon, Ypon, Zpon, gamma, sigma, atten

  !  linear rise
  Tpon = 2*vosc**2*min(1.,t/tpulse) * (sin(omega*t))**2

  sigma = sigma0*sqrt(1.+abs(xd)**2/4./sigma0**2) ! take Rayleigh length 2*sigma0

!  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian
  

  !  Standing wave in vacuum; evanescent in overdense plasma
  if (xd.ge.0) then
     xf = sin(2*omega*xd)
     Xpon = cos(omega*xd)**2
  else
     xf = -2*exp(2*xd)
     Xpon = exp(2*xd)
  endif

  if (abs(yd)<2*sigma) then
     yf = sin(pi*yd/2./sigma)
     Ypon = cos(pi*yd/4./sigma)**2
  else
     yf = 0.
     Ypon = 0.
  endif

  if (abs(zd)<2*sigma) then
     zf = sin(pi*zd/2./sigma)
     Zpon = cos(pi*yd/4./sigma)**2
  else
     zf = 0.
     Zpon = 0.
  endif

  atten = sigma0**2/sigma**2
  phipon = Tpon*Xpon*Ypon*Zpon*atten  ! include attenuation factor

  gamma = sqrt(1.+abs(phipon)/2.)

  Epon_x = Tpon*Ypon*Zpon/4/gamma*xf*atten
  Epon_y = Tpon*Xpon/16./gamma*yf
  Epon_z = Tpon*Xpon/16./gamma*zf


end subroutine fpond
