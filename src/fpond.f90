
! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute relativistic fpond for standing wave field
!
! ==================================================================

subroutine fpond(t,tpulse,sigma0,vosc,omega,x,y,z,epon_x,epon_y,epon_z,phipon)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration
  real, intent(in) :: vosc ! quiver strength
  real, intent(in) :: sigma0 ! pulse width (1/e)
  real, intent(in) :: omega ! laser frequency
  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real :: xf, yf, zf, Tpon, Rpon, Xpon, Ypon, Zpon, gamma, sigma, atten, theta
  real :: r, pi=3.141592654, phi, chi, a02

  !  linear rise
  Tpon = min(1.,t/tpulse) * (sin(omega*t))**2


!  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian


!  Standing wave in vacuum; evanescent in overdense plasma
!  Use standard solution of Helmholtz equation for step profile

  phi = atan(-omega) ! Phase factor given by tan(phi) = -k * l_s
  chi = omega*x + phi  ! Vacuum phase 

  if (x.le.0) then
     xf = omega*sin(2*chi)
     Xpon = sin(chi)    ! laser 'Ez'
     sigma = sigma0*sqrt(1.+abs(x)**2/10./sigma0**2) ! take Rayleigh length 4*sigma0

  else
     xf = -2*sin(phi)**2*exp(-2*x)
     Xpon = sin(phi)*exp(-x)    ! laser Ez inside
     sigma = sigma0  ! Don't expand spot inside target
  endif

  r = sqrt(y**2+z**2)
  theta = pi*r/4./sigma

  if (r < 2*sigma .and. r/= 0.) then
     yf = -pi/4./sigma*y/r*sin(2*theta)
     zf = -pi/4./sigma*z/r*sin(2*theta)
     Rpon = cos(theta)**2
  else if (r >= 2*sigma) then
     yf = 0.
     zf = 0.
     Rpon = 0.
  else
     Rpon = cos(theta)**2

  endif



  atten = sigma0**2/sigma**2
  a02 = vosc**2
  phipon = 4*a02*Xpon**2*Rpon*atten  ! intensity, including attenuation factor

  gamma = sqrt(1.+abs(phipon*Tpon)/2.)  ! relativistic phi_pond

  Epon_x = a02*Tpon*Rpon/gamma*xf*atten
  Epon_y = a02*Tpon*Xpon**2/gamma*yf*atten
  Epon_z = a02*Tpon*Xpon**2/gamma*zf*atten
!  phipon = sqrt(1.+abs(phipon)/2.) ! 
  Epon_y=0.
  Epon_z=0.
end subroutine fpond





