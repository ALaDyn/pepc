
! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute relativistic fpond for standing wave field
!
! ==================================================================

subroutine fpond(t,tpulse,sigma_in,vosc,omega,rho_upper,x,y,z,epon_x,epon_y,epon_z,phipon)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration/rise time
  real, intent(in) :: vosc ! quiver strength
  real, intent(in) :: sigma_in ! pulse width (1/e)
  real, intent(in) :: omega ! laser frequency
  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)
  real, intent(in) :: rho_upper
  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real :: xf, yf, zf, Tpon, Rpon, Xpon, Ypon, Zpon, intensity, gamma, sigma, atten, theta
  real :: r, pi=3.141592654, phi, chi, a02, sigma0
  real :: rho0_up, wp_r, ls_r, gamma_s, f2d

  ! fast oscillations 
 Tpon = (sin(omega*t))**2

  ! intensity envelope
  a02 = vosc**2
  if (t <= 2*tpulse) then 
    intensity = a02*max(0.,sin(pi*t/2./tpulse)**2) 
  else
    Tpon = 0.
  endif

  rho0_up = rho_upper*omega**2   ! Effective density of shelf above xc normalised to rho0

  if (sigma_in <0) then
     sigma0=-sigma_in
     f2d = 0.  ! switch off 2d field components
  else
     sigma0=sigma_in
     f2d=1.
  endif

!  Standing wave in vacuum; evanescent in overdense plasma
!  Use standard solution of Helmholtz equation for step profile

  gamma_s = sqrt(1.+4*intensity/rho_upper)  ! gamma factor for EM solution at x=xc
  wp_r = sqrt(rho0_up/gamma_s)  ! effective plasma frequency of upper shelf
  ls_r = 1./wp_r   ! rel. skin depth

  phi = atan(-omega/wp_r) ! Interface phase factor given by tan(phi) = -k * l_s
  chi = omega*x + phi  ! Vacuum phase 

  if (x.le.0) then
     xf = omega*sin(2*chi)
     Xpon = sin(chi)    ! laser 'Ez'
     sigma = sigma0*sqrt(1.+abs(x)**2/10./sigma0**2) ! take Rayleigh length 4*sigma0

  else
     xf = -2/ls_r*sin(phi)**2*exp(-2*x/ls_r)
     Xpon = sin(phi)*exp(-x/ls_r)    ! laser Ez inside
     sigma = sigma0  ! Don't expand spot inside target
  endif

  r = sqrt(y**2+z**2)

  theta = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot 
!  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian

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
     yf = 0.
     zf = 0.
  endif



  atten = sigma0**2/sigma**2
  phipon = 4*Xpon**2*Rpon*atten*intensity  ! intensity, including attenuation factor

  gamma = sqrt(1.+abs(phipon*Tpon))  ! relativistic phi_pond

  Epon_x = 2*intensity*Tpon*Rpon/gamma*xf*atten
  Epon_y = f2d*2*intensity*Tpon*Xpon**2/gamma*yf*atten
  Epon_z = f2d*2*intensity*Tpon*Xpon**2/gamma*zf*atten
!  Epon_z = 0.   

end subroutine fpond
