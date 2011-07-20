
! ==================================================================
!
!                        TRAVELLING PONDEROMOTIVE FORCE
!
!  Compute relativistic fpond for propagating laser field
!
! ==================================================================

subroutine laser_bullet(tpulse,sigma,vosc,x,y,z,epon_x,epon_y,epon_z,phipon)

  implicit none

  real, intent(in) :: tpulse ! pulse duration
  real, intent(in) :: vosc ! quiver strength
  real, intent(in) :: sigma ! pulse width (1/e)
  real, intent(in) :: x,y,z ! position to evaluate force; distance from laser centre (x0,0,0)

  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real :: xf, yf, zf, Rpon, Xpon, gamma
  real :: pi=3.141592654, phi, a02

  !  sin^2 time envelope

  phi = pi*x/tpulse/2.
  Rpon = exp((-y**2-z**2)/sigma**2)  ! Gaussian radial profile

  if (x.ge.-tpulse .and. x.lt.tpulse) then
     xf = -pi/4.*sin(2*phi)
     Xpon = cos(phi)**2
  else
     xf = 0.
     Xpon = 0.
  endif

  yf = -2*y/sigma**2
  zf = -2*z/sigma**2

  a02 = vosc**2
  phipon = a02*Xpon*Rpon  ! intensity

  gamma = sqrt(1.+abs(phipon)/2.)  ! relativistic phi_pond

  Epon_x = a02*Rpon/gamma*xf
  Epon_y = a02*Xpon*Rpon/gamma*yf
  Epon_z = a02*Xpon*Rpon/gamma*zf

end subroutine laser_bullet




