!     ==================================
!     
!     Electromagnetic fields 
!     
!     - plane wave with finite rise-time and Gaussian spot
!  returns laser fields at particle position x,y,z
!     
!     
!     ==================================


subroutine emplane(t,tpulse,sigma0,a0,w0,x,y,z,ez,by,bx,az,phipon)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration or rise-time
  real, intent(in) :: a0 ! quiver strength
  real, intent(in) :: sigma0 ! pulse width (1/e)
  real, intent(in) :: w0 ! laser frequency
  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

  real, intent(out) :: phipon, ez, by, bx, az ! pond. potential and fields

  real :: tenv, theta_r, gamma0, sigma, atten
  real :: r2, pi=3.141592654, phi, chi, phase, earg
  real :: k0

  !     linear rise

  gamma0 = sqrt(1 + a0**2/2)  ! Sarachik & Schappert gamma
  
  k0=w0
  tenv = min(1.,t/tpulse)
!  tenv = 1.
  phase = w0*t - k0*x
  r2 = y**2+z**2


  if (sigma0.gt.0) then
     !  pulse envelope - Gaussian, 1/e width = sigma
     earg = min(20.,r2/2/sigma0**2)
     theta_r = exp(-earg)
  else 
     !  constant amplitude
     theta_r = 1.
  endif

  ! reconstruct EM fields (s-pol)

  Az = a0*tenv*theta_r*sin(phase)
  Ez = -w0*a0*tenv*theta_r*cos(phase)      ! Ez = -dAz/dt
  By = k0*a0*tenv*theta_r*cos(phase)       ! By = -dAz/dx
  Bx = -y/sigma0**2*a0*tenv*theta_r*sin(phase) ! Bx = dAz/dy


  phipon = Az**2

end subroutine emplane









