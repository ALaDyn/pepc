!     ==================================
!     
!     Electromagnetic fields 
!     
!     - ponderomotive laser model for step-profile
!  returns laser fields at particle position x,y,z relative to critical surface
!     
!     
!     ==================================


subroutine empond(t,tpulse,sigma0,a0,w0,x,y,z,ez,by,bx,az,phipon)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration or rise-time
  real, intent(in) :: a0 ! quiver strength
  real, intent(in) :: sigma0 ! pulse width (1/e)
  real, intent(in) :: w0 ! laser frequency
  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

  real, intent(out) :: phipon, ez, by, bx, az ! pond. potential and fields

  real :: tenv, Rpon, dRpon, gamma, sigma, atten, theta
  real :: r, pi=3.141592654, phi, chi, phase
  real :: k0, lskin, gamma_c, f_helm, g_helm, wp_r, nonc

  !     linear rise

  phase = w0*t
  k0=w0
  nonc = 1./w0**2 ! density normalised to nc
  Tenv = min(1.,t/tpulse) ! time envelope


  !     Standing wave in vacuum; evanescent in overdense plasma
  !    Use standard solution of Helmholtz equation for step profile

  gamma_c = sqrt(1.+4*a0**2/nonc) ! gamma at surface
  wp_r = 1./sqrt(gamma_c)
  lskin = 1./wp_r   ! Rel. skin depth in EM units

  !   Phase factor given by tan(phi) = -k0 * l_s = k0 * c/wp
  phi = atan(-1./wp_r)   
  r = sqrt(y**2+z**2)


  if (x <= 0) then
     !     vacuum solution
     chi = k0*x + phi      ! Vacuum phase 
     f_helm = sin(chi)     ! laser 'Ez'
     g_helm = cos(chi)     ! laser 'By'
     sigma = sigma0*sqrt(1.+abs(x)**2/4./sigma0**2) ! take Rayleigh length 4*sigma0

  else
     !   evanescent wave - need Sudan solution here
     f_helm = sin(phi)*exp(-x/lskin)        ! laser Ez inside
     g_helm = cos(phi)*exp(-x/lskin) ! laser By inside
     sigma = sigma0  ! Don't expand spot inside target
  endif


  !  Radial envelope

  theta = pi*r/4./sigma
  atten = sigma0/sigma

  if (r <= 2.*sigma .and. r/=0. ) then      
     Rpon = atten*cos(theta)
     dRpon = -pi*y/4./sigma/r*atten*sin(theta)

  else
     Rpon = 0.
     dRpon = 0.
  endif

  ! reconstruct EM fields (s-pol)

  az = 2*a0     *           f_helm * sin(phase) * Rpon * tenv

  ez = -2*w0*a0 *           f_helm * cos(phase) * Rpon * tenv ! Ez = -dAz/dt
  by = -2*k0*a0 *           g_helm * sin(phase) * Rpon * tenv ! By = -dAz/dx
  bx =  2*a0 *              f_helm * sin(phase) * dRpon* tenv ! Bx = dAz/dy

  phipon = (2*a0*f_helm*tenv*Rpon)**2

end subroutine empond









