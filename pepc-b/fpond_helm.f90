
! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute ponderomotive fields from Helmholtz solution for EM standing wave 
!
! ==================================================================

subroutine fpond_helm(t,tpulse,sigma_in,vosc,omega, &
                      x,y,z,ux,Az,nxh,xh_start,xh_end,dxh,x_crit, &
                      epon_x,epon_y,epon_z,phipon)

  !   call fpond_helm( tlaser, tpulse,sigma,vosc,omega, &
  !                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
  !	  	      epon_x,epon_y,epon_z,phipon)


  implicit none
  real, intent(in) :: t 	!< time
  real, intent(in) :: tpulse !< pulse duration/rise time
  real, intent(in) :: vosc !< quiver strength
  real, intent(in) :: sigma_in !< pulse width (1/e)
  real, intent(in) :: omega !< laser frequency
  real, intent(in) :: x,y,z !< position to evaluate force;
                            !< x absolute; y,z relative to laser axis
  real, intent(in) :: ux !< forward momentum for gamma factor
  integer, intent(in) :: nxh !< # 1D Helmholtz grid points
  real, intent(in) :: xh_start !< Start point of HH grid
  real, intent(in) :: xh_end !< End point of HH grid
  real, intent(in) :: dxh !< HH grid spacing
  real, intent(in) :: x_crit  !< target surface
  complex, intent(in) :: Az(0:nxh+1) !< Vector pot from Helmholtz solution
  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real ::  yf, zf, xh, Rpon, intensity, gamma, sigma, atten, theta
  real :: Azr_0, Azr_1, Azr_2, Azr_3  ! Vector pot. at control points
  real :: epon1, epon2 ! pond force at control points
  real :: ayi, azi ! vec. pot at particle
  real :: epxi, epyi, epzi ! pond field at particle
  real :: r, pi=3.141592654, sigma0
  real :: pha, xa, b1, b2
  real :: f2d ! 2D switch
  real :: Z_R  ! Rayleigh length
  complex :: yi = (0.,1.)
  integer ::  i1, i2

  ! fast oscillations 
  pha = omega*t

  if (sigma_in <0) then
     sigma0=-sigma_in
     f2d = 0.  ! switch off 2d field components
  else
     sigma0=sigma_in
     f2d=1.
  endif

  ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
  Z_R = omega*sigma0**2


  if (x.le.x_crit) then
     sigma = sigma0*sqrt(1.+abs(x-x_crit)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam
  else
     sigma = sigma0  ! Don't expand spot inside target
  endif

  r = sqrt(y**2+z**2)
  atten = sigma0**2/sigma**2  ! attenuation factor from focal cone

  ! Radial gradients
  theta = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot 

  if (x.lt.xh_start .or. x.gt.xh_end .or. r >= 2*sigma) then
     yf = 0.
     zf = 0.
     Rpon = 0.

  else if (r < 2*sigma .and. r/= 0.) then
     yf = -pi/4./sigma*y/r*sin(2*theta)
     zf = -pi/4./sigma*z/r*sin(2*theta)
     Rpon = cos(theta)**2

  else
     ! special case on axis
     Rpon = cos(theta)**2
     yf = 0.
     zf = 0.
  endif



  xh = x-xh_start  ! particle coord on HH grid
  xa = xh/dxh     ! reduced coord 
  if (xa.ge.1 .and. xa.le.nxh) then
     ! Only compute force for particles inside HH grid
     ! leave one-point buffer at either end for difference dA/dx
     i1 = xa  ! lower NGP
     i2 = i1+1       ! upper NGP
     b2=xa-i1  ! linear weights  W_j = 1-|x_i-x_j|
     b1=1.-b2

     ! Derive fpond from Az - need gradient at both reference points.
     ! include temporal phase factor  
     Azr_0 = Real(Az(i1-1)*cexp(yi*pha))
     Azr_1 = Real(Az(i1)*cexp(yi*pha))
     Azr_2 = Real(Az(i2)*cexp(yi*pha))
     Azr_3 = Real(Az(i2+1)*cexp(yi*pha))

     ! pond force at reference points i1, i2
     ! - 2nd order difference without gamma factor, as in emfield
     Epon1 = .5*omega/dxh*( Azr_2**2-Azr_0**2 )
     Epon2 = .5*omega/dxh*( Azr_3**2-Azr_1**2 )

     !        ayi = b1*ay(i1) + b2*ay(i2)
     ayi = 0.
     azi = b1*azr_1 + b2*azr_2
     gamma = sqrt(1. + ux**2 + ayi**2 + azi**2)

     ! pond field at particle including total particle gamma
     epxi = (b1*epon1 + b2*epon2)/gamma
     epyi = azi**2/gamma
     epzi = azi**2/gamma

     ! potential and fields, correcting for radial dep. and attenuation factor

     phipon = gamma*Rpon*atten 

     Epon_x = epxi*Rpon*atten
     Epon_y = epyi*yf*atten
     Epon_z = epzi*zf*atten
  else
     phipon=0.
     Epon_x=0.
     Epon_y=0.
     Epon_z=0.
  endif


end subroutine fpond_helm
