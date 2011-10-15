!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates laser setup and laser-particle interaction
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_laser
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public setup_laser	!< Laser initialisation	
      public laser		!< Laser control (intensity)
      public force_laser	!< Compute external force from laser
      public laser_hist		!< Track laser history
      public laser_monitor	!< Write laser parameters to stdout
      public emplane		!< Plane wave
      public emplane_lin	!< Plane wave with linear rise-time
      public laser_bullet	!< Travelling fpond
      public empond		!< Standing wave solution
      public emobliq		!< Standing wave, oblique incidence
      public fpond_lin		!< Standing wave, linear rise-time
      public fpond		!< Standing wave
      public fpond_helm		!< Helmholtz solver
      public trisolve		!< Tridiagonal LU decomp
      public em_helmholtz	!< Helmholtz solver 2
      public track_nc		!< Track critical density surface
      public density_helmholtz  !< Density gather for Helmholtz solver

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: pi=3.141592654


      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Setup all local module variables and data that depends on them
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setup_laser()
          use module_physvars
          implicit none

          beam_config=mod(beam_config_in,10)  ! derived config from s,p variations
          navcycle   = 2*pi/dt/omega  ! # timesteps in a laser cycle

        end subroutine setup_laser

subroutine laser_monitor(I_laser)
  use module_physvars
  implicit none
  real :: I_laser
  integer :: ifile

!     write(ifile_cpu,'(//a,i8,(3x,a20,f8.2)/(3x,a,f8.2,a2,f8.2,a4)/a,f9.3,1pe12.3)') 'Timestep ',itime+itime_start &
!          ,' total run time = ',trun &
!          ,' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
!          ,' Laser amplitude, intensity  = ',sqrt(I_laser),I_laser*1.37e18

        do ifile = 6,24,18

              write(ifile,'(//(3x,a,f8.2,a2,f8.2,a4)/4(a20,f9.3/))') &
                    ' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
                   ,' amplitude = ',sqrt(I_laser) &
                   ,' x_crit= ',x_crit &
                   ,' spot size= ',sigma & 
                   ,' theta =  ',theta_beam 

	end do
end subroutine laser_monitor


subroutine laser(I_laser)

  use module_physvars

  implicit none


  real :: I_laser
  real :: amplitude, ampl_min, pha
  real :: Azr(1:nxh)
  complex :: yi=(0.,1.)

  if (beam_config >= 3) tlaser = tlaser + dt  
  pha = tlaser*omega


!  Laser pulse envelope
!  =====================

  laser_model: select case(beam_config_in)

     case(4,34,44,64,6)  	! sin^2 standing wave
	if (tlaser<2*tpulse)  then
           I_laser = vosc**2*max(0.,sin(pi*tlaser/2./tpulse)**2)
	else
	   I_laser=0.
	endif
     case(14,94,54,74,16) 	! standing wave, linear rise
           I_laser = vosc**2*min(1.,tlaser/tpulse)
     case(3) 		! constant
        I_laser = vosc**2
     case default
	I_laser = 0.
  end select laser_model




!    Laser focal position and rezoning
!    ==================================

  laser_focus: select case(beam_config_in)

    case(44,54)  ! Helmholtz solver for vector potential

! Factor-in pulse shape
      amplitude = sqrt(I_laser)  ! amplitude of incoming wave
      call density_helmholtz
      call em_helmholtz(my_rank,itime,nxh,dxh,theta_beam,amplitude,omega,rho_helm,Az_helm)
      Azr(1:nxh) = Real(Az_helm(1:nxh)*cexp(yi*pha))
      ampl_max = maxval(Azr)
      ampl_min = minval(Azr)
      if (ampl_max.lt.abs(ampl_min)) ampl_max=-ampl_min
  
    case(64,74)  ! Helmholtz solver for vector potential circ pol

! Factor-in pulse shape
      amplitude = sqrt(I_laser)  ! amplitude of incoming wave
      call density_helmholtz
      call em_helmholtz(my_rank,itime,nxh,dxh,theta_beam,amplitude,omega,rho_helm,Az_helm)
      Azr(1:nxh) = Real(Az_helm(1:nxh))
      ampl_max = maxval(Azr)
      ampl_min = minval(Azr)
 
    case(5)  ! propagating fpond
     !  Trigger rezoning if laser rezone_frac of the way through plasma
     ! - Only works after restart at present
      if (restart .and. beam_config ==5 .and. focus(1) >= window_min + x_plasma*rezone_frac) then
        if (my_rank==0) then
           write (*,*) 'REZONE'
           !           read (*,*) go
        endif

 !       call rezone
        !        window_min = window_min + dt
      else
        focus(1) = focus(1) + dt  ! propagate forward by c*dt - can include v_g here
        propag_laser=propag_laser + dt
	write (*,*) focus(1)
      endif

    case(4,14,24,94)  ! old fpond model
      if (itime>0) focus(1) = x_crit  ! laser tracks n_c

    case default
     ! leave focal point unchanged

  end select laser_focus

end subroutine laser

!  ===================================================================
!
!                              FORCE_LASER
!
!   $Revision: 1065 $
!
!   Calculate forces due to external fields
! eg: laser, stationary B, E fields
!
!
!
!  ===================================================================


subroutine force_laser(p_start,p_finish)

  use module_physvars
  use module_particle_props

  implicit none
  include 'mpif.h'

  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.
  integer :: p
  real :: xd, yd, zd  ! positions relative to centre of laser spot
  real :: uxd ! x-momentum
  real :: xt, yt, zt, rt  ! positions relative to plasma centre
  real :: Epon_x, Epon_y, Epon_z, Phipon, ez_em, bx_em, by_em, bz_em, az_em

  dxh = (xh_end-xh_start)/nxh  ! HH grid spacing

  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  ! Include force from laser/external field on electrons - ES scheme

  if (beam_config >= 3 .and. beam_config <=7 ) then
     do p = p_start, p_finish
        if (q(p)<0) then
           xd = x(p)-focus(1)
           yd = y(p)-focus(2)
           zd = z(p)-focus(3)

           laser_model: select case(beam_config_in)

           case(3)  ! Uniform sinusoid in x (s-pol)
              Epon_z = 0.
              Epon_y = 0.
              Epon_x = vosc*omega*sin(omega*tlaser)
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(23)  ! Uniform sinusoid in z (s-pol)
              Epon_x = 0.
              Epon_y = 0.
              Epon_z = vosc*omega*sin(omega*tlaser)
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.


           case(4)  ! standing wave fpond, sin2 pulse

              call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon)
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(14)  ! standing wave fpond, linear rise-time, fields reduced
              call fpond_lin( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon)
              epon_y=epon_y/50.
              epon_z=epon_z/50.
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(94)  ! standing wave fpond with transverse fields artificially reduced
              call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper,xd,yd,zd,epon_x,epon_y,epon_z,phipon)
              epon_y=epon_y/10.
              epon_z=epon_z/10.
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(24)  ! oblique incidence standing wave, s-pol
              call emobliq( tlaser, tpulse,sigma,vosc,omega,theta_beam,rho_upper, &
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon,ez_em,bx_em,by_em)
	      Bz_em = 0.

           case(34) ! fpond with fully EM components Gaussian spot, sin^2 pulse
              call empond(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,Epon_z,By_em,Bx_em,az_em,phipon)
	      Epon_x=0.
	      Epon_y=0.
	      Bz_em=0.

           case(44,54)  ! fpond derived from Az_helm; both linear & sin2 pulse forms
	      xd = x(p)
	      uxd = ux(p)
              call fpond_helm( tlaser, sigma, omega, &
                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
		   epon_x,epon_y,epon_z,phipon)

	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(64)  ! fpond derived from Az_helm, c-pol light
	      xd = x(p)
	      uxd = ux(p)
              call fpond_helmc( sigma,omega, &
                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
		   epon_x,epon_y,epon_z,phipon)
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(5)  ! propagating fpond
              call laser_bullet( tpulse,sigma,vosc,xd,yd,zd,epon_x,epon_y,epon_z,phipon)
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(6) ! plane wave with Gaussian spot, sin^2 pulse
              call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,Epon_z,By_em,Bx_em,az_em,phipon)
	      Epon_x=0.
	      Epon_y=0.
	      Bz_em=0.

           case(16) ! plane wave with Gaussian spot, linear rise-time
              call emplane_lin(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,Epon_z,By_em,Bx_em,az_em,phipon)
	      Epon_x=0.
	      Epon_y=0.
	      Bz_em=0.

           case default  ! no laser
              Epon_x=0
              Epon_y=0
              Epon_z=0
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           end select laser_model

 	else    ! (q>0)
           !  ions assumed not to feel laser, so zero fields
           Epon_x=0
           Epon_y=0
           Epon_z=0
           Bx_em = 0.
           By_em = 0.
           Bz_em = 0.

        endif

	if (beam_config_in==7)  then ! Constant B in z-direction - either charge
           Epon_x=0.
           Epon_y=0.
           Epon_z=0.
           Bx_em = 0.
           By_em = 0.
           Bz_em = vosc

	else if (beam_config_in==17)  then !  Z-pinch: circular B in x,y
           xt=x(p) - plasma_centre(1)
           yt=y(p) - plasma_centre(2)
           zt=z(p) - plasma_centre(3)
           rt = sqrt(xt**2+yt**2)
           Epon_x=0.
           Epon_y=0.
           Epon_z=0.
           Bx_em = -vosc*yt/rt
           By_em = vosc*xt/rt
           Bz_em = 0.

	else if (beam_config_in==27)  then !  Tokamak: circular B in theta + const Bz
           xt=x(p) - plasma_centre(1)
           yt=y(p) - plasma_centre(2)
           zt=z(p) - plasma_centre(3)
           rt = sqrt(xt**2+yt**2)
           Epon_x=0.
           Epon_y=0.
           Epon_z=0.
           Bx_em = -vosc/3.*yt/rt
           By_em = vosc/3.*xt/rt
           Bz_em = vosc

	else if (beam_config_in==37)  then !  Mirror in Bz
           xt=x(p) - plasma_centre(1)
           yt=y(p) - plasma_centre(2)
           zt=z(p) - plasma_centre(3)
           rt = sqrt(xt**2+yt**2)
           Epon_x=0.
           Epon_y=0.
           Epon_z=0.          
           Bx_em = 0.
           By_em = 0.
    
           if (zt>-zl/2. .and. zt<zl/2.) then
              Bx_em = -vosc/5.*xt/rt*(2*zt/zl)
              By_em = -vosc/5.*yt/rt*(2*zt/zl)

              Bz_em = vosc*(2*zt/zl)**2+vosc/5.
           else
              Bx_em=0.
              By_em=0.
              Bz_em = 0.
           endif
	endif

        fpon_max = max(fpon_max,abs(Epon_x))
        ! Add external fields to new particle field
        Ex(p) = Ex(p) + Epon_x
        Ey(p) = Ey(p) + Epon_y
        Ez(p) = Ez(p) + Epon_z
        Bx(p) = Bx(p) + Bx_em
        By(p) = By(p) + By_em
        Bz(p) = Bz(p) + Bz_em

     end do
  endif

end subroutine force_laser

!  ================================
!
!         Laser_hist
!
!    time history of laser params
!
!  $Revision 1.6$
!
!  ================================


subroutine laser_hist

  use module_physvars

  implicit none

if (my_rank.eq.0) then
  ! Write out to energy.dat file
  if (itime.eq.1)  write(71,'(a)') '! time  a_L fpond xc'
  write (71,'(f12.5,2(1pe12.3))') &
       tlaser, ampl_max, fpon_max 
  write (*,'(f12.5,2(1pe12.3))') &
       tlaser, ampl_max, fpon_max
endif
end subroutine laser_hist


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

  real :: tenv, theta_r, gamma0
  real :: r2, pi=3.141592654, phase, earg
  real :: k0

  !     linear rise

  gamma0 = sqrt(1 + a0**2/2)  ! Sarachik & Schappert gamma
  
  k0=w0
 if (t <= 2*tpulse) then
   tenv = max(0.,sin(pi*t/2./tpulse)**2)
  else
   tenv= 0.
  endif

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


!     ==================================
!     
!     Electromagnetic fields 
!     
!     - plane wave with finite rise-time and Gaussian spot
!  returns laser fields at particle position x,y,z
!     
!     
!     ==================================


subroutine emplane_lin(t,tpulse,sigma0,a0,w0,x,y,z,ez,by,bx,az,phipon)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration or rise-time
  real, intent(in) :: a0 ! quiver strength
  real, intent(in) :: sigma0 ! pulse width (1/e)
  real, intent(in) :: w0 ! laser frequency
  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

  real, intent(out) :: phipon, ez, by, bx, az ! pond. potential and fields

  real :: tenv, theta_r, gamma0
  real :: r2, phase, earg
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

end subroutine emplane_lin

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

  real :: tenv, Rpon, dRpon, sigma, atten, theta
  real :: r, pi=3.141592654, phi, chi, phase
  real :: k0, lskin, gamma_c, f_helm, g_helm, wp_r, nonc
  real :: Z_R  ! Rayleigh length

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

 ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
  Z_R = w0*sigma0**2

  !   Phase factor given by tan(phi) = -k0 * l_s = k0 * c/wp
  phi = atan(-1./wp_r)   
  r = sqrt(y**2+z**2)


  if (x <= 0) then
     !     vacuum solution
     chi = k0*x + phi      ! Vacuum phase 
     f_helm = sin(chi)     ! laser Ez
     g_helm = cos(chi)     ! laser By
     sigma = sigma0*sqrt(1.+abs(x)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam

  else
     !   evanescent wave - need Sudan solution here
     f_helm = sin(phi)*exp(-x/lskin)        ! laser Ez inside
     g_helm = cos(phi)*exp(-x/lskin) ! laser By inside
     sigma = sigma0  ! Dont expand spot inside target
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



! ==================================================================
!
!              EMOBLIQ
!
!  Compute relativistic fields for oblique-incidence standing wave field
!       s-light
!
! ==================================================================

subroutine emobliq(t,tpulse,sigma_in,vosc,omega,theta,rho_upper,x,y,z,epon_x,epon_y,epon_z,phipon,Ez,Bx,By)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration/rise time
  real, intent(in) :: vosc ! quiver strength
  real, intent(in) :: sigma_in ! pulse width (1/e)
  real, intent(in) :: omega ! laser frequency
  real, intent(in) :: theta ! angle of incidence

  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)
  real, intent(in) :: rho_upper
  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields
  real, intent(out) :: Ez, By, Bx  ! laser fields

  real :: xf, yf, zf, Tpon, Rpon, Xpon, intensity, gamma, sigma, atten, alpha, eps, dXpon
  real :: r, pi=3.141592654, phi_s, chi_s, a02, sigma0
  real :: rho0_up, wp_r, ls_r, gamma_s, f2d
  real :: kx, ky, tphase, k0, thetar

  ! phase factors (normalised to nominal wp)
  thetar=pi*theta/180.
  k0 = omega
  kx = k0*cos(thetar)
  ky = k0*sin(thetar)
  tphase = omega*t - ky*y  ! time phase

  Tpon = sin(tphase)**2

  ! intensity envelope
  a02 = vosc**2
  if (t <= 2*tpulse) then 
!    intensity = a02*max(0.,sin(pi*t/2./tpulse)**2) 
     intensity = a02
  else
    intensity = 0.
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
!  wp_r=1.
  ls_r = 1./wp_r   ! rel. skin depth

  phi_s = atan(-kx*ls_r) ! Interface phase factor given by tan(phi) = -kx * l_r
  chi_s = kx*x + phi_s  ! Vacuum phase (s-pol) 

  if (x.le.0) then
     xf = omega*sin(2*chi_s)
     Xpon = sin(chi_s)    ! laser phase
     dXpon = kx*cos(chi_s) ! derivative
     sigma = sigma0*sqrt(1.+abs(x)**2/10./sigma0**2) ! take Rayleigh length 4*sigma0
     eps = 1.  ! refractive index
  else
     xf = -2/ls_r*sin(phi_s)**2*exp(-2*x/ls_r)
     Xpon = sin(phi_s)*exp(-x/ls_r)    ! laser phase inside
     dXpon = -sin(phi_s)/ls_r*exp(-x/ls_r) ! derivative
     sigma = sigma0  ! Dont expand spot inside target
     eps = 1.-(wp_r/omega)**2  ! refractive index inside target
  endif

  r = sqrt(y**2+z**2)

  alpha = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot 
!  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian

  if (r < 2*sigma .and. r/= 0.) then
     yf = -pi/4./sigma*y/r*sin(2*alpha)
     zf = -pi/4./sigma*z/r*sin(2*alpha)
     Rpon = cos(alpha)**2
  else if (r >= 2*sigma) then
     yf = 0.
     zf = 0.
     Rpon = 0.
  else
     Rpon = cos(alpha)**2
     yf = 0.
     zf = 0.
  endif



  atten = sigma0**2/sigma**2
  phipon = 4*Xpon**2*Rpon*atten*intensity*Tpon  ! intensity, including attenuation factor

  gamma = sqrt(1.+abs(phipon*Tpon))  ! relativistic phi_pond

  Ez = 2*vosc*cos(tphase)*Xpon*sqrt(Rpon*atten)
  Bx = 2*ky*vosc/omega*cos(tphase)*Xpon*sqrt(Rpon*atten)
  By = 2*vosc/omega*sin(tphase)*dXpon*sqrt(Rpon*atten)

  Epon_x = 2*intensity*Tpon*Rpon/gamma*xf*atten
  Epon_y = f2d*2*intensity*Tpon*Xpon**2/gamma*yf*atten
  Epon_z = f2d*2*intensity*Tpon*Xpon**2/gamma*zf*atten
!  Epon_z = 0.   

end subroutine emobliq

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

  real :: xf, yf, zf, Tpon, Rpon, Xpon, intensity, gamma, sigma, atten, theta
  real :: r, pi=3.141592654, phi, chi, a02, sigma0
  real :: rho0_up, wp_r, ls_r, gamma_s, f2d
  real :: Z_R  ! Rayleigh length

  ! fast oscillations 
 Tpon = (sin(omega*t))**2

  ! intensity envelope
  a02 = vosc**2
  if (t <= 2*tpulse) then 
   intensity = a02*max(0.,sin(pi*t/2./tpulse)**2) 
  else 
    intensity = 0.
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

 ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
  Z_R = omega*sigma0**2

  gamma_s = sqrt(1.+4*intensity/rho_upper)  ! gamma factor for EM solution at x=xc
  wp_r = sqrt(rho0_up/gamma_s)  ! effective plasma frequency of upper shelf
  ls_r = 1./wp_r   ! rel. skin depth

  phi = atan(-omega/wp_r) ! Interface phase factor given by tan(phi) = -k * l_s
  chi = omega*x + phi  ! Vacuum phase 

  if (x.le.0) then
     xf = omega*sin(2*chi)
     Xpon = sin(chi)    ! laser 'Ez'
     sigma = sigma0*sqrt(1.+abs(x)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam

  else
     xf = -2/ls_r*sin(phi)**2*exp(-2*x/ls_r)
     Xpon = sin(phi)*exp(-x/ls_r)    ! laser Ez inside
     sigma = sigma0  ! Dont expand spot inside target
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

! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute relativistic fpond for standing wave field, linear rise-time
!
! ==================================================================

subroutine fpond_lin(t,tpulse,sigma_in,vosc,omega,rho_upper,x,y,z,epon_x,epon_y,epon_z,phipon)

  implicit none
  real, intent(in) :: t ! time
  real, intent(in) :: tpulse ! pulse duration/rise time
  real, intent(in) :: vosc ! quiver strength
  real, intent(in) :: sigma_in ! pulse width (1/e)
  real, intent(in) :: omega ! laser frequency
  real, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)
  real, intent(in) :: rho_upper
  real, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

  real :: xf, yf, zf, Tpon, Rpon, Xpon, intensity, gamma, sigma, atten, theta
  real :: r, pi=3.141592654, phi, chi, a02, sigma0
  real :: rho0_up, wp_r, ls_r, gamma_s, f2d
  real :: Xprop
  real :: Z_R  ! Rayleigh length

  ! fast oscillations 
 Tpon = (sin(omega*t))**2

  ! intensity envelope
  a02 = vosc**2
  intensity = a02*min(1.,t/tpulse) 

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

 ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
  Z_R = omega*sigma0**2

  gamma_s = sqrt(1.+4*intensity/rho_upper)  ! gamma factor for EM solution at x=xc
  wp_r = sqrt(rho0_up/gamma_s)  ! effective plasma frequency of upper shelf
  ls_r = 1./wp_r   ! rel. skin depth

  phi = atan(-omega/wp_r) ! Interface phase factor given by tan(phi) = -k * l_s
  chi = omega*x + phi  ! Vacuum phase 

  if (x.le.0) then
     xf = omega*sin(2*chi)
     Xpon = sin(chi)    ! laser 'Ez'
     sigma = sigma0*sqrt(1.+abs(x)**2/Z_R**2) !  Vacuum spot size for Gaussian beam
     Xprop = 1.
  else
     xf = -2/ls_r*sin(phi)**2*exp(-2*x/ls_r)
     Xpon = sin(phi)*exp(-x/ls_r)    ! laser Ez inside
     sigma = sigma0  ! Dont expand spot inside target
     Xprop = exp(-5*x/ls_r)
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
!  phipon = 4*Xpon**2*Rpon*atten*intensity  ! intensity, including attenuation factor
  phipon = Xprop*Rpon*atten*intensity  ! intensity for 'demo' display, including attenuation factor

  gamma = sqrt(1.+abs(phipon*Tpon))  ! relativistic phi_pond

  Epon_x = 2*intensity*Tpon*Rpon/gamma*xf*atten
  Epon_y = f2d*2*intensity*Tpon*Xpon**2/gamma*yf*atten
  Epon_z = f2d*2*intensity*Tpon*Xpon**2/gamma*zf*atten
!  Epon_z = 0.   

end subroutine fpond_lin



! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute ponderomotive fields from Helmholtz solution for EM standing wave 
!
! ==================================================================

subroutine fpond_helm(t, sigma_in, omega, &
                      x,y,z,ux,Az,nxh,xh_start,xh_end,dxh,x_crit, &
                      epon_x,epon_y,epon_z,phipon)

  implicit none
  real, intent(in) :: t 	!< time
!  real, intent(in) :: tpulse !< pulse duration/rise time
!  real, intent(in) :: vosc !< quiver strength
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

  real*8 ::  yf, zf, xh, Rpon, gamma, sigma, atten, theta
  real*8 :: Azr_0, Azr_1, Azr_2, Azr_3  ! Vector pot. at control points
  real*8 :: epon1, epon2 ! pond force at control points
  real*8 :: ayi, azi ! vec. pot at particle
  real*8 :: epxi, epyi, epzi ! pond field at particle
  real*8 :: r, pi=3.141592654, sigma0
  real*8 :: xa, b1, b2
  real :: f2d ! 2D switch
  real :: pha  ! Temporal phase
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
     sigma = sigma0  ! Dont expand spot inside target
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
     i1 = max(int(xa),1)  ! lower NGP
     i2 = min(i1+1,nxh)       ! upper NGP
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
     Epon1 = .5*Azr_1*( Azr_2 - Azr_0 )/dxh
     Epon2 = .5*Azr_2*( Azr_3 - Azr_1 )/dxh

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
!     if (i1==1) write(*,'(a)') 'x,abs(Az(i1)),  azi,  Epon, gamma '
!     write(*,'(i5,5(f12.4))') i1,x,abs(Az(i1)),azi,Epxi,gamma
  else
     phipon=0.
     Epon_x=0.
     Epon_y=0.
     Epon_z=0.
  endif


end subroutine fpond_helm


! ==================================================================
!
!                        PONDEROMOTIVE FORCE
!
!  Compute ponderomotive fields from Helmholtz solution for EM standing wave 
!  c-pol light
!
! ==================================================================

subroutine fpond_helmc(sigma_in,omega, &
                      x,y,z,ux,Az,nxh,xh_start,xh_end,dxh,x_crit, &
                      epon_x,epon_y,epon_z,phipon)

  !   call fpond_helmc( tlaser, tpulse,sigma,vosc,omega, &
  !                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
  !	  	      epon_x,epon_y,epon_z,phipon)


  implicit none
!  real, intent(in) :: t 	!< time
!  real, intent(in) :: tpulse !< pulse duration/rise time
!  real, intent(in) :: vosc !< quiver strength
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

  real*8 ::  yf, zf, xh, Rpon, gamma, sigma, atten, theta
  real*8 :: Azr_0, Azr_1, Azr_2, Azr_3  ! Vector pot. at control points
  real*8 :: epon1, epon2 ! pond force at control points
  real*8 :: ayi, azi ! vec. pot at particle
  real*8 :: epxi, epyi, epzi ! pond field at particle
  real*8 :: r, pi=3.141592654, sigma0
  real*8 :: xa, b1, b2
  real :: f2d ! 2D switch
  real :: Z_R  ! Rayleigh length
  integer ::  i1, i2


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
     sigma = sigma0  ! Dont expand spot inside target
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
     i1 = max(int(xa),1)  ! lower NGP
     i2 = min(i1+1,nxh)       ! upper NGP
     b2=xa-i1  ! linear weights  W_j = 1-|x_i-x_j|
     b1=1.-b2

     ! Derive fpond from Az - need gradient at both reference points.
     Azr_0 = Real(Az(i1-1))
     Azr_1 = Real(Az(i1))
     Azr_2 = Real(Az(i2))
     Azr_3 = Real(Az(i2+1))

     ! pond force at reference points i1, i2
     ! - 2nd order difference without gamma factor, as in emfield
     Epon1 = .5*Azr_1*( Azr_2 - Azr_0 )/dxh
     Epon2 = .5*Azr_2*( Azr_3 - Azr_1 )/dxh

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
!     if (i1==1) write(*,'(a)') 'x,abs(Az(i1)),  azi,  Epon, gamma '
!     write(*,'(i5,5(f12.4))') i1,x,abs(Az(i1)),azi,Epxi,gamma
  else
     phipon=0.
     Epon_x=0.
     Epon_y=0.
     Epon_z=0.
  endif


end subroutine fpond_helmc
!>     ==================================
!>     
!>     Helmholtz solver for electromagnetic fields 
!>     
!>     Based on SPLIM project helmholtz.f90     
!>     ==================================


subroutine em_helmholtz(me,itime,n,dx,theta,a0,w0,rhoe,Az)

  implicit none
  integer, intent(in) :: n !< number of mesh points for 1D Helmholtz fields
  integer, intent(in) :: itime  !< current timestep 
  integer, intent(in) :: me  !< current rank 
  real, intent(in) :: theta  !< angle of incidence
  real, intent(in) ::  a0    !< laser amplitude vosc/c(t)
  real, intent(in) ::  w0    !< laser frequency (w/wp)
  real, intent(in) ::  dx    !< mesh spacing (c/wp)
  real*4, intent(in) :: rhoe(0:n+1)  !< cycle-averaged electron density
!  real*8, intent(out) :: Ezr(0:n+1), Byr(0:n+1), Azr(0:n+1), epond(0:n+1)
  complex, dimension(n) :: alpha,beta,gamma,y !< trisolver work arrays
  complex, dimension(0:n+1) :: Az,Ao	!< Vector potential
  real, dimension(n) :: eps  !< permittivity
  real :: rgam(0:n+1)  !<  relativistic gamma
  real :: err(0:n+1)   !<  error check
  complex :: yi,  carg  !< complex args
  integer :: i,j,n_iter
  real :: pi, s2th, cth, errmax !< constants

  real*8 :: ncrit  !< critical density

  real*8 :: g0, nu_eff
  integer :: iplas, itav

  itav=1
  n_iter = 3
  ncrit = w0**2  ! critical density in norm. units
  nu_eff=0.2  ! Effective collision rate
  pi = asin(1.0)*2
  yi = (0.,1.)
  err=0.
  s2th=sin(pi/180*theta)**2
  cth = cos(pi/180*theta)


! Use gamma from previous step to speed up iteration
  ao=az
  rgam(0)=1.
  rgam(n+1) = 1.
  do i=1,n
     rgam(i) = sqrt(1 + 0.5*abs(ao(i))**2)
     az(i)=(0.,0.) 
  end do


  Az(0)=(0.,0.)
  Az(n+1)=(0.,0.)

  if (a0==0) return

  do j=1,n_iter
     do i=1,n
        !  coefficients as for s-pol light
        ! rhoe normalized to nc defined on own 1D grid along laser axis
        eps(i) = 1.-rhoe(i)/ncrit/rgam(i)/(1.+yi*nu_eff)
        y(i)=(0.,0.)
        alpha(i)=1
        beta(i)=-2 + w0**2*dx**2*(eps(i)-s2th)
        gamma(i)=1
     end do

     !  BCs - note additional factor w0=k0 for phase factors 
     y(1) = 2*yi*a0*sin(w0*dx*cth)
     carg = yi*w0*dx*cth
     beta(1) = beta(1) + cexp(carg)

     call trisolve(alpha,beta,gamma,y,Az(1:n),n-1,n)

     ! relativistic factor - average old and new potentials
     errmax=0.
     do i=1,n
        rgam(i) = sqrt(1 + 0.5*abs(az(i))**2) 
        err(i) = sqrt(abs(az(i)**2-ao(i)**2)/4/a0**2)
     end do

     ! BCs for next iterate (Laplacian of gamma)
     rgam(0) = 2*rgam(1) - rgam(2)
     rgam(n+1) = rgam(n)

     errmax = maxval(err(1:n))
     Ao = Az  ! Store previous iterate

  end do

iplas = n/2
if (me==0) write(*,'(i6,2f12.3)') itime,rhoe(iplas),abs(az(iplas))
if (itime .eq. itav .and. me==0) then
  write (*,'(a20,i2,a10,f12.5,a10,f12.3)') 'Iterate ',j,' error=',errmax,' amplitude',a0
  g0 = sqrt(1+a0**2/2)
  open (40,file='a_error.dat')
  write(40,'(a)') '! x, rho, eps, az/a0, gam/g0, err'
  write(40,'(6(1pe12.3))') (dx*i,rhoe(i),eps(i),abs(az(i))/a0,rgam(i)/g0,err(i),i=1,n)
  close(40)
  
endif

 
  ! Bcs 
  Az(0) = 2*Az(1) - Az(2)  ! gradient continuous
  Az(n+1) = Az(n) ! zero in solid


end subroutine em_helmholtz

! -------------------------------------------------
!
!   Triadiagonal matrix solver.
!
!   Solves equation system:
!
!        alpha_i x_(i-1) + beta_i x_i + gamma_i x_(i+1) = y_i
!
! ------------------------------------------------

subroutine trisolve(alpha,beta,gamma,y,x,n,nmax)

 implicit none
 integer, intent(in) :: n,nmax
 complex, dimension(nmax) :: alpha, beta, gamma,y,x, q
 integer :: i

  q(1) = beta(1)
  x(1) = y(1)/q(1)

  !  forward elimination

  do i = 2,n
     q(i) = beta(i) - alpha(i)*gamma(i-1)/q(i-1)
     x(i) = (y(i) - alpha(i)*x(i-1))/q(i)
  end do

  !  back substitution

  do i=n-1,1,-1
     x(i) = x(i) - gamma(i)/q(i)*x(i+1)
  end do
end subroutine trisolve

!  =================================
!
!    Tracking n_c for laser
!
!   $Revision 1.5$
!
!  =================================

subroutine track_nc

  use module_physvars
  use module_field_grid

  implicit none
  include 'mpif.h'

  real, dimension(0:ngx+1) :: rho1d

  integer :: i, j, k, ng,jfoc, kfoc, icrit
  logical :: found
  real :: xc1, dx, dy, dz

  integer :: icm, icp, ixc, nover
  integer :: ierr

  call densities  ! Compute local ion density rhoi_loc on 3D grid

  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints

  call MPI_ALLREDUCE(rhoi_loc, rhoi, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  dy = yl/ngy
  dz = zl/ngz
  dx = xl/ngx

! density average line-out along laser axis: 5x5 average, converted to n/nc

  jfoc = max(2.0,focus(2)/dy) ! ensure index>=0
  kfoc = max(2.0,focus(3)/dz)
  rho1d(0:ngx+1) = 0.
  do k=kfoc-2,kfoc+2
     do j=jfoc-2,jfoc+2
        rho1d(1:ngx) = rho1d(1:ngx)+rhoi(1:ngx,j,k)/25./omega**2  
     end do
  end do


  ! Determine position of critical density along laser axis
 
! initial leading edge of plasma
  if (target_geometry == 1 .or. target_geometry ==3) then
     xc1 = plasma_centre(1)-r_sphere 
  else
     xc1 = plasma_centre(1)-x_plasma - x_layer(1)  
  endif

! abort tracking if particle number too low
  if (npart_total.le.100) then
	write(*,*) 'fixing xcrit at ',xc1
     x_crit=xc1
     rho_upper=1.0
     return
  endif

 !  start indice
  ixc = xc1/dx
  icm = min(max(ixc,1),ngx-1)
  icp = icm+1
  icrit = icm
  found = .false.
  rho_track = 1.5   ! tracking density

!  sweep down
  i = icm
  do while ( i > 0 .and. .not.found)
     if ( rho1d(i)<= rho_track .and. rho1d(i+1)> rho_track ) then
        found=.true.
        x_crit = i*dx + dx*(1.-rho1d(i))/(rho1d(i+1)-rho1d(i))
        icrit=i
     endif
     i = i-1
  end do

!  sweep up
  i = icm
  do while ( i < ngx .and. .not.found)
     if (rho1d(i)<= rho_track .and. rho1d(i+1)> rho_track) then
        found=.true.
        x_crit = i*dx + dx*(1.-rho1d(i))/(rho1d(i+1)-rho1d(i))
        icrit = i
     endif
     i = i+1
  end do

  if (.not.found .and. itime>0) then 
	beam_config=0 ! switch off laser
	if (my_rank==0) write(6,*) 'Target burnt through - switching off laser'
  endif

  rho_upper=0.
  nover = 5./dx

  do i=icrit+1,min(icrit+nover,ngx)
     rho_upper = rho_upper + rho1d(i)
  end do

  rho_upper = rho_upper/nover  ! ave. upper shelf density

!  if (.not.found) x_crit=xc1  ! original plasma edge 
!  if (.not.found) x_crit=x_crit + dt  ! vacuum propagation
 
if (my_rank==0) then
   write(*,'(/a15,f10.3,a15,f10.3,a15,f10.3)') &
        'plasma edge: ',xc1, ' x_crit: ',x_crit,' n_upper: ',rho_upper
   write(24,*) 'plasma edge: ',xc1, ' x_crit: ',x_crit

endif


end subroutine track_nc

!  =================================
!
!    3D Density gather for rhoi used in Helmholtz solver
!  - 1D grid with 3x3 y-z cells for dummy counts
!
!  =================================

subroutine density_helmholtz

  use module_physvars
  use module_particle_props

  implicit none
  include 'mpif.h'

  real :: rdx, rdy, rdz, dy, dz, cweight
  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za
  real :: yh_start, zh_start ! Start of HH grid in transverse directions
  integer, parameter :: nyh=4 ! # additional points in transverse direction
  integer :: i, ng, i1, i2, j1, j2, k1, k2
  integer :: ierr
  real*4 :: rho_loc(0:nxh+1,0:nyh,0:nyh), rho_glob(0:nxh+1,0:nyh,0:nyh)
  real*4 :: charge_sum, charge_tot

! Helmholtz grid limits
  dxh = (xh_end-xh_start)/nxh
! tranverse resolution defined by particle spacing 
  dy = max(dxh,5*a_ii)
  dz = max(dxh,5*a_ii)
!  dy = dxh
!  dz = dxh
  yh_start = focus(2)-nyh/2*dy
  zh_start = focus(3)-nyh/2*dz

  rdx = 1./dxh
  rdy = 1./dy
  rdz = 1./dz


  !  Any particle outside gets put in ghost cells 0, 2

  cweight = qi*rdx*rdy*rdz       ! charge weighting factor
!  if (me==0)   write(6,'(//a,3f12.3)') 'cw,dx,dy',cweight,dxh,dy

  rho_loc(0:nxh+1,0:nyh,0:nyh) = 0.

  do i=1,np_local

     xa=(x(i) - xh_start)*rdx
     ya=(y(i) - yh_start)*rdy
     za=(z(i) - zh_start)*rdz

     !  indices
     i1=xa+1
     i2=i1+1
     j1=ya
     j2=j1+1
     k1=za
     k2=k1+1

     i1 = min(max(0,i1),nxh+1)
     i2 = min(max(0,i2),nxh+1)
     j1 = min(max(0,j1),nyh)
     j2 = min(max(0,j2),nyh)
     k1 = min(max(0,k1),nyh)
     k2 = min(max(0,k2),nyh)

     !  linear weighting
     fx2=i1-xa  ! Prevent overflow/negative weighting for particles outside box
     fx1=1.-fx2
     fy2=ya-j1
     fy1=1.-fy2
     fz2=za-k1
     fz1=1.-fz2

     !  gather ion charge at nearest grid points
     if (q(i)>0) then

        rho_loc(i1,j1,k1)=rho_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
        rho_loc(i2,j1,k1)=rho_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
        rho_loc(i1,j2,k1)=rho_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
        rho_loc(i2,j2,k1)=rho_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
        rho_loc(i1,j1,k2)=rho_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
        rho_loc(i2,j1,k2)=rho_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
        rho_loc(i1,j2,k2)=rho_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
        rho_loc(i2,j2,k2)=rho_loc(i2,j2,k2) + cweight*fx2*fy2*fz2
     else

     endif
  end do

  ng = (nxh+2)*(nyh+1)**2                         ! total # gridpoints
! gather on all
  call MPI_ALLREDUCE(rho_loc, rho_glob, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

! Store density lineout for Helmholtz solver
  do i=0,nxh+1
    rho_helm(i) = rho_glob(i,nyh/2,nyh/2)
  end do 
 charge_sum = SUM(rho_glob(1:nxh,nyh/2,nyh/2))/cweight
 charge_tot = SUM(rho_glob)/cweight
if (my_rank==0 .and. itime==1) then
  write (6,*) "Charge sum on HH grid lineout/total:",charge_sum,charge_tot
!  write(6,*) "rho20",(rho_glob(i,nyh/2,0),i=1,nxh)
!  write(6,*) "rho21",(rho_glob(i,nyh/2,1),i=1,nxh)
!  write(6,*) "rho22",(rho_glob(i,nyh/2,2),i=1,nxh)
!  write(6,*) "rho23",(rho_glob(i,nyh/2,3),i=1,nxh)
!  write(6,*) "rho24",(rho_glob(i,nyh/2,4),i=1,nxh)
endif
end subroutine density_helmholtz

end module module_laser
