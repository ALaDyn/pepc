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

  use physvars
  use treevars
  use utils
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
              call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon)
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
              call fpond_helm( tlaser, tpulse,sigma,vosc,omega, &
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



