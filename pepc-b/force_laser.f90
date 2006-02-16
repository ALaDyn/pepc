!  ===================================================================
!
!                              FORCE_LASER
!
!   $Revision$
!
!   Calculate forces due to laser fields
!
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
  real :: Epon_x, Epon_y, Epon_z, Phipon, ex_em, ey_em, ez_em, bx_em, by_em, bz_em

  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  ! Include force from laser/external field on electrons - ES scheme

  if (beam_config >= 3 .and. beam_config <=6 ) then
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

           case(14)  ! standing wave fpond, linear rise-time
              call fpond_lin( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                xd,yd,zd,epon_x,epon_y,epon_z,phipon)
              epon_y=epon_y/100.
              epon_z=epon_z/100.
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(94)  ! standing wave fpond with transverse fields artificially reduced
              call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                xd,yd,zd,epon_x,epon_y,epon_z,phipon)
              epon_y=epon_y/100.
              epon_z=epon_z/100.
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

           case(24)  ! oblique incidence standing wave, s-pol
              call emobliq( tlaser, tpulse,sigma,vosc,omega,theta_beam,rho_upper, &
                xd,yd,zd,epon_x,epon_y,epon_z,phipon,ez_em,bx_em,by_em)
	      Bz_em = 0.

           case(5)  ! propagating fpond
              call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, & 
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon)
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = 0.

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

	if (beam_config_in==6)  then ! Constant B in z-direction - either charge
              Epon_x=0.
              Epon_y=0.
              Epon_z=0.
	      Bx_em = 0.
	      By_em = 0.
	      Bz_em = vosc
	endif

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




