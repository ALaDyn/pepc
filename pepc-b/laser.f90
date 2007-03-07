subroutine laser(I_laser)

  use physvars
  implicit none

  integer :: ierr
  real :: I_laser
  real :: amplitude, ampl_min, pha
  real :: Azr(1:nxh)
  complex :: yi=(0.,1.)

  if (beam_config >= 3) tlaser = tlaser + dt  
  pha = tlaser*omega

!  Laser pulse envelope
  laser_model: select case(beam_config_in)

     case(4,34,44,6)  	! sin^2 standing wave
	if (tlaser<2*tpulse)  then
           I_laser = vosc**2*max(0.,sin(pi*tlaser/2./tpulse)**2)
	else
	   I_laser=0.
	endif
     case(14,94,54,16) 	! standing wave, linear rise
           I_laser = vosc**2*min(1.,tlaser/tpulse)
     case(3) 		! constant
        I_laser = vosc**2
     case default
	I_laser = 0.
  end select laser_model

  !  Laser focal position and rezoning

  laser_focus: select case(beam_config)

  case(4)  ! Helmholtz solver for vector potential

! Factor-in pulse shape
     amplitude = sqrt(I_laser)  ! amplitude of incoming wave
     call density_helmholtz
     call em_helmholtz(my_rank,itime,nxh,dxh,theta_beam,amplitude,omega,rho_helm,Az_helm)
     Azr(1:nxh) = Real(Az_helm(1:nxh)*cexp(yi*pha))
     ampl_max = maxval(Azr)
     ampl_min = minval(Azr)
     if (ampl_max.lt.abs(ampl_min)) ampl_max=-ampl_min
 
  case(5)  ! propagating fpond
     !  Trigger rezoning if laser rezone_frac of the way through plasma
     ! - Only works after restart at present
     if (restart .and. beam_config ==5 .and. focus(1) >= window_min + x_plasma*rezone_frac) then
        if (my_rank==0) then
           write (*,*) 'REZONE'
           !           read (*,*) go
        endif

        call rezone
        !        window_min = window_min + dt
     else
        focus(1) = focus(1) + dt  ! propagate forward by c*dt - can include v_g here
        propag_laser=propag_laser + dt
     endif

  case(8)  ! old fpond model
     if (itime>0) focus(1) = x_crit  ! laser tracks n_c

  case default
     ! leave focal point unchanged

  end select laser_focus

end subroutine laser






