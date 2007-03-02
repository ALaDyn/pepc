subroutine laser(I_laser)

  use physvars
  implicit none

  integer :: ierr
  real :: I_laser
  real :: amplitude


  if (beam_config >= 3) tlaser = tlaser + dt  

!  Laser pulse envelope
  laser_model: select case(beam_config_in)

     case(4,34,44)  	! sin^2 standing wave
	if (tlaser<2*tpulse)  then
           I_laser = 4*vosc**2*max(0.,sin(pi*tlaser/2./tpulse)**2)
	else
	   I_laser=0.
	endif
     case(14,94,54) 	! standing wave, linear rise
           I_laser = 4*vosc**2*min(1.,tlaser/tpulse)
     case(3) 		! constant
        I_laser = vosc**2
     case(6) 		! sin^2, propagating
        I_laser = vosc**2*max(0.,sin(pi*tlaser/2./tpulse)**2)
     case(16)		! propag, linear rise
        I_laser = vosc**2*min(1.,tlaser/tpulse)
     case default
	I_laser = 0.
  end select laser_model

  !  Laser focal position and rezoning

  laser_focus: select case(beam_config)

  case(4)  ! Helmholtz solver for vector potential

! Factor-in pulse shape
     amplitude = sqrt(I_laser/4)
     call density_helmholtz
     call em_helmholtz(my_rank,itime,nxh,dxh,theta_beam,amplitude,omega,rho_helm,Az_helm)

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






