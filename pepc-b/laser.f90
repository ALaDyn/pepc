subroutine laser

  use physvars
  implicit none

  integer :: ierr
  real :: amplitude

  !  Laser focal position and rezoning

  laser_focus: select case(beam_config)

  case(4)  ! Helmholtz solver for vector potential

! Factor-in pulse shape
    if (tlaser <= 2*tpulse) then 
      amplitude = vosc*max(0.,sin(pi*tlaser/2./tpulse)) 
    else 
      amplitude = 0.
    endif

     call density_helmholtz
     call em_helmholtz(itime,nxh,dxh,theta_beam,amplitude,omega,rho_helm,Az_helm)

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






