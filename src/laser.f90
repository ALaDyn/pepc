subroutine laser

  use physvars
  use treevars
  implicit none
  integer :: ierr

  !  Laser focal position and rezoning

  laser_focus: select case(beam_config)

  case(4)  ! standing wave fpond
     if (itime>0) focus(1) = x_crit  ! laser tracks n_c

  case(5)  ! propagating fpond
     !  Trigger rezoning if laser rezone_frac of the way through plasma
     ! - Only works after restart at present
     if (restart .and. beam_config ==5 .and. focus(1) >= window_min + x_plasma*rezone_frac) then
        if (me==0) then
           write (*,*) 'REZONE'
           !           read (*,*) go
        endif
        call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

        call rezone
        !        window_min = window_min + dt
     else
        focus(1) = focus(1) + dt  ! propagate forward by c*dt - can include v_g here
        propag_laser=propag_laser + dt
     endif

  case default
     ! leave focal point unchanged

  end select laser_focus

end subroutine laser
