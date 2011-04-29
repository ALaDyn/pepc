! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure
  use physvars
  implicit none

  call special_start(ispecial)

  ! rescale coordinates, since the special_start initializes 1x1x1-boxes
  if ((ispecial .ne. 7) .and. (ispecial .ne. 123)) then
    x = x * x_plasma
    y = y * y_plasma
    z = z * z_plasma
  endif

  !write(*,*) "Particle positions: "
  !do i=1,np_local
  !  write(*,*) my_rank,i,x(i), y(i), z(i), q(i), pelabel(i)
  !end do
  !flush(6)
  !stop

end subroutine configure

