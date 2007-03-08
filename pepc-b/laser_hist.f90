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

  use physvars
  !  use treevars
  use utils

  implicit none

  integer :: ifile
  real*8 :: emag

if (my_rank.eq.0) then
  ! Write out to energy.dat file
  if (itime.eq.1)  write(71,'(a)') '! time  a_L fpond xc'
  write (71,'(f12.5,2(1pe12.3))') &
       tlaser, ampl_max, fpon_max 
  write (*,'(f12.5,2(1pe12.3))') &
       tlaser, ampl_max, fpon_max
endif
end subroutine laser_hist


