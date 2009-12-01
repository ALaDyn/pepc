!  ================================
!
!         ENERGY_CONS
!
!     Find potential, kinetic energies
!
!  $Revision 1.6$
!
!  ================================


subroutine energy_cons(ekine,ekini)

  use physvars
  !  use treevars
  use utils

  implicit none

  integer :: ifile
  real*8 :: epot, ekine, ekini, etot
!  real :: emag = 0.0

  call potenergy(epot)
  call kinenergy(ekine, ekini)

  etot = epot + ekine + ekini


  if (my_rank == 0) then
!  if ( my_rank == 0 .and. db_level.ge.1 ) then
!     do ifile = 6,15,9
!        write (ifile,'(4(a20,1pe12.5/))') &
!	     ' P.E. = ',epot, &
!	     ' Electron K.E. = ',ekine, &
!             ' Ion K.E. = ',ekini, &
!	     ' Total: ',etot
!
!     end do
     ! Write out to energy.dat file
     open(75,file='energy.dat',STATUS='UNKNOWN', POSITION = 'APPEND') 
     if (itime.eq.0)  write(75,'(a)') '! time  Upot  Ukin_e Ukin_i Utot '
     write (75,'(f12.5,5(1pe13.4))') trun, epot, ekine, ekini, ekine+ekini, etot
     close(75)
  endif
end subroutine energy_cons



