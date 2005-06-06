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
  real :: tpon, epot, ekine, ekini, ebeam,  etot, Qplas
  real :: emag

  call potenergy(epot)
  call kinenergy(ekine, ekini)

  etot = epot + emag + ekine + ekini



  if ( my_rank == 0 .and. db_level.ge.1 ) then
     do ifile = 6,15,9
        write (ifile,'(4(a20,1pe12.5/))') &
	     ' P.E. = ',epot, &
	     ' Electron K.E. = ',ekine, &
             ' Ion K.E. = ',ekini, &
	     ' Total: ',etot

     end do
     ! Write out to energy.dat file
     if (itime.eq.0)  write(75,'(a)') '! time  Upot  Ukin_e Ukin_i Utot '
     write (75,'(f12.5,4(1pe13.4))') trun, epot, ekine, ekini, etot
  endif
end subroutine energy_cons



