!  ================================
!
!         ENERGY_CONS
!
!     Find potential, kinetic energies
!
!  $Revision 1.6$
!
!  ================================


subroutine energy_cons(ekine,ekini,emag,ebeam)

  use physvars
  use treevars
  use utils

  implicit none

  real :: tpon, epot, ekine, ekini, ebeam,  etot, Qplas
  real :: emag

  call potenergy(epot,emag)
  call kinenergy(ekine, ekini, ebeam)

  etot = epot + emag + ekine + ekini + ebeam




  laser_energy: select case(beam_config)
  case(4)
     elaser = 3./8.*omega**2*sigma**2*vosc**2*tlaser
  case(5)
     elaser = 3./8.*omega**2*sigma**2*vosc**2*tpulse
  case default
     elaser = 0
  end select laser_energy

  if ( me == 0 ) then
     do ifile = 6,15,9
        write (ifile,'(7(a20,1pe12.5/))') &
	     ' P.E. = ',epot, &
	     ' Magnetic E. = ',emag, &
	     ' Electron K.E. = ',ekine, &
             ' Ion K.E. = ',ekini, &
	     ' Beam K.E.  = ',ebeam, &
	     ' Total: ',etot, &
             ' Laser energy = ',elaser

        write (ifile,'(2(a20,f12.5/))') 'Plasma Te (keV):',convert_kev*ekine,'Ti (keV):',convert_kev*ekini
     end do
! Write out to energy.dat file
if (itime.eq.0)  write(75,'(a)') '! time  Upot  Umag  Ukin_e Ukin_i Ukin_beam Utot Tpon xc'
     write (75,'(f12.5,8(1pe12.3))') trun, convert_kev*epot, convert_kev*emag, convert_kev*ekine, convert_kev*ekini,&
	 convert_kev*ebeam, convert_kev*etot,tpon,x_crit
  endif
end subroutine energy_cons


