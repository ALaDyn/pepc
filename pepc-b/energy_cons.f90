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
  !  use treevars
  use utils

  implicit none

  integer :: ifile
  real*8 :: tpon, epot, ekine, ekini, ebeam,  etot
  real*8 :: emag

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

  if ( my_rank == 0 .and. debug_level.ge.1 ) then
     do ifile = 6,15,9
        write (ifile,'(20x,3a20/6(a20,1pe18.8,8x,1pe12.5,8x,1pe12.5/),a20,1pe18.8)') &
	     'norm  ','keV ','erg ', &
	     ' P.E. = ',epot, epot*convert_keV, epot*convert_erg, &
	     ' Magnetic E. = ',emag, emag*convert_keV, emag*convert_erg, &
	     ' Electron K.E. norm/erg: ',ekine, ekine*convert_keV, ekine*convert_erg, &
             ' Ion K.E. norm/erg: ',ekini, ekini*convert_keV, ekini*convert_erg, &
	     ' Beam K.E.  = ',ebeam, ebeam*convert_keV, ebeam*convert_erg, &
	     ' Total: ',etot, etot*convert_keV, etot*convert_erg, &
             ' Laser energy = ',elaser

        write (ifile,'(2(a20,f12.5/))') 'Plasma Te (keV):',convert_kev*ekine/max(1,ne),'Ti (keV):',convert_kev*ekini/max(1,ni)
     end do
     ! Write out to energy.dat file
     if (current_step.eq.1)  write(75,'(a)') '! time  Upot  Umag  Ukin_e Ukin_i Ukin_beam Utot Tpon xc'
     write (75,'(f12.5,5(1pe12.3),1pe13.5)') &
          trun, convert_kev*epot, convert_kev*emag, convert_kev*ekine, convert_kev*ekini,&
          convert_kev*ebeam, convert_kev*etot
  endif
end subroutine energy_cons


