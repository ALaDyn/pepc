!  ================================
!
!         ENERGY_CONS
!
!     Find potential, kinetic energies
!
!  ================================


subroutine energy_cons

  use treevars
  use utils

  implicit none

  real :: epot, ekine, ekini, ebeam,  etot, Qplas, conv_kev

  call potenergy(epot)
  call kinenergy(ekine, ekini, ebeam)

  etot = epot + ekine + ekini + ebeam

  if (ne>0) then
     Qplas = abs(qe)*ne
  else
     Qplas = abs(qi)*ni
  endif

  conv_kev = 2./3./Qplas*511

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
        write (ifile,'(7(a20,1pe12.5/))') 'P.E. = ',epot,' Electron K.E. = ',ekine, &
             ' Ion K.E. = ',ekini,' Beam K.E.  = ',ebeam,' Total: ',etot, &
             ' Laser energy = ',elaser, &
             ' Laser focus = ',focus(1)
        write (ifile,'(2(a20,f12.5/))') 'Plasma Te (keV):',conv_kev*ekine,'Ti (keV):',conv_kev*ekini
     end do
     write (75,'(f12.5,6(1pe12.3))') (itime+itime_start)*dt, conv_kev*epot, conv_kev*ekine, conv_kev*ekini, conv_kev*ebeam, conv_kev*etot,x_crit
  endif
end subroutine energy_cons


