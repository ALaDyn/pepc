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
  elaser = 3./8.*omega**2*sigma**2*vosc**2*tlaser

  if ( me == 0 ) then
     do ifile = 6,15,9
        write (ifile,'(5(a,1pe12.5/))') 'P.E. = ',epot,' elec K.E. = ',ekine, &
             ' Ions ',ekini,' Beam = ',ebeam,' Total: ',etot, &
             ' Deposited laser energy = ',elaser 
        write (ifile,'(a10,2(1pe12.5),a4)') 'Plasma Te, Ti ',conv_kev*ekine,conv_kev*ekini,' keV'
     end do
     write (75,'(f12.5,6(1pe12.3))') (itime+itime_start)*dt, conv_kev*epot, conv_kev*ekine, conv_kev*ekini, conv_kev*ebeam, conv_kev*etot,x_crit
  endif
end subroutine energy_cons


