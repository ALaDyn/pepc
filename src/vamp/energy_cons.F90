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
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none

  real :: epot, ekine, ekini, ebeam,  etot, Qplas, conv_kev

!VAMPINST subroutine_start
       CALL VTENTER(IF_energy_cons,VTNOSCL,VTIERR)
!      write(*,*) 'VT: energy_cons S>',VTIERR,
!     *    IF_energy_cons,ICLASSH
!
  call potenergy(epot)
  call kinenergy(ekine, ekini, ebeam)

  etot = epot + ekine + ekini + ebeam
  Qplas = abs(qe)*ne
  conv_kev = 2./3./Qplas*511

  if ( me == 0 ) then
     do ifile = 6,15,9
        write (ifile,'(5(a,1pe12.5))') 'P.E. = ',epot,' elec K.E. = ',ekine,' Ions ',ekini,' Beam = ',ebeam,' Total: ',etot
        write (ifile,'(a10,2(1pe12.5),a4)') 'Plasma Te, Ti ',conv_kev*ekine,conv_kev*ekini,' keV'
     end do
     write (75,'(f12.5,6(1pe12.3))') (itime+itime_start)*dt, conv_kev*epot, conv_kev*ekine, conv_kev*ekini, conv_kev*ebeam, conv_kev*etot,x_crit
  endif
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: energy_cons S<',VTIERR,ICLASSH
!
end subroutine energy_cons
