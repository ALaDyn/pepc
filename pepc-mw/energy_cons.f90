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
  use module_fmm_framework
  use module_laser
  use module_pusher
  use module_io
  use module_units

  implicit none

  real*8 :: epot, ekine, ekini, etot, tempe, tempi
  integer :: ifile

  call potenergy(epot)
  call kinenergy(ekine, ekini, tempe, tempi)

  ! rescale total potential and kinetic energy to energy per particle
  epot = epot / npart_total
  ekine = ekine / ne
  ekini = ekini / ni
  potnearfield = potnearfield / npart_total
  potfarfield  = potfarfield  / npart_total

  etot = epot + ekine + ekini

  laser_energy: select case(beam_config)
  case(4)
     elaser = 3./8.*omega**2*sigma**2*vosc**2*tlaser
  case(5)
     elaser = 3./8.*omega**2*sigma**2*vosc**2*tpulse
  case default
     elaser = 0
  end select laser_energy

  if (my_rank == 0) then
     do ifile = 6,15,9
        write (ifile,'(/"-- ENERGY --"/20x,2a20/7(a20,1pe18.8,8x,1pe18.8/))') &
	     'norm  ','eV  ', &
         ' P.E. (total)      = ', epot,                 epot*unit_Ryd_in_eV, &
         ' P.E. (near field) = ', potnearfield, potnearfield*unit_Ryd_in_eV, &
         ' P.E. (far field)  = ', potfarfield,   potfarfield*unit_Ryd_in_eV, &
	     ' K.E. (electrons)  = ', ekine,               ekine*unit_Ryd_in_eV, &
         ' K.E. (ions)       = ', ekini,               ekini*unit_Ryd_in_eV, &
	     ' K.E. (total)      = ', etot,                 etot*unit_Ryd_in_eV, &
         ' Laser energy      = ', elaser,             elaser*unit_Ryd_in_eV

        write (ifile,'(2(a20,e18.8,8x,e18.8/))') &
                 'w/wo drift Te (eV):',unit_Ryd_in_eV*ekine/unit_kB*2./3., unit_Ryd_in_eV*tempe, &
                            'Ti (eV):',unit_Ryd_in_eV*ekini/unit_kB*2./3., unit_Ryd_in_eV*tempi
     end do
     ! Write out to energy.dat file
       if (itime.eq.1)  write(file_energy_dat,'(a)') '# time  Upot(total)  Upot(near field) Upot(far field)  Ukin_e Ukin_i Ukin_tot Ukine(wo drift) Ukini(wo drift) Utot Te0  Te_uncor  chie  delta_Te  Ti0  Ti_uncor  chii  delta_Ti'
       write (file_energy_dat,'(f12.5,17(1pe13.4))') trun*unit_t0_in_fs, epot, potnearfield, potfarfield, ekine, ekini, ekine+ekini, 3./2.*unit_kB*tempe, 3./2.*unit_kB*tempi, etot, &
             Te0, Te_uncor, chie, delta_Te, Ti0, Ti_uncor, chii, delta_Ti
  endif
end subroutine energy_cons



