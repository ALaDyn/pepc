!  ================================
!
!         ENERGY_CONS
!
!     Find potential, kinetic energies
!
!  ================================


subroutine energy_cons

  use treevars
  use physvars
  use utils

  implicit none

  real*8 :: epot_dust, epot_star, ekin_dust, ekin_star,  etot, conv_e
  integer :: ifile

  call potenergy(epot_dust,epot_star)
  call kinenergy(ekin_dust, ekin_star)

  etot = epot_dust + epot_star + ekin_dust + ekin_star 

!  energy conversion factor for output
  conv_e = 1


  if ( me == 0 .and. mod(itime,iprot)==0 ) then
     do ifile = 6,15,9
        write (ifile,'(6(a20,1pe12.5/))') &
        ' dust P.E. = ',epot_dust, &
        ' star P.E. = ',epot_star, &
        ' dust K.E. = ',ekin_dust, &
        ' star K.E. = ',ekin_star, &
        ' Total: ',etot
     end do
     if (itime==0) write (75,'(a)') '! time pot_dust, pot_stars, kin_dust, kin_stars, total'
     write (75,'(f12.5,6(1pe13.5))') (itime+itime_start)*dt, epot_dust, epot_star, ekin_dust, ekin_star, etot

  endif
end subroutine energy_cons




