!  ===================================================================
!
!                              VELOCITIES
!
!   $Revision$
!
!   Calculate velocities from accelerations
! 
!   apply thermodynamic constraint according to ensemble
!
!
!
!  ===================================================================


subroutine velocities(p_start,p_finish,delta_t)


  use treevars
  use utils
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none
  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer p, i, ne_loc
  real, dimension(nppm) :: uhx, uhy, uhz
  real :: sum_vxe, sum_vye, sum_vze, sum_v2e, sum_2ve, Te0, Te_uncor, Ti0, Ti_uncor, chie, chii
  real :: sum_vxi, sum_vyi, sum_vzi, sum_v2i, sum_2vi, mass_eqm
  real :: global_v2e, global_v2i, gammah, delta_Te, delta_Ti, Te_loc

!  Available ensemble modes
!      1 = NVE - total energy conserved
!      2 = NVT - global Te, Ti conserved
!      3 = global NVT electron Te conserved; ions frozen
!      4 = local NVT: each PE keeps Te clamped; ions frozen
!      5 = local NVT, ions only; electrons added at end of run



!VAMPINST subroutine_start
       CALL VTENTER(IF_velocities,VTNOSCL,VTIERR)
!      write(*,*) 'VT: velocities S>',VTIERR,
!     *    IF_velocities,ICLASSH
!
  if (ensemble == 2) then
     ! Conserve kinetic energies of electrons and ions (initial Te const)
     ! adapted from
     !  Allen and Tildesley p230, Brown & Clark, Mol. Phys. 51, 1243 (1984)

     !  Definitions:
     !
     !   unconstrained velocities    uh(x,y,z) = v'(x,y,z)

     !  1)  Unconstrained half-step for electrons

     sum_vxe=0.0  ! partial sums
     sum_vye=0.0
     sum_vze=0.0
     sum_v2e=0.0
     sum_vxi=0.0  ! partial sums (ions)
     sum_vyi=0.0
     sum_vzi=0.0
     sum_v2i=0.0

     do p=1,npp
           uhx(p) = ux(p) + 0.5*delta_t*ax(p)
           uhy(p) = uy(p) + 0.5*delta_t*ay(p)
           uhz(p) = uz(p) + 0.5*delta_t*az(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
        if (pelabel(p)<=ne) then
           ! electrons
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        else if (pelabel(p)<=ne+ni) then
           ! ions
           sum_vxi  = sum_vxi  + uhx(p)/gammah
           sum_vyi  = sum_vyi  + uhy(p)/gammah
           sum_vzi  = sum_vzi  + uhz(p)/gammah
           sum_v2i = sum_v2i + gammah-1.
        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2
     sum_2vi = sum_vxi**2 + sum_vyi**2 + sum_vzi**2

     ! Find global KE sums
     call MPI_ALLREDUCE(sum_v2e, global_v2e, one, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(sum_v2i, global_v2i, one, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*global_v2e/ne  ! This should equal 3/2 kT for 3v Maxwellian
     Te0 = Te_keV  ! normalised electron temp

     Ti_uncor = 511*mass_i/mass_e*2./3.*global_v2i/ni  ! This should equal 3/2 kT for 3v Maxwellian
     Ti0 = Ti_keV  ! normalised electron temp


     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
     chii = sqrt(abs(Ti0/Ti_uncor)) 

     if (me==0)	write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie
     if (me==0)	write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii

     !  3)  Complete full step

     do p=1,npp
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*ax(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*ay(p)
           uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*az(p)

        elseif (pelabel(p)<=ne+ ni) then
           ux(p) = (2*chii-1.)*ux(p) + chii*delta_t*ax(p)
           uy(p) = (2*chii-1.)*uy(p) + chii*delta_t*ay(p)
           uz(p) = (2*chii-1.)*uz(p) + chii*delta_t*az(p)
        endif
     end do




     delta_Ti=0.


     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating
     !     heate = heate + delta_Te


  else if (ensemble ==3) then

! electrons clamped, ions frozen

     sum_vxe=0.0  ! partial sums
     sum_vye=0.0
     sum_vze=0.0
     sum_v2e=0.0
     ne_loc = 0  ! # local electrons

     do p=1,npp

        if (pelabel(p)<=ne) then
           ! electrons
           ne_loc = ne_loc + 1
           uhx(p) = ux(p) + 0.5*delta_t*ax(p)
           uhy(p) = uy(p) + 0.5*delta_t*ay(p)
           uhz(p) = uz(p) + 0.5*delta_t*az(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2

     ! Find global KE sums
     call MPI_ALLREDUCE(sum_v2e, global_v2e, one, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*global_v2e/ne  ! This should equal 3/2 kT for 3v Maxwellian
     Te_loc = 511*2./3.*sum_v2e/ne_loc
     Te0 = Te_keV  ! normalised electron temp
     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
     chie = min(1.,max(chie,0.5))
     if (me==0)	write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie
     write(*,*) 'Tloc ',Te_loc

     !  3)  Complete full step

     do p=1,npp
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*ax(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*ay(p)
           uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*az(p)
        endif
     end do

     delta_Ti=0.
     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating

  else if (ensemble ==4) then


! electrons clamped locally, ions frozen
! - require T=T_e on each PE to avoid local drifts

     sum_vxe=0.0  ! partial sums
     sum_vye=0.0
     sum_vze=0.0
     sum_v2e=0.0
     ne_loc = 0  ! # local electrons

     do p=1,npp

        if (pelabel(p)<=ne) then
           ! electrons
           ne_loc = ne_loc+1  
           uhx(p) = ux(p) + 0.5*delta_t*ax(p)
           uhy(p) = uy(p) + 0.5*delta_t*ay(p)
           uhz(p) = uz(p) + 0.5*delta_t*az(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2

     ! Find local KE temp

    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*sum_v2e/ne_loc  ! This should equal 3/2 kT for 3v Maxwellian
     Te0 = Te_keV  ! normalised electron temp
     ! exponent of chie should be 1/2 - take 1/3 to soften oscillations
     chie = (abs(Te0/Te_uncor))**0.5     ! multipliers from Temperature ratio - reset once every cycle

     chie = min(1.25,max(chie,0.75))  ! Set bounds of +- 50%

     write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

     !  3)  Complete full step

     do p=1,npp
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*ax(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*ay(p)
           uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*az(p)
        endif
     end do

     delta_Ti=0.
     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


  else if (ensemble == 5) then
     ! Conserve kinetic energy of ions only (initial Ti const)
     mass_eqm = 20.  ! artificial ion mass for eqm stage
     sum_vxi=0.0  ! partial sums (ions)
     sum_vyi=0.0
     sum_vzi=0.0
     sum_v2i=0.0

     do p=1,npp
           uhx(p) = ux(p) + 0.5*delta_t*ax(p)
           uhy(p) = uy(p) + 0.5*delta_t*ay(p)
           uhz(p) = uz(p) + 0.5*delta_t*az(p)

        if (pelabel(p)>=ne) then
           ! ions
           sum_vxi  = sum_vxi  + uhx(p)
           sum_vyi  = sum_vyi  + uhy(p)
           sum_vzi  = sum_vzi  + uhz(p)
           sum_v2i = sum_v2i + 0.5*(uhx(p)**2+uhy(p)**2+uhz(p)**2)
        endif
     end do
     sum_2vi = sum_vxi**2 + sum_vyi**2 + sum_vzi**2

     ! Find global KE sums
!     call MPI_ALLREDUCE(sum_v2i, global_v2i, one, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     Ti_uncor = 511*mass_eqm*2./3.*sum_v2i/npp  ! This should equal 3/2 kT for 3v Maxwellian
     Ti0 = Ti_keV  ! normalised electron temp


     chii = sqrt(abs(Ti0/Ti_uncor)) 
     chii = min(1.2,max(chii,0.8))  ! Set bounds of +- 50%


     !  3)  Complete full step

     do p=1,npp

        if (pelabel(p)>=ne) then
       ! make ions lighter for eqm phase
           ux(p) = (2*chii-1.)*ux(p) + chii*delta_t*mass_i/mass_eqm*ax(p)
           uy(p) = (2*chii-1.)*uy(p) + chii*delta_t*mass_i/mass_eqm*ay(p)
           uz(p) = (2*chii-1.)*uz(p) + chii*delta_t*mass_i/mass_eqm*az(p)
        endif
     end do
     delta_Ti = 2*Ti0*(1.0/chii**2-1.0)       !  heating
     if (me==0)	write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii,' heating:',delta_Ti



  else

     ! unconstrained motion by default

     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * ax(p)
	uy(p) = uy(p) + delta_t * ay(p)
	uz(p) = uz(p) + delta_t * az(p)
     end do

  endif

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: velocities S<',VTIERR,ICLASSH
!
end subroutine velocities
