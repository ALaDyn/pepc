! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2014 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!  ===================================================================
!
!                              VELOCITIES
!
!   $Revision: 1068 $
!
!   Calculate velocities from accelerations
! 
!   apply thermodynamic constraint according to ensemble
!
!
!
!  ===================================================================


subroutine velocities(p_start,p_finish,delta_t)


  use physvars
  use treevars
  use utils
  implicit none
  include 'mpif.h'

  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer p, i, ne_loc, ierr
  real*8, dimension(nppm) :: uhx, uhy, uhz, accx, accy, accz
  real*8 :: sum_vxe, sum_vye, sum_vze, sum_v2e, sum_2ve, acmax
  real :: Te0, Te_uncor, Ti0, Ti_uncor, chie, chii, delta_u
  real*8 :: sum_vxi, sum_vyi, sum_vzi, sum_v2i, sum_2vi, mass_eqm
  real*8 :: global_v2e, global_v2i, gammah, delta_Te, delta_Ti, Te_local

!  Available ensemble modes
!      1 = NVE - total energy conserved
!      2 = NVT - global Te, Ti conserved
!      3 = global NVT electron Te conserved; ions frozen
!      4 = local NVT: each PE keeps Te clamped; ions frozen
!      5 = local NVT, ions only; electrons added at end of run

! Accelerations
  acmax=0.
  do i=1,npp
     accx(i) = q(i)*ex(i)/m(i)
     accy(i) = q(i)*ey(i)/m(i)
     accz(i) = q(i)*ez(i)/m(i)
     acmax = max(abs(accx(i)),abs(accy(i)),abs(accz(i)),acmax)
  end do

  delta_u = acmax*delta_t

 pusher: select case(scheme)

 case(2)
     ! Conserve kinetic energies of electrons and ions (initial Te const)
     ! adapted from
     !  Allen and Tildesley p230, Brown & Clark, Mol. Phys. 51, 1243 (1984)

     !  Definitions:
     !
     !   unconstrained velocities    uh(x,y,z) = v*(x,y,z)

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
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)
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
     call MPI_ALLREDUCE(sum_v2e, global_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(sum_v2i, global_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


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
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*accx(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*accz(p)

        elseif (pelabel(p)<=ne+ ni) then
           ux(p) = (2*chii-1.)*ux(p) + chii*delta_t*accx(p)
           uy(p) = (2*chii-1.)*uy(p) + chii*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chii-1.)*uz(p) + chii*delta_t*accz(p)
        endif
     end do




     delta_Ti=0.


     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating
     !     heate = heate + delta_Te


  case(3)

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
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2

     ! Find global KE sums
     call MPI_ALLREDUCE(sum_v2e, global_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*global_v2e/ne  ! This should equal 3/2 kT for 3v Maxwellian
     Te_local = 511*2./3.*sum_v2e/ne_loc
     Te0 = Te_keV  ! normalised electron temp
     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
     chie = min(1.25,max(chie,0.75))  ! Set bounds of +- 50%

     if (me==0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

     !  3)  Complete full step

     do p=1,npp
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*accx(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*accz(p)
        endif
     end do

     delta_Ti=0.
     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


  case(4)

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
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)
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

     if (me.eq.0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

     !  3)  Complete full step

     do p=1,npp
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*accx(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*accz(p)
        endif
     end do

     delta_Ti=0.
     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


  case(5)

     ! Conserve kinetic energy of ions only (initial Ti const)
     mass_eqm = 20.  ! artificial ion mass for eqm stage
     sum_vxi=0.0  ! partial sums (ions)
     sum_vyi=0.0
     sum_vzi=0.0
     sum_v2i=0.0
!  Scale down accelerations if too big
     do p=1,npp
	accx(p)=accx(p)/max(1.d0,acmax)
	accy(p)=accy(p)/max(1.d0,acmax)
	accz(p)=accz(p)/max(1.d0,acmax)
     end do

     do p=1,npp
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)

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
!     call MPI_ALLREDUCE(sum_v2i, global_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     Ti_uncor = 511*2./3.*sum_v2i/npp  ! This should equal 3/2 kT for 3v Maxwellian
     Ti0 = Ti_keV  ! normalised electron temp


     chii = sqrt(abs(Ti0/Ti_uncor)) 
     chii = min(1.2,max(chii,0.6))  ! Set bounds of +- 50%


     !  3)  Complete full step

     do p=1,npp

        if (pelabel(p)>=ne) then
       ! make ions lighter for eqm phase
           ux(p) = (2*chii-1.)*ux(p) + chii*delta_t*accx(p)
           uy(p) = (2*chii-1.)*uy(p) + chii*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chii-1.)*uz(p) + chii*delta_t*accz(p)
!           ux(p) = ux(p)*sqrt(Ti0/Ti_uncor)
!           uy(p) = uy(p)*sqrt(Ti0/Ti_uncor)
!           uz(p) = uz(p)*sqrt(Ti0/Ti_uncor)

        endif
     end do
     delta_Ti = 2*Ti0*(1.0/chii**2-1.0)       !  heating
     if (me==0)	then
	write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii,' heating:',delta_Ti,' bond',bond_const
        write (*,*) 'Max delta-ux: ',delta_u
     endif


  case default
  if (me==0) then
	write(6,*) 'PEPC | velocities scheme',scheme
	write(24,*) 'PEPC | velocities scheme',scheme
  endif
     ! unconstrained motion by default (scheme=1,7)
   if (idim==3) then
     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * accx(p)
	uy(p) = uy(p) + delta_t * accy(p)
	uz(p) = uz(p) + delta_t * accz(p)
     end do
   else if (idim==2) then
     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * accx(p)
	uy(p) = uy(p) + delta_t * accy(p)
     end do
   else
     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * accx(p)
     end do
   endif
 
  end select pusher

end subroutine velocities








