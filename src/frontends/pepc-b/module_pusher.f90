! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2017 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything concerning the actual particle movement (integrator, pusher, boundaries, constraints, etc)
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_pusher
      use module_units
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, public :: integrator_scheme = 1
      real*8, public :: Te0 = 0., Te_uncor = 0., chie = 0., delta_Te = 0.
      real*8, public :: Ti0 = 0., Ti_uncor = 0., chii = 0., delta_Ti = 0.
      logical, public :: enable_drift_elimination = .false. !< if .true., the global particle drift is included during velocity rescaling for constant temperature regime (i.e. it is eliminated)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public integrator
      public push_em
      public push_nonrel
      public push_full3v

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Velocity and position update - integration of the equation of motion
		!> Boundary conditions are also applied here
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine integrator(p_start,p_finish,scheme)

		  use physvars
		  implicit none
          integer, intent(in) :: p_start,p_finish,scheme

		 if (my_rank==0) then
            write( 6,'(/"-- PUSHER --"/a20,i8)') 'scheme = ',scheme
            write(24,'(/"-- PUSHER --"/a20,i8)') 'scheme = ',scheme
		 endif

		  pusher: select case(scheme)

		  case(1,2,3,4,5)
		     call velocities(p_start,p_finish,scheme)  ! pure ES, NVT ensembles
		     call push_x(p_start,p_finish,dt)  ! update positions

		  case(6)
		     call push_full3v(p_start,p_finish,dt)  ! full EM pusher (all E, B components)
		     call push_x(p_start,p_finish,dt)  ! update positions

		  case(7)
		     call velocities(p_start,p_finish,scheme)  ! nonrelativistic push
		     call push_nonrel(p_start,p_finish,dt)
		  case default
		     ! do nothing!

		  end select pusher


		end subroutine integrator



		!  ===============================================================
		!
		!                           PMOVE
		!
		!   Update particle positions - used with leap-frog scheme
		!
		!  ===============================================================

		subroutine push_x(ips,ipf,delt)
		  use physvars
		  use module_units
		  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
		  real*8, intent(in) :: delt
		  integer :: p
		  real*8 :: gam

		  !  relativistic particle push in space
		  do p=ips,ipf
		       gam  = sqrt(1.0 + (ux(p)**2 + uy(p)**2 + uz(p)**2)/unit_c2)
		                      x(p) = x(p) + ux(p)/gam*delt
		       if (idim > 1)  y(p) = y(p) + uy(p)/gam*delt
		       if (idim == 3) z(p) = z(p) + uz(p)/gam*delt
		  end do

		end subroutine push_x


		!  ===================================================================
		!
		!                              VELOCITIES
		!
		!   $Revision: 1264 $
		!
		!   Calculate velocities from accelerations
		!
		!   apply thermodynamic constraint according to ensemble
		!
		!
		!
		!  ===================================================================


		subroutine velocities(p_start,p_finish,scheme)


		  use physvars
		  use utils
		  implicit none
		  include 'mpif.h'

		  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.
		  integer, intent(in) :: scheme

		  integer p, i, ne_loc, ierr
		  real*8, dimension(nppm) :: uhx, uhy, uhz, accx, accy, accz
		  real*8 :: sum_vxe, sum_vye, sum_vze, sum_v2e, sum_2ve, acmax
		  real*8 :: delta_u
		  real*8 :: sum_vxi, sum_vyi, sum_vzi, sum_v2i, sum_2vi, mass_eqm
		  real*8 :: global_v2e, gammah, Te_local
		  real*8 :: uprime(1:3), uprime2, sum_ve(1:3), sum_vi(1:3)

		!  Available ensemble modes
		!      1 = NVE - total energy conserved
		!      2 = NVT - global Te, Ti conserved
		!      3 = global NVT electron Te conserved; ions frozen
		!      4 = local NVT: each PE keeps Te clamped; ions frozen
		!      5 = local NVT, ions only; electrons added at end of run

		! Accelerations
		  acmax=0.
		  do i=p_start,p_finish
		     accx(i) = q(i)*ex(i)/m(i)
		     accy(i) = q(i)*ey(i)/m(i)
		     accz(i) = q(i)*ez(i)/m(i)
		     acmax = max(abs(accx(i)),abs(accy(i)),abs(accz(i)),acmax)
		  end do

		  delta_u = acmax*dt

		 pusher: select case(scheme)

		 case(2)
		     ! Conserve kinetic energies of electrons and ions (initial Te const)
		     ! adapted from
		     !  Allen and Tildesley p230, Brown & Clark, Mol. Phys. 51, 1243 (1984)

		     !  Definitions:
		     !
		     !   unconstrained velocities    uh(x,y,z) = v'(x,y,z)

		     !  1)  Unconstrained half-step for electrons and ions

             sum_v2e=0.0
             sum_v2i=0.0
             sum_ve =0.0
             sum_vi =0.0

		     do p=p_start,p_finish
		           uprime(1) = ux(p) + 0.5*dt*accx(p)
		           uprime(2) = uy(p) + 0.5*dt*accy(p)
		           uprime(3) = uz(p) + 0.5*dt*accz(p)
		           uprime2   = dot_product(uprime,uprime)
		           gammah = sqrt(1.0 + uprime2/unit_c2)

		        if (pelabel(p)<=ne) then
		           ! electrons
		           sum_v2e = sum_v2e + uprime2/gammah**2.
		           sum_ve  = sum_ve  + uprime/gammah
		        else if (pelabel(p)<=ne+ni) then
		           ! ions
                   sum_v2i = sum_v2i + uprime2/gammah**2.
                   sum_vi  = sum_vi  + uprime/gammah
		        endif
		     end do

		     ! Find global KE sums
             call MPI_ALLREDUCE(MPI_IN_PLACE,sum_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE,sum_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE, sum_ve, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
             call MPI_ALLREDUCE(MPI_IN_PLACE, sum_vi, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

             sum_v2e = sum_v2e / ne
             sum_ve  = sum_ve  / ne
             sum_v2i = sum_v2i / ni
             sum_vi  = sum_vi  / ni

             ! 2) unconstrained temperatures (particle velocities corrected by drift)
             if (.not. enable_drift_elimination) then
		       Te_uncor = mass_e/(3.*unit_kB)*(sum_v2e - dot_product(sum_ve,sum_ve))  ! This should equal kT for the unconstrained half-step velocities (see Allen & Tildesley, (2.49) and p231)
		       Ti_uncor = mass_i/(3.*unit_kB)*(sum_v2i - dot_product(sum_vi,sum_vi))
		     else
               Te_uncor = mass_e/(3.*unit_kB)*sum_v2e  ! This should equal kT for the unconstrained half-step velocities (see Allen & Tildesley, (2.49) and p231)
               Ti_uncor = mass_i/(3.*unit_kB)*sum_v2i
		     endif
		     Te0 = Te  ! normalised (desired) electron temp
		     Ti0 = Ti  ! normalised (desired) ion temp
		     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
		     chii = sqrt(abs(Ti0/Ti_uncor))

             ! reduce artificial friction factor by 1./2. for both species
             chie = 1./( (1./2.)*(1./chie - 1) + 1.)
             chii = 1./( (1./2.)*(1./chii - 1) + 1.)

		     !  3)  Complete full step
		     do p=p_start,p_finish
		        if (pelabel(p)<=ne) then
		           ux(p) = (2.*chie-1.)*ux(p) + chie*dt*accx(p)
		           uy(p) = (2.*chie-1.)*uy(p) + chie*dt*accy(p)
		           if (idim==3) uz(p) = (2.*chie-1.)*uz(p) + chie*dt*accz(p)

		        elseif (pelabel(p)<=ne+ ni) then
		           ux(p) = (2.*chii-1.)*ux(p) + chii*dt*accx(p)
		           uy(p) = (2.*chii-1.)*uy(p) + chii*dt*accy(p)
		           if (idim==3) uz(p) = (2.*chii-1.)*uz(p) + chii*dt*accz(p)
		        endif
		     end do

		     delta_Ti = 2*Ti0*(1.0/chii**2-1.0)
		     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating

             if (my_rank==0) then
               write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie, 'delta_Te', delta_Te
               write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii, 'delta_Ti', delta_Ti

               if (enable_drift_elimination) then
                 write(*,'(a20)') 'Drift elim. ENABLED'
               else
                 write(*,'(a20)') 'Drift elim. DISABLED'
               endif

               write (*,'(a20,3(e12.2,8x))') 'e-drift = ', sum_ve
               write (*,'(a20,3(e12.2,8x))') 'i-drift = ', sum_vi
               write (*,*) ''
             endif

		  case(3)

		! electrons clamped, ions frozen

		     sum_vxe=0.0  ! partial sums
		     sum_vye=0.0
		     sum_vze=0.0
		     sum_v2e=0.0
		     ne_loc = 0  ! # local electrons

		     do p=p_start,p_finish

		        if (pelabel(p)<=ne) then
		           ! electrons
		           ne_loc = ne_loc + 1
		           uhx(p) = ux(p) + 0.5*dt*accx(p)
		           uhy(p) = uy(p) + 0.5*dt*accy(p)
		           uhz(p) = uz(p) + 0.5*dt*accz(p)
		           gammah = sqrt(1.0 +(uhx(p)**2 + uhy(p)**2 + uhz(p)**2)/unit_c2)
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
		     Te0 = Te_eV/1000.  ! normalised electron temp
		     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
		     chie = min(1.25_8,max(chie,0.75_8))  ! Set bounds of +- 50%

		     if (my_rank==0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

		     !  3)  Complete full step

		     do p=p_start,p_finish
		        if (pelabel(p)<=ne) then
		           ux(p) = (2*chie-1.)*ux(p) + chie*dt*accx(p)
		           uy(p) = (2*chie-1.)*uy(p) + chie*dt*accy(p)
		           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*dt*accz(p)
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

		     do p=p_start,p_finish

		        if (pelabel(p)<=ne) then
		           ! electrons
		           ne_loc = ne_loc+1
		           uhx(p) = ux(p) + 0.5*dt*accx(p)
		           uhy(p) = uy(p) + 0.5*dt*accy(p)
		           uhz(p) = uz(p) + 0.5*dt*accz(p)
		           gammah = sqrt(1.0 +(uhx(p)**2 + uhy(p)**2 + uhz(p)**2)/unit_c2)
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
		     Te0 = Te_eV/1000.  ! normalised electron temp
		     ! exponent of chie should be 1/2 - take 1/3 to soften oscillations
		     chie = (abs(Te0/Te_uncor))**0.5     ! multipliers from Temperature ratio - reset once every cycle

		     chie = min(1.25_8,max(chie,0.75_8))  ! Set bounds of +- 50%

		     if (my_rank.eq.0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

		     !  3)  Complete full step

		     do p=p_start,p_finish
		        if (pelabel(p)<=ne) then
		           ux(p) = (2*chie-1.)*ux(p) + chie*dt*accx(p)
		           uy(p) = (2*chie-1.)*uy(p) + chie*dt*accy(p)
		           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*dt*accz(p)
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
		     do p=p_start,p_finish
		       accx(p)=accx(p)/max(1.d0,acmax)
		       accy(p)=accy(p)/max(1.d0,acmax)
		       accz(p)=accz(p)/max(1.d0,acmax)
		     end do

		     do p=p_start,p_finish
		           uhx(p) = ux(p) + 0.5*dt*accx(p)
		           uhy(p) = uy(p) + 0.5*dt*accy(p)
		           uhz(p) = uz(p) + 0.5*dt*accz(p)

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
		     Ti_uncor = 511*2./3.*sum_v2i/(p_finish-p_start)  ! This should equal 3/2 kT for 3v Maxwellian
		     Ti0 = Ti_eV/1000.  ! normalised electron temp


		     chii = sqrt(abs(Ti0/Ti_uncor))
		     chii = min(1.2_8,max(chii,0.6_8))  ! Set bounds of +- 50%


		     !  3)  Complete full step

		     do p=p_start,p_finish

		        if (pelabel(p)>=ne) then
		       ! make ions lighter for eqm phase
		           ux(p) = (2*chii-1.)*ux(p) + chii*dt*accx(p)
		           uy(p) = (2*chii-1.)*uy(p) + chii*dt*accy(p)
		           if (idim==3) uz(p) = (2*chii-1.)*uz(p) + chii*dt*accz(p)
		!           ux(p) = ux(p)*sqrt(Ti0/Ti_uncor)
		!           uy(p) = uy(p)*sqrt(Ti0/Ti_uncor)
		!           uz(p) = uz(p)*sqrt(Ti0/Ti_uncor)

		        endif
		     end do
		     delta_Ti = 2*Ti0*(1.0/chii**2-1.0)       !  heating
		     if (my_rank==0) then
		        write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii,' heating:',delta_Ti
		        write (*,*) 'Max delta-ux: ',delta_u
		     endif


		  case default
		     ! unconstrained motion by default (scheme=1,7)
		   if (idim==3) then
		     do p = p_start, p_finish
		       ux(p) = ux(p) + dt * accx(p)
		       uy(p) = uy(p) + dt * accy(p)
		       uz(p) = uz(p) + dt * accz(p)
		     end do
		   else if (idim==2) then
		     do p = p_start, p_finish
		       ux(p) = ux(p) + dt * accx(p)
		       uy(p) = uy(p) + dt * accy(p)
		     end do
		   else
		     do p = p_start, p_finish
		       ux(p) = ux(p) + dt * accx(p)
		     end do
		   endif

		  end select pusher

		end subroutine velocities




		!  ==============================================
		!
		!     3v particle pusher
		!
		!  TM fields only (s-pol):  (0,0,Ez)  (Bx,By,0)
		!  TE fields to follow (Ex,Ey,0) (0,0,Bz)
		!  ==============================================


		subroutine push_em(p_start,p_finish,dts)
		  use physvars
		  use module_laser
		  implicit none
		  integer, intent(in) :: p_start, p_finish
		  real*8, intent(in) :: dts

		  integer :: p
		  real*8 :: beta, gam1,tt, sy, sz, tz, ty
		  real*8 :: uxd,uyp,uzp,uxp,uxm,uym,uzm
		  real*8 :: exi, eyi, ezi, bxi, byi, bzi
		  real*8 :: phipon, ez_em, bx_em, by_em, az_em
		  real*8 :: xd, yd, zd


		  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

		  do p = p_start, p_finish
		     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

		     xd = x(p)-focus(1)
		     yd = y(p)-focus(2)
		     zd = z(p)-focus(3)

		     ! evaluate external fields at particle positions

		     if (beam_config.eq.4) then
		! pond. standing wave on step-profile
		        call empond(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipon)

		     else if (beam_config.eq.6) then
		! plane wave with Gaussian spot
		        call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipon)
		     endif

		     !  Sum internal and external fields
		     exi = ex(p)
		     eyi = ey(p)
		     ezi = ez(p)+ez_em
		     bxi = bx_em
		     byi = by_em
		     bzi = 0.

		     ! transverse momentum from pz=az including thermal motion
		     !  uzi = -q(p)/m(p)*az_em + uz(p)


		     !   first half-accn
		     uxm = ux(p) + beta*exi
		     uym = uy(p) + beta*eyi
		     uzm = uz(p) + beta*ezi
		     !     uzm = uz(p) - q(p)/m(p)*az_em
		     !   rotation
		     gam1=dsqrt(1.d0 + (uxm**2 + uym**2 + uzm**2)/unit_c2)
		     ty = beta*byi/gam1
		     tz = beta*bzi/gam1
		     tt = 1.0 + ty**2+tz**2
		     sy = 2.0*ty/tt
		     sz = 2.0*tz/tt

		     uxd = uxm + uym*tz - uzm*ty
		     uyp = uym - uxd*sz
		     uzp = uzm + uxd*sy
		     uxp = uxd + uyp*tz - uzp*ty

		     !   second half-accn
		     ux(p) = uxp + beta*exi
		     uy(p) = uyp + beta*eyi
		     uz(p) = uzp + beta*ezi

		  end do


		end subroutine push_em



		!  ==============================================
		!
		!     Full 3v particle pusher (Boris rotation scheme)
		!
		!  TM fields only (s-pol):  (Ex,0,Ez)  (0,By,0)
		!  TE fields  (0,Ey,0) (Bx,0,Bz)
		!  ==============================================


		subroutine push_full3v(p_start,p_finish,dts)
		  use physvars
		  use module_laser
		  implicit none
		  integer, intent(in) :: p_start, p_finish
		  real*8, intent(in) :: dts

		  integer :: p
		  real*8 :: beta, gam1, xd, yd, zd
		  real*8 :: tt, sx, sy, sz, tz, ty, tx
		  real*8 :: uxd, uyd, uzd, uyp, uzp, uxp, uxm, uym, uzm
		  real*8 :: exi, eyi, ezi, bxi, byi, bzi

		  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

		  do p = p_start, p_finish
		     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

		     xd = x(p)-focus(1)
		     yd = y(p)-focus(2)
		     zd = z(p)-focus(3)

		     exi = ex(p)
		     eyi = ey(p)
		     ezi = ez(p)
		     ! magnetic fields are currently not implemented in the code
		     ! TODO: use the missing stuff from sum_bfield.f90, fields_p.f90
		     ! and put it into calc_force_per_interaction()
		     ! then reactivate this stuff here and some lines below
		     bxi = Bx(p)
		     byi = By(p)
		     bzi = Bz(p)

		     !   first half-accn
		     uxm = ux(p) + beta*exi
		     uym = uy(p) + beta*eyi
		     uzm = uz(p) + beta*ezi

		     !   rotation
		     gam1=dsqrt(1.d0 + (uxm**2 + uym**2 + uzm**2)/unit_c2)

		     tx = beta*bxi/gam1
		     ty = beta*byi/gam1
		     tz = beta*bzi/gam1
		     tt = 1.0 + tx**2 + ty**2 + tz**2

		     sx = 2.0*tx/tt
		     sy = 2.0*ty/tt
		     sz = 2.0*tz/tt

		     uxd = uxm + uym*tz - uzm*ty
		     uyd = uym + uzm*tx - uxm*tz
		     uzd = uzm + uxm*ty - uym*tx

		     uxp = uxm + uyd*sz - uzd*sy
		     uyp = uym + uzd*sx - uxd*sz
		     uzp = uzm + uxd*sy - uyd*sx

		     !   second half-accn
		     ux(p) = uxp + beta*exi
		     uy(p) = uyp + beta*eyi
		     uz(p) = uzp + beta*ezi

		  end do


		end subroutine push_full3v



		!  ===============================================================
		!
		!                           PUSH_NONREL
		!
		!   Nonrelativistic particle position update - used with leap-frog scheme
		!
		!  ===============================================================

		subroutine push_nonrel(ips,ipf,delt)

		  use physvars
		  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
		  real*8, intent(in) :: delt
		  integer :: p


		  do p=ips,ipf

		     x(p)=x(p)+ux(p)*delt
		     y(p)=y(p)+uy(p)*delt
		     z(p)=z(p)+uz(p)*delt

		  end do

		end subroutine push_nonrel



end module module_pusher
