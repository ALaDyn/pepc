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

      !> possible values for integrator_scheme
      integer, public, parameter :: INTEGRATOR_SCHEME_NVE = 1 !      1 = NVE - total energy conserved
      integer, public, parameter :: INTEGRATOR_SCHEME_NVT = 2 !      2 = NVT - global Te, Ti conserved
      integer, public, parameter :: INTEGRATOR_SCHEME_NVT_IONS_FROZEN       = 3 ! 3 = global NVT electron Te conserved; ions frozen
      integer, public, parameter :: INTEGRATOR_SCHEME_LOCAL_NVT_IONS_FROZEN = 4 ! 4 = local NVT: each PE keeps Te clamped; ions frozen
      integer, public, parameter :: INTEGRATOR_SCHEME_LOCAL_NVT_IONS_ONLY   = 5 ! 5 = local NVT, ions only; electrons added at end of run
      integer, public, parameter :: INTEGRATOR_SCHEME_FULL_EM  = 6 ! full EM pusher (all E, B components)
      integer, public, parameter :: INTEGRATOR_SCHEME_NONREL   = 7 ! nonrelativistic push
      integer, public, parameter :: INTEGRATOR_SCHEME_NVE_IONS_FROZEN = 8 ! NVE - total electron energy conserved, ions frozen

      integer, public :: integrator_scheme = INTEGRATOR_SCHEME_NVE
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
      public reorder_particles

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>
	!> reorders particles for more efficient force computation (e.g. if only forces
        !> for some particles are needed)
	!>
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine reorder_particles(np_local, particles, num_force_particles)
          use module_pepc_types
          implicit none
          integer, intent(in) :: np_local
          type(t_particle), intent(inout) :: particles(1:np_local)
          integer, intent(out) :: num_force_particles

          type(t_particle) :: temp(1:np_local)
          integer :: num_electrons, num_ions, i

          select case (integrator_scheme)
            case(INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
                temp(1:np_local) = particles(1:np_local)
                num_electrons = 0
                num_ions      = 0

                do i=1,np_local
                   if (temp(i)%data%q < 0.) then
                     num_electrons = num_electrons + 1
                     particles(num_electrons) = temp(i)
                   else
                     num_ions = num_ions + 1
                     particles(np_local-num_ions+1) = temp(i)
                   end if
                end do

                num_force_particles = num_electrons

            case default
              ! nothing to do
              num_force_particles = np_local
          end select

        end subroutine


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

		  case(INTEGRATOR_SCHEME_NVE,INTEGRATOR_SCHEME_NVT,INTEGRATOR_SCHEME_NVT_IONS_FROZEN,INTEGRATOR_SCHEME_LOCAL_NVT_IONS_FROZEN,INTEGRATOR_SCHEME_LOCAL_NVT_IONS_ONLY,INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
		     call velocities(p_start,p_finish,scheme)  ! pure ES, NVT ensembles
		     call push_x(p_start,p_finish,dt)  ! update positions

		  case(INTEGRATOR_SCHEME_FULL_EM)
!		     call push_full3v(p_start,p_finish,dt)  ! full EM pusher (all E, B components)
!		     call push_x(p_start,p_finish,dt)  ! update positions

		  case(INTEGRATOR_SCHEME_NONREL)
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
		  integer :: p,i
		  real*8 :: gam

		  !  relativistic particle push in space
		  do p=ips,ipf
		       gam  = sqrt(1.0 + (dot_product(particles(p)%data%v,particles(p)%data%v))/unit_c2)

		       do i=1,idim
                  particles(p)%x(i) = particles(p)%x(i) + particles(p)%data%v(i)/gam*delt
               end do
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
		  implicit none
		  include 'mpif.h'

		  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.
		  integer, intent(in) :: scheme

		  integer p, i, ne_loc, ierr
		  real*8, dimension(np_local) :: uhx, uhy, uhz, accx, accy, accz
		  real*8 :: sum_vxe, sum_vye, sum_vze, sum_v2e, sum_2ve, acmax
		  real*8 :: delta_u
		  real*8 :: sum_vxi, sum_vyi, sum_vzi, sum_v2i, sum_2vi, mass_eqm
		  real*8 :: global_v2e, gammah, Te_local
		  real*8 :: uprime(1:3), uprime2, sum_ve(1:3), sum_vi(1:3)

		  real*8 :: dimfac(3)

		!  Available ensemble modes
		!      1 = NVE - total energy conserved
		!      2 = NVT - global Te, Ti conserved
		!      3 = global NVT electron Te conserved; ions frozen
		!      4 = local NVT: each PE keeps Te clamped; ions frozen
		!      5 = local NVT, ions only; electrons added at end of run

		! Accelerations
		  acmax=0.
		  do i=p_start,p_finish
		     accx(i) = particles(i)%data%q*particles(i)%results%e(1)/particles(i)%data%m
		     accy(i) = particles(i)%data%q*particles(i)%results%e(2)/particles(i)%data%m
		     accz(i) = particles(i)%data%q*particles(i)%results%e(3)/particles(i)%data%m
		     acmax = max(abs(accx(i)),abs(accy(i)),abs(accz(i)),acmax)
		  end do

		  dimfac         = 0._8
		  dimfac(1:idim) = 1._8

		  delta_u = acmax*dt

		 pusher: select case(scheme)

		 case(INTEGRATOR_SCHEME_NVT)
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
		           uprime(1) = particles(p)%data%v(1) + 0.5*dt*accx(p)
		           uprime(2) = particles(p)%data%v(2) + 0.5*dt*accy(p)
		           uprime(3) = particles(p)%data%v(3) + 0.5*dt*accz(p)
		           uprime2   = dot_product(uprime,uprime)
		           gammah = sqrt(1.0 + uprime2/unit_c2)

		        if (particles(p)%label<=ne) then
		           ! electrons
		           sum_v2e = sum_v2e + uprime2/gammah**2.
		           sum_ve  = sum_ve  + uprime/gammah
		        else if (particles(p)%label<=ne+ni) then
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
		        if (particles(p)%label<=ne) then
		           particles(p)%data%v(1) = (2.*chie-1.)*particles(p)%data%v(1) + chie*dt*accx(p)
		           particles(p)%data%v(2) = (2.*chie-1.)*particles(p)%data%v(2) + chie*dt*accy(p)
		           if (idim==3) particles(p)%data%v(3) = (2.*chie-1.)*particles(p)%data%v(3) + chie*dt*accz(p)

		        elseif (particles(p)%label<=ne+ ni) then
		           particles(p)%data%v(1) = (2.*chii-1.)*particles(p)%data%v(1) + chii*dt*accx(p)
		           particles(p)%data%v(2) = (2.*chii-1.)*particles(p)%data%v(2) + chii*dt*accy(p)
		           if (idim==3) particles(p)%data%v(3) = (2.*chii-1.)*particles(p)%data%v(3) + chii*dt*accz(p)
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

		  case(INTEGRATOR_SCHEME_NVT_IONS_FROZEN)

		! electrons clamped, ions frozen

		     sum_vxe=0.0  ! partial sums
		     sum_vye=0.0
		     sum_vze=0.0
		     sum_v2e=0.0
		     ne_loc = 0  ! # local electrons

		     do p=p_start,p_finish

		        if (particles(p)%label<=ne) then
		           ! electrons
		           ne_loc = ne_loc + 1
		           uhx(p) = particles(p)%data%v(1) + 0.5*dt*accx(p)
		           uhy(p) = particles(p)%data%v(2) + 0.5*dt*accy(p)
		           uhz(p) = particles(p)%data%v(3) + 0.5*dt*accz(p)
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
		        if (particles(p)%label<=ne) then
		           particles(p)%data%v(1) = (2*chie-1.)*particles(p)%data%v(1) + chie*dt*accx(p)
		           particles(p)%data%v(2) = (2*chie-1.)*particles(p)%data%v(2) + chie*dt*accy(p)
		           if (idim==3) particles(p)%data%v(3) = (2*chie-1.)*particles(p)%data%v(3) + chie*dt*accz(p)
		        endif
		     end do

		     delta_Ti=0.
		     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


		  case(INTEGRATOR_SCHEME_LOCAL_NVT_IONS_FROZEN)

		! electrons clamped locally, ions frozen
		! - require T=T_e on each PE to avoid local drifts

		     sum_vxe=0.0  ! partial sums
		     sum_vye=0.0
		     sum_vze=0.0
		     sum_v2e=0.0
		     ne_loc = 0  ! # local electrons

		     do p=p_start,p_finish

		        if (particles(p)%label<=ne) then
		           ! electrons
		           ne_loc = ne_loc+1
		           uhx(p) = particles(p)%data%v(1) + 0.5*dt*accx(p)
		           uhy(p) = particles(p)%data%v(2) + 0.5*dt*accy(p)
		           uhz(p) = particles(p)%data%v(3) + 0.5*dt*accz(p)
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
		        if (particles(p)%label<=ne) then
		           particles(p)%data%v(1) = (2*chie-1.)*particles(p)%data%v(1) + chie*dt*accx(p)
		           particles(p)%data%v(2) = (2*chie-1.)*particles(p)%data%v(2) + chie*dt*accy(p)
		           if (idim==3) particles(p)%data%v(3) = (2*chie-1.)*particles(p)%data%v(3) + chie*dt*accz(p)
		        endif
		     end do

		     delta_Ti=0.
		     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


		  case(INTEGRATOR_SCHEME_LOCAL_NVT_IONS_ONLY)

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
		           uhx(p) = particles(p)%data%v(1) + 0.5*dt*accx(p)
		           uhy(p) = particles(p)%data%v(2) + 0.5*dt*accy(p)
		           uhz(p) = particles(p)%data%v(3) + 0.5*dt*accz(p)

		        if (particles(p)%label>=ne) then
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

		        if (particles(p)%label>=ne) then
		       ! make ions lighter for eqm phase
		           particles(p)%data%v(1) = (2*chii-1.)*particles(p)%data%v(1) + chii*dt*accx(p)
		           particles(p)%data%v(2) = (2*chii-1.)*particles(p)%data%v(2) + chii*dt*accy(p)
		           if (idim==3) particles(p)%data%v(3) = (2*chii-1.)*particles(p)%data%v(3) + chii*dt*accz(p)
		!           particles(p)%data%v(1) = particles(p)%data%v(1)*sqrt(Ti0/Ti_uncor)
		!           particles(p)%data%v(2) = particles(p)%data%v(2)*sqrt(Ti0/Ti_uncor)
		!           particles(p)%data%v(3) = particles(p)%data%v(3)*sqrt(Ti0/Ti_uncor)

		        endif
		     end do
		     delta_Ti = 2*Ti0*(1.0/chii**2-1.0)       !  heating
		     if (my_rank==0) then
		        write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii,' heating:',delta_Ti
		        write (*,*) 'Max delta-ux: ',delta_u
		     endif

          case (INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
           ! unconstrained motion for negatively charged particles, frozen positively charged particles
           do p = p_start, p_finish
             if (particles(p)%data%q<0.) then
               particles(p)%data%v(1) = particles(p)%data%v(1) + dt * accx(p) * dimfac(1)
               particles(p)%data%v(2) = particles(p)%data%v(2) + dt * accy(p) * dimfac(2)
               particles(p)%data%v(3) = particles(p)%data%v(3) + dt * accz(p) * dimfac(3)
             else
               particles(p)%data%v(1) = 0.
               particles(p)%data%v(2) = 0.
               particles(p)%data%v(3) = 0.
             endif
           end do

		  case default
		     ! unconstrained motion by default (scheme=INTEGRATOR_SCHEME_NVE,INTEGRATOR_SCHEME_NONREL)
		   if (idim==3) then
		     do p = p_start, p_finish
		       particles(p)%data%v(1) = particles(p)%data%v(1) + dt * accx(p)
		       particles(p)%data%v(2) = particles(p)%data%v(2) + dt * accy(p)
		       particles(p)%data%v(3) = particles(p)%data%v(3) + dt * accz(p)
		     end do
		   else if (idim==2) then
		     do p = p_start, p_finish
		       particles(p)%data%v(1) = particles(p)%data%v(1) + dt * accx(p)
		       particles(p)%data%v(2) = particles(p)%data%v(2) + dt * accy(p)
		     end do
		   else
		     do p = p_start, p_finish
		       particles(p)%data%v(1) = particles(p)%data%v(1) + dt * accx(p)
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


		  do p = p_start, p_finish
		     beta=particles(p)%data%q/particles(p)%data%m*dts*0.5  ! charge/mass constant

		     xd = particles(p)%x(1)-focus(1)
		     yd = particles(p)%x(2)-focus(2)
		     zd = particles(p)%x(3)-focus(3)

		     ! evaluate external fields at particle positions

		     if (modulo(beam_config_in,10).eq.4) then
		! pond. standing wave on step-profile
		        call empond(t_laser,t_pulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipon)

		     else if (modulo(beam_config_in,10).eq.6) then
		! plane wave with Gaussian spot
		        call emplane(t_laser,t_pulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipon)
		     endif

		     !  Sum internal and external fields
		     exi = particles(p)%results%e(1)
		     eyi = particles(p)%results%e(2)
		     ezi = particles(p)%results%e(3)+ez_em
		     bxi = bx_em
		     byi = by_em
		     bzi = 0.

		     ! transverse momentum from pz=az including thermal motion
		     !  uzi = -particles(p)%data%q/particles(p)%data%m*az_em + particles(p)%data%v(3)


		     !   first half-accn
		     uxm = particles(p)%data%v(1) + beta*exi
		     uym = particles(p)%data%v(2) + beta*eyi
		     uzm = particles(p)%data%v(3) + beta*ezi
		     !     uzm = particles(p)%data%v(3) - particles(p)%data%q/particles(p)%data%m*az_em
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
		     particles(p)%data%v(1) = uxp + beta*exi
		     particles(p)%data%v(2) = uyp + beta*eyi
		     particles(p)%data%v(3) = uzp + beta*ezi

		  end do


		end subroutine push_em


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

		     particles(p)%x(1)=particles(p)%x(1)+particles(p)%data%v(1)*delt
		     particles(p)%x(2)=particles(p)%x(2)+particles(p)%data%v(2)*delt
		     particles(p)%x(3)=particles(p)%x(3)+particles(p)%data%v(3)*delt

		  end do

		end subroutine push_nonrel



end module module_pusher
