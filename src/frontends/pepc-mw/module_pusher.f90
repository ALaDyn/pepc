! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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
    public reorder_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        use module_io
        implicit none
        integer, intent(in) :: p_start,p_finish,scheme

        if (my_rank==0) then
            write(file_stdout  ,'(/"-- PUSHER --"/a20,i8)') 'scheme = ',scheme
            write(file_pepc_out,'(/"-- PUSHER --"/a20,i8)') 'scheme = ',scheme
        endif

        pusher: select case(scheme)

            case(INTEGRATOR_SCHEME_NVE,INTEGRATOR_SCHEME_NVT,INTEGRATOR_SCHEME_NVT_IONS_FROZEN,INTEGRATOR_SCHEME_LOCAL_NVT_IONS_FROZEN,INTEGRATOR_SCHEME_LOCAL_NVT_IONS_ONLY,INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
                call velocities(p_start,p_finish,scheme)  ! pure ES, NVT ensembles
                call push_rel(p_start,p_finish,dt)  ! update positions

            case(INTEGRATOR_SCHEME_FULL_EM)
                call velocities_boris(p_start,p_finish,dt)  ! full EM pusher (all E, B components)
                call push_rel(p_start,p_finish,dt)  ! update positions

            case(INTEGRATOR_SCHEME_NONREL)
                call velocities(p_start,p_finish,scheme)
                call push_nonrel(p_start,p_finish,dt)       ! nonrelativistic push
            case default
           ! do nothing!

        end select pusher


    end subroutine integrator


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate velocities from accelerations due to electric field
    !> respect certain thermodynamic constraints due to using specific ensembles
    !> (does not (!) apply to bors pusher)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine velocities(p_start,p_finish,scheme)

        use physvars
        implicit none
        include 'mpif.h'

        integer, intent(in) :: p_start,p_finish  ! min, max particle nos.
        integer, intent(in) :: scheme

        integer p, i, ne_loc, ierr
        real*8, dimension(3, np_local) :: uh, acc
        real*8 :: sum_vxe, sum_vye, sum_vze, acmax
        real*8 :: delta_u
        real*8 :: sum_vxi, sum_vyi, sum_vzi, mass_eqm
        real*8 :: global_v2e, gammah, Te_local
        real*8 :: uprime(1:3), uprime2, sum_ve(1:3), sum_vi(1:3)

        integer, parameter :: V2E = 1
        integer, parameter :: VEX = 2
        integer, parameter :: VEY = 3
        integer, parameter :: VEZ = 4
        integer, parameter :: V2I = 5
        integer, parameter :: VIX = 6
        integer, parameter :: VIY = 7
        integer, parameter :: VIZ = 8

        real*8 :: sums(1:8)

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
            acc(1:3,i) = particles(i)%data%q * particles(i)%results%e(1:3) / particles(i)%data%m
        end do

        acmax = maxval(abs(acc(:,:)))

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

                sums = 0.0

                !  1)  Unconstrained half-step for electrons and ions
                do p=p_start,p_finish
                    uprime(1:3) = particles(p)%data%v(1:3) + 0.5*dt*acc(1:3, p)
                    uprime2     = dot_product(uprime,uprime)
                    gammah      = sqrt(1.0 + uprime2/unit_c2)

                    if (particles(p)%data%q < 0.) then
                        ! electrons
                        sums(V2E)     = sums(V2E)     + uprime2 / gammah**2.
                        sums(VEX:VEZ) = sums(VEX:VEZ) + uprime  / gammah
                    else if (particles(p)%data%q > 0.) then
                        ! ions
                        sums(V2I)     = sums(V2I)     + uprime2 / gammah**2.
                        sums(VIX:VIZ) = sums(VIX:VIZ) + uprime  / gammah
                    endif
                end do

                ! Find global KE sums
                call MPI_ALLREDUCE(MPI_IN_PLACE, sums, size(sums), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

                sums(V2E:VEZ) = sums(V2E:VEZ) / ne
                sums(V2I:VIZ) = sums(V2I:VIZ) / ni

                ! 2) unconstrained temperatures (particle velocities corrected by drift)
                if (.not. enable_drift_elimination) then
                    Te_uncor = mass_e/(3.*unit_kB)*(sums(V2E) - dot_product(sums(VEX:VEZ),sums(VEX:VEZ)))  ! This should equal kT for the unconstrained half-step velocities (see Allen & Tildesley, (2.49) and p231)
                    Ti_uncor = mass_i/(3.*unit_kB)*(sums(V2I) - dot_product(sums(VIX:VIZ),sums(VIX:VIZ)))
                else
                    Te_uncor = mass_e/(3.*unit_kB)*sums(V2E)  ! This should equal kT for the unconstrained half-step velocities (see Allen & Tildesley, (2.49) and p231)
                    Ti_uncor = mass_i/(3.*unit_kB)*sums(V2I)
                endif

                Te0 = Te  ! normalised (desired) electron temp
                Ti0 = Ti  ! normalised (desired) ion temp

                chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
                chii = sqrt(abs(Ti0/Ti_uncor))

                ! reduce artificial friction factor by 1./2. for both species - this does not comply with Allen & Tildesley but works fine :-S
                chie = 1./( (1./2.)*(1./chie - 1) + 1.)
                chii = 1./( (1./2.)*(1./chii - 1) + 1.)

                !  3)  Complete full step
                do p=p_start,p_finish
                    if (particles(p)%data%q < 0.) then
                        particles(p)%data%v(1:3) = (2.*chie-1.)*particles(p)%data%v(1:3) + chie*dt*acc(1:3,p)
                    elseif (particles(p)%data%q > 0.) then
                        particles(p)%data%v(1:3) = (2.*chii-1.)*particles(p)%data%v(1:3) + chii*dt*acc(1:3,p)
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

                    write (*,'(a20,3(e12.2,8x))') 'e-drift = ', sums(VEX:VEZ)
                    write (*,'(a20,3(e12.2,8x))') 'i-drift = ', sums(VIX:VIZ)
                    write (*,*) ''
                endif

!            case(INTEGRATOR_SCHEME_NVT_IONS_FROZEN)
!
!                ! electrons clamped, ions frozen
!
!                sum_vxe=0.0  ! partial sums
!                sum_vye=0.0
!                sum_vze=0.0
!                sum_v2e=0.0
!                ne_loc = 0  ! # local electrons
!
!                do p=p_start,p_finish
!
!                    if (particles(p)%label<=ne) then
!                        ! electrons
!                        ne_loc = ne_loc + 1
!                        uhx(p) = particles(p)%data%v(1) + 0.5*dt*accx(p)
!                        uhy(p) = particles(p)%data%v(2) + 0.5*dt*accy(p)
!                        uhz(p) = particles(p)%data%v(3) + 0.5*dt*accz(p)
!                        gammah = sqrt(1.0 +(uhx(p)**2 + uhy(p)**2 + uhz(p)**2)/unit_c2)
!                        sum_vxe  = sum_vxe  + uhx(p)/gammah
!                        sum_vye  = sum_vye  + uhy(p)/gammah
!                        sum_vze  = sum_vze  + uhz(p)/gammah
!                        sum_v2e = sum_v2e + gammah-1.
!
!                    endif
!                end do
!                sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2
!
!                ! Find global KE sums
!                call MPI_ALLREDUCE(sum_v2e, global_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
!
!
!                ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
!                Te_uncor = 511*2./3.*global_v2e/ne  ! This should equal 3/2 kT for 3v Maxwellian
!                Te_local = 511*2./3.*sum_v2e/ne_loc
!                Te0 = Te_eV/1000.  ! normalised electron temp
!                chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
!                chie = min(1.25_8,max(chie,0.75_8))  ! Set bounds of +- 50%
!
!                if (my_rank==0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie
!
!                !  3)  Complete full step
!
!                do p=p_start,p_finish
!                    if (particles(p)%label<=ne) then
!                        particles(p)%data%v(1) = (2*chie-1.)*particles(p)%data%v(1) + chie*dt*accx(p)
!                        particles(p)%data%v(2) = (2*chie-1.)*particles(p)%data%v(2) + chie*dt*accy(p)
!                        if (idim==3) particles(p)%data%v(3) = (2*chie-1.)*particles(p)%data%v(3) + chie*dt*accz(p)
!                    endif
!                end do
!
!                delta_Ti=0.
!                delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating
!
!
!            case(INTEGRATOR_SCHEME_LOCAL_NVT_IONS_FROZEN)
!
!                ! electrons clamped locally, ions frozen
!                ! - require T=T_e on each PE to avoid local drifts
!
!                sum_vxe=0.0  ! partial sums
!                sum_vye=0.0
!                sum_vze=0.0
!                sum_v2e=0.0
!                ne_loc = 0  ! # local electrons
!
!                do p=p_start,p_finish
!
!                    if (particles(p)%label<=ne) then
!                        ! electrons
!                        ne_loc = ne_loc+1
!                        uhx(p) = particles(p)%data%v(1) + 0.5*dt*accx(p)
!                        uhy(p) = particles(p)%data%v(2) + 0.5*dt*accy(p)
!                        uhz(p) = particles(p)%data%v(3) + 0.5*dt*accz(p)
!                        gammah = sqrt(1.0 +(uhx(p)**2 + uhy(p)**2 + uhz(p)**2)/unit_c2)
!                        sum_vxe  = sum_vxe  + uhx(p)/gammah
!                        sum_vye  = sum_vye  + uhy(p)/gammah
!                        sum_vze  = sum_vze  + uhz(p)/gammah
!                        sum_v2e = sum_v2e + gammah-1.
!
!                    endif
!                end do
!                sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2
!
!                ! Find local KE temp
!
!                ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
!                Te_uncor = 511*2./3.*sum_v2e/ne_loc  ! This should equal 3/2 kT for 3v Maxwellian
!                Te0 = Te_eV/1000.  ! normalised electron temp
!                ! exponent of chie should be 1/2 - take 1/3 to soften oscillations
!                chie = (abs(Te0/Te_uncor))**0.5     ! multipliers from Temperature ratio - reset once every cycle
!
!                chie = min(1.25_8,max(chie,0.75_8))  ! Set bounds of +- 50%
!
!                if (my_rank.eq.0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie
!
!                !  3)  Complete full step
!
!                do p=p_start,p_finish
!                    if (particles(p)%label<=ne) then
!                        particles(p)%data%v(1) = (2*chie-1.)*particles(p)%data%v(1) + chie*dt*accx(p)
!                        particles(p)%data%v(2) = (2*chie-1.)*particles(p)%data%v(2) + chie*dt*accy(p)
!                        if (idim==3) particles(p)%data%v(3) = (2*chie-1.)*particles(p)%data%v(3) + chie*dt*accz(p)
!                    endif
!                end do
!
!                delta_Ti=0.
!                delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating
!
!
!            case(INTEGRATOR_SCHEME_LOCAL_NVT_IONS_ONLY)
!
!                ! Conserve kinetic energy of ions only (initial Ti const)
!                mass_eqm = 20.  ! artificial ion mass for eqm stage
!                sum_vxi=0.0  ! partial sums (ions)
!                sum_vyi=0.0
!                sum_vzi=0.0
!                sum_v2i=0.0
!                !  Scale down accelerations if too big
!                do p=p_start,p_finish
!                    accx(p)=accx(p)/max(1.d0,acmax)
!                    accy(p)=accy(p)/max(1.d0,acmax)
!                    accz(p)=accz(p)/max(1.d0,acmax)
!                end do
!
!                do p=p_start,p_finish
!                    uhx(p) = particles(p)%data%v(1) + 0.5*dt*accx(p)
!                    uhy(p) = particles(p)%data%v(2) + 0.5*dt*accy(p)
!                    uhz(p) = particles(p)%data%v(3) + 0.5*dt*accz(p)
!
!                    if (particles(p)%label>=ne) then
!                        ! ions
!                        sum_vxi  = sum_vxi  + uhx(p)
!                        sum_vyi  = sum_vyi  + uhy(p)
!                        sum_vzi  = sum_vzi  + uhz(p)
!                        sum_v2i = sum_v2i + 0.5*(uhx(p)**2+uhy(p)**2+uhz(p)**2)
!                    endif
!                end do
!                sum_2vi = sum_vxi**2 + sum_vyi**2 + sum_vzi**2
!
!                ! Find global KE sums
!                !     call MPI_ALLREDUCE(sum_v2i, global_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
!                Ti_uncor = 511*2./3.*sum_v2i/(p_finish-p_start)  ! This should equal 3/2 kT for 3v Maxwellian
!                Ti0 = Ti_eV/1000.  ! normalised electron temp
!
!
!                chii = sqrt(abs(Ti0/Ti_uncor))
!                chii = min(1.2_8,max(chii,0.6_8))  ! Set bounds of +- 50%
!
!
!                !  3)  Complete full step
!
!                do p=p_start,p_finish
!
!                    if (particles(p)%label>=ne) then
!                        ! make ions lighter for eqm phase
!                        particles(p)%data%v(1) = (2*chii-1.)*particles(p)%data%v(1) + chii*dt*accx(p)
!                        particles(p)%data%v(2) = (2*chii-1.)*particles(p)%data%v(2) + chii*dt*accy(p)
!                        if (idim==3) particles(p)%data%v(3) = (2*chii-1.)*particles(p)%data%v(3) + chii*dt*accz(p)
!                    !           particles(p)%data%v(1) = particles(p)%data%v(1)*sqrt(Ti0/Ti_uncor)
!                    !           particles(p)%data%v(2) = particles(p)%data%v(2)*sqrt(Ti0/Ti_uncor)
!                    !           particles(p)%data%v(3) = particles(p)%data%v(3)*sqrt(Ti0/Ti_uncor)
!
!                    endif
!                end do
!                delta_Ti = 2*Ti0*(1.0/chii**2-1.0)       !  heating
!                if (my_rank==0) then
!                    write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii,' heating:',delta_Ti
!                    write (*,*) 'Max delta-ux: ',delta_u
!                endif

            case (INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
                ! unconstrained motion for negatively charged particles, frozen positively charged particles
                do p = p_start, p_finish
                    if (particles(p)%data%q<0.) then
                        particles(p)%data%v(1:3) = particles(p)%data%v(1:3) + dt * acc(1:3,p) * dimfac(1:3)
                    else
                        particles(p)%data%v(1:3) = 0.
                    endif
                end do

            case default
                ! unconstrained motion by default (scheme=INTEGRATOR_SCHEME_NVE,INTEGRATOR_SCHEME_NONREL)
                do p = p_start, p_finish
                    particles(p)%data%v(1:3) = particles(p)%data%v(1:3) + dt * acc(1:3,p) * dimfac(1:3)
                end do

        end select pusher

    end subroutine velocities



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> 3v particle velocity update (Boris rotation scheme)
    !> Notation follows Birdsall and Langdon:
    !>   "Plasma Physics via Computer Simulations", Chapter 15-4
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine velocities_boris(ips,ipf,delt)
        use physvars
        use module_math_tools
        implicit none
        integer, intent(in) :: ips,ipf
        real*8, intent(in) :: delt

        integer :: p
        real*8 :: beta, gam
        real*8, dimension(3) :: uminus, uprime, uplus, t, s


        do p=ips,ipf
            ! charge/mass*time-constant
            beta   = particles(p)%data%q / (2. * particles(p)%data%m) * delt
            ! first half step with electric field
            uminus = particles(p)%data%v + beta * particles(p)%results%e
            ! gamma factor
            gam    = sqrt( 1.0 + ( dot_product(uminus, uminus) ) / unit_c2 )
            ! rotation with magnetic field
            t      = beta/gam * particles(p)%data%b
            uprime = uminus + cross_product(uminus, t)
            s      = 2. * t / (1 + dot_product(t, t))
            uplus  = uminus + cross_product(uprime, s)
            ! second half step with electric field
            particles(p)%data%v = uplus + beta * particles(p)%results%e
        end do

    end subroutine velocities_boris


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Nonrelativistic particle position update - used with leap-frog scheme
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push_nonrel(ips,ipf,delt)

        use physvars
        integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
        real*8, intent(in) :: delt
        integer :: p

        do p=ips,ipf
            particles(p)%x = particles(p)%x + particles(p)%data%v * delt
        end do

    end subroutine push_nonrel


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Relativistic particle position update - used with leap-frog scheme
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push_rel(ips,ipf,delt)

        use physvars
        integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
        real*8, intent(in) :: delt
        integer :: p
        real*8 :: gam

        do p=ips,ipf
            gam  = sqrt( 1.0 + ( dot_product(particles(p)%data%v,particles(p)%data%v) ) / unit_c2 )

            particles(p)%x = particles(p)%x + particles(p)%data%v/gam * delt
        end do

    end subroutine push_rel


end module module_pusher
