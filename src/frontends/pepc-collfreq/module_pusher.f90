! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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
    use module_histogram
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
    integer, public, parameter :: INTEGRATOR_SCHEME_NVT_NOSE_HOOVER = 9 ! NVT - global Te, Ti conserved, Nose-Hover (extended Lagrangian) approach following [J.Chem.Phys 107, 9514]
    integer, public, parameter :: INTEGRATOR_SCHEME_NVT_BERENDSEN   = 10! NVT - global Te, Ti (damped) velocity rescaling by Berendsen

    integer, public :: integrator_scheme = INTEGRATOR_SCHEME_NVE
    real*8, public :: Te0 = 0., Te_uncor = 0., chie = 0., delta_Te = 0.
    real*8, public :: Ti0 = 0., Ti_uncor = 0., chii = 0., delta_Ti = 0.
    logical, public :: enable_drift_elimination = .false. !< if .true., the global particle drift is included during velocity rescaling for constant temperature regime (i.e. it is eliminated)
    
    real*8, public :: nose_hoover_e(3) = [0., 0., 0.]
    real*8, public :: nose_hoover_i(3) = [0., 0., 0.]
    real*8, public :: nose_hoover_Q_e = -1.0
    real*8, public :: nose_hoover_Q_i = -1.0
    integer, public, parameter :: NH_ETA  = 1
    integer, public, parameter :: NH_VETA = 2
    integer, public, parameter :: NH_AETA = 3
    
    real*8, public :: tau_temp_relaxation = 5.0 !< temperature relaxation time in units of the (finally determined) simulation timestep dt
    
    logical :: nh_firstcall = .true.
    logical :: nvt_firstcall = .true.
    logical :: nve_firstcall = .true.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public integrator

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type(t_threedhist) :: velhist
    type(t_threedhist) :: velhist_laseravg

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Velocity and position update - integration of the equation of motion
    !> Boundary conditions are also applied here
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine integrator(particles,scheme)
        use module_pepc_types
        use module_io
        use physvars, only: dt, my_rank
        implicit none
        type(t_particle), intent(inout) :: particles(:)
        integer :: scheme

        if (my_rank==0) then
            write(file_stdout  ,'(/"-- PUSHER --"/a20,i8)') 'scheme = ',scheme
            write(file_pepc_out,'(/"-- PUSHER --"/a20,i8)') 'scheme = ',scheme
        endif

        pusher: select case(scheme)

            case(INTEGRATOR_SCHEME_NVE,                   &
                 INTEGRATOR_SCHEME_NVT,                   &
                 INTEGRATOR_SCHEME_NVT_IONS_FROZEN,       &
                 INTEGRATOR_SCHEME_LOCAL_NVT_IONS_FROZEN, &
                 INTEGRATOR_SCHEME_LOCAL_NVT_IONS_ONLY,   &
                 INTEGRATOR_SCHEME_NVE_IONS_FROZEN,       &
                 INTEGRATOR_SCHEME_NVT_NOSE_HOOVER,       &
                 INTEGRATOR_SCHEME_NVT_BERENDSEN           )
                call velocities(particles,scheme)  ! pure ES, NVT ensembles
                call push_rel(particles,dt)  ! update positions

            case(INTEGRATOR_SCHEME_FULL_EM)
                call velocities_boris(particles,dt)  ! full EM pusher (all E, B components)
                call push_rel(particles,dt)  ! update positions

            case(INTEGRATOR_SCHEME_NONREL)
                call velocities(particles,scheme)
                call push_nonrel(particles,dt)       ! nonrelativistic push
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
    subroutine velocities(particles,scheme)
        use module_pepc_types
        use physvars
        implicit none
        include 'mpif.h'

        type(t_particle), intent(inout) :: particles(:)
        integer, intent(in) :: scheme

        integer(kind_particle) ::  i, p
        integer(kind_default) :: ierr
        real*8 :: acmax, delta_u, gammah
        real*8 :: uprime(1:3), uprime2
        real*8 :: ex_e, ex_i, Ue_uncor, Ui_uncor

        integer, parameter :: V2E = 1
        integer, parameter :: VEX = 2
        integer, parameter :: VEY = 3
        integer, parameter :: VEZ = 4
        integer, parameter :: V2I = 5
        integer, parameter :: VIX = 6
        integer, parameter :: VIY = 7
        integer, parameter :: VIZ = 8

        real*8 :: sums(1:8), dimfac(3)
        real*8, allocatable, dimension(:,:) :: uh, acc

        allocate(uh(3, size(particles, kind=kind_particle)))
        allocate(acc(3, size(particles, kind=kind_particle)))

        !  Available ensemble modes
        !      1 = NVE - total energy conserved
        !      2 = NVT - global Te, Ti conserved
        !      3 = global NVT electron Te conserved; ions frozen
        !      4 = local NVT: each PE keeps Te clamped; ions frozen
        !      5 = local NVT, ions only; electrons added at end of run

        ! Accelerations
        acmax=0.
        do i=1,size(particles, kind=kind(p))
            acc(1:3,i) = particles(i)%data%q * particles(i)%results%e(1:3) / particles(i)%data%m
        end do

        acmax = maxval(abs(acc(:,:)))

        dimfac         = 0._8
        dimfac(1:idim) = 1._8

        delta_u = acmax*dt

        pusher: select case(scheme)  ! be careful: there is another select-case block in module_energies::energy_kin()

            case(INTEGRATOR_SCHEME_NVT, &
                 INTEGRATOR_SCHEME_NVT_BERENDSEN )
                ! Conserve kinetic energies of electrons and ions (initial Te const)
                ! adapted from
                !  Allen and Tildesley p230, Brown & Clark, Mol. Phys. 51, 1243 (1984)

                !  Definitions:
                !
                !   unconstrained velocities    uh(x,y,z) = v'(x,y,z)
                !
                !  Berendsen: smooth velocity scaling with certain time constant

                sums = 0.0

                !  1)  Unconstrained half-step for electrons and ions
                do p=1,size(particles,kind=kind(p))
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

                ! drift is averaged cumulative velocity
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

                pusher_inner: select case (scheme)
                  case (INTEGRATOR_SCHEME_NVT) ! allen & Tildesley resp. Brown & Clarke
                    chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
                    chii = sqrt(abs(Ti0/Ti_uncor))

                    ! reduce artificial friction factor by 1./2. for both species - this does not comply with Allen & Tildesley but works fine :-S
                    chie = 1./( (1./2.)*(1./chie - 1) + 1.)
                    chii = 1./( (1./2.)*(1./chii - 1) + 1.)
                  case (INTEGRATOR_SCHEME_NVT_BERENDSEN)
                    chie = sqrt(1+dt/tau_temp_relaxation*(Te0/Te_uncor - 1.))
                    chii = sqrt(1+dt/tau_temp_relaxation*(Ti0/Ti_uncor - 1.))
                end select pusher_inner

                !  3)  Complete full step
                do p=1,size(particles,kind=kind(p))
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

                call NVT_diagnostics(particles, Te_uncor, Ti_uncor)
                
            case (INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
                ! unconstrained motion for negatively charged particles, frozen positively charged particles
                do p = 1, size(particles, kind=kind(p))
                    if (particles(p)%data%q<0.) then
                        particles(p)%data%v(1:3) = particles(p)%data%v(1:3) + dt * acc(1:3,p) * dimfac(1:3)
                    else
                        particles(p)%data%v(1:3) = 0.
                    endif
                end do
                
            case (INTEGRATOR_SCHEME_NVT_NOSE_HOOVER)
                ! NVT - global Te, Ti conserved, Nose-Hoover (extended Lagrangian) approach following [J.Chem.Phys 107, 9514]
                
                ! determine current temperatures
                sums = 0.0

                do p=1, size(particles,kind=kind(p))
                    uprime(1:3) = particles(p)%data%v(1:3)
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
                call MPI_ALLREDUCE(MPI_IN_PLACE, sums, size(sums), MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
                
                ! drift is averaged cumulative velocity
                sums(V2E:VEZ) = sums(V2E:VEZ) / ne
                sums(V2I:VIZ) = sums(V2I:VIZ) / ni

                if (.not. enable_drift_elimination) then
                    Ue_uncor = mass_e/2.*ne*(sums(V2E) - dot_product(sums(VEX:VEZ),sums(VEX:VEZ))) ! This should equal 3/2NkT
                    Ui_uncor = mass_i/2.*ni*(sums(V2I) - dot_product(sums(VIX:VIZ),sums(VIX:VIZ)))
                else
                    Ue_uncor = mass_e/2.*ne*sums(V2E)  ! This should equal 3/2NkT
                    Ui_uncor = mass_i/2.*ni*sums(V2I)
                endif


                Te0 = Te  ! normalised (desired) electron temp
                nose_hoover_e(NH_AETA) = ( Ue_uncor - 3./2.*ne*unit_kB*Te0 )/nose_hoover_Q_e
                nose_hoover_e(NH_VETA) = nose_hoover_e(NH_VETA) + nose_hoover_e(NH_AETA) * dt
                nose_hoover_e(NH_ETA ) = nose_hoover_e(NH_ETA ) + nose_hoover_e(NH_VETA) * dt
                ex_e                   = exp(-dt/2.*nose_hoover_e(NH_VETA))
                
                Ti0 = Ti  ! normalised (desired) ion temp
                nose_hoover_i(NH_AETA) = ( Ui_uncor - 3./2.*ni*unit_kB*Ti0 )/nose_hoover_Q_i
                nose_hoover_i(NH_VETA) = nose_hoover_i(NH_VETA) + nose_hoover_i(NH_AETA) * dt
                nose_hoover_i(NH_ETA ) = nose_hoover_i(NH_ETA ) + nose_hoover_i(NH_VETA) * dt
                ex_i                   = exp(-dt/2.*nose_hoover_i(NH_VETA))

                do p = 1, size(particles, kind=kind(p))
                    if (particles(p)%data%q<0.) then
                        particles(p)%data%v(1:3) = particles(p)%data%v(1:3)*ex_e*ex_e + dt * acc(1:3,p) * ex_e
                    else
                        particles(p)%data%v(1:3) = particles(p)%data%v(1:3)*ex_i*ex_i + dt * acc(1:3,p) * ex_i
                    endif
                end do
                
                call nose_hoover_diagnostics(particles)

            case default
                ! unconstrained motion by default (scheme=INTEGRATOR_SCHEME_NVE,INTEGRATOR_SCHEME_NONREL)
                do p = 1, size(particles, kind=kind(p))
                    particles(p)%data%v(1:3) = particles(p)%data%v(1:3) + dt * acc(1:3,p) * dimfac(1:3)
                end do

                call NVE_diagnostics(particles)
                
        end select pusher

        deallocate(uh, acc)

    end subroutine velocities



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> 3v particle velocity update (Boris rotation scheme)
    !> Notation follows Birdsall and Langdon:
    !>   "Plasma Physics via Computer Simulations", Chapter 15-4
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine velocities_boris(particles, delt)
        use module_math_tools
        use module_pepc_types
        implicit none
        type(t_particle), intent(inout) :: particles(:)
        real*8, intent(in) :: delt

        integer :: p
        real*8 :: beta, gam
        real*8, dimension(3) :: uminus, uprime, uplus, t, s


        do p=1, size(particles, kind=kind(p))
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
    subroutine push_nonrel(particles,delt)
        use module_pepc_types
        type(t_particle), intent(inout) :: particles(:)
        real*8, intent(in) :: delt
        integer :: p

        do p=1,size(particles, kind=kind(p))
          particles(p)%x = particles(p)%x + particles(p)%data%v * delt
        end do

    end subroutine push_nonrel


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Relativistic particle position update - used with leap-frog scheme
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push_rel(particles,delt)
        use module_pepc_types
        type(t_particle), intent(inout) :: particles(:)
        real*8, intent(in) :: delt
        integer :: p
        real*8 :: gam

        do p=1, size(particles, kind=kind(p))
            gam  = sqrt( 1.0 + ( dot_product(particles(p)%data%v,particles(p)%data%v) ) / unit_c2 )

            particles(p)%x = particles(p)%x + particles(p)%data%v/gam * delt
        end do

    end subroutine push_rel


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Computes value of Nose_Hoover Hamiltonian and prints it
    !> (and some other diagnostics) to a file
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nose_hoover_diagnostics(particles)
      use module_pepc_types
      use module_units
      use physvars
      implicit none
      include 'mpif.h'
      type(t_particle), intent(in) :: particles(:)
      real*8 :: ex_e, ex_i
      real*8 :: epot_e, epot_i, H_e, H_i
      real*8 :: eta_e, veta_e, eta_i, veta_i, Ue_uncor, Ui_uncor
      real*8 :: uprime(1:3), uprime2, gammah, acc(1:3)
      integer(kind_particle) :: p
      integer(kind_default) :: ierr
      integer, parameter :: file_nose_hoover_dat = 93
      
      integer, parameter :: V2E = 1
      integer, parameter :: VEX = 2
      integer, parameter :: VEY = 3
      integer, parameter :: VEZ = 4
      integer, parameter :: V2I = 5
      integer, parameter :: VIX = 6
      integer, parameter :: VIY = 7
      integer, parameter :: VIZ = 8

      real*8 :: sums(1:8)
      
      ! initialize velocity histogram
      !call histograms_init(p_finish-p_start)
      
      ! determine current temperatures
      sums = 0.0
      
      ex_e = exp(+dt/2.*nose_hoover_e(NH_VETA))
      ex_i = exp(+dt/2.*nose_hoover_i(NH_VETA))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Kinetic energy
      ! Find local KE sums
      do p=1, size(particles, kind=kind(p))
      
        acc(1:3) = particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m
        

        if (particles(p)%data%q < 0.) then
          ! electrons

          ! perform one half step backwards
          uprime(1:3) = particles(p)%data%v(1:3) * (1 + ex_e*ex_e)/2. - acc(1:3)*dt/2.*ex_e
          uprime2     = dot_product(uprime,uprime)
          gammah      = sqrt(1.0 + uprime2/unit_c2)

          ! add values to velocity histogram (only for electrons)
          !call histograms_add(uprime(1:3))

          sums(V2E)     = sums(V2E)     + uprime2 / gammah**2.
          sums(VEX:VEZ) = sums(VEX:VEZ) + uprime  / gammah
        else if (particles(p)%data%q > 0.) then
          ! ions

          ! perform one half step backwards
          uprime(1:3) = particles(p)%data%v(1:3) * (1 + ex_i*ex_i)/2. - acc(1:3)*dt/2.*ex_i
          uprime2     = dot_product(uprime,uprime)
          gammah      = sqrt(1.0 + uprime2/unit_c2)

          sums(V2I)     = sums(V2I)     + uprime2 / gammah**2.
          sums(VIX:VIZ) = sums(VIX:VIZ) + uprime  / gammah
        endif
      end do
      
      ! Find global KE sums
      call MPI_ALLREDUCE(MPI_IN_PLACE, sums, size(sums), MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

      ! drift is averaged cumulative velocity
      sums(V2E:VEZ) = sums(V2E:VEZ) / ne
      sums(V2I:VIZ) = sums(V2I:VIZ) / ni

      if (.not. enable_drift_elimination) then
        Ue_uncor = mass_e/2.*ne*(sums(V2E) - dot_product(sums(VEX:VEZ),sums(VEX:VEZ))) ! This should equal 3/2NkT
        Ui_uncor = mass_i/2.*ni*(sums(V2I) - dot_product(sums(VIX:VIZ),sums(VIX:VIZ)))
      else
        Ue_uncor = mass_e/2.*ne*sums(V2E)  ! This should equal 3/2NkT
        Ui_uncor = mass_i/2.*ni*sums(V2I)
      endif
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Potential energy
      epot_e = 0.
      epot_i = 0.
      
      ! Find local PE sums
      do p=1, size(particles, kind=kind(p))
      
        if (particles(p)%data%q < 0.) then
          ! electrons
          epot_e = epot_e + 0.5 * particles(p)%data%q * particles(p)%results%pot
        else if (particles(p)%data%q > 0.) then
          ! ions
          epot_i = epot_i + 0.5 * particles(p)%data%q * particles(p)%results%pot
        endif
      end do
      
      ! Find global PE sums
      call MPI_ALLREDUCE(MPI_IN_PLACE, epot_e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, epot_i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! extended coordinate and velocity
      eta_e  = nose_hoover_e(NH_ETA) - nose_hoover_e(NH_VETA) * dt/2.
      veta_e = nose_hoover_e(NH_VETA)

      eta_i  = nose_hoover_i(NH_ETA) - nose_hoover_i(NH_VETA) * dt/2.
      veta_i = nose_hoover_i(NH_VETA)
      
      Te0 = Te  ! normalised (desired) electron temp
      Ti0 = Ti  ! normalised (desired) ion temp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! extended Hamiltonian    
      
      H_e = Ue_uncor + epot_e + nose_hoover_Q_e * veta_e*veta_e/2. + 3.*unit_kB*ne*Te0*eta_e
      H_i = Ui_uncor + epot_i + nose_hoover_Q_i * veta_i*veta_i/2. + 3.*unit_kB*ni*Ti0*eta_i
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! output
      if (my_rank == 0) then
        if (nh_firstcall .and. .not. restart) then
          nh_firstcall = .false.
          open(file_nose_hoover_dat, FILE='nose_hoover.dat',STATUS='UNKNOWN', POSITION = 'REWIND')
          write(file_nose_hoover_dat,'("#",16(1x,a20))') "time", &
                                                         "epot_e",   "Te_uncor",   "Te0", "eta_e", "v_eta_e", "nose_hoover_Q_e", "H_e",  &
                                                         "epot_i",   "Ti_uncor",   "Ti0", "eta_i", "v_eta_i", "nose_hoover_Q_i", "H_i",  &
                                                         "H_e + H_i"
        else
          open(file_nose_hoover_dat, FILE='nose_hoover.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
        endif

        write(file_nose_hoover_dat,'(" ",1(1x,f20.6),15(1x,1pe20.10))') trun*unit_t0_in_fs, &
                                                         epot_e,   Te_uncor,   Te0, eta_e, veta_e, nose_hoover_Q_e, H_e, &
                                                         epot_i,   Ti_uncor,   Ti0, eta_i, veta_i, nose_hoover_Q_i, H_i, &
                                                         H_e + H_i
        close(file_nose_hoover_dat)
      endif

    
      ! write velocity histogram
      !call histograms_finalize(Ue_uncor, 'NVT Nose-Hoover thermostat')
    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Computes removed energy etc. and prints it
    !> (and some other diagnostics) to a file
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine NVT_diagnostics(particles,Te_before,Ti_before)
      use module_units
      use module_pepc_types
      use physvars
      implicit none
      include 'mpif.h'
      type(t_particle), intent(in) :: particles(:)
      real*8, intent(in) :: Te_before, Ti_before
      real*8 :: betae, betai
      real*8 :: epot_e, epot_i, H_e, H_i
      real*8 :: Ue_uncor, Ui_uncor, delta_Ue, delta_Ui
      real*8 :: uprime(1:3), uprime2, gammah, acc(1:3)
      integer(kind_particle) :: p
      integer(kind_default) :: ierr
      integer, parameter :: file_thermostat_dat = 93
      
      real*8, save :: delta_Ue_cum = 0.
      real*8, save :: delta_Ui_cum = 0.
      
      integer, parameter :: V2E = 1
      integer, parameter :: VEX = 2
      integer, parameter :: VEY = 3
      integer, parameter :: VEZ = 4
      integer, parameter :: V2I = 5
      integer, parameter :: VIX = 6
      integer, parameter :: VIY = 7
      integer, parameter :: VIZ = 8

      real*8 :: sums(1:8)
      
      ! initialize velocity histogram
      !call histograms_init(p_finish-p_start)
      
      ! determine current temperatures
      sums = 0.0
      
      betae = chie/(2*chie-1.)
      betai = chii/(2*chii-1.)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Kinetic energy
      ! Find local KE sums
      do p=1, size(particles, kind=kind(p))
      
        acc(1:3) = particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m
        
        if (particles(p)%data%q < 0.) then
          ! electrons

          ! perform one half step backwards
          uprime(1:3) = betae*(particles(p)%data%v(1:3) - acc(1:3) * dt/2.)
          uprime2     = dot_product(uprime,uprime)
          gammah      = sqrt(1.0 + uprime2/unit_c2)
          
          ! add values to velocity histogram (only for electrons)
          !call histograms_add(uprime(1:3))

          sums(V2E)     = sums(V2E)     + uprime2 / gammah**2.
          sums(VEX:VEZ) = sums(VEX:VEZ) + uprime  / gammah
        else if (particles(p)%data%q > 0.) then
          ! ions

          ! perform one half step backwards
          uprime(1:3) = betai*(particles(p)%data%v(1:3) - acc(1:3) * dt/2.)
          uprime2     = dot_product(uprime,uprime)
          gammah      = sqrt(1.0 + uprime2/unit_c2)

          sums(V2I)     = sums(V2I)     + uprime2 / gammah**2.
          sums(VIX:VIZ) = sums(VIX:VIZ) + uprime  / gammah
        endif
      end do
      
      ! Find global KE sums
      call MPI_ALLREDUCE(MPI_IN_PLACE, sums, size(sums), MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

      ! drift is averaged cumulative velocity
      sums(V2E:VEZ) = sums(V2E:VEZ) / ne
      sums(V2I:VIZ) = sums(V2I:VIZ) / ni

      if (.not. enable_drift_elimination) then
        Ue_uncor = mass_e/2.*ne*(sums(V2E) - dot_product(sums(VEX:VEZ),sums(VEX:VEZ))) ! This should equal 3/2NkT
        Ui_uncor = mass_i/2.*ni*(sums(V2I) - dot_product(sums(VIX:VIZ),sums(VIX:VIZ)))
      else
        Ue_uncor = mass_e/2.*ne*sums(V2E)  ! This should equal 3/2NkT
        Ui_uncor = mass_i/2.*ni*sums(V2I)
      endif
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Potential energy
      epot_e = 0.
      epot_i = 0.
      
      ! Find local PE sums
      do p=1, size(particles, kind=kind(p))
      
        if (particles(p)%data%q < 0.) then
          ! electrons
          epot_e = epot_e + 0.5 * particles(p)%data%q * particles(p)%results%pot
        else if (particles(p)%data%q > 0.) then
          ! ions
          epot_i = epot_i + 0.5 * particles(p)%data%q * particles(p)%results%pot
        endif
      end do
      
      ! Find global PE sums
      call MPI_ALLREDUCE(MPI_IN_PLACE, epot_e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, epot_i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Hamiltonian    
      
      H_e = Ue_uncor + epot_e
      H_i = Ui_uncor + epot_i
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! removed energy
      
      delta_Ue     = 3./2.*ne*unit_kB*Te_before - Ue_uncor
      delta_Ue_cum = delta_Ue_cum + delta_Ue

      delta_Ui     = 3./2.*ni*unit_kB*Ti_before - Ui_uncor
      delta_Ui_cum = delta_Ui_cum + delta_Ui

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! output
      if (my_rank == 0) then
        if (nvt_firstcall .and. .not. restart) then
          nvt_firstcall = .false.
          open(file_thermostat_dat, FILE='thermostat.dat',STATUS='UNKNOWN', POSITION = 'REWIND')
          write(file_thermostat_dat,'("#",18(1x,a20))') "time", &
                                                         "epot_e",   "Ue_before", "Ue_now", "Ue_target", "delta_Ue", "delta_Ue_cum", "H_e",  &
                                                         "epot_i",   "Ui_before", "Ui_now", "Ui_target", "delta_Ui", "delta_Ui_cum", "H_i",  &
                                                         "delta_Ue+i", "delta_Ue+i_cum", "H_e + H_i"
        else
          open(file_thermostat_dat, FILE='thermostat.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
        endif

        write(file_thermostat_dat,'(" ",1(1x,f20.6),17(1x,1pe20.10))') trun*unit_t0_in_fs, &
                                                         epot_e, 3./2.*ne*unit_kB*Te_before, Ue_uncor, 3./2.*ne*unit_kB*Te, delta_Ue, delta_Ue_cum, H_e, &
                                                         epot_i, 3./2.*ni*unit_kB*Ti_before, Ui_uncor, 3./2.*ni*unit_kB*Ti, delta_Ui, delta_Ui_cum, H_i, &
                                                         delta_Ue + delta_Ui, delta_Ue_cum+delta_Ui_cum, H_e + H_i
        close(file_thermostat_dat)
      endif

    

      ! write velocity histogram
      !call histograms_finalize(Ue_uncor, 'NVT thermostat')
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Computes removed energy etc. and prints it
    !> (and some other diagnostics) to a file
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine NVE_diagnostics(particles)
      use module_units
      use module_pepc_types
      use physvars
      implicit none
      include 'mpif.h'
      type(t_particle), intent(in) :: particles(:)
      real*8 :: epot_e, epot_i, H_e, H_i
      real*8 :: Ue_uncor, Ui_uncor
      real*8 :: uprime(1:3), uprime2, gammah, acc(1:3)
      integer(kind_particle) :: p
      integer(kind_default) :: ierr
      integer, parameter :: file_nve_dat = 93
      
      integer, parameter :: V2E = 1
      integer, parameter :: VEX = 2
      integer, parameter :: VEY = 3
      integer, parameter :: VEZ = 4
      integer, parameter :: V2I = 5
      integer, parameter :: VIX = 6
      integer, parameter :: VIY = 7
      integer, parameter :: VIZ = 8

      real*8 :: sums(1:8)
      
      ! initialize velocity histograms
      !call histograms_init(p_finish-p_start)
      
      ! determine current temperatures
      sums = 0.0
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Kinetic energy
      ! Find local KE sums
      do p=1, size(particles, kind=kind(p))
      
        acc(1:3) = particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m
        
        if (particles(p)%data%q < 0.) then
          ! electrons

          ! perform one half step backwards
          uprime(1:3) = particles(p)%data%v(1:3) - acc(1:3) * dt/2.
          uprime2     = dot_product(uprime,uprime)
          gammah      = sqrt(1.0 + uprime2/unit_c2)
          
          ! add values to velocity histogram (only for electrons)
          !call histograms_add(uprime(1:3))

          sums(V2E)     = sums(V2E)     + uprime2 / gammah**2.
          sums(VEX:VEZ) = sums(VEX:VEZ) + uprime  / gammah
        else if (particles(p)%data%q > 0.) then
          ! ions

          ! perform one half step backwards
          uprime(1:3) = particles(p)%data%v(1:3) - acc(1:3) * dt/2.
          uprime2     = dot_product(uprime,uprime)
          gammah      = sqrt(1.0 + uprime2/unit_c2)

          sums(V2I)     = sums(V2I)     + uprime2 / gammah**2.
          sums(VIX:VIZ) = sums(VIX:VIZ) + uprime  / gammah
        endif
      end do
      
      ! Find global KE sums
      call MPI_ALLREDUCE(MPI_IN_PLACE, sums, size(sums), MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

      ! drift is averaged cumulative velocity
      sums(V2E:VEZ) = sums(V2E:VEZ) / ne
      sums(V2I:VIZ) = sums(V2I:VIZ) / ni

      if (.not. enable_drift_elimination) then
        Ue_uncor = mass_e/2.*ne*(sums(V2E) - dot_product(sums(VEX:VEZ),sums(VEX:VEZ))) ! This should equal 3/2NkT
        Ui_uncor = mass_i/2.*ni*(sums(V2I) - dot_product(sums(VIX:VIZ),sums(VIX:VIZ)))
      else
        Ue_uncor = mass_e/2.*ne*sums(V2E)  ! This should equal 3/2NkT
        Ui_uncor = mass_i/2.*ni*sums(V2I)
      endif
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Potential energy
      epot_e = 0.
      epot_i = 0.
      
      ! Find local PE sums
      do p=1, size(particles, kind=kind(p))
      
        if (particles(p)%data%q < 0.) then
          ! electrons
          epot_e = epot_e + 0.5 * particles(p)%data%q * particles(p)%results%pot
        else if (particles(p)%data%q > 0.) then
          ! ions
          epot_i = epot_i + 0.5 * particles(p)%data%q * particles(p)%results%pot
        endif
      end do
      
      ! Find global PE sums
      call MPI_ALLREDUCE(MPI_IN_PLACE, epot_e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, epot_i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Hamiltonian    
      
      H_e = Ue_uncor + epot_e
      H_i = Ui_uncor + epot_i
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! output
      if (my_rank == 0) then
        if (nve_firstcall .and. .not. restart) then
          nve_firstcall = .false.
          open(file_nve_dat, FILE='nve.dat',STATUS='UNKNOWN', POSITION = 'REWIND')
          write(file_nve_dat,'("#",8(1x,a20))')   "time", &
                                                         "epot_e",   "Ue_now", "H_e",  &
                                                         "epot_i",   "Ui_now", "H_i",  &
                                                         "H_e + H_i"
        else
          open(file_nve_dat, FILE='nve.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
        endif

        write(file_nve_dat,'(" ",1(1x,f20.6),7(1x,1pe20.10))') trun*unit_t0_in_fs, &
                                                         epot_e, Ue_uncor, H_e, &
                                                         epot_i, Ui_uncor, H_i, &
                                                         H_e + H_i
        close(file_nve_dat)
      endif

    

      ! write velocity histogram
      !call histograms_finalize(Ue_uncor, 'NVE dynamics')
    end subroutine
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine histograms_init(nvals)
      use module_laser
      use physvars
      implicit none
      integer, intent(in) :: nvals
      integer :: oldcycle, currcycle
      
      return
      
      call velhist%init(nvals)
      
      currcycle = floor((t_laser/dt   )/navcycle)
      oldcycle  = floor((t_laser/dt-1.)/navcycle)
      if ((t_laser > 0.) .and. (currcycle.ne.oldcycle)) then
        call velhist_laseravg%init(ceiling(navcycle)*nvals) ! FIXME: this is in priciple just an estimate, since particle number per processor could vary during the laser cycle
      endif
    end subroutine
    

    subroutine histograms_add(vel)
      implicit none
      real*8, intent(in) :: vel(3)
      
      return
      
      call velhist%add(vel)
      
      if (velhist_laseravg%initialized) then
        call velhist_laseravg%add(vel)
      endif
    end subroutine
    

    subroutine histograms_finalize(Ue_uncor, status)
      use module_units
      use module_laser
      use module_utils
      use physvars
      implicit none
      real*8, intent(in) :: Ue_uncor
      character(*), intent(in) :: status
      character(50) :: histfile
      integer :: newcycle, currcycle
           
      ! these may not change during simulation run since it affects bin position and widthin histogram
      real*8, save :: sigmabins     = -1.
      real*8, save :: sigmabins_avg = -1.      
      
      return

      call create_directory('hist')
      
      write(histfile,'("hist/histogram_ve.",I6.6)') itime
        if (sigmabins < 0.) sigmabins = sqrt(2./3.*Ue_uncor/mass_e/ne)
      call velhist%dump(sqrt(2./3.*Ue_uncor/mass_e/ne), sigmabins, trim(histfile), my_rank, MPI_COMM_PEPC, status)
      call velhist%finalize()
      
      currcycle = floor((t_laser/dt   )/navcycle)
      newcycle  = floor((t_laser/dt+1.)/navcycle)
      if ((currcycle > 0) .and. (currcycle.ne.newcycle)) then
        write(histfile,'("hist/histogram_ve_laseravg.",I6.6)') itime
        if (sigmabins_avg < 0.) sigmabins_avg = sqrt(2./3.*Ue_uncor/mass_e/ne)
        call velhist_laseravg%dump(sqrt(2./3.*Ue_uncor/mass_e/ne), sigmabins_avg, trim(histfile), my_rank, MPI_COMM_PEPC, status)
        call velhist_laseravg%finalize()
      endif
    
    end subroutine

end module module_pusher