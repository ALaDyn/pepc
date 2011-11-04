!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates programmatic time-dependent changes of the simulation setup
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_workflow
      use treevars
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, public :: workflow_setup = 0 !< time-dependent setup (0 = no time dependence of configuration, other values: see workflow()-routine)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public workflow

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Sets up/modifies configuration according to current simulation time
		!> by calling different
		!> workflows are defined via workflow_setup variable
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine workflow(my_rank, itime, trun, dt)
		  implicit none
		  integer, intent(in) :: my_rank, itime
		  real*8, intent(in) :: trun, dt
		  character(150) :: setup_name

		  select case(workflow_setup)
		    case(0) ! no time variation
              setup_name = 'no time variation'

            case(1) ! temperature rescaling incl drift after every full laser cycle
              setup_name = 'temp. rescaling after full laser cycle'
              call workflow_temp_rescaling_after_full_laser_cycle(itime, trun, dt)

            case(2) ! [PRE 71, 056408 (2005)] P. Hilse et al: Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.
              setup_name = '[PRE 71, 056408 (2005)] P. Hilse et al: Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.'
              call workflow_PRE_71_056408(itime, trun, dt)

            case(3) ! [J. Phys. A: Math. Theor 42 (2009), 214048] Th. Raitza et al: Collision frequency of electrons in laser excited small clusters
              setup_name = '[J. Phys. A: Math. Theor 42 (2009), 214048] Th. Raitza et al: Collision frequency of electrons in laser excited small clusters'
              call workflow_JPhysA_42_214048(itime, trun, dt)

		  end select

          if (my_rank == 0) write( *,'(/"-- WORKFLOW --"/a20,i8," ", a)') 'workflow_setup = ', workflow_setup, setup_name

		end subroutine workflow


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> rescaling of electron and ion temperature incl. drift elimination
        !> after every full laser cycle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_temp_rescaling_after_full_laser_cycle(itime, trun, dt)
          use module_laser
          use module_pusher
          implicit none
          integer, intent(in) :: itime
          real*8, intent(in) :: trun, dt

          real*8 :: remainder

          integrator_scheme = 2
          remainder = mod(1._8*trun, navcycle)
          enable_drift_elimination = ( (remainder >= -0.5_8) .and. (remainder < 0.5_8) )

        end subroutine workflow_temp_rescaling_after_full_laser_cycle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Workflow from [PRE 71, 056408 (2005)]
        !> Hilse et al.:
        !> Collisional absorption of dense plasmas in strong laser
        !> fields: Quantum statistical results and simulation
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_PRE_71_056408(itime, trun, dt)
          use module_laser
          use module_pusher
          use module_units
          implicit none
          integer, intent(in) :: itime
          real*8, intent(in) :: trun, dt

          real*8 :: time_fs

          time_fs = trun*unit_t0_in_fs

          if      (time_fs <= 1.25) then        ! phase of 'establishment of correlations':
            beam_config_in = 0                   ! relaxation (not into equilibrium due to large ion mass)
            call laser_setup()                   ! ==> 2-temp. plasma
            integrator_scheme = 1

          elseif (time_fs <= 2.25) then        ! 'thermalization':
            beam_config_in = 0                   ! velocity rescaling ==> coupling to heat bath
            call laser_setup()
            integrator_scheme = 2
            enable_drift_elimination = .true.

          elseif (time_fs <= 3.50) then        ! 'heat bath off', equilibrium
            beam_config_in = 0                   ! (pot. and kin. energy constant)
            call laser_setup()
            integrator_scheme = 1
            enable_drift_elimination = .false.

          else                                 ! laser switched on
            beam_config_in = 3
            call laser_setup()
            integrator_scheme = 1

          endif

        end subroutine workflow_PRE_71_056408


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Workflow from [J. Phys. A: Math. Theor 42 (2009), 214048]
        !> Th. Raitza et al.:
        !> Collision frequency of electrons in laser excited small clusters
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_JPhysA_42_214048(itime, trun, dt)
          use module_laser
          use module_pusher
          use module_units
          implicit none
          integer, intent(in) :: itime
          real*8, intent(in) :: trun, dt
          real*8, parameter :: pulse_length_fs = 100.0
          real*8 :: time_fs

          integer, save :: origbeamconfig

          if (itime == 1) origbeamconfig = beam_config_in

          time_fs = trun*unit_t0_in_fs

          if      (time_fs <= pulse_length_fs) then       ! initial laser pulse
            beam_config_in = origbeamconfig
            t_pulse = pulse_length_fs
            call laser_setup()
            integrator_scheme = 1

          else                                 ! laser switched off
            beam_config_in = 0
            call laser_setup()
            integrator_scheme = 1

          endif

        end subroutine workflow_JPhysA_42_214048


end module module_workflow
