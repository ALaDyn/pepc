! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2016 Juelich Supercomputing Centre, 
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

      integer, public, parameter :: WORKFLOW_STEP_PRE  = 0
      integer, public, parameter :: WORKFLOW_STEP_POST = 1

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
	subroutine workflow(my_rank, itime, trun, dt, workflow_step)
	  use physvars, only: workflow_setup
          use module_pusher
	  implicit none
	  integer, intent(in) :: my_rank, itime, workflow_step
	  real*8, intent(in) :: trun
      real*8, intent(inout) :: dt
	  character(200) :: setup_name
	  integer :: stage = -1

	  select case(workflow_setup)
	    case(0) ! no time variation
              setup_name = 'no time variation'

            case(1) ! temperature rescaling incl drift after every full laser cycle
              setup_name = 'temp. rescaling after full laser cycle'
              call workflow_temp_rescaling_after_full_laser_cycle(itime, trun, dt, workflow_step, stage)

            case(2) ! [PRE 71, 056408 (2005)] P. Hilse et al: Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.
              setup_name = '[PRE 71, 056408 (2005)] P. Hilse et al: Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.'
              call workflow_PRE_71_056408(itime, trun, dt, workflow_step, stage)

            case(4) ! [PRE 71, 056408 (2005)] P. Hilse et al: Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.
              setup_name = '[PRE 71, 056408 (2005)] P. Hilse et al: Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation; Fixed vo/vte'
              call workflow_PRE_71_056408_fixed_v0_vte(itime, trun, dt, workflow_step, stage)

            case(5) ! [PRE 57, 4698 (1998)] Pfalzner & Gibbon: Direct calculation of inverse-bremsstrahlung absoprtion...
                    ! similar to [PRE 71, 056408 (2005)] but with permanent temperature rescaling
              setup_name = '[PRE 57, 4698 (1998)] Pfalzner & Gibbon: Direct calculation of inverse-bremsstrahlung absoprtion... with Nose-Hoover thermostat'
              call workflow_PRE_57_4698_generic(itime, trun, dt, workflow_step, stage, INTEGRATOR_SCHEME_NVT_NOSE_HOOVER)

            case(6) ! [PRE 57, 4698 (1998)] Pfalzner & Gibbon: Direct calculation of inverse-bremsstrahlung absoprtion...
                    ! similar to [PRE 71, 056408 (2005)] but with permanent temperature rescaling
              setup_name = '[PRE 57, 4698 (1998)] Pfalzner & Gibbon: Direct calculation of inverse-bremsstrahlung absoprtion... with Berendsen thermostat'
              call workflow_PRE_57_4698_generic(itime, trun, dt, workflow_step, stage, INTEGRATOR_SCHEME_NVT_BERENDSEN)

            case(7) ! [PRE 57, 4698 (1998)] Pfalzner & Gibbon: Direct calculation of inverse-bremsstrahlung absoprtion...
                    ! similar to [PRE 71, 056408 (2005)] but with permanent temperature rescaling
              setup_name = '[PRE 57, 4698 (1998)] Pfalzner & Gibbon: Direct calculation of inverse-bremsstrahlung absoprtion... with brute-force velocity scaling thermostat'
              call workflow_PRE_57_4698_generic(itime, trun, dt, workflow_step, stage, INTEGRATOR_SCHEME_NVT)

	  end select

          if (my_rank == 0) write( *,'(/"-- WORKFLOW --"/a20,i8," ", a/a20,i8/)') 'workflow_setup = ', workflow_setup, setup_name, 'stage = ', stage

		end subroutine workflow


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> rescaling of electron and ion temperature incl. drift elimination
        !> after every full laser cycle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_temp_rescaling_after_full_laser_cycle(itime, trun, dt, workflow_step, stage)
          use module_laser
          use module_pusher
          implicit none
          integer, intent(in) :: itime, workflow_step
          real*8, intent(in) :: trun, dt
          integer, intent(out) :: stage

          real*8 :: remainder

          integrator_scheme = INTEGRATOR_SCHEME_NVT
          remainder = mod(1._8*trun, navcycle)
          enable_drift_elimination = ( (remainder >= -0.5_8) .and. (remainder < 0.5_8) )
          stage = 0

        end subroutine workflow_temp_rescaling_after_full_laser_cycle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Workflow from [PRE 71, 056408 (2005)]
        !> Hilse et al.:
        !> Collisional absorption of dense plasmas in strong laser
        !> fields: Quantum statistical results and simulation
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_PRE_71_056408(itime, trun, dt, workflow_step, stage)
          use module_laser
          use module_pusher
          use module_units
          use module_interaction_specific, only : kelbg_invsqrttemp
          use physvars, only: Te, maxdt, dt_ori, &
            MAXDT_OMEGA_PLASMA, MAXDT_VELUPDATE_VTE, MAXDT_POSUPDATE_VTE, &
            MAXDT_OMEGA_LASER,  MAXDT_VELUPDATE_VTO, MAXDT_POSUPDATE_VTO
          implicit none
          integer, intent(in) :: itime, workflow_step
          real*8, intent(in) :: trun
          real*8, intent(inout) :: dt
          integer, intent(out) :: stage
          logical, save :: firstcall = .true.
          integer, save :: origbeamconfig

          real*8 :: time_fs

          if (firstcall) then
            firstcall = .false.
            origbeamconfig = beam_config_in
          endif

          time_fs = trun*unit_t0_in_fs

          select case (workflow_step)

            case (WORKFLOW_STEP_PRE)

              if      (time_fs <= 1.25) then        ! phase of 'establishment of correlations':
                beam_config_in = 0                   ! relaxation (not into equilibrium due to large ion mass)
                call laser_setup()                   ! ==> 2-temp. plasma
                integrator_scheme = INTEGRATOR_SCHEME_NVE
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTE)))
                stage = 1

                ! we set the temperature of the kelbg interaction to the desired electron temperature here instead of actual temperature
                ! to avoid problems due to invalid rescaling in next stage
                kelbg_invsqrttemp = 1._8/sqrt(Te)

              elseif (time_fs <= 2.25) then        ! 'thermalization':
                beam_config_in = 0                   ! velocity rescaling ==> coupling to heat bath
                call laser_setup()
                integrator_scheme = INTEGRATOR_SCHEME_NVT
                enable_drift_elimination = .true.
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTE)))
                stage = 2

                ! we set the temperature of the kelbg interaction to the desired electron temperature here instead of actual temperature
                ! to avoid problems due to invalid rescaling in this stage
                kelbg_invsqrttemp = 1._8/sqrt(Te)

              elseif (time_fs <= 3.50) then        ! 'heat bath off', equilibrium
                beam_config_in = 0                   ! (pot. and kin. energy constant)
                call laser_setup()
                integrator_scheme = INTEGRATOR_SCHEME_NVE
                enable_drift_elimination = .false.
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTE)))
                stage = 3

              else                                 ! laser switched on
                beam_config_in = origbeamconfig
                call laser_setup()
                integrator_scheme = INTEGRATOR_SCHEME_NVE
                enable_drift_elimination = .false.
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTO)))
                stage = 4

              endif

          case (WORKFLOW_STEP_POST)
            ! no special diagnostics here
          end select


        end subroutine workflow_PRE_71_056408



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Workflow from [PRE 71, 056408 (2005)]
        !> Hilse et al.:
        !> Collisional absorption of dense plasmas in strong laser
        !> fields: Quantum statistical results and simulation
        !>
        !> Additionally, v0/vte is kept at its initial value by adjusting 
        !> the laser intensity accordingly
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_PRE_71_056408_fixed_v0_vte(itime, trun, dt, workflow_step, stage)
          use module_laser
          use module_pusher
          use module_units
          use module_interaction_specific, only : kelbg_invsqrttemp
          use physvars, only : mass_e, qe
          use physvars, only : vte
          use module_prepare, only : adjust_maxdt
          implicit none
          integer, intent(in) :: itime, workflow_step
          real*8, intent(in) :: trun
          real*8, intent(inout) :: dt
          integer, intent(out) :: stage

          call workflow_PRE_71_056408(itime, trun, dt, workflow_step, stage)
          
          if (stage == 4) then
            vosc = vosc_vte * vte
            E0   = vosc*mass_e*omega/abs(qe)
            I0_Wpercm2 = (unit_epsilon0 * unit_c * E0**2 / 2. ) * unit_P0_in_W / (100*unit_abohr_in_m)**2

            call laser_setup()

            call adjust_maxdt()
          endif

        end subroutine workflow_PRE_71_056408_fixed_v0_vte


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Workflow from [PRE 57, 4698 (1998)]
        !> Pfalzner & Gibbon:
        !> Direct calculation of inverse-bremsstrahlung absorption in ...
        !>
        !> with possibility to choose thermostat during final simulation stage
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine workflow_PRE_57_4698_generic(itime, trun, dt, workflow_step, stage, thermostat)
          use module_laser
          use module_pusher
          use module_units
          use module_interaction_specific, only : kelbg_invsqrttemp
          use physvars, only: Te, maxdt, dt_ori, &
            MAXDT_OMEGA_PLASMA, MAXDT_VELUPDATE_VTE, MAXDT_POSUPDATE_VTE, &
            MAXDT_OMEGA_LASER,  MAXDT_VELUPDATE_VTO, MAXDT_POSUPDATE_VTO
          implicit none
          integer, intent(in) :: itime, workflow_step
          real*8, intent(in) :: trun
          real*8, intent(inout) :: dt
          integer, intent(out) :: stage
          integer, intent(in) :: thermostat
          logical, save :: firstcall = .true.
          integer, save :: origbeamconfig

          real*8 :: time_fs

          if (firstcall) then
            firstcall = .false.
            origbeamconfig = beam_config_in
          endif

          time_fs = trun*unit_t0_in_fs

          select case (workflow_step)

            case (WORKFLOW_STEP_PRE)

              if      (time_fs <= 1.25) then        ! phase of 'establishment of correlations':
                beam_config_in = 0                   ! relaxation (not into equilibrium due to large ion mass)
                call laser_setup()                   ! ==> 2-temp. plasma
                integrator_scheme = INTEGRATOR_SCHEME_NVE
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTE)))
                stage = 1
                

                ! we set the temperature of the kelbg interaction to the desired electron temperature here instead of actual temperature
                ! to avoid problems due to invalid rescaling in next stage
                kelbg_invsqrttemp = 1._8/sqrt(Te)

              elseif (time_fs <= 2.25) then        ! 'thermalization':
                beam_config_in = 0                   ! velocity rescaling ==> coupling to heat bath
                call laser_setup()
                integrator_scheme = INTEGRATOR_SCHEME_NVT
                enable_drift_elimination = .true.
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTE)))
                stage = 2

                ! we set the temperature of the kelbg interaction to the desired electron temperature here instead of actual temperature
                ! to avoid problems due to invalid rescaling in this stage
                kelbg_invsqrttemp = 1._8/sqrt(Te)

              elseif (time_fs <= 3.50) then        ! we still keep the heat bath on and will later measure the amount of removed energy during laser heating
                beam_config_in = 0                   ! (pot. and kin. energy constant)
                call laser_setup()
                integrator_scheme = INTEGRATOR_SCHEME_NVE
                enable_drift_elimination = .false.
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTE)))
                stage = 3

              else                                 ! laser switched on
                beam_config_in = origbeamconfig
                call laser_setup()
                integrator_scheme = thermostat
                dt = min(dt_ori, minval(maxdt(MAXDT_OMEGA_PLASMA:MAXDT_VELUPDATE_VTO))) !adjust timestep to resolve laser oscillation
                stage = 4

              endif

          case (WORKFLOW_STEP_POST)
            ! no special diagnostics here
          end select


        end subroutine workflow_PRE_57_4698_generic

end module module_workflow
