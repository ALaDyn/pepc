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

program pepc
  ! pepc modules
  use module_pepc
  use module_pepc_types
  use module_timings
  use module_debug
  use module_vtk
  ! frontend helper routines
  use pepcboris_helper
  use pepcboris_integrator
  use pepcboris_diagnostics

  use pf_mod_verlet, only: pf_verlet_create, pf_verlet_destroy
  use pfm_helper
  use pfm_encap
  use pfm_feval
  use pfm_hooks
  use pf_mod_comm_mpi, only: pf_mpi_create, pf_mpi_destroy, pf_mpi_setup
  use pf_mod_pfasst, only: pf_pfasst_create, pf_pfasst_setup, pf_pfasst_destroy
  use pf_mod_parallel, only: pf_pfasst_run
  use pf_mod_hooks
  use pf_mod_options

  implicit none

  include 'mpif.h'

  ! variables for pfasst
  integer(kind_default):: MPI_COMM_SPACE, MPI_COMM_TIME, mpi_err
  type(pf_pfasst_t) :: pf
  type(pf_nml_t) :: pf_nml
  type(pf_comm_t) :: tcomm
  type(pf_encap_t), target :: encap
  type(pf_sweeper_t), target :: sweeper
  type(level_params_t), pointer :: level_params(:)
  type(app_data_t), pointer :: y0, yend

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable, target :: particles(:)
  integer :: step

  ! Take care of communication stuff
  call pfm_init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)

  ! store communicators and ranks/sizes
  pepcboris_nml%comm_space = MPI_COMM_SPACE
  pepcboris_nml%comm_time  = MPI_COMM_TIME
  call get_mpi_rank(MPI_COMM_WORLD,           pepcboris_nml%rank_world, pepcboris_nml%nrank_world)
  call get_mpi_rank(pepcboris_nml%comm_space, pepcboris_nml%rank_space, pepcboris_nml%nrank_space)
  call get_mpi_rank(pepcboris_nml%comm_time,  pepcboris_nml%rank_time,  pepcboris_nml%nrank_time)

  ! initialize pepc library and MPI
  call pepc_initialize('pepc-boris-pfasst', init_mpi=.false., db_level_in=DBG_STATUS, comm=pepcboris_nml%comm_space)
  ! this is not an MPI-parallel application
  if (pepcboris_nml%nrank_space > 1 .and. pepcboris_nml%rank_world == 0) then
    DEBUG_WARNING(*, "You are requesting spatial parallelism. This is still considered untested. Maybe you wanted to set num_space_instances=num_ranks to only use temporal parallelism?")
  endif
  ! frontend parameter initialization, particle configuration etc.
  call pepcboris_init(particles, dt=pf_nml%tend/pf_nml%nsteps, nt=pf_nml%nsteps)  ! we use the finest (i.e. highest) level here
  ! commit all internal pepc variables
  call pepc_prepare(dim)
  ! prepare table with level-dependent parameters                          # FIXME: v-- this should be pepcboris_nml%rank_time==pepcboris_nml%nrank_time-1, but for that rank the hook is not called apparently
  call pfm_setup_solver_level_params(particles, level_params, pf_nml, dim)
  ! initial potential will be needed for energy computation - using finest level here
  call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, comm=MPI_COMM_SPACE, clearresults=.true.) ! again, use parameters of finest level

  if (.not. pepcboris_nml%setup_params(PARAMS_OMEGAB)**2 - 4._8*pepcboris_nml%setup_params(PARAMS_OMEGAE)**2 > 0) then
    DEBUG_WARNING(*, 'Trapping condition is not fulfilled due to inappropriate choice of PARAMS_OMEGAB and PARAMS_OMEGAE.')
  endif
  ! initial particle dump
  call dump_particles(VTK_STEP_FIRST, 0, 0._8, particles, MPI_COMM_SPACE, do_average=.false.)
  call dump_energy(0._8, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)

  select case (pepcboris_nml%workingmode)
    case (WM_BORIS_SDC, WM_BORIS_MLSDC, WM_BORIS_SDC_REF)
      ! Set up PFASST object
      call pf_mpi_create(tcomm, MPI_COMM_TIME)
      call pfm_encap_init(encap, particles)
      call pf_verlet_create(sweeper, calc_Efield, build_rhs, impl_solver)

      ! if we were requested to perform only single-level SDC, we use finest-level params
      if (pepcboris_nml%workingmode == WM_BORIS_SDC) then
        level_params(1) = level_params(pf_nml%nlevels)
        pf_nml%nlevels = 1
      endif

      call pf_pfasst_create(pf, tcomm, pf_nml%nlevels)
      call pfm_fill_pfasst_object(pf, encap, sweeper, pf_nml, level_params)

      ! Initial conditions for pfasst
      call feval_init(y0, yend, pf_nml%nlevels, pf%levels(pf_nml%nlevels)%levelctx, encap%encapctx, particles)

      call pf_mpi_setup(tcomm, pf)
      call pf_pfasst_setup(pf)

      ! Add user-defined calls, e.g. diagnostics, here
      call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, dump_particles_hook)
      !call pf_add_hook(pf, pf_nml%nlevels, PF_PRE_STEP, dump_particles_hook)

      ! some informative output about what we are actually doing
      if (pepcboris_nml%rank_world == 0) call pf_print_options(pf)

      ! Here we go       pfasst-object, initial value, dt, t_end, number of steps, final solution
      call pf_pfasst_run(pf, c_loc(y0), pf_nml%tend/pf_nml%nsteps, pf_nml%tend, nsteps=pf_nml%nsteps, qend=c_loc(yend))

      ! Remove everything (hopefully)
      call pf_mpi_destroy(tcomm)
      call feval_finalize(y0,yend)
      call pf_pfasst_destroy(pf)
      call pf_verlet_destroy(sweeper)
      call pfm_finalize_solver_level_params(level_params, pf_nml%nlevels)

    case (WM_BORIS)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)
        if (.not. pepcboris_nml%setup_params(PARAMS_OMEGAB)*dt < 1.) then
          DEBUG_WARNING(*, 'Gyrofrequency too high or timestep too large. The gyroradius will increase linearly during simulation. Compare J. comp. Phys. 116, 386 (1995)')
        endif

        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles_velocity_verlet_boris(particles, dt)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call update_velocities_velocity_verlet_boris(particles, dt)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

     case (WM_BORIS_TANALPHA)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)
        if (.not. pepcboris_nml%setup_params(PARAMS_OMEGAB)*dt < 1.) then
          DEBUG_WARNING(*, 'Gyrofrequency too high or timestep too large. The gyroradius will increase linearly during simulation. Compare J. comp. Phys. 116, 386 (1995)')
        endif

        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles_velocity_verlet_boris_tanalpha(particles, dt)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call update_velocities_velocity_verlet_boris_tanalpha(particles, dt)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

    case (WM_BORIS_LEAP_FROG)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)

        call update_velocities_boris(particles,dt/2._8)

        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles(particles, dt)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call backup_velocities(particles) ! for particle and energy dumping we will average over old and new velocities to get velocities on integer timesteps
          call update_velocities_boris(particles, dt)
          ! ATTENTION: here, velocities are defined on timestep step+1/2, that is why we have to average over old and new velocities
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.true.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.true.)
        end do
      end associate

    case (WM_TAJIMA_LEAP_FROG_IMPLICIT)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)

        call update_velocities_tajima_implicit(particles,dt/2._8)

        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles(particles, dt)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call backup_velocities(particles) ! for particle and energy dumping we will average over old and new velocities to get velocities on integer timesteps
          call update_velocities_tajima_implicit(particles, dt)
          ! ATTENTION: here, velocities are defined on timestep step+1/2, that is why we have to average over old and new velocities
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.true.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.true.)
        end do
      end associate

    case (WM_TAJIMA_LEAP_FROG_EXPLICIT)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)

        call update_velocities_tajima_explicit(particles,dt/2._8)

        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles(particles, dt)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call backup_velocities(particles) ! for particle and energy dumping we will average over old and new velocities to get velocities on integer timesteps
          call update_velocities_tajima_explicit(particles, dt)
          ! ATTENTION: here, velocities are defined on timestep step+1/2, that is why we have to average over old and new velocities
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.true.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.true.)
        end do
      end associate

    case (WM_CYCLOTRONIC)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)
        do step=1,nt
          call print_timestep(step, nt, dt)
          call drift_cyclotronic(particles, dt/2._8)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call kick_cyclotronic(particles, dt)
          call drift_cyclotronic(particles, dt/2._8)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

    case (WM_CYCLOTRONIC_NOTAN)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)
        do step=1,nt
          call print_timestep(step, nt, dt)
          call drift_cyclotronic(particles, dt/2._8, no_tan_trans=.true.)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call kick_cyclotronic(particles, dt)
          call drift_cyclotronic(particles, dt/2._8, no_tan_trans=.true.)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

    case (WM_BORIS_PATACCHINI)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt, &
                 params => pepcboris_nml%setup_params)
        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles(particles, dt/2._8)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call kick_boris_patacchini(particles, dt)
          call push_particles(particles, dt/2._8)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

    case (WM_BORIS_PATACCHINI_NOTAN)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt, &
                 params => pepcboris_nml%setup_params)
        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles(particles, dt/2._8)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, MPI_COMM_SPACE, clearresults=.true.)
          call kick_boris_patacchini(particles, dt, no_tan_trans=.true.)
          call push_particles(particles, dt/2._8)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

    case (WM_MATRIX_VERLET)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt, &
                 params => pepcboris_nml%setup_params)
        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_matrix_verlet(particles, dt)
          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate

    case (WM_ANALYTIC)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt, &
                 params => pepcboris_nml%setup_params)
        block
          complex*16 :: u, udot, Rp, Rm
          real*8 :: Omegap, Omegam, Rscrm, Rscrp, Iscrm, Iscrp, Omegasq
          complex*16, parameter :: ic = (0._8,1._8)
          real*8, parameter :: sqrttwo = sqrt(2._8)

          Omegasq = sqrt((params(PARAMS_OMEGAB)**2)/4._8 - params(PARAMS_OMEGAE)**2)
          Omegap  = params(PARAMS_OMEGAB)/2._8 + Omegasq
          Omegam  = params(PARAMS_OMEGAB)/2._8 - Omegasq

          Rscrm = (Omegap*params(PARAMS_X0) + params(PARAMS_VY0)) / (Omegap - Omegam)
          Iscrm = (Omegap*params(PARAMS_Y0) - params(PARAMS_VX0)) / (Omegap - Omegam)
          Rscrp = params(PARAMS_X0) - Rscrm
          Iscrp = params(PARAMS_Y0) - Iscrm

          Rp = Rscrp + ic*Iscrp
          Rm = Rscrm + ic*Iscrm

          do step=1,nt
            call print_timestep(step, nt, dt)

            u    =   Rp * exp(-ic*Omegap*step*dt) &
                   + Rm * exp(-ic*Omegam*step*dt)
            udot = -ic*Omegap*Rp * exp(-ic*Omegap*step*dt) &
                   -ic*Omegam*Rm * exp(-ic*Omegam*step*dt)

            particles(1)%x(1) = real(u)
            particles(1)%x(2) = aimag(u)
            particles(1)%x(3) = params(PARAMS_Z0) * cos(sqrttwo*params(PARAMS_OMEGAE)*step*dt) + &
              params(PARAMS_VZ0)/(sqrttwo*params(PARAMS_OMEGAE)) * sin(sqrttwo*params(PARAMS_OMEGAE)*step*dt)

            particles(1)%data%v(1) = real(udot)
            particles(1)%data%v(2) = aimag(udot)
            particles(1)%data%v(3) = -sqrttwo*params(PARAMS_OMEGAE)*params(PARAMS_Z0) * sin(sqrttwo*params(PARAMS_OMEGAE)*step*dt) + &
              params(PARAMS_VZ0) * cos(sqrttwo*params(PARAMS_OMEGAE)*step*dt)

            call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
            call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
          end do
        end block
      end associate

    case default
        DEBUG_ERROR(*,'Invalid working mode:', pepcboris_nml%workingmode)
  end select

  call dump_particles(VTK_STEP_LAST, pepcboris_nml%nt, pepcboris_nml%dt, particles, MPI_COMM_SPACE, do_average=.false.)
  call dump_nfeval(pepcboris_nml%rank_world, MPI_COMM_WORLD)

  call pepc_finalize()
  call MPI_COMM_FREE(MPI_COMM_TIME, mpi_err)
  call MPI_COMM_FREE(MPI_COMM_SPACE, mpi_err)
  call MPI_FINALIZE( mpi_err )

  write(*,*)

  contains
    subroutine print_timestep(step, nt, dt)
      implicit none
      real*8, intent(in) :: dt
      integer, intent(in) :: step, nt

      if(pepcboris_nml%rank_world==0) then
        if (step == 0) write(*,*)
        write(*,'(a1, " == computing step",i12,"/",i0, " == simulation time: ", f 12.4)', advance='no')  char(13), step, nt, step*dt
        if (step == nt) write(*,*)
      end if

    end subroutine

end program pepc

