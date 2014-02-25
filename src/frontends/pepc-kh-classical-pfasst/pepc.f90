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
   use module_vtk

   use encap
   use pepc_helper
   use time_helper
   use field_helper
   use physics_helper
   use checkpoint_helper
   use module_rng
   use pfm_feval

   use pf_mod_verlet, only: pf_verlet_create, pf_verlet_destroy
   use pfm_helper
   use pfm_encap
   use pfm_hooks
   use pf_mod_comm_mpi, only: pf_mpi_create, pf_mpi_destroy, pf_mpi_setup
   use pf_mod_pfasst, only: pf_pfasst_create, pf_pfasst_setup, pf_pfasst_destroy
   use pf_mod_parallel, only: pf_pfasst_run
   use pf_mod_hooks
   use pf_mod_options
   use pepcboris_paralleldump

   implicit none

  include 'mpif.h'

  ! variables for pfasst
  integer(kind_default):: MPI_COMM_SPACE, MPI_COMM_TIME, mpi_err
  type(pf_pfasst_t) :: pf
  type(pf_nml_t) :: pf_nml
  type(pf_comm_t) :: tcomm
  type(pf_encap_t), target :: enc
  type(pf_sweeper_t), target :: sweeper
  type(level_params_t), pointer :: level_params(:)
  type(app_data_t), pointer :: y0, yend


   type(pepc_pars_t) :: pepc_pars
   type(time_pars_t) :: time_pars
   type(physics_pars_t) :: physics_pars
   type(field_grid_t) :: field_grid
  ! particle data (position, velocity, mass, charge)
   type(t_particle), dimension(:), allocatable, target :: particles

   integer(kind = 4) :: step
   real(kind = 8) :: timer_total, timer_init, timer_step, timer_pcomp, &
      timer_pio, timer_fcomp, timer_fio, timer_dynamics, timer_chkpt

   logical :: do_pdump, do_fdump, do_cdump

   call t_start(timer_total)
   call t_start(timer_init)

  ! Take care of communication stuff
  call pfm_init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)

  ! store communicators and ranks/sizes
   pepc_pars%pepc_comm%comm_space = MPI_COMM_SPACE
   pepc_pars%pepc_comm%comm_time  = MPI_COMM_TIME
   call get_mpi_rank(MPI_COMM_WORLD,           pepc_pars%pepc_comm%rank_world, pepc_pars%pepc_comm%nrank_world)
   call get_mpi_rank(pepc_pars%pepc_comm%comm_space, pepc_pars%pepc_comm%rank_space, pepc_pars%pepc_comm%nrank_space)
   call get_mpi_rank(pepc_pars%pepc_comm%comm_time,  pepc_pars%pepc_comm%rank_time,  pepc_pars%pepc_comm%nrank_time)
   pepc_pars%pepc_comm%root_stdio  = (pepc_pars%pepc_comm%rank_space == 0) .and. (pepc_pars%pepc_comm%rank_time  == pepc_pars%pepc_comm%nrank_time-1)
   pepc_pars%pepc_comm%root_file   = (pepc_pars%pepc_comm%rank_space == 0)

   ! initialize pepc library and MPI
   call pepc_initialize('pepc-kh-classical-pfasst', init_mpi=.false., db_level_in=DBG_STATUS, comm=pepc_pars%pepc_comm%comm_space)
   ! has to be called after pepc_initialize because there MPI_BUFFER_ATTACH is called (and paralleldump needs this buffer)
   call paralleldump_init()
   ! this is not an MPI-parallel application
   if (pepc_pars%pepc_comm%nrank_space > 1 .and. pepc_pars%pepc_comm%rank_world == 0) then
     DEBUG_WARNING(*, "You are requesting spatial parallelism. This is still considered untested. Maybe you wanted to set num_space_instances=num_ranks to only use temporal parallelism?")
   endif
  ! frontend parameter initialization, particle configuration etc.
   call pepc_setup(pepc_pars)
   call setup_time(time_pars, pepc_pars%pepc_comm)
   pf_nml%tend = time_pars%te
   pf_nml%nsteps = time_pars%nsteps
   call rng_init(pepc_pars%pepc_comm%rank_space + 1)
   call setup_physics(physics_pars, time_pars, particles, pepc_pars)
   call setup_field_grid(field_grid, pepc_pars%pepc_comm)

  ! commit all internal pepc variables
  call pepc_prepare(2_kind_dim)
  ! prepare table with level-dependent parameters                          # FIXME: v-- this should be pepc_pars%pepc_comm%rank_time==pepc_pars%pepc_comm%nrank_time-1, but for that rank the hook is not called apparently
  call pfm_setup_solver_level_params(particles, level_params, pf_nml, 2_kind_dim, pepc_pars, physics_pars)
  ! initial potential will be needed for energy computation - using finest level here
  call eval_force(particles, level_params(pf_nml%nlevels), comm=MPI_COMM_SPACE, clearresults=.true.) ! again, use parameters of finest level

  ! initial particle dump
  if (pepc_pars%pepc_comm%root_stdio) then
!TODO    call dump_particles(VTK_STEP_FIRST, 0, 0._8, particles, MPI_COMM_SPACE, do_average=.false.)
!TODO    call dump_energy(0._8, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
  endif

  call t_stop(timer_init)

  if (pepc_pars%pepc_comm%root_stdio) then
    print *, "== [pepc-kh-classical-pfasst]"
    print *, "   running on", pepc_pars%pepc_comm%nrank_space, " MPI ranks."
    print *, "   pdump   = ", pepc_pars%pdump
    print *, "   fdump   = ", pepc_pars%fdump
    print *, "   cdump   = ", pepc_pars%cdump
    print *, ""
    write(*,'(a,es12.4)') " == time in setup (s)                            : ", timer_init
  end if

  do_pdump = .false.
  do_fdump = .false.
  do_cdump = .false.

  select case (pepc_pars%workingmode)
    case (WM_BORIS_SDC, WM_BORIS_MLSDC)
      ! Set up PFASST object
      call pf_mpi_create(tcomm, MPI_COMM_TIME)
      call pfm_encap_init(enc, particles)
      call pf_verlet_create(sweeper, calc_Efield, build_rhs, impl_solver)

      ! if we were requested to perform only single-level SDC, we use finest-level params
      if (pepc_pars%workingmode == WM_BORIS_SDC) then
        level_params(1) = level_params(pf_nml%nlevels)
        pf_nml%nlevels = 1
      endif

      call pf_pfasst_create(pf, tcomm, pf_nml%nlevels)
      call pfm_fill_pfasst_object(pf, enc, sweeper, pf_nml, level_params, pepc_pars%pepc_comm%root_file)

      ! Initial conditions for pfasst
      call feval_init(y0, yend, pf_nml%nlevels, pf%levels(pf_nml%nlevels)%levelctx, enc%encapctx, particles)

      call pf_mpi_setup(tcomm, pf)
      call pf_pfasst_setup(pf)

      ! Add user-defined calls, e.g. diagnostics, here
!TODO      call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, dump_particles_hook)
      !call pf_add_hook(pf, pf_nml%nlevels, PF_PRE_STEP, dump_particles_hook)

      ! some informative output about what we are actually doing
      if (pepc_pars%pepc_comm%rank_world == 0) call pf_print_options(pf)

      ! Here we go       pfasst-object, initial value, dt, t_end, number of steps, final solution
      call pf_pfasst_run(pf, c_loc(y0), pf_nml%tend/pf_nml%nsteps, pf_nml%tend, nsteps=pf_nml%nsteps, qend=c_loc(yend))

      ! Remove everything (hopefully)
      call pf_mpi_destroy(tcomm)
      call feval_finalize(y0,yend)
      call pf_pfasst_destroy(pf)
      call pf_verlet_destroy(sweeper)
      call pfm_finalize_solver_level_params(level_params, pf_nml%nlevels)

    case (WM_BORIS)
      associate (dt => time_pars%dt, &
                 nt => time_pars%nsteps)
        do step=1,nt
          call print_timestep(step, nt, dt)
          call push_particles_velocity_verlet_boris(particles, dt, [0.0_8, 0.0_8, physics_pars%B0])
          call eval_force(particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, clearresults=.true.)
          call update_velocities_velocity_verlet_boris(particles, dt, [0.0_8, 0.0_8, physics_pars%B0])
!TODO          call dump_particles(VTK_STEP_NORMAL, step, dt, particles, MPI_COMM_SPACE, do_average=.false.)
!TODO          call dump_energy(step*dt, particles, level_params(pf_nml%nlevels), MPI_COMM_SPACE, do_average=.false.)
        end do
      end associate
  case (WM_BENEDIKT)
       do step = time_pars%nresume, time_pars%nsteps
        if(pepc_pars%pepc_comm%root_stdio) then
          write(*,*) " "
          write(*,'(a,i12)')    " ====== computing step  :", step
          write(*,'(a,es12.4)') " ====== simulation time :", step * time_pars%dt
          write(*,*) ""
        end if

        call t_start(timer_step)

        if (pepc_pars%pdump .ne. 0) do_pdump = mod(step, pepc_pars%pdump) .eq. 0
        if (pepc_pars%fdump .ne. 0) do_fdump = mod(step, pepc_pars%fdump) .eq. 0
        if (pepc_pars%cdump .ne. 0) do_cdump = mod(step, pepc_pars%cdump) .eq. 0

        if (do_cdump) then
          call t_start(timer_chkpt)
          call write_checkpoint(pepc_pars, time_pars, step, physics_pars, field_grid, particles)
          call t_stop(timer_chkpt)
        end if

        call t_start(timer_pcomp)

        call pepc_particleresults_clear(particles)

        call pepc_grow_tree(particles)
        call pepc_traverse_tree(particles)
        particles(:)%results%e(1) = particles(:)%results%e(1) * force_const
        particles(:)%results%e(2) = particles(:)%results%e(2) * force_const
        particles(:)%results%pot  = particles(:)%results%pot  * force_const

        call apply_external_field(particles)

        call t_stop(timer_pcomp)

        if(do_fdump) then
          call t_start(timer_fcomp)
          call compute_field(pepc_pars, field_grid, particles)
          call t_stop(timer_fcomp)

          call t_start(timer_fio)
          call write_field_on_grid(pepc_pars%pepc_comm, time_pars, step, &
              physics_pars, field_grid)
          call t_stop(timer_fio)
        end if

        if (do_pdump) then
           call t_start(timer_pio)
           call write_domain(time_pars, step, particles)
           call write_particles(pepc_pars, time_pars, step, particles)
           call t_stop(timer_pio)
        end if

        !call pepc_restore_particles(p)
        call pepc_timber_tree()

        call physics_dump(pepc_pars, physics_pars, time_pars, step, particles)

        call t_start(timer_dynamics)
        call push_particles(time_pars, physics_pars, particles)
        call constrain_particles(physics_pars, particles)
        call t_stop(timer_dynamics)

        call timings_GatherAndOutput(step, 0, step == 0)

        call t_stop(timer_step)

        if (pepc_pars%pepc_comm%root_stdio) then
           write(*,'(a,es12.4)') " == time in step (s)                              : ", timer_step
           if (do_cdump) &
              write(*,'(a,es12.4)') " == time in checkpoint I/O (s)                    : ", timer_chkpt
           write(*,'(a,es12.4)') " == time in force computation (s)                 : ", timer_pcomp
           if (do_pdump) &
              write(*,'(a,es12.4)') " == time in particle I/O (s)                      : ", timer_pio
           if (do_fdump) then
              write(*,'(a,es12.4)') " == time in field computation (s)                 : ", timer_fcomp
              write(*,'(a,es12.4)') " == time in field I/O (s)                         : ", timer_fio
           end if
           write(*,'(a,es12.4)') " == time in pusher (s)                            : ", timer_dynamics
        end if

      end do
    case default
        DEBUG_ERROR(*,'Invalid working mode:', pepc_pars%workingmode)
  end select

  deallocate(particles)

  call t_stop(timer_total)

!TODO  if (pepc_pars%pepc_comm%root_stdio) call dump_particles(VTK_STEP_LAST, time_pars%nsteps, time_pars%dt, particles, MPI_COMM_SPACE, do_average=.false.)
  call dump_nfeval(pepc_pars%pepc_comm%rank_world, MPI_COMM_WORLD, pepc_pars%workingmode + IFILE_SUMMAND_NFEVAL)

  if(pepc_pars%pepc_comm%root_stdio) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time without setup (s): ", timer_total - timer_init
    write(*,'(a,es12.4)') " ===== total run time with setup (s):    ", timer_total
  end if

  call paralleldump_cleanup()

  !!! cleanup pepc and MPI
  call pepc_finalize()
  call MPI_COMM_FREE(MPI_COMM_TIME, mpi_err)
  call MPI_COMM_FREE(MPI_COMM_SPACE, mpi_err)
  call MPI_FINALIZE( mpi_err )


contains

  subroutine t_start(t)
    implicit none
    real(kind = 8), intent(out) :: t
    t = -get_time()
  end subroutine t_start

  subroutine t_stop(t)
    implicit none
    real(kind = 8), intent(inout) :: t
    t = t + get_time()
  end subroutine t_stop

  subroutine t_resume(t)
    implicit none
    real(kind = 8), intent(inout) :: t
    t = t - get_time()
  end subroutine t_resume

  subroutine print_timestep(step, nt, dt)
    implicit none
    real*8, intent(in) :: dt
    integer, intent(in) :: step, nt

    if(pepc_pars%pepc_comm%rank_world==0) then
      if (step == 0) write(*,*)
      write(*,'(a1, " == computing step",i12,"/",i0, " == simulation time: ", f 12.4)', advance='no')  char(13), step, nt, step*dt
      if (step == nt) write(*,*)
    end if

  end subroutine
end program pepc

