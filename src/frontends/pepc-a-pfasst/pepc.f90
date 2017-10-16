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

program pepc
  ! pepc modules
  use module_pepc
  use module_pepc_kinds
  use module_pepc_types
  use module_timings
  use module_debug
  use module_checkpoint
  ! frontend helper routines
  use pepca_helper
  use pepca_integrator
  use pepca_diagnostics
  use pepca_units

  use pf_mod_verlet, only: pf_verlet_create, pf_verlet_destroy
  use pfm_helper
  use pfm_encap
  use pfm_feval
  use pfm_hooks
  use pf_mod_comm_mpi, only: pf_mpi_create, pf_mpi_destroy, pf_mpi_setup
  use pf_mod_pfasst, only: pf_pfasst_create, pf_pfasst_setup, pf_pfasst_destroy
  use pf_mod_parallel, only: pf_pfasst_run
  use pf_mod_hooks
  use module_mirror_boxes, only: periodicity

  implicit none

  integer :: step, dumpstep
  character(100) :: filename

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
  real*8 :: energies(E_MAXIDX)

  ! Take care of communication stuff
  call pfm_init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)
  ! initialize pepc library and MPI
  pepca_nml%comm=MPI_COMM_SPACE
  call pepc_initialize('pepc-a-pfasst', pepca_nml%rank, pepca_nml%nrank, .false., db_level_in=DBG_STATUS, comm=pepca_nml%comm)
  ! frontend parameter initialization, particle configuration etc.
  call pepca_init(pepca_nml, particles, dt=pf_nml%tend/pf_nml%nsteps/(pf_nml%nnodes(pf_nml%nlevels)-1), nt=pf_nml%nsteps*(pf_nml%nnodes(pf_nml%nlevels)-1))  ! we use the finest (i.e. highest) level here
  ! commit all internal pepc variables
  call pepc_pre(dim)
  ! prepare table with level-dependent parameters
  call pfm_setup_solver_level_params(particles, level_params, pf_nml, dim, pepca_nml%rank, MPI_COMM_SPACE)
  ! initial potential will be needed for energy computation
  call eval_force(particles, level_params(pf_nml%nlevels), step=0, comm=MPI_COMM_SPACE, clearresults=.true.) ! again, use parameters of finest level

  if (pepca_nml%use_pfasst) then
      ! Set up PFASST object
      call pf_mpi_create(tcomm, MPI_COMM_TIME)
      call pf_pfasst_create(pf, tcomm, pf_nml%nlevels)

      call pfm_encap_init(encap, particles)
      call pf_verlet_create(sweeper, eval_acceleration)
      call pfm_fill_pfasst_object(pf, encap, sweeper, pf_nml, level_params)

      ! Initial conditions for pfasst
      call feval_init(y0, yend, pf_nml%nlevels, pf%levels(pf_nml%nlevels)%levelctx, encap%encapctx, particles)

      call pf_mpi_setup(tcomm, pf)
      call pf_pfasst_setup(pf)

      ! Add user-defined calls, e.g. diagnostics, here
      call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, track_energy_hook)
      call pf_add_hook(pf, pf_nml%nlevels, PF_POST_ITERATION, compare_checkpoint_hook)
      call pf_add_hook(pf, pf_nml%nlevels, PF_PRE_STEP, compare_checkpoint_hook)
      if (any(periodicity)) call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, constrain_particles_hook)
      call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, dump_particles_hook)

      !if (pf_nml%echo_errors) then
      !    call pf_add_hook(pf, pf_nml%nlevels, PF_POST_ITERATION, echo_stats)
      !    call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, gather_stats)
      !    do l = 1, pf_nml%nlevels
      !        !call pf_add_hook(pf, l, PF_POST_SWEEP, echo_residual)
      !    end do
      !end if

      ! call pf_logger_attach(pf)

      ! Here we go       pfasst-object, initial value, dt, t_end, number of steps, final solution
      call pf_pfasst_run(pf, c_loc(y0), pf_nml%tend/pf_nml%nsteps, pf_nml%tend, nsteps=pf_nml%nsteps, qend=c_loc(yend))

    !  ! FIXME: Call dumping routine
    !  call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
    !  if (pf_nml%echo_errors) then
    !      pf%state%iter = 0
    !      if (pf%rank == pf%comm%nproc-1) then
    !            call dump_grid(yend,pf%state%iter)
    !      end if
    !  end if

      ! Remove everything (hopefully)
      call pf_mpi_destroy(tcomm)
      call feval_finalize(y0,yend)
      call pf_pfasst_destroy(pf)
      call pf_verlet_destroy(sweeper)
      call pfm_finalize_solver_level_params(level_params, pf_nml%nlevels)
   else ! use PEPC
      associate (Ngrid => pepca_nml%Ngrid, &
                    dt => pepca_nml%dt,    &
                    nt => pepca_nml%nt)
        do step=0,nt
          if(pepca_nml%rank==0) then
            write(*,*) " "
            write(*,'(a,i12,"/",i0)')   ' ====== computing step       :', step, nt
            write(*,'(a,f12.4)')        ' ====== simulation time (fs) :', step*dt*unit_time_fs_per_simunit
          end if

          call eval_force(particles, level_params(pf_nml%nlevels), step, MPI_COMM_SPACE, clearresults=.true.)

          if (step > 0) then
            ! first half step to synchronize velocities with positions
            call update_velocities(particles, dt/2.)
          end if

          ! do diagnostics etc here
          if (dumpnow(pepca_nml%output_interval(OI_PARTICLES_VTK), step, nt)) call write_particles_vtk(particles, step, nt, dt*step*unit_time_fs_per_simunit, MPI_COMM_SPACE)
          if (dumpnow(pepca_nml%output_interval(OI_DOMAIN_VTK   ), step, nt)) call write_domain(particles, step, nt, dt*step*unit_time_fs_per_simunit)

          dumpstep = get_checkpoint_id(dt*step)
          if (dumpnow(pepca_nml%output_interval(OI_VERIFY_PARTICLES), step, nt)) then
            block
              real*8 :: xerr, verr
              call compare_particles_to_checkpoint(particles, dumpstep, MPI_COMM_SPACE, xerr, verr)
              if(pepca_nml%rank==0) write(*,'(" step: ",i5," t=", es10.3, " Ex: ",es14.7," Ev: ",es14.7)') step, step*dt, xerr, verr
            end block
          end if
          if (dumpnow(pepca_nml%output_interval(OI_PARTICLES_MPI), step, nt)) call write_particles_mpiio(MPI_COMM_SPACE, dumpstep, 2*pepca_nml%numparts_total, particles, filename)

          if (dumpnow(pepca_nml%output_interval(OI_PARTICLES_ASC), step, nt)) call write_particles_ascii(pepca_nml%rank, dumpstep, particles, filename)
          if (dumpnow(pepca_nml%output_interval(OI_DENSITIES_VTK), step, nt)) call gather_and_write_densities(particles, Ngrid, dumpstep, nt, dt*step*unit_time_fs_per_simunit, pepca_nml%rank)

          call diagnose_energy(particles, energies, step, dt*step*unit_time_fs_per_simunit, MPI_COMM_SPACE, pepca_nml%rank==0)

          ! second half step: velocities are one half step ahead again
          call update_velocities(particles, dt/2.)
          ! full step for positions: now positions are one half step ahead
          call push_particles(particles, dt)

          call timings_GatherAndOutput(step, 0, 0==step)

        end do
      end associate
  endif

  call pepc_finalize()
  call MPI_COMM_FREE(MPI_COMM_TIME, mpi_err)
  call MPI_COMM_FREE(MPI_COMM_SPACE, mpi_err)
  call MPI_FINALIZE( mpi_err )

end program pepc

