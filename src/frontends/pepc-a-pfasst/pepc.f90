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

program pepc
  ! pepc modules
  use module_pepc
  use module_pepc_types
  use module_timings
  use module_debug
  ! frontend helper routines
  use pepca_helper
  use pepca_integrator
  use pepca_diagnostics
  use pepca_globals

  use pfasst
  use pf_mod_imex
  use pfm_helper
  use pfm_encap
  use pfm_feval
  
  implicit none

  integer :: step
  
  ! variables for pfasst
  integer(kind_default):: MPI_COMM_SPACE, MPI_COMM_TIME, mpi_err
  type(pf_pfasst_t) :: pf
  type(pf_nml_t) :: pf_nml
  type(pf_comm_t) :: tcomm
  type(pf_encap_t), target :: encap
  type(pf_sweeper_t), target :: sweeper
  type(app_params_t), pointer :: level_params(:)
  type(app_data_t), pointer :: y0, yend
  

  !integer(kind_particle), parameter :: numparts = 296568
  integer(kind_particle), parameter :: numparts = 2500
  
  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)

  debug_level = DBG_STATUS

  ! Take care of communication stuff, set up PFASST and PMG
  call pfm_init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)

  ! initialize pepc library and MPI
  call pepc_initialize('pepc-a-pfasst', rank_space, nrank_space, .false., db_level_in=DBG_STATUS, comm=MPI_COMM_SPACE)
  root_space = rank_space.eq.0

  ! Set up PFASST object
  call pf_mpi_create(tcomm, MPI_COMM_TIME)
  call pf_pfasst_create(pf, tcomm, pf_nml%nlevels)

  call pfm_encap_create(encap)
  call pf_verlet_create(sweeper, eval_acceleration)
  call pfm_setup_solver_level_params(level_params, pf_nml%nlevels, numparts) ! numparts is per species, so total number of particles will be 2*numparts
  call pfm_fill_pfasst_object(pf, encap, sweeper, pf_nml, level_params)

  call pf_mpi_setup(tcomm, pf)
  call pf_pfasst_setup(pf)

  !FIXME: Add user-defined calls, e.g. diagnostics, here
  !if (pf_nml%echo_errors) then
  !    call pf_add_hook(pf, pf_nml%nlevels, PF_POST_ITERATION, echo_stats)
  !    call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, gather_stats)
  !    do l = 1, pf_nml%nlevels
  !        !call pf_add_hook(pf, l, PF_POST_SWEEP, echo_residual)
  !    end do
  !end if

  call set_parameter() ! FIXME: this function sets parameters (timestep, etc) which should be used in subsequent calls
  ! Initial conditions
  call feval_init(y0, yend, pf_nml%nlevels, pf%levels(pf_nml%nlevels)%ctx, encap%ctx)
  ! call pf_logger_attach(pf)

  ! Here we go       pfasst-object, initial value, dt, t_end, number of steps, final solution
  call pf_pfasst_run(pf, c_loc(y0), pf_nml%te/pf_nml%nsteps, pf_nml%te, c_loc(yend))

!  ! FIXME: Call PMG's dumping routine
!  call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
!  if (pf_nml%echo_errors) then
!      pf%state%iter = 0
!      if (pf%rank == pf%comm%nproc-1) then
!            call dump_grid(yend,pf%state%iter)
!      end if
!  end if

  ! Remove everything (hopefully)
  call feval_finalize(y0,yend)
  call pf_imex_destroy(sweeper)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(tcomm)
  call pfm_finalize_solver_level_params(level_params, pf_nml%nlevels)
  call pepc_finalize()
  call MPI_Finalize( mpi_err )


  stop 0
!!!!!!!!!!!!!!!! original PEPC-A !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  associate (my_rank => rank_space, root => root_space)

    do step=0, nt - 1
      if(root) then
        write(*,*) " "
        write(*,'(a,i12,"/",i0)')   ' ====== computing step       :', step, nt-1
        write(*,'(a,f12.4)') ' ====== simulation time (fs) :', step*dt*unit_time_fs_per_simunit
      end if

      call timer_start(t_user_step)

      call pepc_particleresults_clear(particles)
      call pepc_grow_tree(particles)
      if(root) write(*,'(a,f12.2," s")') ' ====== tree grow time   :', timer_read(t_fields_tree)
      call pepc_traverse_tree(particles)
      if(root) write(*,'(a,f12.2," s")') ' ====== tree walk time   :', timer_read(t_fields_passes)

      if (dbg(DBG_STATS)) call pepc_statistics(step)
      call pepc_timber_tree()

      if (step > 0) then
        ! first half step to synchronize velocities with positions
        call update_velocities(particles, dt/2.)
      end if

      ! do diagnostics etc here
      call timer_start(t_user_diag)
      if ((particle_output_interval>0) .and. ((mod(step, particle_output_interval)==0) .or. (step==nt-1))) then
        call write_particles_vtk(particles, step, dt*step*unit_time_fs_per_simunit)
        call write_particles_ascii(my_rank, step, particles)
        call gather_and_write_densities(particles, step, dt*step*unit_time_fs_per_simunit)
      endif
      if ((domain_output_interval  >0) .and. ((mod(step, domain_output_interval)  ==0) .or. (step==nt-1))) call write_domain(   particles, step, dt*step*unit_time_fs_per_simunit)
      call diagnose_energy(particles, step, dt*step*unit_time_fs_per_simunit)
      call timer_stop(t_user_diag)
      if(root) write(*,'(a,f12.2," s")') ' ====== diagnostics time :', timer_read(t_user_diag)

      ! second half step: velocities are one half step ahead again
      call update_velocities(particles, dt/2.)
      ! full step for positions: now positions are one half step ahead
      call push_particles(particles, dt)    

      call timer_stop(t_user_step)
      if(root) write(*,'(a,f12.2," s")') ' == time in step: ', timer_read(t_user_step)

      call timings_GatherAndOutput(step, 0, 0 == step)

    end do 

  end associate
 
end program pepc

