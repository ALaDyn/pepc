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


! ==============================================================
!
!
!                  PEPC-MW
!
!    Parallel Efficient Parallel Coulomb-solver
!
!   Driver code for Coulomb-solver library lpepc
!
!  ==============================================================

program pepc

  use module_prepare
  use module_namelist
  use module_pepc_types
  use physvars
  use module_timings
  use module_pepc
  use module_mirror_boxes, only : do_periodic, constrain_periodic
  use module_laser
  use module_pusher
  use module_io
  use module_diagnostics
  use module_workflow
  use module_units
  use module_param_dump
  use module_energy
  use module_particle_setup
  use module_debug
  use module_interaction_specific, only : kelbg_invsqrttemp
  implicit none
  include 'mpif.h'

  logical :: para_file_available
  character(255) :: para_file_name

  integer :: ifile

! particle arrays
  type(t_particle), allocatable :: particles(:)     ! position
    
  ! Allocate array space for tree
  call pepc_initialize("pepc-collfreq", my_rank, n_cpu, .true., comm=MPI_COMM_PEPC)
  call pepc_read_parameters_from_first_argument(para_file_available, para_file_name)

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(file_stdout,1)
  if (my_rank==0) call stamp(file_pepc_out,1)

  ! Each CPU gets copy of initial data
  if (para_file_available) call read_frontend_parameters_from_file(para_file_name)

  call frontend_prepare()

  ! Set up particles - in case if ispecial ==-1, particles and parameters are automatically read from a file again
  call particle_setup(particles, ispecial)


  if (my_rank == 0) write(*,*) "Starting PEPC-COLLFREQ with",n_cpu," Processors, simulating",size(particles,kind=kind_particle), " Particles on each Processor in",nt,"timesteps..."

  ! initial particle output
  ! no initial checkpoint since this would override the current checkpoint if in resume-mode
  call write_particles(particles, .false.)

  ! time-dependent setup stuff
  call workflow(my_rank, 0, 0._8, dt, WORKFLOW_STEP_PRE)

  call pepc_prepare(idim)

  ! parameter output
  if (my_rank == 0) then
    call PrintPhysicalParameters(file_stdout)
    call PrintPhysicalParameters(file_pepc_out)
  endif

  ! Loop over all timesteps
  do while (trun < tend)
   itime = itime + 1
   trun  = trun  + dt

   if (my_rank==0 ) then
      ifile=6
         write(ifile,'(//"===================================================================")')
         write(ifile,'(//"Timestep  ",i8,3x," total run time = ",f12.3,"  (",f12.3," fs)")') &
                  itime , trun, trun*unit_t0_in_fs
         write(ifile,'(//"Steps left (approx.) ",i8)') int((tend-trun)/dt)
   endif

   ! set temperature in each timestep for kelbg interaction
   kelbg_invsqrttemp = 1._8/sqrt(tempe)

   ! time-dependent setup stuff
   call workflow(my_rank, itime, trun, dt, WORKFLOW_STEP_PRE)

   call timer_start(t_tot)

   ! laser propagation according to beam_config
   call laser_update()
   call PrintLaserParameters()

   call pepc_particleresults_clear(particles)

   if (.not. directforce) then
    call pepc_grow_tree(particles)
   endif

   if (.not. directforce) then
    call pepc_traverse_tree(particles(1:size(particles,kind=kind_particle)))
    call debug_barrier()
    if (dbg(DBG_STATS)) call pepc_statistics(itime)

    call debug_barrier()
    !call pepc_restore_particles(particles)
    call pepc_timber_tree()

   else

    call compute_force_direct(particles)

   endif

   particles(1:size(particles,kind=kind_particle))%results%e(1) = particles(1:size(particles,kind=kind_particle))%results%e(1) * force_const
   particles(1:size(particles,kind=kind_particle))%results%e(2) = particles(1:size(particles,kind=kind_particle))%results%e(2) * force_const
   particles(1:size(particles,kind=kind_particle))%results%e(3) = particles(1:size(particles,kind=kind_particle))%results%e(3) * force_const
   particles(1:size(particles,kind=kind_particle))%results%pot  = particles(1:size(particles,kind=kind_particle))%results%pot  * force_const

   ! add any external forces (laser field etc)
   call force_laser(particles)

   ! periodic particle dump
   call write_particles(particles, .true.)

   ! Velocity and position update - explicit schemes only
   call integrator(particles, integrator_scheme)

   ! periodic systems demand periodic boundary conditions
   if (do_periodic) call constrain_periodic(particles)

   call energies(particles, Ukine, Ukini)

   ! time-dependent diagnostics stuff
   call workflow(my_rank, itime, trun, dt, WORKFLOW_STEP_POST)

   ! timings dump
   call timer_stop(t_tot) ! total loop time without diags

   call timings_LocalOutput(itime)
   call timings_GatherAndOutput(itime)

   call flushfiles()

  end do

  ! final particle dump
  call write_particles(particles, .true.)

  deallocate (particles)
  
  ! Time stamp
  if (my_rank==0) call stamp(file_stdout,2)
  if (my_rank==0) call stamp(file_pepc_out,2)

  ! Tidy up O/P files
  call closefiles

 ! call benchmark_end

  ! cleanup of lpepc static data
  call pepc_finalize(MPI_COMM_PEPC)

end program pepc
