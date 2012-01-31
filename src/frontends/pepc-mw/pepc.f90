
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
  use module_fields
  use module_acf
  use module_diagnostics
  use module_workflow
  use module_units
  use module_param_dump
  use module_treediags
  use module_vtk
  use module_energy
  use module_particle_setup
  use module_debug
  implicit none
  include 'mpif.h'

  integer :: vtk_step, ierr, num_force_particles
  logical :: para_file_available
  character(255) :: para_file_name

  integer :: ifile
  type(acf) :: momentum_acf

  ! Initialization of signal handler - deactivated atm since it outputs the call stack for every mpi rank which results in a very messy output
  !call InitSignalHandler()

  ! Allocate array space for tree
  call pepc_initialize("pepc-mw", my_rank, n_cpu, .true.)
  call pepc_read_parameters_from_first_argument(para_file_available, para_file_name)

  ! prepare a copy of the MPI-communicator
  call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_PEPC, ierr)

!  call benchmark_pre

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(file_stdout,1)
  if (my_rank==0) call stamp(file_pepc_out,1)

  ! Each CPU gets copy of initial data
  if (para_file_available) call read_frontend_parameters_from_file(para_file_name)

  call pepcmw_prepare()

  if (my_rank == 0) write(*,*) "Starting PEPC-MW with",n_cpu," Processors, simulating",np_local, " Particles on each Processor in",nt,"timesteps..."

  ! Set up particles - in case if ispecial ==-1, particles and parameters are automatically read from a file again
  call particle_setup(ispecial)

  ! parameter output
  if (my_rank == 0) then
    call PrintPhysicalParameters(file_stdout)
    call PrintPhysicalParameters(file_pepc_out)
  endif

  ! initial particle output
  ! no initial checkpoint since this would override the current checkpoint if in resume-mode
  call write_particles(.false.)
  if (( idump .gt. 0 ) .and. ((ispecial==9).or.(ispecial==10).or.(ispecial==11))) call sum_radial(itime)

  ! time-dependent setup stuff
  call workflow(my_rank, 0, 0._8, dt)

  call momentum_acf%initialize(nt-momentum_acf_from_timestep, dt*unit_t0_in_fs, my_rank, n_cpu, MPI_COMM_PEPC)
  if (restart) call momentum_acf%from_file("momentum_electrons_Kt.dat")
!  call benchmark_inner

  call pepc_prepare(idim)

  ! Loop over all timesteps
  do while (itime < nt)
    itime = itime + 1
    trun  = trun  + dt

     if (my_rank==0 ) then
        ifile=6
           write(ifile,'(//a)') "==================================================================="
           write(ifile,'(//a,i8,3x,a,f12.3,"  (",f12.3," fs)")') &
                ' Timestep ',itime &
                ,' total run time = ',trun, trun*unit_t0_in_fs
     endif
     
     ! time-dependent setup stuff
     call workflow(my_rank, itime, trun, dt)

     ! dump trajectory
 !    if (my_rank == 0 .and. itime == nt) call dump_trajectory()

     call timer_start(t_tot)

     ! laser propagation according to beam_config
     call laser_update()
     call PrintLaserParameters()

     call pepc_particleresults_clear(particles, np_local)

     if (.not. directforce) then
       call pepc_grow_tree(np_local, npart_total, particles)
     endif

     ! if necessary, reorder particles here: particles(1:nep) - electrons, particles(nep+1:nep*nip)) - ions
     ! e.g. if we only traverse the tree for electrons but not for ions
     call reorder_particles(np_local, particles, num_force_particles)

     if (.not. directforce) then
       call pepc_traverse_tree(num_force_particles, particles)
       if (dbg(DBG_STATS)) call pepc_statistics(itime)

       !call verifydirect(particles, np_local, [1, 2, np_local-1, np_local], 3, my_rank, n_cpu, MPI_COMM_PEPC)

       ! output of tree diagnostics
       if (treediags) then
         if (itime == 1) then
           vtk_step = VTK_STEP_FIRST
         else if (itime == nt) then
           vtk_step = VTK_STEP_LAST
         else
           vtk_step = VTK_STEP_NORMAL
         endif

         call write_branches_to_vtk(itime,   trun*unit_t0_in_fs, vtk_step)
         call write_spacecurve_to_vtk(itime, trun*unit_t0_in_fs, vtk_step, particles)
       endif

       call pepc_restore_particles(np_local, particles)
       call pepc_timber_tree()

     else

       call compute_force_direct(num_force_particles, particles)

     endif

     particles(1:np_local)%results%e(1) = particles(1:np_local)%results%e(1) * force_const
     particles(1:np_local)%results%e(2) = particles(1:np_local)%results%e(2) * force_const
     particles(1:np_local)%results%e(3) = particles(1:np_local)%results%e(3) * force_const
     particles(1:np_local)%results%pot  = particles(1:np_local)%results%pot  * force_const

     ! add any external forces (laser field etc)
     call force_laser(1, np_local)

!     if (itime == nt) call gather_particle_diag()

     ! Velocity and position update - explicit schemes only
     call integrator(1, np_local, integrator_scheme)

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(particles(1:np_local)%x(1), particles(1:np_local)%x(2), particles(1:np_local)%x(3),np_local)

     call energies(Ukine,Ukini)

     ! periodic particle dump
     call write_particles(.true.)

     if ( idump .gt. 0 ) then
       if ( mod(itime, idump ) .eq. 0) then
         if ((ispecial==9).or.(ispecial==10).or.(ispecial==11)) call sum_radial(itime)

         call field_dump(itime)
       end if
     endif

     !> special diagnostics for [J. Phys. A: Math. Theor 42 (2009), 214048] Th. Raitza et al: Collision frequency of electrons in laser excited small clusters
     if (workflow_setup == 3) call cluster_diagnostics(itime, trun*unit_t0_in_fs, momentum_acf)

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     call flushfiles()

  end do

 ! call benchmark_post

  ! final particle dump
  call write_particles(.true.)

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(file_stdout,2)
  if (my_rank==0) call stamp(file_pepc_out,2)

  ! Tidy up O/P files
  call closefiles

 ! call benchmark_end

  ! cleanup of lpepc static data
  call pepc_finalize()

end program pepc
