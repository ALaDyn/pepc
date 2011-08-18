
! ==============================================================
!
!
!                  PEPC-MW
!
!    Parallel Efficient Parallel Coulomb-solver
!
!   $ Revision $
!
!   Driver code for Coulomb-solver library lpepc
!
!  ==============================================================

program pepc

  use treetypes
  use physvars
  use benchmarking
  use timings
  use module_pepcfields
  use module_fmm_framework
  use module_laser
  use module_pusher
  use module_io
  use module_fields
  use module_acf
  use module_diagnostics
  use module_workflow
  use module_units
  use module_setup
  use module_param_dump
  use module_treediags
  use module_vtk
  use module_directsum
  implicit none
  include 'mpif.h'

  integer :: vtk_step

  integer :: ierr, ifile, provided
  integer, parameter :: MPI_THREAD_LEVEL = MPI_THREAD_FUNNELED ! `The process may be multi-threaded, but the application
                                                                  !  must ensure that only the main thread makes MPI calls.`
  type(acf) :: momentum_acf
   real*8 :: mom(4)

  ! Initialization of signal handler - deactivated atm since it outputs the call stack for every mpi rank which results in a very messy output
  !call InitSignalHandler()

  ! Initialize the MPI system (thread safe version, will fallback automatically if thread safety cannot be guaranteed)
  call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, ierr)

  ! prepare a copy of the MPI-communicator
  call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_PEPC, ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_PEPC, my_rank, ierr)

  ! inform the user about possible issues concerning MPI thread safety
  if ((my_rank == 0) .and. (provided < MPI_THREAD_LEVEL)) then
    write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
         MPI_THREAD_LEVEL, provided
    write(*,*) "Initializing with provided level of multithreading. Stability is possibly not guaranteed."
  end if

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_PEPC, n_cpu, ierr)

  call benchmark_pre

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup()

  ! Allocate array space for tree
  call libpepc_setup(my_rank,n_cpu,db_level)

  ! Set up particles
  call particle_setup(ispecial)

  ! parameter output
  if (my_rank == 0) then
    call PrintPhysicalParameters(6)
    call PrintPhysicalParameters(24)
  endif

  ! initial particle output
  ! no initial checkpoint since this would override the current checkpoint if in resume-mode
  call write_particles(.false.)
  if (( idump .gt. 0 ) .and. ((ispecial==9).or.(ispecial==10).or.(ispecial==11))) call sum_radial(itime)

  call momentum_acf%initialize(nt, dt*unit_t0_in_fs, my_rank, n_cpu, MPI_COMM_PEPC)

  call benchmark_inner

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
     if (my_rank == 0 .and. itime == nt) call dump_trajectory()

     call timer_start(t_tot)

     ! laser propagation according to beam_config
     call laser()

     call pepc_fields(np_local,npart_total,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult,mac, theta, calc_force_params(eps, force_const, 3), &
                      itime, weighted, curve_type, &
                      num_neighbour_boxes, neighbour_boxes, treediags)

     call verifydirect(x, y, z, q, ex, ey, ez, pot, np_local, [1, 2, np_local-1, np_local], &
                       calc_force_params(eps, force_const, 3), 1, my_rank, n_cpu, MPI_COMM_PEPC)

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
       call write_spacecurve_to_vtk(itime, trun*unit_t0_in_fs, vtk_step)
     endif

     ! add any external forces (laser field etc)
     call force_laser(1, np_local)

     if (itime == nt) call gather_particle_diag()

     ! Velocity and position update - explicit schemes only
     call integrator(1, np_local, integrator_scheme)

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(x(1:np_local),y(1:np_local),z(1:np_local),np_local)

     call energy_cons(Ukine,Ukini)

     ! periodic particle dump
     call write_particles(.true.)

     if ( idump .gt. 0 ) then
       if ( mod(itime, idump ) .eq. 0) then
         if ((ispecial==9).or.(ispecial==10).or.(ispecial==11)) call sum_radial(itime)

         call field_dump(itime)
       end if
     endif

     ! output total momentum of all negatively charged particles
     call write_total_momentum('momentum_electrons.dat', itime, trun, q(1:np_local) < 0., mom)
     call momentum_acf%addval(mom(1:3))
     call momentum_acf%to_file("momentum_electrons_Kt.dat")

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     call flushfiles()

  end do

  call benchmark_post

  ! cleanup of lpepc static data
  call libpepc_finalize()

  ! final particle dump
  call write_particles(.true.)

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  call benchmark_end

  ! End the MPI run
  call MPI_FINALIZE(ierr)

end program pepc
