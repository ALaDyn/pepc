
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

  use physvars
  use benchmarking
  use timings
  use module_fmm_framework
  use module_laser
  use module_pusher
  use module_io
  use module_fields
  use module_acf
  use module_diagnostics
  use module_workflow
  use module_units
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile, nppm_ori, provided
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
  call setup()

  ! Allocate array space for tree
  call pepc_setup(my_rank,n_cpu,npart_total,db_level,np_mult,nppm_ori)

  ! Set up particles
  call configure

  ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
  call fmm_framework_init(my_rank, wellsep = 1)

  ! initial particle output
  ! no initial checkpoint since this would override the current checkpoint if in resume-mode
  call write_particles(.false.)
  if (( idump .gt. 0 ) .and. ((ispecial==9).or.(ispecial==10).or.(ispecial==11))) call sum_radial(itime)

  call momentum_acf%initialize(nt)

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

     call pepc_fields(np_local,npart_total,nppm_ori,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult,mac, theta, eps, force_const, &
                      itime, weighted, &
                      num_neighbour_boxes, neighbour_boxes)
      

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
         call momentum_acf%to_file("momentum_Kt.dat")
       end if
     endif

     call write_total_momentum(itime, trun, mom)
     call momentum_acf%addval(mom(1:3))

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     call flushfiles()

  end do

  call benchmark_post

  ! final particle dump
  call write_particles(.true.)

  ! deallocate array space for tree
  call pepc_cleanup(my_rank,n_cpu)

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