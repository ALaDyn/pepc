
! ==============================================================
!
!
!                  PEPC-V
!
!    Parallel Efficient Parallel Coulomb-solver: Vortex particles
!
!  ==============================================================

program pepcv

  use module_pepcfields
  use physvars
  use manipulate_particles
  use timings
  use files
  use diagnostics
  implicit none
  include 'mpif.h'

  integer :: ierr, provided
  real :: trun                     ! total run time including restarts and offset
  integer :: itime, stage
  type(t_calc_force_params) ::cf_par
  integer, parameter :: MPI_THREAD_LEVEL = MPI_THREAD_FUNNELED ! "The process may be multi-threaded, but the application
                                                                  !  must ensure that only the main thread makes MPI calls."
  ! Initialize the MPI system (thread safe version, will fallback automatically if thread safety cannot be guaranteed)
  call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! inform the user about possible issues concerning MPI thread safety
  if ((my_rank == 0) .and. (provided < MPI_THREAD_LEVEL)) then
    write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
         MPI_THREAD_LEVEL, provided
    write(*,*) "Initializing with provided level of multithreading. Stability is possibly not guaranteed."
  end if

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

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
  call special_start()

  ! initialize calc force params
  cf_par%theta       = theta
  cf_par%mac         = mac
  cf_par%eps         = eps
  cf_par%force_const = force_const
  cf_par%force_law   = 2

  itime = 0
  trun = ts
  ! Loop over all timesteps
  do while (itime < nt)

     if (my_rank==0 ) write(*,'(a5,i8,a3,i8,a7,i8,a)') 'Step',itime,' of',nt,', using',n,' particles -------------'
     
     ! Runge-Kutta time integration
     do stage = 1,rk_stages

        if (my_rank == 0) write(*,*) 'Time:',trun,'/',te

        call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
        call timer_start(t_tot)

        call pepc_fields(np, n, vortex_particles, np_mult, cf_par, itime, weighted, curve_type, &
                         1, [0, 0, 0], .false., .false.)
        call push_rk2(stage)

        call timer_stop(t_tot)   ! total loop time without diags

        call timings_LocalOutput(itime,stage)
        call timings_GatherAndOutput(itime,stage)

        flush(6)

        trun  = trun  + dt/rk_stages

     end do

     call linear_diagnostics(itime,trun)

     itime = itime + 1

  end do

  ! cleanup of lpepc static data
  call libpepc_finalize()

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  ! End the MPI run
  call MPI_FINALIZE(ierr)

end program pepcv
