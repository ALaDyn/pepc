
! ==============================================================
!
!
!                  PEPC-MINI
!
!    Parallel Efficient Parallel Coulomb-solver: Electrostatics 
!
!   $ Revision $
!
!   Driver code for Coulomb-solver library lpepc
!
!  Major changes:
!   June 2005: Development begun
!
!   See README.compile for summary of program units
!  ==============================================================

program pepce

  use treetypes, only: &
       t_calc_force_params

  use physvars, only: &
       particle_results, &
       trun, &
       particles, &
       weighted, &
       nt, &
       npart_total, &
       np_mult, &
       np_local, &
       n_cpu, &
       my_rank, &
       itime, &
       ispecial, &
       dt, &
       db_level, &
       curve_type

  use module_neighbour_test, only: &
       validate_n_nearest_neighbour_list, &
       draw_neighbours
  

  use timings
  use module_mirror_boxes
  use module_pepcfields
  use files
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile, provided, i
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
  call special_start(ispecial)

  ! initialize calc force params
  cf_par%mac         = 1 ! NN MAC
  cf_par%force_law   = 5 ! NN "Interaction"

  ! Loop over all timesteps
  do while (itime < nt)
     itime = itime + 1
     trun  = trun  + dt

     if (my_rank==0 ) then
        ifile=6
           write(ifile,'(//a,i8,(3x,a,f12.3))') &
                ' Timestep ',itime &
                ,' total run time = ',trun 
     endif
     
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call timer_start(t_tot)

    call pepc_fields(np_local, npart_total, particles, particle_results, &
        np_mult, cf_par, itime, weighted, curve_type, num_neighbour_boxes, neighbour_boxes, .true., .true.)

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     call draw_neighbours(np_local, npart_total, particles, particle_results, itime)


     do i=1, np_local
       write(37+my_rank,*) i, "|", particle_results(i)%neighbour_nodes(:)
       flush(6)
     end do

     call validate_n_nearest_neighbour_list(np_local, npart_total, particles, particle_results, &
       itime, num_neighbour_boxes, neighbour_boxes)

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

end program pepce
