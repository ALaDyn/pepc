
! ==============================================================
!
!
!                  PEPC-E
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

  use physvars
  use particle_pusher
  use benchmarking
  use timings
  use module_fmm_framework
  use files
  use energies
  use module_diagnostics
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile, nppm_ori, provided
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

  call benchmark_inner

  flush(6)

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

     call pepc_fields(np_local,npart_total,nppm_ori,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult, mac, theta, eps, force_const, &
                      itime, weighted, &
                      num_neighbour_boxes, neighbour_boxes)

     if (itime == nt) then
        call gather_particle_diag()
        if (my_rank == 0) call benchmarking_dump_diagnostics()
     end if

     ! Integrator
     call velocities(1,np_local,dt)  ! pure ES, NVT ensembles
     call push(1,np_local,dt)  ! update positions

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(x(1:np_local),y(1:np_local),z(1:np_local),np_local)

     call energy_cons(Ukine,Ukini)

     ! periodic particle dump
     call write_particles(.true.)

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     flush(6)

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

end program pepce