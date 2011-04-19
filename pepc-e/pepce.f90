
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
  use module_laser
  use module_fields
  use module_acf
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile, nppm_ori, provided
  integer, parameter :: MPI_THREAD_LEVEL = MPI_THREAD_FUNNELED ! "The process may be multi-threaded, but the application
                                                                  !  must ensure that only the main thread makes MPI calls."
  type(acf) :: momentum_acf
  real*8 :: mom(4)

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

  call OutputMemUsage(1, "[pepce startup]", (my_rank==0), 59)

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
  if( idump .gt. 0 ) then
    call write_particles(0)
    if ((ispecial==9).or.(ispecial==10).or.(ispecial==11)) call sum_radial(itime)
  end if

  call benchmark_inner

  flush(6)

  if (experiment) call momentum_acf%initialize(nt)

  ! Loop over all timesteps
  do itime = 1,nt
     trun = trun + dt

     if (my_rank==0 ) then
        ifile=6
           write(ifile,'(//a,i8,(3x,a,f12.3))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 
     endif
     
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call timer_start(t_tot)

     ! laser propagation according to beam_config
     if (experiment) call laser()

     call OutputMemUsage(2, "[pepce before fields]", (db_level==7) .and. (my_rank==0), 59)

     call pepc_fields(np_local,nppm_ori,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult, mac, theta, eps, force_const, err_f, &
                      itime, choose_sort,weighted, &
                      num_neighbour_boxes, neighbour_boxes)

     ! add any external forces (laser field etc)
     if (experiment) then
       call force_laser(1, np_local)
     end if

     if (itime == nt) then
        call gather_particle_diag()
        if (my_rank == 0) call benchmarking_dump_diagnostics()
     end if

     call OutputMemUsage(7, "[pepce after fields]", (db_level==7) .and. (my_rank==0), 59)

     ! Integrator
     call velocities(1,np_local,dt)  ! pure ES, NVT ensembles
     call push(1,np_local,dt)  ! update positions

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(x(1:np_local),y(1:np_local),z(1:np_local),np_local)

     call energy_cons(Ukine,Ukini)

     ! periodic particle dump
     if ( idump .gt. 0 ) then
       if ( mod(itime, idump ) .eq. 0) then
         call write_particles(itime)
         if ((ispecial==9).or.(ispecial==10).or.(ispecial==11)) call sum_radial(itime)

         if (experiment) call field_dump(itime)
       end if
     endif

     if (experiment) then
       call momentum_dump(itime, trun, mom)
       call momentum_acf%addval(mom(1:3))
     endif

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     flush(6)

  end do

  if (experiment) call momentum_acf%finalize()

  call benchmark_post

  ! final particle dump
  if ( idump .gt. 0 ) then
    if ( mod(nt,idump) .ne. 0 ) call write_particles(nt)
  endif

  ! deallocate array space for tree
  call pepc_cleanup(my_rank,n_cpu)

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  call OutputMemUsage(8, "[pepce end of program]", (db_level==7) .and. (my_rank==0), 59)

  ! Tidy up O/P files
  call closefiles

  call benchmark_end

  ! End the MPI run
  call MPI_FINALIZE(ierr)

end program pepce
