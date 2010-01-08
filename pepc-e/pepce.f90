
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
!
!  ==============================================================

program pepce

  use physvars
  use utils
  use benchmarking
  use timings
  implicit none
  include 'mpif.h'

  ! timing stuff
  real*8 :: t0, t1, t2, t3, t4, ttot

  integer :: ierr, lvisit_active, ifile, debug, i, k, init_mb, nppm_ori
  integer :: p
  character(50) :: cme, cfile
  character*5 :: ci

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  call benchmark_pre

  ! Time stamp
!  if (my_rank==0) call stamp(6,1)
!  if (my_rank==0) call stamp(15,1)

  ! Set up O/P files
  call openfiles

  ! Each CPU gets copy of initial data
  call setup(init_mb)

  ! Allocate array space for tree
  call pepc_setup(my_rank,n_cpu,npart_total,theta,db_level,np_mult,fetch_mult,init_mb,nppm_ori)  

  ! Set up particles
  call configure

  ! initial particle output
  if( idump .gt. 0 ) call write_particles(0)

  call benchmark_inner

  do itime = 1,nt
     trun = trun + dt

     if (my_rank==0 ) then
        ifile=6
           write(ifile,'(//a,i8,(3x,a,f12.3))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 
     endif
     
     ! dump trajectory
     if (my_rank == 0 .and. itime == nt) call dump_trajectory()

     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     t0 = MPI_WTIME()
     t1 = MPI_WTIME()

     call pepc_fields(np_local,nppm_ori,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult,fetch_mult,mac, theta, eps, force_const, err_f, xl, yl, zl, &
                      itime, scheme, choose_sort,weighted,choose_build,init_mb)
      

    ! Integrator
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     t2 = MPI_WTIME()

     call velocities(1,np_local,dt)  ! pure ES, NVT ensembles
     call push_x(1,np_local,dt)  ! update positions

     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     t3 = MPI_WTIME()

     call energy_cons(Ukine,Ukini)

     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     t4 = MPI_WTIME()

     ! periodic particle dump
     if ( idump .gt. 0 ) then
       if ( mod(itime, idump ) .eq. 0) call write_particles(itime)
     endif

     if (my_rank==0) then
        ttot = t3-t0 ! total loop time without diags
        open(112,file = 'timing.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
        if (choose_build == 0) then
           if (itime==1) then
              write(112,*) "# trun,t0_domains,t0_allocate,t0_build,t0_branches,t0_fill,t0_properties,t0_walk,t0_walkc,t0_force,t0_restore,t0_deallocate,t0_all,ttot"
           endif
           write(112,*) trun,t0_domains,t0_allocate,t0_build,t0_branches,t0_fill,t0_properties,t0_walk,t0_walkc,t0_force,t0_restore,t0_deallocate,t0_all,ttot
        else
           if (itime==1) then
              write(112,*) "# trun,t0_domains,t0_allocate,t0_local,t0_exchange,t0_global,t0_walk,t0_walkc,t0_force,t0_restore,t0_deallocate,t0_all,ttot"
           endif
           write(112,*) trun,t0_domains,t0_allocate,t0_local,t0_exchange,t0_global,t0_walk,t0_walkc,t0_force,t0_restore,t0_deallocate,t0_all,ttot
        end if
        close(112)
        write(*,*) "t_all ", t0_all
        write(*,*) "ttot ", ttot
        write(*,*) "ttot-t_all ", ttot-t0_all
     endif

  end do

  call benchmark_post

  ! final particle dump
  if ( idump .gt. 0 ) then
    if ( mod(nt,idump) .ne. 0 ) call write_particles(nt)
  endif

  call pepc_cleanup(my_rank,n_cpu)

  call closefiles      ! Tidy up O/P files


  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  call benchmark_end

  ! End the MPI run
  call MPI_FINALIZE(ierr)

end program pepce
