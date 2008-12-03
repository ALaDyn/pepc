
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
  implicit none
  include 'mpif.h'

  ! timing stuff
  real :: t0, t1, t_key, t_domain=0., t_build=0., t_branches=0., t_fill=0., t_properties=0., t_restore=0.
  real :: t_begin=0.,t_walk=0., t_walkc=0., t_force=0.
  real :: t_stuff_1=0., t_stuff_2=0., t_en, t_mpi=0., t_end=0.
  real :: t_push, t_diag, t_start_push, t_prefetch=0., Tpon, ttot, t_laser, t_all = 0.
  integer :: tremain ! remaining wall_clock seconds
  integer :: ierr, lvisit_active, ifile, debug, i, k, init_mb, nppm_ori
  character*5 :: ci

!  integer, allocatable :: indxl(:),irnkl(:), label_tmp(:)
!  integer, allocatable :: islen(:),irlen(:)
!  integer, allocatable :: fposts(:),gposts(:)

!  integer :: npnew

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)


  call openfiles       ! Set up O/P files

  call setup(init_mb)           ! Each CPU gets copy of initial data

!  debug=2

  call pepc_setup(my_rank,n_cpu,npart_total,theta,db_level,np_mult,fetch_mult,init_mb,nppm_ori)  ! Allocate array space for tree

!  allocate(indxl(nppm),irnkl(nppm),islen(n_cpu),irlen(n_cpu),fposts(n_cpu+1),gposts(n_cpu+1),label_tmp(nppm))

!  call param_dump      ! Dump initial data
  call configure       ! Set up particles

!  ------------------ VISIT ------------------
#ifdef VISIT_NBODY
  if (my_rank ==0 .and. vis_on) then
     write(*,*) "Initialzing VISIT..."
     call flvisit_nbody2_init ! Start up VISIT interface to xnbody
     call flvisit_nbody2_check_connection(lvisit_active)
  end if
#else
#endif
!  ------------------ VISIT ------------------

!  npnew = np_local	

  do itime = 1,nt
     trun = trun + dt

     if (my_rank==0 ) then
        ifile=6
!        do ifile = 6 !,15,9
           write(ifile,'(//a,i8,(3x,a,f12.3))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 

 !       end do

     endif

!  ------------------ VISIT ------------------
#ifdef VISIT_NBODY
     call vis_parts_nbody(itime)
#else
#endif
!  ------------------ VISIT ------------------

     call cputime(t0)

     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

!    do i = 1,np_local
!    	label_tmp(i) = pelabel(i) 
!    end do	

!     call pepc_fields_p1(np_local,x(1:np_local),y(1:np_local),z(1:np_local), &
!                 q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
!		 ex,ey,ez,pot,npnew,indxl,irnkl,islen,irlen,fposts,gposts,np_mult,fetch_mult, &
!	         mac, theta, eps, force_const, err_f, xl, yl, zl, itime, &
!	         t_begin,t_domain,t_build,t_branches,t_fill,t_properties,t_prefetch, &
!	         t_stuff_1,t_stuff_2,t_walk,t_walkc,t_force,t_restore,t_mpi,t_end,t_all)


     call cputime(t1)

!     call restore_p1(npnew,np_local,indxl(1:np_local),irnkl(1:npnew),islen(1:n_cpu),irlen(1:n_cpu),fposts(1:n_cpu+1),gposts(1:n_cpu+1))
    
!     np_local = npnew

     call pepc_fields(np_local,nppm_ori,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult,fetch_mult,mac, theta, eps, force_const, err_f, xl, yl, zl, itime, &
	              t_begin,t_domain,t_build,t_branches,t_fill,t_properties,t_prefetch, &
		      t_stuff_1,t_stuff_2,t_walk,t_walkc,t_force,t_restore,t_mpi,t_end,t_all,init_mb)

! TODO: need proper mac selection instead of beam_config

!     call error_test(npp)
!     call pepc_cleanup(my_rank,n_cpu)

!     stop
	
 
    ! Integrator
     call cputime(t_start_push)

     call velocities(1,np_local,dt)  ! pure ES, NVT ensembles

     call push_x(1,np_local,dt)  ! update positions

     call cputime(t_push)

     call energy_cons(Ukine,Ukini)

     call cputime(t_diag)

!     t_restore = t_start_push - t1    

     if (my_rank==0 .and. db_level .ge.1) then
        ttot = t_push-t0 ! total loop time without diags

!        write(*,*) t_start_push-t0,t_diag-t_start_push,ttot,t_diag-t0
        write(112,*) trun,t_domain,t_build,t_branches,t_fill,t_properties,t_walk,t_walkc,t_force,t_restore,t_all,ttot
        write(*,*) t_all,ttot,ttot-t_all
!        write(112,'(7f12.3)') trun,t_domain,t_build,t_walk+t_walkc,t_force,t_restore,ttot      

!        if (itime ==1 .or. mod(itime,iprot).eq.0) then
!           ifile = 6
!        else 
           ifile = 15
!        endif
!        write(ifile,'(//a/)') 'Timing:  Routine   time(s)  percentage'
!        write(ifile,'(a20,2f12.3,a1)') 'Domains: ',t_domain,100*t_domain/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Build: ',t_build,100*t_build/ttot

!        write(ifile,'(a20,2f12.3,a1)') 'Prefetch: ',t_prefetch,100*t_prefetch/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Walk serial: ',t_walk,100*t_walk/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Walk comm: ',t_walkc,100*t_walkc/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Restore: ',t_restore,100*t_restore/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push-t_start_push,100*(t_push-t_start_push)/ttot
!        write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag-t_push,100*(t_diag-t_push)/ttot

!        write(ifile,'(a20,2f12.3,a1)') 'Total: ',ttot,100.
!        write(ifile,'(a20,i4,7f12.3)') 'Timing format: ',n_cpu,t_domain,t_build,t_prefetch,t_walk+t_walkc,t_force,t_restore,ttot
     endif
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

!    k = 0
!    do i = 1,np_local
!    	if (label_tmp(i) <> pelabel(i)) k = k+1
!    end do	
!    write(*,*) my_rank, k

  end do

  call pepc_cleanup(my_rank,n_cpu)

  call closefiles      ! Tidy up O/P files


!  ------------------ VISIT ------------------
#ifdef VISIT_NBODY
  if (my_rank==0 .and. vis_on) then 
     call flvisit_nbody2_close ! Tidy up VISIT interface to xnbody
  endif
#else
#endif
!  ------------------ VISIT ------------------


  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)
  ! End the MPI run
  call MPI_FINALIZE(ierr)


end program pepce




