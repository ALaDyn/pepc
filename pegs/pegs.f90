
! ==============================================================
!
!
!                  PEGS
!
!    Pretty Efficient Gravity Solver
!
!    Parallel Gravity Tree code 
!    Uses tree-code library lpepc to compute gravitational forces
!
!
!  See README.compile for summary of program units
!
!   $ Revision $

!  ==============================================================

program pegs

  use utils
  use physvars

  implicit none
  include 'mpif.h'

  integer :: ierr, ifile
  real :: t0, t_key, t_domain, t_build, t_branch, t_fill, t_props, t_walk, t_en, t_force
  real :: t_push, t_diag, t_start_push, t_prefetch, t_walkc, go,tend, ttot

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

!  if (me ==0 .and. vis_on) call flvisit_spk_init() ! Start up VISIT

  call openfiles       ! Set up O/P files

  call setup           ! Initialise simulation parameters from parts_info.in

! Allocate array space for tree
  call pepc_setup(my_rank,n_cpu,npart_total,theta,debug_tree,np_mult,fetch_mult,debug_rank) 

! Read/setup disc particles
  call setup_particles

! Read/setup stars
  call setup_stars

  call dump_inputs
  call cputime(t0)

  do itime = 1,nt
     trun = trun + dt
     if (my_rank==0 .and. mod(itime,iprot)==0) then
        do ifile=6,15,9
           write(ifile,'(//a,i8,a,f10.5)') 'Timestep ',itime+itime_start,' t=',(itime_start+itime)*dt
        end do
     endif

     ! Compute fields (accelerations) and potential of dust particles
     ! Uses internal particle arrays from library (defined in pepc_setup)
     ! # particles and labelling on CPU may change due to re-sort

     call forces(np_local, walk_scheme, mac, theta, ifreeze, eps, force_tolerance, balance, gamma, bond_const, &
          dt, xl, yl, zl, itime, &
          t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force, iprot,work_tot) 


!     if (my_rank==0) write(*,*) 'Computing star forces ...'
     call stars(dt)
!     if (my_rank==0) write(*,*) '... done'
     call cputime(t_start_push)
     call velocities(1,np_local,dt)
     call push(1,np_local,dt)

     call cputime(t_push)

     call diagnostics
     call cputime(t_diag)
     call cputime(tend)
     t_diag=t_diag-t_push
     t_push=t_push-t_start_push
     ttot=t_domain+t_build+t_walk+t_walkc+t_prefetch+t_force+t_push+t_diag

     if (my_rank==0 .and. mod(itime,iprot)==0) then
        do ifile=6,15,9
           write(ifile,'(//a/)') 'Timing:  Routine   time(s)  percentage'
           write(ifile,'(a20,2f12.3,a1)') 'Domains: ',t_domain,100*(t_domain)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Build: ',t_build,100*t_build/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Walk local: ',t_walk,100*t_walk/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Walk comm: ',t_walkc,100*t_walkc/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Prefetch: ',t_prefetch,100*t_prefetch/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push,100*t_push/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag,100*t_diag/ttot

           write(ifile,'(a20,2f12.4,a1)') 'Total: ',ttot,100.
           write(ifile,'(a20,i4,7f12.3)') 'Timing format: ',n_cpu,t_domain,t_build,t_walk,t_walkc,t_prefetch,t_force,ttot
        end do
     endif
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  end do

  if (mod(nt,idump).ne.0) call dump(nt+itime_start)

  call closefiles      ! Tidy up O/P files

  ! End the MPI run
  call MPI_FINALIZE(ierr)

  stop
end program pegs
