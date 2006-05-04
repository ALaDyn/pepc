
! ==============================================================
!
!
!                  PEGS
!
!    Pretty Efficient Gravity Solver
!
!    Parallel Gravity Tree code 
!    Uses tree-code library lpepc to solve for gravitational forces
!
!
!  See README.compile for summary of program units
!
!   $ Revision $

!  ==============================================================

program pegs

  use utils
  use physvars

  real :: t0, t_key, t_domain, t_build, t_branch, t_fill, t_props, t_walk, t_en, t_force
  real :: t_push, t_diag, t_start_push, t_prefetch, Tpon, go,tend

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

!  if (me ==0 .and. vis_on) call flvisit_spk_init() ! Start up VISIT

  call openfiles       ! Set up O/P files


  call setup           ! Each PE gets copy of initial data
  call configure       ! Set up particles

  do itime = 1,nt
     trun = trun + dt
     if (me==0 .and. mod(itime,iprot)==0) then
        do ifile=6,15,9
           write(ifile,'(//a,i8,a,f10.5)') 'Timestep ',itime+itime_start,' t=',(itime_start+itime)*dt
        end do
     endif
     write(ipefile,'(//a,i8,a,f10.5)') 'Timestep ',itime+itime_start,' t=',(itime+itime_start)*dt

     call cputime(t0)
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call tree_domains    ! Domain decomposition: allocate particle keys to PEs

     call cputime(t_domain)
     call tree_build      ! Build trees from local particle lists
     call cputime(t_build)
     call tree_branches   ! Determine and concatenate branch nodes
     call cputime(t_branch)
     call tree_fill       ! Fill in remainder of local tree
     call cputime(t_fill)
     call tree_properties ! Compute multipole moments for local tree
     call cputime(t_props)
     call tree_prefetch
     call cputime(t_prefetch)

     call forces(1,npp,dt,t_walk,t_force)   ! Compute force with balanced shortlists

     call cputime(t_start_push)
     call velocities(1,npp,dt)
     call push(1,npp,dt)
     call stars(dt)

     call cputime(t_push)

     call diagnostics
     call cputime(t_diag)
     call cputime(tend)
     ttot=tend-t0

     if (me==0 .and. mod(itime,iprot)==0) then
        do ifile=6,15,9
           write(ifile,'(//a/)') 'Timing:  Routine   time(s)  percentage'
           write(ifile,'(a20,2f12.3,a1)') 'Domains: ',t_domain-t0,100*(t_domain-t0)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Build: ',t_build-t_domain,100*(t_build-t_domain)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Branches: ',t_branch-t_build,100*(t_branch-t_build)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Fill: ',t_fill-t_branch,100*(t_fill-t_branch)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Props: ',t_props-t_fill,100*(t_props-t_fill)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Prefetch: ',t_prefetch-t_props,100*(t_prefetch-t_props)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Walk: ',t_walk,100*t_walk/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Keys+domains: ',t_domain-t0,100*(t_domain-t0)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Build+props: ',t_prefetch-t_domain,100*(t_prefetch-t_domain)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push-t_start_push,100*(t_push-t_start_push)/ttot
           write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag-t_push,100*(t_diag-t_push)/ttot

           write(ifile,'(a20,2f12.4,a1)') 'Total: ',ttot,100.
           write(ifile,'(a20,i4,5f12.3)') 'Timing format: ',num_pe,t_domain-t0,t_props-t_domain,t_walk,t_force,ttot
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
