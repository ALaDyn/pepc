
! ==============================================================
!
!
!                  PEPC
!
!    Parallel Electrostatic Plasma Tree Code 
!
!   $Revision $
!   Based on Hashed Oct Tree method of Warren & Salmon
!    plus chunks of vectorised f90 SPK code.
!
!  Major changes:
!   Juelich 27 September 2001: Development begun
!   October 2002:  Completed asynchronous, latency-hiding tree traversal
!   November 2002: Incorporated real-time VISIT interface for beam-plasma system (Wolgang Frings)
!   February 2003: Parallel sort with load-balancing
!   March 2003:    Ported to IBM p690 cluster
!   April 2003:    tree_walk improved by collating multipole info before shipping 
!   July 2003:     Prefetch added to cut barrier-time in tree_walk
!
!  See README.compile for summary of program units
!
!  ==============================================================

program treemp

  use treevars
  use utils

  real :: t0, t_key, t_domain, t_build, t_branch, t_fill, t_props, t_walk, t_en, t_force
  real :: t_push, t_diag, t_start_push, t_prefetch, Tpon

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, num_pe, ierr)

  if (me ==0 .and. vis_on) call flvisit_spk_init() ! Start up VISIT

  call openfiles       ! Set up O/P files


  call setup           ! Each PE gets copy of initial data
  call configure       ! Set up particles

  do itime = 1,nt
     trun = trun + dt
     if (beam_config == 4) tlaser = tlaser + dt  
     write(ipefile,'(//a,i8,a,f10.5)') 'Timestep ',itime,' t=',itime*dt

     if (me==0) then
        Tpon = 2*vosc**2*min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
        write(6,*) 'Laser intensity: ',Tpon, ' Amplitude: ',vosc, ' Spot size: ',sigma
        do ifile = 6,15,9
           write(ifile,'(//a,i8,(3x,a,f8.2)/(3x,a,f8.2,a2,f8.2,a4)/a,f9.3)') 'Timestep ',itime+itime_start &
                ,' total run time = ',trun &
                ,' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
                ,' intensity= ',Tpon
           write(ifile,*) 'new npp: ',npp,' new npart: ',npart
           write (ifile,'(a,i8,a3,i8)') 'Max length of all interaction lists: ',max_list_length,' / ',nintmax
        end do
     endif

     call cputime(t0)
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call make_domains    ! Domain decomposition: allocate particle keys to PEs

     call cputime(t_domain)
     call tree_build      ! Build trees from local particle lists
     call cputime(t_build)
     call make_branches   ! Determine and concatenate branch nodes
     call cputime(t_branch)
     call tree_fill       ! Fill in remainder of local tree
     call cputime(t_fill)
     call tree_properties ! Compute multipole moments for local tree
     call cputime(t_props)
     call tree_prefetch
     call cputime(t_prefetch)


     if (walk_balance) then
     	call forces_bal(1,npp,dt,t_walk,t_force)   ! Compute force with balanced shortlists
     else    
     	call forces(1,npp,dt,t_walk,t_force)   ! Compute forces with partial or no balancing
     endif

     call cputime(t_start_push)
     call velocities(1,npp,dt)
     call push(1,npp,dt)
     if ( particle_bcs == 2 ) then 
	call constrain   ! relective particle bcs for temperature-clamp mode
     else if (particle_bcs == 3) then
	call earth_plate  ! special bcs for grounded target end 
     endif

     call cputime(t_push)


     call diagnostics
     call cputime(t_diag)



     if (me==0) then
        ttot = t_push-t0 ! total loop time without diags

        if (itime ==1 .or. mod(itime,iprot).eq.0) then
           ifile = 6
        else 
           ifile = 15
        endif
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

        write(ifile,'(a20,2f12.3,a1)') 'Total: ',ttot,100.
        write(ifile,'(a20,i4,5f12.3)') 'Timing format: ',num_pe,t_domain-t0,t_props-t_domain,t_walk,t_force,ttot

     endif
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  end do

  if (ensemble ==5 ) then
     !  ion eqm mode: add electrons before dumping particle positions
     call add_electrons
     call dump(nt+itime_start)
  endif

  call closefiles      ! Tidy up O/P files

  if (me ==0 .and. vis_on) call flvisit_spk_close()  ! Tidy up VISIT

  ! End the MPI run
  call MPI_FINALIZE(ierr)

  stop
end program treemp
