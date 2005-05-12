
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

  use physvars
  use treevars
  use utils
  implicit none

  ! timing stuff
  real :: t0, t_key, t_domain, t_build, t_branch, t_fill, t_props, t_walk, t_walkc, t_en, t_force
  real :: t_push, t_diag, t_start_push, t_prefetch, Tpon, ttot, t_laser
  integer :: tremain ! remaining wall_clock seconds
  integer :: llwrem  ! function to enquire remaining wall_clock time
  integer :: ierr

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, num_pe, ierr)

  ! Time stamp
  if (me==0) call stamp(6,1)
  if (me==0) call stamp(15,1)

!  if (me ==0 .and. vis_on) call flvisit_spk_init() ! Start up VISIT
  if (vis_on) call flvisit_nbody2_init ! Start up VISIT interface to xnbody

  call openfiles       ! Set up O/P files


  call setup           ! Each PE gets copy of initial data
  call setup_arrays    ! Allocate array space
  call param_dump      ! Dump initial data
  call configure       ! Set up particles

  do itime = 1,nt
     trun = trun + dt
     if (beam_config >= 3) tlaser = tlaser + dt  
     write(ipefile,'(//a,i8,a,f10.5)') 'Timestep ',itime,' t=',trun

     if (beam_config==4 .and. tlaser<= 2*tpulse) then
        Tpon = 2*vosc**2*max(0.,sin(3.14*tlaser/2./tpulse)**2) 
        ! * (sin(omega*tlaser))**2
     else if (beam_config==3) then
        Tpon = vosc**2
     else
	Tpon = 0.
     endif
     write(ipefile,*) 'Laser intensity: ',Tpon
     write(ipefile,'(//a,i8,(3x,a20,f8.2)/(3x,a,f8.2,a2,f8.2,a4)/a,f9.3)') 'Timestep ',itime+itime_start &
          ,' total run time = ',trun &
          ,' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
          ,' Laser intensity  = ',Tpon


     !     tremain=llwrem(0)
     if (me==0) then
        do ifile = 6,15,9
           write(ifile,'(//a,i8,(3x,a,f8.2)/(3x,a,f8.2,a2,f8.2,a4)/4(a20,f9.3/))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun &
                ,' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
                ,' intensity= ',Tpon &
                ,' x_crit= ',x_crit &
                ,' spot size= ',sigma & 
                ,' theta =  ',theta_beam 
           !                ,' remaining wall-clock time (s)= ',tremain 
        end do
        if (beam_config==5) then 
           write(6,'(4(a,f8.2/))') 'Laser amplitude =',vosc &
                , 'Pulse length',tpulse &
                , 'Pulse width', sigma &
                , 'Focal position',focus(1)
        endif
     endif

     call cputime(t0)
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call make_domains(xl,yl,zl)    ! Domain decomposition: allocate particle keys to PEs

     call cputime(t_domain)
     call tree_build      ! Build trees from local particle lists
     call cputime(t_build)
     call make_branches   ! Determine and concatenate branch nodes
     call cputime(t_branch)
     call tree_fill       ! Fill in remainder of local tree
     call cputime(t_fill)
     call tree_properties ! Compute multipole moments for local tree
     call cputime(t_props)
!     call diagnose_tree
     if (num_pe>1 .and. prefetch) call tree_prefetch(itime)
     call cputime(t_prefetch)
!     call MPI_FINALIZE(ierr)
!     call closefiles
!     stop
     if (coulomb .or. lenjones) then
        call forces(1,npp,dt,t_walk,t_walkc,t_force)   ! Compute forces with load balancing
     endif
     call cputime(t_start_push)

     ! Integrator
     if (scheme == 6 ) then
        call push_full3v(1,npp,dt)  ! full EM pusher (all E, B components)
     else
        call velocities(1,npp,dt)  ! pure ES, NVT ensembles
     endif

     call push_x(1,npp,dt)  ! update positions


     boundaries: select case(particle_bcs)

     case(2)
	call constrain   ! relective particle bcs for temperature-clamp mode

     case(3)
	call earth_plate  ! special bcs for grounded target end 

     case default
        ! do nothing

     end select boundaries

     call cputime(t_push)

     call laser            ! laser propagation according to beam_config
     call cputime(t_laser)

     if (.not. perf_anal) call diagnostics
      if (me.eq.0 .and. .not. perf_anal) then
        do ifile = 6,15,9
	   write(ifile,'(/a)') 'Tree stats:'
           write(ifile,'(a50,2i8,a3,i8,a1)') 'new npp, npart, (max): ',npp,npart,'(',nppm,')'
           write(ifile,'(a50,2i8)') 'local # leaves, twigs: ',nleaf_me,ntwig_me
           write(ifile,'(a50,3i8,a3,i8,a1)') 'final # leaves, twigs, keys, (size_tree): ',nleaf,ntwig,nleaf+ntwig,'(',size_tree,')'
           write(ifile,'(a50,2i8,a3,i8,a1)') 'local, global # branches, (max): ',nbranch,nbranch_sum,'(',nbranch_max,')'
           write (ifile,'(a50,i8,a3,i8)') 'Max length of all interaction lists: ',max_list_length,' / ',nintmax
           write (ifile,'(a50,i8)') 'Max # traversals ',maxtraverse
           write (ifile,'(a50,i8,a3,i8,a1)') 'Max # multipole ships/traversal, (size_fetch):',maxships,'(',size_fetch,')'
           write (ifile,'(a50,i8)') 'Total # multipole ships/iteration ',sumships
           write (ifile,'(a50,i8,a3,i8,a1)') 'Total # multipole ships/prefetch, (numpe*size_fetch): ',sumprefetches,'(',num_pe*size_fetch,')'
        end do
      endif
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
        write(ifile,'(a20,2f12.3,a1)') 'Walk serial: ',t_walk,100*t_walk/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk comm: ',t_walkc,100*t_walkc/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Keys+domains: ',t_domain-t0,100*(t_domain-t0)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Build+props: ',t_props-t_domain,100*(t_props-t_domain)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push-t_start_push,100*(t_push-t_start_push)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag-t_push,100*(t_diag-t_push)/ttot

        write(ifile,'(a20,2f12.3,a1)') 'Total: ',ttot,100.
        write(ifile,'(a20,i4,6f12.3)') 'Timing format: ',num_pe,t_domain-t0,t_props-t_domain,t_prefetch-t_props,t_walk+t_walkc,t_force,ttot

     endif
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  end do



  if (scheme ==5 ) then

     if (ramp) call add_ramp  ! add exponential ramp to front of target
     !  ion eqm mode: add electrons before dumping particle positions
     call add_electrons
     call dump(nt+itime_start)
  endif


  call closefiles      ! Tidy up O/P files

!  if (me ==0 .and. vis_on) call flvisit_spk_close()  ! Tidy up VISIT
  if (me==0 .and. vis_on) call flvisit_nbody2_close ! Tidy up VISIT interface to xnbody

  ! Time stamp
  if (me==0) call stamp(6,2)
  if (me==0) call stamp(15,2)
  ! End the MPI run
  call MPI_FINALIZE(ierr)


end program treemp




