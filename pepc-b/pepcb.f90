
! ==============================================================
!
!
!                  PEPC-B
!
!    Pretty Efficient Plasma Coulomb-solver + B (magnetic fields)
!
!    - Magnetoinductive version of PEPC  
!
!   $ Revision $
!
!   Based on Hashed Oct Tree method of Warren & Salmon
!    plus chunks of vectorised f90 SPK code.
!
!  Major changes:
!   Juelich 27 September 2001: Development begun
!   October 2002:  Completed asynchronous, latency-hiding tree traversal
!   November 2002: Incorporated real-time VISIT interface for beam-plasma system (Wolfgang Frings)
!   February 2003: Parallel sort with load-balancing
!   March 2003:    Ported to IBM p690 cluster
!   April 2003:    tree_walk improved by collating multipole info before shipping 
!   July 2003:     Prefetch added to cut barrier-time in tree_walk
!   June 2005:     Separation of tree routines (now in lpepc) from physics applications
!                  currently include: pepc-e (electrostatics), pepc-b and pepc-g (gravitation)
!
!  See README.compile for summary of program units
!
!  ==============================================================

program pepcb

  use physvars
  use utils
  implicit none
  include 'mpif.h'

  ! timing stuff
  real :: t0, t_key, t_domain, t_build, t_branch, t_fill, t_props, t_walk, t_walkc, t_en, t_force
  real :: t_push, t_diag, t_start_push, t_prefetch, Tpon, ttot, t_laser
  real :: t_loop, t_start_loop, t_end_loop, t_start_prog, t_end_prog
  real :: t_record(10000)
  integer :: tremain ! remaining wall_clock seconds
  integer :: llwrem  ! function to enquire remaining wall_clock time
  integer :: ierr, lvisit_active, ifile, incdf
  integer :: vbufcols = 22, irecord=0

!POMP$ INST INIT
!POMP$ INST BEGIN(pepcb)

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  call cputime(t_start_prog)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)


  call openfiles       ! Set up O/P files


  call setup           ! Read input deck, setup plasma config
  call setup_arrays    ! Set up field arrays
  call pepc_setup(my_rank,n_cpu,npart_total,theta,debug_tree,np_mult,fetch_mult)  ! Allocate array space for tree
  call param_dump      ! Dump initial data
  !  if (my_rank ==0 .and. vis_on) call flvisit_spk_init() ! Start up VISIT
  if (my_rank ==0 .and. vis_on) then
!     call flvisit_nbody2_init ! Start up VISIT interface to xnbody
!     call flvisit_nbody2_check_connection(lvisit_active)
     ndom_max = 1000  ! Should be set to nbranch_amx
     nbuf_max = npart_total+ndom_max
!     call ncnbody_open(nbuf_max,vbufcols,ncid,incdf)
  endif

  call configure       ! Set up particles

  call cputime(t_start_loop)

  do itime = 1,nt
     trun = trun + dt
     if (my_rank==0 ) then
        do ifile = 6,15,9
           write(ifile,'(/a)') '========================================================================'
        end do
           write(ifile_cpu,'(/a)') '========================================================================'
     endif


     if (beam_config >= 3) tlaser = tlaser + dt  

     if (beam_config==4 .and. tlaser<= 2*tpulse) then
        Tpon = 2*vosc**2*max(0.,sin(3.14*tlaser/2./tpulse)**2) 
        ! * (sin(omega*tlaser))**2
     else if (beam_config==3) then
        Tpon = vosc**2
     else
	Tpon = 0.
     endif

     write(ifile_cpu,'(//a,i8,(3x,a20,f8.2)/(3x,a,f8.2,a2,f8.2,a4)/a,f9.3)') 'Timestep ',itime+itime_start &
          ,' total run time = ',trun &
          ,' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
          ,' Laser intensity  = ',Tpon


     !     tremain=llwrem(0)
     if (my_rank==0 ) then
        do ifile = 6,15,9
           write(ifile,'(//a,i8,(3x,a,f8.2))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 
           if (debug_level >= 2)  then
              write(ifile,'(//(3x,a,f8.2,a2,f8.2,a4)/4(a20,f9.3/))') &
                    ' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
                   ,' intensity= ',Tpon &
                   ,' x_crit= ',x_crit &
                   ,' spot size= ',sigma & 
                   ,' theta =  ',theta_beam 
              !                ,' remaining wall-clock time (s)= ',tremain 
           endif
        end do
        if (beam_config==5) then 
           write(6,'(4(a,f8.2/))') 'Laser amplitude =',vosc &
                , 'Pulse length',tpulse &
                , 'Pulse width', sigma &
                , 'Focal position',focus(1)
        endif
     endif

     ! Compute E, B-fields and pot using lpepc
     ! Uses internal particle arrays from library (setup up in configure step)
     ! # particles on CPU may change due to resort

!POMP$ INST BEGIN(fields)
     call pepc_fields_p(np_local, mac, theta, ifreeze, eps, force_tolerance, balance, force_const, bond_const, &
          dt, xl, yl, zl, itime, &
          coulomb, bfields, bonds, lenjones, &
          t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force, iprot)   
!POMP$ INST END(fields)
     call cputime(t_start_push)
     call force_laser(1,np_local)


!POMP$ INST BEGIN(integ)
     call integrator
!POMP$ INST END(integ)
     call cputime(t_push)
     call laser            ! laser propagation according to beam_config
     call cputime(t_laser)
!POMP$ INST BEGIN(diagno)
     if (debug_level>=1) then
        call diagnostics
     endif
!POMP$ INST END(diagno)
     call cputime(t_diag)

     t_diag = t_diag - t_laser
     t_laser = t_laser - t_push
     t_push = t_push - t_start_push
     ttot = t_domain+t_build+t_prefetch+t_walkc+t_walk+t_force + t_push + t_laser

     if (my_rank==0 .and. debug_level .ge.0 .and. mod(itime,iprot).eq.0) then
        irecord = irecord+1
        ifile = 6
        write(ifile,'(//a/)') 'Timing:  Routine   time(s)  percentage'
        write(ifile,'(a20,2f12.3,a1)') 'Domains: ',t_domain,100*t_domain/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Build: ',t_build,100*t_build/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Prefetch: ',t_prefetch,100*t_prefetch/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk local: ',t_walk,100*t_walk/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk comm: ',t_walkc,100*t_walkc/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push,100*t_push/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag,100*t_diag/ttot

        write(ifile,'(a20,2f12.3,a1)') 'Total: ',ttot,100.
        write(ifile,'(a20/a5,7a12/i5,7f12.3)') 'Timing format:', &
             ' #CPU','domains','build',' prefetch','walk-local','walk-comm','force','tot' &
	  ,n_cpu,t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force,ttot
        t_record(irecord) = ttot
     endif

  end do
  
  call cputime(t_end_loop)
  t_loop=SUM(t_record(1:irecord))


  if (scheme ==5 ) then

     if (ramp) call add_ramp  ! add exponential ramp to front of target
     !  ion eqm mode: add electrons before dumping particle positions
     call add_electrons
     call dump(nt+itime_start)
  endif


  call closefiles      ! Tidy up O/P files

  !  if (my_rank ==0 .and. vis_on) call flvisit_spk_close()  ! Tidy up VISIT
  if (my_rank==0 .and. vis_on) then 
 !    call flvisit_nbody2_close ! Tidy up VISIT interface to xnbody
 !    call ncnbody_close(ncid,incdf)
  endif

  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  call cputime(t_end_prog)

  if (my_rank==0) then
	
     write(*,'(a20,2f12.3)') 'Loop total, average: ',t_loop, t_loop/irecord
     write(*,'(a20,2f12.3)') 'Loop+diags: ',t_end_loop-t_start_loop
     write (*,'(a20,f14.4)') 'Total run time: ',t_end_prog-t_start_prog
  endif 
  ! End the MPI run
  call MPI_FINALIZE(ierr)
!POMP$ INST END(pepcb)

end program pepcb




