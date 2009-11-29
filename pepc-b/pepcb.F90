
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
!                  currently include: pepc-e (serial electrostatic version), pepc-b and pepc-g (gravitation)
!
!  See README.compile for summary of program units
!
!  ==============================================================

program pepcb
  use treevars
  use physvars
  use utils
  implicit none
  include 'mpif.h'

  ! timing stuff
  real*8 :: t0, t_key, t_domain, t_build(4), t_branch, t_fill, t_props, t_walk, t_walkc, t_en, t_force
  real*8 :: t_push, t_diag, t_start_push, t_prefetch, I_laser, ttot, t_laser
  real*8 :: t_loop, t_start_loop, t_end_loop, t_start_prog, t_end_prog
  real*8 :: t_record(10000)
  integer :: tremain ! remaining wall_clock seconds
  integer :: llwrem  ! function to enquire remaining wall_clock time
  integer :: ierr, lvisit_active, ifile, incdf,i, init_mb, nppm_ori
  integer :: vbufcols = 22, irecord=0, ico

!POMP$ INST INIT
!POMP$ INST BEGIN(pepcb)

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  t_start_prog=MPI_WTIME()

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  call setup           ! Read input deck, setup plasma config
  call setup_arrays(init_mb)    ! Set up field arrays
  call openfiles       ! Set up O/P files

! Allocate array space for tree
  call pepc_setup(my_rank,n_cpu,npart_total,theta,debug_tree,np_mult,fetch_mult,nint_max,init_mb,nppm_ori) 

  if (.not.dynamic_memalloc) call tree_allocate(theta,init_mb)

! call closefiles
!  call MPI_FINALIZE(ierr)
! stop 
! ---- Preprocess VISIT setup -----------
 
#ifdef VISIT_NBODY
  if (my_rank ==0 .and. vis_on) then
     call flvisit_nbody2_init ! Start up VISIT interface to xnbody
     call flvisit_nbody2_check_connection(lvisit_active)
     write(*,*) 'Opening visit'
     if (lvisit_active.ne.0) write(*,*) 'Visit connection OK'
  endif
#else
!  ---- No VISIT installed ---------
#endif

#ifdef NETCDFLIB
  if (my_rank ==0 .and. netcdf) then
    write(*,*) 'Opening NETCDF file'
    call ncnbody_open(nbuf_max,vbufcols,ngx, ngy, ngz,ncid,incdf)
  endif
#else
!  ---- No NETCDF installed ---------
#endif


! ---- end of preprocess -------------


  if (launch) then
    call configure(nppm_ori,init_mb)       ! Set up particles

  else
    ico=1
    do while (.not. launch)
	call configure(nppm_ori,init_mb)
#ifdef VISIT_NBODY
	call vis_config
     	if ( mod(ico,ivis) ==0 ) call vis_parts_nbody(ico)       
#endif
     	if ( mod(ico,ivis_fields)==0 ) then
        !     call pot_grid
          call densities
          call sum_fields
#ifdef VISIT_NBODY
          call vis_fields_nbody(0)
#endif
       endif
       ico=ico+1
    end do
  endif
     write (ipefile,'(/a/a/(z21,i6,5f12.4))') 'Particle list after configure:', &
          '  key,   , label  coords, ux, ex', &
          (pekey(i),pelabel(i),x(i),y(i),z(i),ux(i),ex(i),i=1,npp) 


  if (debug_level >0) then
	call diagnostics     ! Initial config

  endif

  t_start_loop=MPI_WTIME()

  if (dynamic_memalloc) call tree_deallocate(nppm_ori)
  
  do itime = 1,nt
     t0=MPI_WTIME()
     trun = trun + dt
     if (my_rank==0 ) then
        do ifile = 6,24,9
           write(ifile,'(/a)') '========================================================================'
        end do
           write(ifile_cpu,'(/a)') '========================================================================'
     endif



     call laser(I_laser)         ! laser propagation according to beam_config
     t_laser=MPI_WTIME()

!     write(ifile_cpu,'(//a,i8,(3x,a20,f8.2)/(3x,a,f8.2,a2,f8.2,a4)/a,f9.3,1pe12.3)') 'Timestep ',itime+itime_start &
!          ,' total run time = ',trun &
!          ,' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
!          ,' Laser amplitude, intensity  = ',sqrt(I_laser),I_laser*1.37e18


     !     tremain=llwrem(0)
     if (my_rank==0 ) then
        do ifile = 6,24,18
           write(ifile,'(//a,i8,(3x,a,f8.2))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 
           if (debug_level >= 2)  then
              write(ifile,'(//(3x,a,f8.2,a2,f8.2,a4)/4(a20,f9.3/))') &
                    ' tlaser = ',tlaser,' (',tlaser*convert_fs,' fs)' &
                   ,' amplitude = ',sqrt(I_laser) &
                   ,' x_crit= ',x_crit &
                   ,' spot size= ',sigma & 
                   ,' theta =  ',theta_beam 
              !                ,' remaining wall-clock time (s)= ',tremain 
           endif
        end do
        if (beam_config==5 .or. beam_config==6) then 
           write(ifile,'(5(a,f8.2/))') 'Laser amplitude =',sqrt(I_laser) &
                , 'Pulse length',tpulse &
                , 'Pulse width', sigma &
                , 'Focal position',focus(1) &
                , 'Elapsed',tlaser 
        endif
    endif

     ! Compute internal E-, B-fields and pot using tree algorithm
     ! Uses internal particle arrays from library (setup up in configure step)
     ! # particles on CPU may change due to re-sort

!POMP$ INST BEGIN(fields)


     call pepc_fields_p(np_local, nppm_ori, walk_scheme, mac, theta, ifreeze, eps, force_tolerance, balance, &
	  force_const, bond_const, &
          dt, xl, yl, zl, itime+itime_start, &
          coulomb, bfields, bonds, lenjones, &
          t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force, iprot,work_tot, init_mb) 

  
!POMP$ INST END(fields)


     t_start_push=MPI_WTIME()

!  Compute external fields (laser, stationary E,B fields)
     fpon_max=0.
     call force_laser(1,np_local)


!POMP$ INST BEGIN(integ)

!  Velocity and position update - explicit schemes only
     call integrator


!POMP$ INST END(integ)

     t_push=MPI_WTIME()

!POMP$ INST BEGIN(diagno)

!  Particle and field diagnostics; visualisation; tree stats

     if (debug_level>=1) then
        call diagnostics
     endif

!POMP$ INST END(diagno)

     t_diag=MPI_WTIME()

     t_diag = t_diag - t_push
     t_laser = t_laser - t0
     t_push = t_push - t_start_push
     ttot = t_domain+SUM(t_build)+t_prefetch+t_walkc+t_walk+t_force + t_push + t_laser
     ops_per_sec = work_tot/ttot

     if (my_rank==0 .and. debug_level .ge.0 .and. mod(itime,iprot).eq.0) then
        irecord = irecord+1
!        do ifile = 6,15,9
        ifile=15
        write(ifile,'(//a/)') 'Timing:  Routine   time(s)  percentage'
        write(ifile,'(a20,2f12.3,a1)') 'Domains: ',t_domain,100*t_domain/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Tree-build: ',t_build(1),100*t_build(1)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Branches: ',t_build(2),100*t_build(2)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Fill: ',t_build(3),100*t_build(3)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Props: ',t_build(4),100*t_build(4)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Prefetch: ',t_prefetch,100*t_prefetch/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk local: ',t_walk,100*t_walk/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk comm: ',t_walkc,100*t_walkc/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push,100*t_push/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag,100*t_diag/ttot

        write(ifile,'(a20,2f12.3,a1)') 'Total: ',ttot,100.
        write(ifile,'(a20/a5,7a12/i5,7f12.3)') 'Timing format:', &
             ' #CPU','domains','build',' prefetch','walk-local','walk-comm','force','tot' &
	  ,n_cpu,t_domain,SUM(t_build),t_prefetch,t_walk,t_walkc,t_force,ttot
!	end do
        t_record(irecord) = ttot
     endif

     if (dynamic_memalloc) call tree_deallocate(nppm_ori)

  end do
  
  t_end_loop=MPI_WTIME()
  t_loop=SUM(t_record(1:irecord))


  if (scheme ==5 ) then

     if (ramp) call add_ramp(x_plasma)  ! add exponential ramp to front of target
     !  ion eqm mode: add electrons before dumping particle positions
     call add_electrons
     call dump(nt+itime_start)
  endif

  if (.not.dynamic_memalloc) call tree_deallocate(nppm_ori)

  call closefiles      ! Tidy up O/P files


! ---- Preprocess VISIT setup -----------
 
#ifdef VISIT_NBODY
if (my_rank==0 .and. vis_on) then 
  call flvisit_nbody2_close ! Tidy up VISIT interface to xnbody
  write(*,*) 'Closing visit'
endif
#else

!  ---- No VISIT installed ---------
#endif

#ifdef NETCDFLIB
  if (my_rank==0 .and. netcdf) then
	call ncnbody_close(ncid,incdf)
	write(*,*) 'Closed netcdf file'
  endif
#endif


  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  t_end_prog=MPI_WTIME()

  if (my_rank==0) then
    do ifile=6,24,18	
     write(ifile,'(a20,2f12.3)') 'Loop total, average: ',t_loop, t_loop/irecord
     write(ifile,'(a20,2f12.3)') 'Loop+diags, average: ',t_end_loop-t_start_loop,(t_end_loop-t_start_loop)/nt
     write(ifile,'(a20,i8)') '# timesteps: ',nt
     write (ifile,'(a20,f14.4)') 'Total run time: ',t_end_prog-t_start_prog
    end do
  endif 

  ! End the MPI run
  call MPI_FINALIZE(ierr)

!POMP$ INST END(pepcb)

end program pepcb
