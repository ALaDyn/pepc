
! ==============================================================
!
!
!                  PEPC-B
!
!    Pretty Efficient Plasma Coulomb-solver + B (magnetic fields)
!
!    - Laser-plasma version of PEPC  
!
!   $ Revision $
!
!
!  History:
!   Juelich 27 September 2001: Development begun: Hashed Oct Tree method of Warren & Salmon
!    plus chunks of vectorised f90 SPK code.
!   October 2002:  Completed asynchronous, latency-hiding tree traversal
!   November 2002: Incorporated real-time VISIT interface for beam-plasma system (Wolfgang Frings)
!   February 2003: Parallel sort with load-balancing
!   March 2003:    Ported to IBM p690 cluster
!   April 2003:    tree_walk improved by collating multipole info before shipping 
!   July 2003:     Prefetch added to cut barrier-time in tree_walk
!   June 2005:     Separation of tree routines (now in lpepc) from physics applications
!                  currently include: pepc-e (electrostatic version), pepc-b and pepc-g (gravitation)
!   January 2010   Modularization of pepc-b; clean separation from lpepc kernel
!
!
!  ==============================================================

program pepcb

  use module_particle_props
  use module_physvars
  use module_utilities
  use module_geometry
  use module_laser
  use module_field_grid
  use module_io
  use module_diagnostics
  use module_fmm_framework
  use module_pepcfields

  implicit none
  include 'mpif.h'

  ! timing stuff
  real*8 :: t0, t_domain, t_build(4), t_walk, t_walkc, t_force
  real*8 :: t_push, t_diag, t_start_push, t_prefetch, ttot, t_laser
  real*8 :: t_loop, t_start_loop, t_end_loop, t_start_prog, t_end_prog
  real*8 :: t_record(10000)
  real*4 :: I_laser  ! current laser amplitude


  integer :: ierr, ifile, i, init_mb, nppm_ori
  integer :: irecord=0

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
  call setup_particle_arrays(np_alloc)    ! Allocate particle memory

  call openfiles       ! Set up O/P files



  ! Allocate array space for tree
  call libpepc_setup(my_rank, n_cpu, debug_level)

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



    call configure       ! Set up particles


 
    call param_dump

    ! Compute initial field values

    if (my_rank==0) write(*,*) 'Computing initial fields'
   
      
    call pepc_fields(np_local,npart_total,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult, mac, theta, eps, force_const, &
                      itime, weighted, curve_type, &
                      num_neighbour_boxes, neighbour_boxes, .false.)
     
! TODO - shift to lpepcsrc
!     if (np_error>0) then
!        call error_test(np_error)
!        call diagnose_tree   ! Printed tree info (htable etc)
!        call draw_tree2d(xl)     ! Draw PE-trees
!        call draw_lists      ! Draw interaction lists
!        call draw_domains()   ! Domains
!        stop
!     endif

  if (debug_level > 1) then
	  write (ipefile,'(/a/a/(i6,5f12.4))') 'Particle list after configure:', &
          '  key,   , label  coords, ux, ex', &
          (pelabel(i),x(i),y(i),z(i),ux(i),ex(i),i=1,np_local) 

  else if (debug_level >0) then
    call diagnostics     ! Initial config
  endif

  t_start_loop=MPI_WTIME()

  
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

     !     tremain=llwrem(0)
     if (my_rank==0 ) then
        do ifile = 6,24,9
           write(ifile,'(//a,i8,(3x,a,f8.2))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 
        end do
! laser protocol
        if (debug_level >= 2 .and. beam_config.ne.0)  call laser_monitor(I_laser)
    endif

     ! Compute internal E-, B-fields and pot using tree algorithm
     ! Uses internal particle arrays from library (setup up in configure step)
     ! # particles on CPU may change due to re-sort

!POMP$ INST BEGIN(fields)


    call pepc_fields(np_local,npart_total,x(1:np_local),y(1:np_local),z(1:np_local), &
	              q(1:np_local),work(1:np_local),pelabel(1:np_local), &
        	      ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
              	      np_mult, mac, theta, eps, force_const, &
                      itime, weighted, curve_type, &
                      num_neighbour_boxes, neighbour_boxes, .false.)
  
!POMP$ INST END(fields)


     t_start_push=MPI_WTIME()

!  Compute external fields (laser, stationary E,B fields)
     fpon_max=0.
     call force_laser(1,np_local)


!POMP$ INST BEGIN(integ)

!  Velocity and position update - explicit schemes only
     call integrator

     ! periodic systems demand periodic boundary conditions
 !    if (do_periodic) call constrain_periodic(x(1:np_local),y(1:np_local),z(1:np_local),np_local)

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



  end do
  
  t_end_loop=MPI_WTIME()
  t_loop=SUM(t_record(1:irecord))


  if (scheme ==5 ) then
     if (ramp) call add_ramp(x_plasma)  ! add exponential ramp to front of target
     !  ion eqm mode: add electrons before dumping particle positions
     call add_electrons
     call dump(nt+itime_start)
  endif

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