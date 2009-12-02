!  ===================================================================
!
!                              FIELDS
!
!
!   Calculate fields and potential from coordinates x,y,z:
!
!
!   ** Returns fields Ex, Ey, Ez and potential pot excluding external terms **
!
!
!  ===================================================================


subroutine pepc_fields(np_local,nppm_ori,p_x, p_y, p_z, p_q, p_m, p_w, p_label, &
     p_Ex, p_Ey, p_Ez, p_pot, t_np_mult,t_fetch_mult, &
     mac, theta, eps, force_const, err_f, xl, yl, zl, itime, walk_scheme, choose_sort,weighted, init_mb)

  use treevars
  use utils
  use timings
  implicit none
  include 'mpif.h'

  integer, intent(in) :: np_local, t_fetch_mult,nppm_ori  ! # particles on this CPU
  real, intent(in) :: theta, t_np_mult       ! multipole opening angle
  real, intent(in) :: err_f       ! max tolerated force error (rms)
  real, intent(in) :: force_const       ! scaling factor for fields & potential
  real, intent(in) :: eps         ! potential softening distance
  real, intent(in) :: xl, yl, zl         ! box dimensions
  integer, intent(in) :: itime  ! timestep
  integer, intent(in) :: mac, choose_sort, weighted
  real*8, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc 
!  real*8, intent(in),  dimension(np_local) :: p_vx, p_vy, p_vz  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc 
  real*8, intent(in), dimension(np_local) :: p_q, p_m ! charges, masses
  integer, intent(in), dimension(np_local) :: p_label  ! particle label 
  real*8, intent(out), dimension(np_local) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return

  real*8, dimension(nppm_ori) :: ex_tmp,ey_tmp,ez_tmp,pot_tmp,w_tmp
  real*8, dimension(np_local) :: p_w ! work loads
  integer, intent(in) :: init_mb, walk_scheme
  
  integer :: npnew,npold

  integer :: indxl(nppm_ori),irnkl(nppm_ori)
  integer :: islen(num_pe),irlen(num_pe)
  integer :: fposts(num_pe+1),gposts(num_pe+1)

  integer, parameter :: npassm=100000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps,  max_npass,nshort_list, ipe, k
  real*8 :: ttrav, tfetch ! timing integrals
  integer :: pshortlist(nshortm),nshort(npassm),pstart(npassm) ! work balance arrays
  integer :: hashaddr ! Key address 

  real*8 :: ts1b, ts1e, ts2b, ts2e, ta1b, ta1e, ta2b, ta2e

  integer :: max_local,  timestamp
  integer :: ierr
  integer :: iprot = 50  ! frequency for load balance dump

  real*8 :: phi_coul, ex_coul, ey_coul, ez_coul ! partial forces/pot
  real :: ax_ind, ay_ind, az_ind, bx_ind, by_ind, bz_ind
  real :: Epon_x, Epon_y, Epon_z, Phipon, ex_em, ey_em, ez_em, bx_em, by_em, bz_em
  real :: xd, yd, zd  ! positions relative to centre of laser spot
  real :: load_average, load_integral, total_work, average_work
  integer :: total_parts
  character(30) :: cfile, ccol1, ccol2
  character(4) :: cme
  integer :: key2addr        ! Mapping function to get hash table address from key

!  real*8 :: p_ex_nps(nshortm),p_ey_nps(nshortm),p_ez_nps(nshortm)

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  np_mult = t_np_mult
  fetch_mult = t_fetch_mult

  npp = np_local

  if (force_debug) then
     write (*,'(a7,a50/2i5,4f15.2)') 'PEPC | ','Params: itime, mac, theta, eps, force_const, err:', &
		itime, mac, theta, eps, force_const, err_f
     write (*,'(a7,a20/(i16,4f15.3,i8))') 'PEPC | ','Initial buffers: ',(p_label(i), p_x(i), p_y(i), p_z(i), p_q(i), &
		p_label(i),i=1,npp) 
  endif

 ! Copy particle buffers to tree arrays
  do i=1,npp
     x(i) = p_x(i)
     y(i) = p_y(i)
     z(i) = p_z(i)

     ux(i) = 0.  ! No B-fields for now
     uy(i) = 0.
     uz(i) = 0.
     q(i) = p_q(i)
     m(i) = p_m(i)
     if (p_label(i) <=0) then
        ! Trap bad particle labels
        write (*,*) '*** Error: particle labels must be positive integers (1,2,3,...)! '
        write (*,*) p_label(1:20)
        call MPI_ABORT(ierr)
        stop
     else
        pelabel(i) = p_label(i)     
     endif
     pepid(i) = me
     if (num_pe==1 .or. p_w(i)==0) then
        work(i) = 1.
     else
       work(i) = p_w(i)
     endif
     ax(i) = 0.
     ay(i) = 0.
     az(i) = 0.
  end do

  ta1e = MPI_WTIME()
  t_fields_begin = ta1e-ta1b
  ta1b = MPI_WTIME()

  ! Domain decomposition: allocate particle keys to PEs
  call tree_domains(xl,yl,zl,indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold, choose_sort, weighted)  

  call tree_allocate(theta,init_mb)

  call tree_build      ! Build trees from local particle lists

  call tree_branches   ! Determine and concatenate branch nodes

  call tree_fill       ! Fill in remainder of local tree

  call tree_properties ! Compute multipole moments for local tree

  ta1e = MPI_WTIME()
  t_fields_tree = ta1e-ta1b
  ta1b = MPI_WTIME()

  t_walk=0.
  t_walkc=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  work_local = 0  ! total workload
  maxtraverse=0   ! max # traversals
  sum_fetches=0      ! total # multipole fetches/iteration
  sum_ships=0      ! total # multipole shipments/iteration

  nfetch_total=0     ! Zero key fetch counters if prefetch mode off
  nreqs_total=0

  !  # passes needed to process all particles
  nshort_list =nshortm/4 
  npass = max(1,npart/num_pe/nshort_list)   ! ave(npp)/nshort_list   - make nshort_list a power of 2
  load_average = SUM(work(1:npp))/npass   ! Ave. workload per pass: same for all PEs if already load balanced
  nshort(1:npass+1) = 0

  load_integral = 0.
  jpass = 1
  pstart(jpass) = 1
  do i=1,npp
     load_integral = load_integral + work(i)   ! integrate workload

     if (i-pstart(jpass) + 1 == nshortm) then ! Need to check that nshort < nshortm
        write(*,*) 'Warning from PE: ',me,' # parts on pass ',jpass,' in shortlist exceeds array limit ',nshortm
        write(*,*) 'Putting spill-over into following pass'
        nshort(jpass) = nshortm
        jpass = jpass + 1
        pstart(jpass) = i+1

     else if (load_integral >= load_average * jpass .or. i==npp ) then
        nshort(jpass) = i-pstart(jpass) + 1
        jpass = jpass + 1
        pstart(jpass) = i+1

     endif
  end do

  if (me==0 .and. walk_summary) write (*,*) 'LPEPC | Passes:',npass
  if (jpass-1 > npass ) then
     write(*,*) 'Step',itime,', PE',me,' missed some:',nshort(npass+1)
     if (nshort(npass) + nshort(npass+1) <= nshortm) then
        nshort(npass) = nshort(npass) + nshort(npass+1)
     else
        npass = npass+1
     endif
  endif

  if (force_debug)   write (ipefile,*) 'Shortlists: ',(nshort(j),j=1,npass+1)

  max_npass = npass

  ip1 = 1

  ta1e = MPI_WTIME()
  t_fields_nshort = ta1e-ta1b
  ta1b = MPI_WTIME()

  do jpass = 1,max_npass
     !  make short-list
     nps = nshort(jpass)
     ip1 = pstart(jpass)
     pshortlist(1:nps) = (/ (ip1+i-1, i=1,nps) /)

     if (force_debug) then
       	write(*,*) 'pass ',jpass,' of ',max_npass,': # parts ',ip1,' to ',ip1+nps-1
        write(ipefile,*) 'pass ',jpass,' # parts ',ip1,' to ',ip1+nps-1
     endif
     
     !  build interaction list: 
     ! tree walk creates intlist(1:nps), nodelist(1:nps) for particles on short list
     
     if (walk_scheme == 1) then
        call tree_walkc(pshortlist,nps,jpass,theta,itime,mac,ttrav,tfetch)
     else
        call tree_walk(pshortlist,nps,jpass,theta,eps,itime,mac,ttrav,tfetch)
     end if

     t_walk = t_walk + ttrav  ! traversal time (serial)
     t_walkc = t_walkc + tfetch  ! multipole swaps

     ta2b = MPI_WTIME()
     do i = 1, nps

       p = pshortlist(i)    ! local particle index
       
       pot_tmp(p) = 0.
       ex_tmp(p) = 0.
       ey_tmp(p) = 0.
       ez_tmp(p) = 0.

       !  compute Coulomb fields and potential of particle p from its interaction list
       call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul, work(p))

       pot_tmp(p) = pot_tmp(p)+force_const * phi_coul
       ex_tmp(p) = ex_tmp(p)+force_const * ex_coul
       ey_tmp(p) = ey_tmp(p)+force_const * ey_coul
       ez_tmp(p) = ez_tmp(p)+force_const * ez_coul
       w_tmp(p) = work(p)  ! send back work load for next iteration
       work_local = work_local+nterm(i)

    end do

    ta2e= MPI_WTIME()
    t_force = t_force + ta2e-ta2b

    max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list

!     if (dump_tree) call diagnose_tree
     if ((me == 0).and. tree_debug .and. (mod(jpass,max_npass/10+1)==0)) &
          write(*,'(a26,a10,f12.4,a2)') ' LPEPC | TREE WALK (AS) --','Completed',100.0*jpass/max_npass,' %'
  end do

! restore initial particle order specified by calling routine to reassign computed forces
! notice the swapped order of the index-fields -> less changes in restore.f90 compared to tree_domains.f90 

  ta1e = MPI_WTIME()
  t_fields_passes = ta1e-ta1b
  ta1b = MPI_WTIME()
  ts2b = MPI_WTIME()

  call restore(npnew,npold,nppm_ori,irnkl,indxl,irlen,islen,gposts,fposts, &
       pot_tmp(1:npnew),ex_tmp(1:npnew),ey_tmp(1:npnew),ez_tmp(1:npnew),w_tmp(1:npnew),p_pot,p_ex,p_ey,p_ez,p_w)    

  ta1e = MPI_WTIME()
  t_restore_async = ta1e-ta1b
  ts2e = MPI_WTIME()
  t_restore = ts2e-ts2b 
  ta1b = MPI_WTIME()

  if (tree_debug .and. mod(itime,iprot)==0) then
     call MPI_GATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
     call MPI_GATHER(npp, 1, MPI_INTEGER, npps, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  ! Gather particle distn
  end if


  timestamp = itime

  nkeys_total = nleaf+ntwig  

  if (me ==0 .and. mod(itime,iprot)==0 .and. tree_debug) then
     total_work = SUM(work_loads)
     average_work = total_work/num_pe
     cme = achar(timestamp/1000+48) // achar(mod(timestamp/100,10)+48) &
          // achar(mod(timestamp/10,10)+48) // achar(mod(timestamp,10)+48) 
     cfile="load_"//cme//".dat"
     total_parts=SUM(npps)
     open(60, file=cfile)
!     write(60,'(a/a,i8,2(a,1pe15.6))')  '! Full balancing','Parts: ',total_parts,' Work: ',total_work, &
!          ' Ave. work:',average_work        
!     write(60,'(2i8,f12.3)')  (i-1,npps(i),work_loads(i)/average_work,i=1,num_pe)
     close(60)
  endif

  if (force_debug) then
     write (ipefile,101)
     write (*,101)
     write (ipefile,101) force_const

     do i=1,npp
        write (ipefile,102) pelabel(i), & 
             q(i), m(i), ux(i), p_pot(i), p_ex(i)
        write (*,102) pelabel(i), x(i), & 
             q(i), m(i), ux(i), p_pot(i)
     end do

101  format('Tree forces:'/'   p    q   m   ux   pot  ',f8.2)
102  format(1x,i7,5(1pe14.5))

  endif

  if (tree_debug) call tree_stats(itime)

  ta1e = MPI_WTIME()
  t_fields_stats = ta1e-ta1b
  ta1b = MPI_WTIME()

  call tree_deallocate(nppm_ori)

  ts1e = MPI_WTIME()
  t_deallocate = ts1e-ts1b
  t_all = ts1e-ts1b

  call MPI_REDUCE(t_domains,t0_domains,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_build,t0_build,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_allocate,t0_allocate,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_branches,t0_branches,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_fill,t0_fill,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_properties,t0_properties,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  t0_walk = t_walk
  t0_walkc = t_walkc
  t0_force = t_force
  call MPI_REDUCE(t_restore,t0_restore,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_deallocate,t0_deallocate,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(t_all,t0_all,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)


  write(cfile,'(a,i6.6,a)') "load_", me, ".dat"  
  open(60, file=cfile,STATUS='UNKNOWN', POSITION = 'APPEND')
  write(60,*) itime,' ',work_local,' ',npp
  close(60)   

  write(cfile,'(a,i6.6,a)') "timing_", me, ".dat"  
  open(60, file=cfile,STATUS='UNKNOWN', POSITION = 'APPEND')
  write(60,*) itime,' ',t_fields_begin,' ',&
       t_domains_keys,' ',t_domains_sort,' ',t_domains_ship,' ',t_domains_bound,' ',&
       t_allocate_async,' ',&
       t_build_neigh,' ',t_build_part,' ',t_build_byte,' ',&
       t_branches_find,' ',t_branches_exchange,' ',t_branches_integrate,' ',&
       t_fill_local,' ',t_fill_global,' ',&
       t_props_leafs,' ',t_props_twigs,' ',t_props_branches,' ',t_props_global,' ',&
       t_fields_nshort,' ',t_fields_passes,' ',t_restore_async,' ',t_fields_stats,' ',&
       t_domains_sort_pure,' ',t_walk,' ',t_walkc,' ',t_force
  close(60)   

end subroutine pepc_fields









