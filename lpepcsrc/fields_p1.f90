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

subroutine pepc_fields_p1(np_local, p_x, p_y, p_z, p_q, p_m, p_w, p_label, &
     p_Ex, p_Ey, p_Ez, p_pot, npnew, indxl, irnkl, islen, irlen, fposts, gposts, &
     mac, theta, eps, force_const, err_f, xl, yl, zl, itime, &
     t_begin,t_domain,t_build,t_branches,t_fill,t_properties,t_prefetch,t_stuff_1,t_stuff_2,t_walk,t_walkc,t_force,t_restore,t_mpi,t_end,t_all)

  use treevars
  use utils
  implicit none
  include 'mpif.h'

  integer, intent(in) :: np_local  ! # particles on this CPU
  real, intent(in) :: theta       ! multipole opening angle
  real, intent(in) :: err_f       ! max tolerated force error (rms)
  real, intent(in) :: force_const       ! scaling factor for fields & potential
  real, intent(in) :: eps         ! potential softening distance
  real, intent(in) :: xl, yl, zl         ! box dimensions
  integer, intent(in) :: itime  ! timestep
  integer, intent(in) :: mac  ! choice of mac
  real*8, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc 
  real*8, intent(in), dimension(np_local) :: p_q, p_m ! charges, masses
  integer, intent(in), dimension(np_local) :: p_label  ! particle label 
  real*8, intent(out), dimension(nppm) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return
  real*8, dimension(nppm) :: p_w ! work loads

  integer :: npnew, npold

  integer, intent(out) :: indxl(nppm),irnkl(nppm)
  integer, intent(out) :: islen(num_pe),irlen(num_pe)
  integer, intent(out) :: fposts(num_pe+1),gposts(num_pe+1)

  integer, parameter :: npassm=10000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps,  max_npass,nshort_list, ipe, k
  real :: t_begin,t_domain, t_build, t_fill, t_branches, t_properties,t_prefetch, t_walk, t_walkc, t_force, t_restore, t_all
  real :: t_stuff_1, t_stuff_2, t_end, t_mpi,ttrav, tfetch, t1, t2, t3, t4  ! timing integrals
  real :: tb1, tb2, td1, td2, tb3, td3, tb4, td4, td5, tb5, tb6, td6, tb7, td7
  integer :: pshortlist(nshortm),nshort(npassm),pstart(npassm) ! work balance arrays
  real*8 :: ex_sl(nshortm),ey_sl(nshortm),ez_sl(nshortm)  ! Fields in particle short list sent to tree_walk
  integer :: hashaddr ! Key address 

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


!!!!!!!!!!!!!!!!!!!!!merge mit trunk

  real :: ww_list
  integer :: a
  

!  force_debug=.true.
!  tree_debug=.false.
!  build_debug=.false.
!  domain_debug = .false.
!  branch_debug=.false.
!  prefetch_debug=.false.
!  walk_debug=.false.
!  walk_summary=.true.
!  dump_tree=.true.

  call cputime(tb1)

  npp = np_local
  npnew = np_local

  if (force_debug) then
     write (*,'(a7,a50/2i5,4f15.2)') 'PEPC | ','Params: itime, mac, theta, eps, force_const, err:',itime, mac, theta, eps, force_const, err_f
     write (*,'(a7,a20/(i16,4f15.3,i8))') 'PEPC | ','Initial buffers: ',(p_label(i), p_x(i), p_y(i), p_z(i), p_q(i), p_label(i),i=1,npp) 
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

  call cputime(tb2)

!  if (me == 0) write(*,*) "tree_domains"
! Domain decomposition: allocate particle keys to PEs
  call tree_domains(xl,yl,zl,indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold)    



! particles now sorted according to keys assigned in tree_domains.
! Serial mode: incoming labels can serve to restore order - now kept in pelabel 
! Parallel mode: have to retain sort vectors from tree_domains

!  if (me == 0) write(*,*) "tree_build"
  call cputime(tb3)
  call tree_build      ! Build trees from local particle lists
!  if (me == 0) write(*,*) "tree_branch"
  call cputime(tb4)
  call tree_branches   ! Determine and concatenate branch nodes
!  if (me == 0) write(*,*) "tree_fill"
  call cputime(tb5)
  call tree_fill       ! Fill in remainder of local tree
!  if (me == 0) write(*,*) "tree_props"
  call cputime(tb6)
  call tree_properties ! Compute multipole moments for local tree
  call cputime(tb7)
!  if (num_pe>1) call tree_prefetch(itime)
  call cputime(td7)

  t_walk=0.
  t_walkc=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  max_list_length = 0
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

  call cputime(td6)

  do jpass = 1,max_npass
     !  make short-list
     nps = nshort(jpass)
     ip1 = pstart(jpass)
     pshortlist(1:nps) = (/ (ip1+i-1, i=1,nps) /)
! copy fields to shortlist
     ex_sl(1:nps) = p_ex(ip1:ip1+nps-1)
     ey_sl(1:nps) = p_ey(ip1:ip1+nps-1)
     ez_sl(1:nps) = p_ez(ip1:ip1+nps-1)

     

     if (force_debug) then
       	write(*,*) 'pass ',jpass,' of ',max_npass,': # parts ',ip1,' to ',ip1+nps-1
        write(ipefile,*) 'pass ',jpass,' # parts ',ip1,' to ',ip1+nps-1
     endif
     
     !  build interaction list: 
     ! tree walk creates intlist(1:nps), nodelist(1:nps) for particles on short list
   
     
     call tree_walk(pshortlist,nps,jpass,theta,eps,itime,mac,ttrav,tfetch,ex_sl,ey_sl,ez_sl,np_local)
  
    ! call tree_walk(pshortlist,nps,jpass,theta,eps,itime,mac,ttrav,tfetch)   
    
     
     t_walk = t_walk + ttrav  ! traversal time (serial)
     t_walkc = t_walkc + tfetch  ! multipole swaps
    
     call cputime(t2)   ! timing
     do i = 1, nps

       p = pshortlist(i)    ! local particle index
       
       !  compute Coulomb fields and potential of particle p from its interaction list
       call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul, work(p))

       p_pot(p) = force_const * phi_coul
       p_ex(p) = force_const * ex_coul
       p_ey(p) = force_const * ey_coul
       p_ez(p) = force_const * ez_coul
       p_w(p) = work(p)  ! send back work load for next iteration
       work_local = work_local+nterm(i)
       
    end do

    call cputime(t3)   ! timing
     t_force = t_force + t3-t2

     max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list


!     do a=1,nps
!        ww_list = ww_list + nterm(a)
        !      write(*,*) nterm(a)
!     end do
     

     if (dump_tree) call diagnose_tree
  
  end do
  
!  ww_list = ww_list / nps
  
!  write(*,*)"Durchschnittliche Laenge der WW-Liste:",ww_list



! restore initial particle order specified by calling routine to reassign computed forces
! notice the swapped order of the index-fields -> less changes in restore.f90 compared to tree_domains.f90 

!  if (me == 0) write(*,*) "tree restore"

  call cputime(td4)

!  call restore(npnew,npold,irnkl,indxl,irlen,islen,gposts,fposts, &
!       pot_tmp(1:npnew),ex_tmp(1:npnew),ey_tmp(1:npnew),ez_tmp(1:npnew),w_tmp(1:npnew),p_pot,p_ex,p_ey,p_ez,p_w)    


  call cputime(td3)


  call MPI_ALLREDUCE(max_local, max_list_length, 1, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
  call MPI_GATHER(npp, 1, MPI_INTEGER, npps, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  ! Gather particle distn

!  timestamp = itime + itime_start
  timestamp = itime

  call cputime(td2)

  if (me ==0 .and. mod(itime,iprot)==0) then
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

  call cputime(td1)

  t_all = td1-tb1
  t_begin = tb2-tb1
  t_domain = tb3-tb2
  t_build = tb4-tb3
  t_branches = tb5-tb4
  t_fill = tb6-tb5
  t_properties = tb7-tb6
  t_prefetch = td7-tb7
  t_stuff_1 = td6-td7
  t_stuff_2 = td4-td6
  t_restore = td3-td4
  t_mpi = td2-td3  
  t_end = td1-td2  


end subroutine pepc_fields_p1








