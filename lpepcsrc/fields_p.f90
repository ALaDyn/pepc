!  ===================================================================
!
!                              FIELDS - internal parallel version
!
!   $Revision$
!
!   Calculate fields and potential from INTERNAL coordinates x,y,z, ux, uy, uz in treevars module:
!
!
!   ** Returns fields Ex, Ey, Ez and potential pot excluding external terms **
!
!
!  ===================================================================


subroutine pepc_fields_p(np_local, nppm_ori, walk_scheme, mac, theta, ifreeze, eps, err_f, balance, force_const, bond_const, &
     delta_t,  xl, yl, zl, itime, &
     coulomb, bfield_on, bonds, lenjones, &
     t_domain,t_build,t_prefetch, t_walk, t_walkc, t_force, iprot,total_work, init_mb)

  use treevars
  use utils
  implicit none
  include 'mpif.h'

  integer :: np_local        ! # particles on CPU - can be changed in tree_domains
  integer, intent(in) :: nppm_ori
  integer :: init_mb  ! max # particles, initial MB
  real, intent(in) :: theta       ! multipole opening angle
  real, intent(in) :: err_f       ! max tolerated force error (rms)
  real, intent(in) :: delta_t       ! timestep 
  real, intent(in) :: force_const, bond_const       ! scaling factors for Coulomb/Len-Jones fields & potential
  real, intent(in) :: eps         ! potential softening distance
  logical, intent(in) :: coulomb  ! Compute Coulomb fields
  logical, intent(in) :: bfield_on  ! Switch for including B-fields
  logical, intent(in) :: bonds  ! Include bond forces
  logical, intent(in) :: lenjones ! Include Lennard-Jones potential
  real, intent(in) :: xl, yl, zl         ! box dimensions
  integer, intent(in) :: itime  ! timestep
  integer, intent(in) :: walk_scheme  ! choice of tree walk 
  integer, intent(in) :: mac  ! choice of tree walk 
  integer, intent(in) :: balance  ! balancing
  integer :: ifreeze



  integer, parameter :: npassm=100000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps,  max_npass,nshort_list, ipe
  real :: t_domain, t_build(4), t_prefetch, t_walk, t_walkc, t_force, ttrav, tfetch, t1, t2, t3  ! timing integrals
  real :: tb1, tb2, tb3, tb4, td1, td2, tp1, tp2
  integer :: pshortlist(nshortm),nshort(npassm),pstart(npassm) ! work balance arrays
  real*8 :: ex_sl(nshortm),ey_sl(nshortm),ez_sl(nshortm)  ! Fields in shortlist from previous timestep
  integer :: hashaddr ! Key address 

  integer :: max_local,  timestamp
  integer :: ierr
  integer :: iprot  ! frequency for load balance dump

  real*8 :: fsx, fsy, fsz, phi, phi_coul, ex_coul, ey_coul, ez_coul
  real*8 :: fljmax
  real*8 :: ax_ind, ay_ind, az_ind, bx_ind, by_ind, bz_ind
  real ::  load_average, load_integral, total_work, average_work
  integer :: total_parts
  character(30) :: cfile, ccol1, ccol2
  character(4) :: cme
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer :: npnew,npold
  integer :: indxl(nppm),irnkl(nppm)      ! Merge-Sort arrays
  integer :: islen(num_pe),irlen(num_pe)
  integer :: fposts(num_pe+1),gposts(num_pe+1)

  !  force_debug=.true.
  !  tree_debug=.true.
  !  build_debug=.false.
  !  domain_debug = .false.
  !  branch_debug=.false.
  !  prefetch_debug=.false.
  !  walk_debug=.true.
  !  walk_summary=.true.
  !  dump_tree=.true.
  !  npp = np_local  ! assumed lists matched for now

   

  if (walk_scheme /= 3) ifreeze=1

  load_balance=balance


  if (force_debug) then
     if (me==0) write (*,*)
     if (me==0) write (*,'(a8,a60/a7,2i5,6f11.2)') 'LPEPC | ','Params itime, walk_scheme, theta, eps, force_const, bond_const, err, delta_t:', &
          'LPEPC | ',itime, walk_scheme, theta, eps, force_const, bond_const, err_f, delta_t
     if (me==0) write (*,'(a8,a17,4l4)') 'LPEPC | ','Force switches: ',coulomb,bfield_on,lenjones,bonds
     write (ipefile,'(a8,a20/(i16,4f15.3))') 'LPEPC | ','Initial buffers: ',(pelabel(i), x(i), y(i), z(i), q(i),i=1,npp) 
  endif

  if (mod(itime-1,ifreeze)==0) then
     if (me==0) write (*,'(a23)') 'LPEPC | REBUILDING TREE'
 !    stop
     call cputime(td1)
     !POMP$ INST BEGIN(domains)
!     call tree_domains(xl,yl,zl)    ! Domain decomposition: allocate particle keys to PEs
     call tree_domains(xl,yl,zl,indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold)    
     !POMP$ INST END(domains)
     ! particles now sorted according to keys assigned in tree_domains.


     if (dynamic_memalloc) call tree_allocate(theta,init_mb)

     call cputime(tb1)

     !POMP$ INST BEGIN(build)
     call tree_build      ! Build trees from local particle lists

     call cputime(tb2)
     call tree_branches   ! Determine and concatenate branch nodes
!call closefiles
!call MPI_FINALIZE(ierr)
!stop
     call cputime(tb3)

     call tree_fill       ! Fill in remainder of local tree
     !POMP$ INST END(build)
     call cputime(tb4)
     t_domain = tb1-td1
     t_build(1) = tb2-tb1
     t_build(2) = tb3-tb2
     t_build(3) = tb4-tb3

  else 
     if (me==0) write (*,'(a19)') 'LPEPC | FREEZE MODE'
     t_domain=0.
     t_build=0.
  endif

     call cputime(tp1)
  !POMP$ INST BEGIN(properties)
  call tree_properties ! Compute multipole moments for local tree
!call closefiles
!call MPI_FINALIZE(ierr)
!stop
  !POMP$ INST END(properties)
     call cputime(tp2)
	t_build(4) = tp2-tp1

  call cputime(tp1)
  if (walk_scheme==3) then
     if (mod(itime-1,ifreeze) /= 0) then
        ! freeze mode - re-fetch nonlocal multipole info
        !POMP$ INST BEGIN(update)
        call tree_update(itime)
        !POMP$ INST END(update)

     else
        ! tree just rebuilt so check for missing nodes before re-fetching
        !POMP$ INST BEGIN(prefetch)
!        call tree_prefetch(itime)
        !POMP$ INST END(prefetch)
        nfetch_total=0     ! Zero key fetch/request counters if fresh tree walk needed
        nreqs_total=0
        sum_fetches=0      ! total current # multipole fetches (per tree update)
        sum_ships=0      ! total current # multipole shipments

     endif

  else if (walk_scheme==2 .and. num_pe>1) then
     !POMP$ INST BEGIN(prefetch)
     call tree_prefetch(itime)
     !POMP$ INST END(prefetch)

  else 
     ! fresh walk in asynch. (walk_scheme=0)  or collective mode (walk_scheme=1)
     nfetch_total=0     ! Zero key fetch/request counters if fresh tree walk needed
     nreqs_total=0
     sum_fetches=0      ! total current # multipole fetches (per tree update)
     sum_ships=0      ! total current # multipole shipments
  endif

  call cputime(tp2)


  t_prefetch = tp2-tp1
  t_walk=0.
  t_walkc=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  max_list_length = 0
  work_local = 0  ! total workload
  maxtraverse=0   ! max # traversals


  !  # passes needed to process all particles
  nshort_list =nshortm/5 
  npass = max(1,npart/num_pe/nshort_list)   ! ave(npp)/nshort_list   - make nshort_list a power of 2
  load_average = SUM(work(1:npp))/npass   ! Ave. workload per pass: same for all PEs if already load balanced
  nshort(1:npass+1) = 0

  load_integral = 0.
  jpass = 1
  pstart(jpass) = 1

! Debugging for load balance
!  if (domain_debug) then
!	write(*,*) 'PE ',me,': npp',npp
!	write(*,'((f12.3))') (work(i),i=1,npp)
!  endif

  do i=1,npp
     load_integral = load_integral + work(i)   ! integrate workload
     
     if (i-pstart(jpass) + 1 == nshortm) then ! Need to check that nshort < nshortm
        write(*,*) 'Warning from PE: ',me,' # parts ',i-pstart(jpass)+1,' on pass ',jpass,' in shortlist exceeds array limit ',nshortm
        write(*,*) 'Load=',load_integral,' Average ',load_average
	write(*,*) 'npp=',npp
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

! TODO:  This bit is buggy - can cause freeze-up
! Something wrong with book-keeping here !

  if (jpass-1 > npass ) then
     write(*,*) 'LPEPC | PE',me,' missed some:',nshort(npass+1)
     write (*,*) 'LPEPC | npass',npass,(nshort(i),i=1,npass)
     if (nshort(npass) + nshort(npass+1) <= nshortm) then
        nshort(npass) = nshort(npass) + nshort(npass+1)
     else
        npass = npass+1
     endif
  endif

  if (force_debug)   write (ipefile,*) 'Shortlists: ',(nshort(j),j=1,npass+1)

  max_npass = npass

  ip1 = 1
  fljmax = 0.

  do jpass = 1,max_npass
     !  make short-list
     nps = nshort(jpass)
     ip1 = pstart(jpass)
     pshortlist(1:nps) = (/ (ip1+i-1, i=1,nps) /)
! Copy fields for tree walk (WS mac=3)
     ex_sl(1:nps) = ex(ip1:ip1+nps-1)
     ey_sl(1:nps) = ey(ip1:ip1+nps-1)
     ez_sl(1:nps) = ez(ip1:ip1+nps-1)
     
     
     if (force_debug) then
       	if (me==0) write(*,'(a14,i4,a4,i4,a20,i8,a4,i8)') &
             'LPEPC |  pass ',jpass,' of ',max_npass,': # particles ',ip1,' to ',ip1+nps-1
        write(ipefile,*) 'LPEPC |  pass ',jpass,' # particles ',ip1,' to ',ip1+nps-1
     endif

     !  build interaction list: 
     ! tree walk creates intlist(1:nps), nodelist(1:nps) for particles on short list

     if (walk_scheme==2 .or. walk_scheme==1) then
   ! collective walk
        call tree_walkc(pshortlist,nps,jpass,theta,itime,mac,ttrav,tfetch)
    else
   ! asynchronous walk  (0,3)
!       call tree_walk(pshortlist,nps,jpass,theta,eps,itime,mac,ttrav,tfetch,ex_sl(1:nps),ey_sl(1:nps),ez_sl(1:nps),np_local)
       call tree_walk(pshortlist,nps,jpass,theta,eps,itime,mac,ttrav,tfetch)
    endif

     t_walk = t_walk + ttrav  ! traversal time (serial)
     t_walkc = t_walkc + tfetch  ! multipole swaps

     !POMP$ INST BEGIN(force)
     call cputime(t2)   ! timing
     do i = 1, nps

        p = pshortlist(i)    ! local particle index

        ! zero field sums
        Ex(p) = 0.
        Ey(p) = 0.
        Ez(p) = 0.
        pot(p) = 0.
        Axo(p) = Ax(p)
        Ayo(p) = Ay(p)
        Azo(p) = Az(p)
        Ax(p) = 0.
        Ay(p) = 0.
        Az(p) = 0.
        Bx(p) = 0.
        By(p) = 0.
        Bz(p) = 0.

        if (coulomb) then
           !  compute Coulomb fields and potential of particle p from its interaction list
           call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul, work(p))
!           call sum_force_split(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul, work(p))

           pot(p) = pot(p) + force_const * phi_coul
           Ex(p) = Ex(p) + force_const * ex_coul
           Ey(p) = Ey(p) + force_const * ey_coul
           Ez(p) = Ez(p) + force_const * ez_coul
        endif

        if (bfield_on) then
           !  compute magnetic fields and vector potential 
           call sum_bfield(p, nterm(i), nodelist( 1:nterm(i),i), eps, bx_ind, by_ind, bz_ind, ax_ind, ay_ind, az_ind )

           Ax(p) =  force_const * ax_ind
           Ay(p) =  force_const * ay_ind
           Az(p) =  force_const * az_ind
           Bx(p) =  force_const * bx_ind
           By(p) =  force_const * by_ind
           Bz(p) =  force_const * bz_ind
	endif

        if (coulomb .and. bfield_on) then
           ! Adjust E-field with inductive term if ES fields on
	   Ex(p) = Ex(p) - (Ax(p)-Axo(p))/delta_t
	   Ey(p) = Ey(p) - (Ay(p)-Ayo(p))/delta_t
   	   Ez(p) = Ez(p) - (Az(p)-Azo(p))/delta_t
        endif

        if (bonds) then
           !  compute short-range forces and potential of particle p from its interaction list
           !          call sum_bond(p, nterm(i), nodelist( 1:nterm(i),i ), fsx, fsy, fsz, phi )

           !          pot(p) = pot(p) + bond_const * phi
           !          fx(p) = fx(p) + bond_const * fsx
           !          fy(p) = fy(p) + bond_const * fsy
           !          fz(p) = fz(p) + bond_const * fsz
        endif

        if (lenjones) then
           !  compute short-range Lennard-Jones forces and potential of particle p from its interaction list
           call sum_lennardjones(p, nterm(i), nodelist( 1:nterm(i),i ), eps, fsx, fsy, fsz, phi )
	   fljmax=max(fljmax,sqrt(fsx**2+fsy**2+fsz**2))
           pot(p) = pot(p) + bond_const * phi
           Ex(p) = Ex(p) + bond_const * fsx
           Ey(p) = Ey(p) + bond_const * fsy
           Ez(p) = Ez(p) + bond_const * fsz
        endif
	
        if (num_pe.gt.1) then
	  work(p) = nterm(i)        ! Should really compute this in sum_force to allow for leaf/twig terms
        else
          work(p)=1
        endif
        work_local = work_local+nterm(i)
     end do

     call cputime(t3)   ! timing
     !POMP$ INST END(force)

     t_force = t_force + t3-t2

     max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TODO: Put all stats in tree_stats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Retain total # keys stored locally for stats and load balancing
  nkeys_total = nleaf+ntwig  

  !  timestamp = itime + itime_start
  timestamp = itime
  if (mod(itime,iprot)==0) then
     call MPI_ALLREDUCE(max_local, max_list_length, 1, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )
     call MPI_GATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
     call MPI_GATHER(npp, 1, MPI_INTEGER, npps, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  ! Gather particle distn
     part_imbal_max = MAXVAL(npps) 
     part_imbal_min = MINVAL(npps)
     part_imbal = (part_imbal_max-part_imbal_min)/1.0/npart*num_pe	
     total_work = SUM(work_loads)
     average_work = total_work/num_pe
     work_imbal_max = MAXVAL(work_loads)/average_work
     work_imbal_min = MINVAL(work_loads)/average_work
     work_imbal = 0.
     do i=1,num_pe
	work_imbal = work_imbal + abs(work_loads(i) - average_work)/average_work/num_pe
     end do
     if (me ==0 ) then

        cme = achar(timestamp/1000+48) // achar(mod(timestamp/100,10)+48) &
             // achar(mod(timestamp/10,10)+48) // achar(mod(timestamp,10)+48) 
        cfile="log/load_"//cme//".dat"
        total_parts=SUM(npps)
        open(60, file=cfile)
        write(60,'(a/a,i8,2(a,1pe15.6))')  '! Full balancing','Parts: ',total_parts,' Work: ',total_work, &
             ' Ave. work:',average_work        
        write(60,'(2i8,f12.3)')  (i-1,npps(i),work_loads(i)/average_work,i=1,num_pe)
        close(60)
     endif
  endif

  if (lenjones) then
	write(*,*) 'Max force: ',fljmax
  endif
  if (force_debug) then
!     if (me==0) write (*,101) force_const
     write (ipefile,101) 'LPEPC | Tree forces:','   p    q   m   ux   pot  ',force_const

     do i=1,npp
        write (ipefile,102) pelabel(i), & 
             q(i), m(i), ux(i), pot(i), ex(i)
        !        if (me==0) write (*,102) pelabel(i), x(i), & 
        !             q(i), m(i), ux(i), pot(i)
     end do

101  format(a/a,f8.2)
102  format(1x,i7,5(1pe14.5))

  endif

  np_local = npp   ! reset local # particles for calling routine

end subroutine pepc_fields_p




