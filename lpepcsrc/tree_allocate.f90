subroutine tree_allocate(theta)

  use treevars
  use timings
  implicit none
  include 'mpif.h'  

  real, intent(in) :: theta

  real*8 :: ts1b=0., ts1e=0., ta1b=0., ta1e=0.

  integer :: k, nintest

  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  nppm=npp
  ! Estimate of interaction list length - Hernquist expression
  if (theta >0.01 ) then
     nintest = int(35.*log(1.*npartm)/(theta**2))
  else
     nintest = npartm
  endif

  nintmax=min(nintest,npartm)   

  !  Space for # table and tree arrays
  !  TODO: need good estimate for max # branches
  size_tree = max(30*nintmax+4*npp,10000)

  if (np_mult>0) then
!     nbaddr = max(log(1.*size_tree)/log(2.) + 1,15.)
!     maxaddress = 2**nbaddr
      nbaddr = int(max(log(1.*size_tree)/log(2.),15.))
      maxaddress = size_tree
  else
     maxaddress = int(abs(np_mult)*10000)
     nbaddr = int(max(log(1.*maxaddress)/log(2.) ,15.))
  endif
  
  if (num_pe > 1) then
     size_fetch = fetch_mult*size_tree/6
     nbranch_max = int(.75*maxaddress)         ! Global max # branches
     nbranch_local_max = 2*nbranch_max/num_pe  ! Local max # branches
  else
     size_fetch=npp
     nbranch_max = 5*nintmax
     nbranch_local_max = 5*nintmax
  endif

  if (num_pe==1) then
    maxleaf = npart 
    maxtwig = maxaddress/2
  else if (num_pe.lt.1024) then
    maxleaf = maxaddress/3
    maxtwig = 2*maxleaf
  else
!  required # twigs increases with P because of branches
    maxleaf = int(maxaddress/(1.+log(1.*num_pe)/3.))
    maxtwig = maxaddress-maxleaf
  endif 

  hashconst = 2**nbaddr-1

  if (me==0 .and. tree_debug) then
!  if (me==0) then
     write(*,'(//a/)') 'Allocating new multipole fields'

    if (dynamic_memalloc) write(*,'(/a/)') 'Dynamic memory management switched on!'	  
    write(*,*) '# procs',num_pe
    write(*,*) 'npart=',npart
    write(*,*) 'N/P=',npart/num_pe
    write(*,*) 'npp= ',npp
    write(*,*) 'nppm= ',nppm
    write(*,*) 'nshortm= ',nshortm
    write(*,*) 'nintest/max=',nintest,nintmax
    write(*,*) 'size_tree= ',size_tree
    write(*,*) 'max address = ',maxaddress
    write(*,*) 'address bits = ',nbaddr
    write(*,*) '# const = ',hashconst
    write(*,*) 'max leaf = ',maxleaf
    write(*,*) 'max twig = ',maxtwig
    write(*,*) 'max branches = ',nbranch_max
    write(*,*) 'size_fetch= ',size_fetch
    write(*,*) 'np_mult= ',np_mult
    write(*,*) 'fetch_mult= ',fetch_mult
    write(*,'(a/)') '... done'
  endif

  allocate ( nterm(nshortm), intlist(nintmax,nshortm), nodelist(nintmax,nshortm) )! interaction key-, node-lists

  allocate ( htable(0:maxaddress), all_addr(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
       treekey(maxaddress), branch_key(nbranch_max), branch_owner(nbranch_max), &
       pebranch(nbranch_max), leaf_key(maxaddress), twig_key(maxaddress), &
       fetched_owner(size_fetch), fetched_keys(size_fetch), requested_owner(size_fetch), requested_keys(size_fetch) )

  all_addr = (/ (k,k=0,maxaddress) /)      ! List of all possible # table addresses
  free_addr = 0
  htable%node = 0
  htable%key = 0
  htable%link = -1
  htable%leaves = 0
  htable%childcode = 0

  ! Allocate memory for tree node properties

  allocate ( first_child(-maxtwig:maxleaf), n_children(-maxtwig:maxleaf), node_level(-maxtwig:maxleaf) )

  allocate ( charge(-maxtwig:maxleaf), &                    ! charge
       abs_charge(-maxtwig:maxleaf), &                ! absolute charge
       xcoc(-maxtwig:maxleaf), ycoc(-maxtwig:maxleaf), zcoc(-maxtwig:maxleaf), &    ! centre of charge 
       xshift(-maxtwig:maxleaf), yshift(-maxtwig:maxleaf), zshift(-maxtwig:maxleaf), &    ! shift vector
       size_node(-maxtwig:maxleaf), & ! rms node size
       xdip(-maxtwig:maxleaf), ydip(-maxtwig:maxleaf), zdip(-maxtwig:maxleaf), &          ! dipole moment
       xxquad(-maxtwig:maxleaf), yyquad(-maxtwig:maxleaf), zzquad(-maxtwig:maxleaf), &       ! quadrupole moment
       xyquad(-maxtwig:maxleaf), yzquad(-maxtwig:maxleaf), zxquad(-maxtwig:maxleaf), &
       jx(-maxtwig:maxleaf), jy(-maxtwig:maxleaf), jz(-maxtwig:maxleaf), &      ! current
       magmx(-maxtwig:maxleaf), magmy(-maxtwig:maxleaf), magmz(-maxtwig:maxleaf)) ! magnetic moment 

  allocate ( pack_child(size_fetch),get_child(size_fetch) )    ! Multipole shipping buffers

  ta1e = MPI_WTIME()
  t_allocate_async = ta1e-ta1b
  ts1e = MPI_WTIME()
  t_allocate = ts1e-ts1b


end subroutine tree_allocate

  
