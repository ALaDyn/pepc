subroutine tree_allocate(theta,init_mb)

  use treevars
  implicit none
  
  real, intent(in) :: theta
  integer, intent(in) :: init_mb
  integer :: k, npsize
  integer :: mem_parts, mem_multipoles, mem_fields, mem_tree, mem_prefetch, mem_tot, mem_lists
  real, parameter :: mb=2.**20

  !nppm_old = nppm
  nppm = npp
  
  ! Estimate of interaction list length - Hernquist expression
  if (theta >0 ) then
     nintmax = max(1.*24*log(2.*npartm)/max(theta**2,0.5),2200.)
  else
     nintmax = npartm
  endif
  max_list_length = 0 ! current max length of all interaction lists
  
  
  !  Space for # table and tree arrays
  !  TODO: need good estimate for max # branches
  npsize=nppm
  size_tree = max(4*nintmax+2*npsize,10000)
  if (np_mult>0) then
     nbaddr = max(log(1.*size_tree)/log(2.) + 1,15.)
     maxaddress = 2**nbaddr
  else
     maxaddress = abs(np_mult)*10000
     nbaddr = max(log(1.*maxaddress)/log(2.) ,15.)
  endif
  
  if (num_pe > 1) then
     size_fetch = fetch_mult*size_tree/4
     nbranch_max = .75*maxaddress              ! Global max # branches
     nbranch_local_max = 4*nbranch_max/num_pe  ! Local max # branches
  else
     size_fetch=nppm
     nbranch_max = 5*nintmax
     nbranch_local_max = 5*nintmax
  endif

  if (num_pe==1) then
    maxleaf = npart 
    maxtwig = maxaddress/2
  else if (num_pe.lt.2048) then
    maxleaf = maxaddress/3
    maxtwig = 2*maxleaf
  else
!  required # twigs increases with P because of branches
    maxleaf = maxaddress/4
    maxtwig = 3*maxleaf
  endif 

  hashconst = 2**nbaddr-1

! Memory estimate
  mem_lists = nshortm + nshortm*nintmax*(8+4)
  mem_tree =  maxaddress * (36 + 4 + 4 + 4) & ! # htable stuff
                      + num_pe * (4+4+4+4)  & ! request stuff
                      + maxaddress * (3*8) & ! keys
                      + size_fetch * 2*(22*8+2*4) & ! get_child, pack_child buffers
                      + nbranch_max * (8 + 4 + 8 + 8*23)  ! branches
  mem_multipoles = maxaddress * (8+3*4 + 23*8 + 8) 
  mem_prefetch = size_fetch*(8 + 4) + num_pe*4 *11 + maxaddress*8*2*2
  mem_tot = init_mb+mem_tree+mem_prefetch+mem_multipoles+mem_lists

  if (me==0 .and. tree_debug) then
     write(*,'(//a/)') 'Actual memory allocation:'
     write(*,'(6(a15,f14.3,a3/))') 'Inital Alloc:',init_mb/mb,' MB', &
                               'Tree:',mem_tree/mb,' MB', &
                               'Lists:',mem_lists/mb,' MB', &
                               'Prefetch:',mem_prefetch/mb,' MB', &
                               'Multipoles:',mem_multipoles/mb,' MB', &
                               'TOTAL: ',mem_tot/mb,' MB'
    write(*,*) '# procs',num_pe
    write(*,*) 'npart=',npart
    write(*,*) 'N/P=',npart/num_pe
    write(*,*) 'nppm= ',nppm
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

end subroutine tree_allocate

  
