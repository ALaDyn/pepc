subroutine pepc_setup(my_rank,n_cpu,npart_total,theta,db_level,np_mult,fetch_mult)
  use treevars
  implicit none
  include 'mpif.h'

  real, intent(in) :: theta     ! Multipole clumping parameter
  integer, intent(in) :: my_rank  ! MPI cpu rank
  integer, intent(in) :: n_cpu  ! MPI # CPUs
  integer, intent(in) :: npart_total  ! total (max) # simulation particles
  integer, intent(in) :: db_level
  real, intent(in) :: np_mult ! particle array size multiplication factor (1.5)
  integer, intent(in) :: fetch_mult ! fetch array size multiplication factor (10)

  integer :: ibig, machinebits,k
  integer :: ierr,npsize
  integer :: mem_parts, mem_multipoles, mem_fields, mem_tree, mem_prefetch, mem_tot, mem_lists
  character(3) :: cme
  character(30) :: cfile
  real, parameter :: mb=2.**20

! copy call parameters to treevars module

  me = my_rank
  num_pe = n_cpu
  npart = npart_total
  ipefile = 20

  if (db_level==1) then
!      domain_debug = .true.
	
  else if (db_level==2) then
      tree_debug=.true.
      force_debug=.false.
      walk_summary=.true.
      prefetch_debug=.false. 
      domain_debug = .true.

  else if (db_level==3) then
     tree_debug=.true.
     domain_debug = .true.
     build_debug=.true.
     branch_debug=.true.
     prefetch_debug=.true.
     walk_debug=.true.
     walk_summary=.true.
     force_debug=.true.

  else if (db_level==4) then
     dump_tree=.true.

  else if (db_level==5) then
     tree_debug=.true.
     domain_debug = .true.
     build_debug=.true.
     branch_debug=.true.
     prefetch_debug=.true.
     walk_debug=.true.
     walk_summary=.true.
     force_debug=.true.
     dump_tree=.true.
  else
! all off by default

  endif

  npartm = npart
  if (np_mult>0) then 
    nppm = abs(np_mult)*5*max(npartm/num_pe,1000) ! allow 50% fluctuation
  else
    nppm = 5*max(npartm/num_pe,1000) ! allow 50% fluctuation
  endif
  
  nshortm = 1000    ! Max shortlist length: leave safety factor for nshort_list in FORCES

  ! Estimate of interaction list length - Hernquist expression
  if (theta >0 ) then
     nintmax = max(1.*24*log(2.*npartm)/theta**2,2200.)
  else
     nintmax = npartm
  endif
  max_list_length = 0 ! current max length of all interaction lists

  ! tree stuff

  nlev = 20                     ! max refinement level
  iplace = 2_8**(3*nlev)           ! place holder bit

  !  Space for # table and tree arrays
  !  TODO: need good estimate for max # branches
   npsize=nppm
   size_tree = max(4*nintmax+npsize,10000)
   if (np_mult>0) then
     nbaddr = max(log(1.*size_tree)/log(2.) + 1,15.)
     maxaddress = 2**nbaddr
   else
     nbaddr = 17   ! fixed address range
     maxaddress = abs(np_mult)*10000
   endif
   size_fetch = fetch_mult*size_tree/2
   nbranch_max = maxaddress
!   nbranch_max = 4*nintmax*max(1.,log(1.*num_pe))
   if (num_pe==1) size_fetch=size_tree
!  maxaddress = 512
   hashconst = 2**nbaddr-1
  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)
  if (me==0 .and. db_level>0) then
    write(*,*) 'size_tree= ',size_tree
    write(*,*) 'max address = ',maxaddress
    write(*,*) 'max branches = ',nbranch_max
    write(*,*) 'size_fetch= ',size_fetch
    write(*,*) 'np_mult= ',np_mult
    write(*,*) 'fetch_mult= ',fetch_mult
  endif 

  work_local = 1 ! Initial value for local load

! Memory estimate
  mem_parts = nppm*(22*8 + 2*4 + 8)
  mem_lists = nshortm + nshortm*nintmax*(8+4)
  mem_tree =  maxaddress * (36 + 4 + 4 + 4) & ! # htable stuff
                      + num_pe * (4+4+4+4)  & ! request stuff
                      + maxaddress * (3*8) & ! keys
                      + size_fetch * 2*(22*8+2*4) & ! get_child, pack_child buffers
                      + nbranch_max * (8 + 4 + 8 + 8*23)  ! branches
  mem_multipoles = maxaddress * (8+3*4 + 23*8 + 8) 
  mem_prefetch = size_fetch*(8 + 4) + num_pe*4 *11 + maxaddress*8*2*2
  mem_tot = mem_parts+mem_tree+mem_prefetch+mem_multipoles+mem_lists

  if (me==0 .and. db_level>0) then
     write(*,'(//a/)') 'Initial memory allocation:'
     write(*,'(6(a15,f14.3,a3/)/)') 'Particles: ',mem_parts/mb,' MB', &
                               'Tree:',mem_tree/mb,' MB', &
                               'Lists:',mem_lists/mb,' MB', &
                               'Prefetch:',mem_prefetch/mb,' MB', &
                               'Multipoles:',mem_multipoles/mb,' MB', &
                               'TOTAL: ',mem_tot/mb,' MB'
     write(*,'(a)') 'Allocating particle and tree arrays ...'
  endif

  ! array allocation

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), work(nppm), &
       Ex(nppm), Ey(nppm), Ez(nppm), pot(nppm), &
       Ax(nppm), Ay(nppm), Az(nppm), &
       Bx(nppm), By(nppm), Bz(nppm),  &
       Axo(nppm), Ayo(nppm), Azo(nppm), &
       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE



  allocate ( nterm(nshortm), intlist(nintmax,nshortm), nodelist(nintmax,nshortm) )! interaction key-, node-lists

  allocate ( htable(0:maxaddress), all_addr(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
       nbranches(num_pe+2), igap(num_pe+3), &
       treekey(maxaddress), branch_key(nbranch_max), branch_owner(nbranch_max), &
       pebranch(nbranch_max), leaf_key(maxaddress), twig_key(maxaddress), &
       fetched_owner(size_fetch), fetched_keys(size_fetch), requested_owner(size_fetch), requested_keys(size_fetch), &
       nreqs_total(0:num_pe-1), nfetch_total(0:num_pe-1) )



  all_addr = (/ (k,k=0,maxaddress) /)      ! List of all possible # table addresses
  htable%node = 0
  htable%key = 0
  htable%link = -1
  htable%leaves = 0
  htable%childcode = 0

  ! Allocate memory for tree node properties

  maxleaf = maxaddress/3
  maxtwig = 2*maxleaf

  nreqs_total(0:num_pe-1) = 0   ! Zero cumulative fetch/ship counters for non-local nodes
  nfetch_total(0:num_pe-1) = 0  

  allocate ( first_child(-maxtwig:maxleaf), n_children(-maxtwig:maxleaf), node_level(-maxtwig:maxleaf) )

  allocate ( charge(-maxtwig:maxleaf), &                    ! charge
       abs_charge(-maxtwig:maxleaf), &                ! absolute charge
       xcoc(-maxtwig:maxleaf), ycoc(-maxtwig:maxleaf), zcoc(-maxtwig:maxleaf), &    ! centre of charge 
       xshift(-maxtwig:maxleaf), yshift(-maxtwig:maxleaf), zshift(-maxtwig:maxleaf), &    ! shift vector
       xdip(-maxtwig:maxleaf), ydip(-maxtwig:maxleaf), zdip(-maxtwig:maxleaf), &          ! dipole moment
       xxquad(-maxtwig:maxleaf), yyquad(-maxtwig:maxleaf), zzquad(-maxtwig:maxleaf), &       ! quadrupole moment
       xyquad(-maxtwig:maxleaf), yzquad(-maxtwig:maxleaf), zxquad(-maxtwig:maxleaf), &
       jx(-maxtwig:maxleaf), jy(-maxtwig:maxleaf), jz(-maxtwig:maxleaf), &      ! current
       magmx(-maxtwig:maxleaf), magmy(-maxtwig:maxleaf), magmz(-maxtwig:maxleaf) ) ! magnetic moment 

  allocate ( pack_child(size_fetch),get_child(size_fetch) )    ! Multipole shipping buffers

 

! work balance arrays
  allocate  (work_loads(num_pe),npps(num_pe),pivots(num_pe+1))  ! Work load & Particle distrib amoung PEs


  ! Create new contiguous datatype for shipping particle properties (15 arrays)

  blocklengths(1:nprops_particle) = 1   


  types(1:12) = MPI_REAL8
  types(13) = MPI_INTEGER8
  types(14:15) = MPI_INTEGER

  call MPI_ADDRESS( get_props%x, receive_base, ierr )  ! Base address for receive buffer
  call MPI_ADDRESS( ship_props%x, send_base, ierr )  ! Base address for send buffer



  call MPI_ADDRESS( ship_props%x, address(1), ierr )
  call MPI_ADDRESS( ship_props%y, address(2), ierr )
  call MPI_ADDRESS( ship_props%z, address(3), ierr )
  call MPI_ADDRESS( ship_props%ux, address(4), ierr )
  call MPI_ADDRESS( ship_props%uy, address(5), ierr )
  call MPI_ADDRESS( ship_props%uz, address(6), ierr )
  call MPI_ADDRESS( ship_props%q, address(7), ierr )
  call MPI_ADDRESS( ship_props%m, address(8), ierr )
  call MPI_ADDRESS( ship_props%work, address(9), ierr )
  call MPI_ADDRESS( ship_props%ax, address(10), ierr )
  call MPI_ADDRESS( ship_props%ay, address(11), ierr )
  call MPI_ADDRESS( ship_props%az, address(12), ierr )
  call MPI_ADDRESS( ship_props%key, address(13), ierr )
  call MPI_ADDRESS( ship_props%label, address(14), ierr )
  call MPI_ADDRESS( ship_props%pid, address(15), ierr )

  displacements(1:nprops_particle) = address(1:nprops_particle) - send_base  !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, mpi_type_particle, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_particle, ierr)

  ! Create new contiguous datatype for shipping multipole properties (25 arrays)

  blocklengths(1:nprops_multipole) = 1   


  types(1) = MPI_INTEGER8
  types(2:4) = MPI_INTEGER
  types(5) = MPI_INTEGER8
  types(6:25) = MPI_REAL8

  call MPI_ADDRESS( node_dummy%key, send_base, ierr )  ! Base address for send buffer



  call MPI_ADDRESS( node_dummy%key, address(1), ierr )
  call MPI_ADDRESS( node_dummy%byte, address(2), ierr )
  call MPI_ADDRESS( node_dummy%leaves, address(3), ierr )
  call MPI_ADDRESS( node_dummy%owner, address(4), ierr )
  call MPI_ADDRESS( node_dummy%next, address(5), ierr )
  call MPI_ADDRESS( node_dummy%q, address(6), ierr )
  call MPI_ADDRESS( node_dummy%absq, address(7), ierr )
  call MPI_ADDRESS( node_dummy%xcoc, address(8), ierr )
  call MPI_ADDRESS( node_dummy%ycoc, address(9), ierr )
  call MPI_ADDRESS( node_dummy%zcoc, address(10), ierr )
  call MPI_ADDRESS( node_dummy%xdip, address(11), ierr )
  call MPI_ADDRESS( node_dummy%ydip, address(12), ierr )
  call MPI_ADDRESS( node_dummy%zdip, address(13), ierr )
  call MPI_ADDRESS( node_dummy%xxquad, address(14), ierr )
  call MPI_ADDRESS( node_dummy%yyquad, address(15), ierr )
  call MPI_ADDRESS( node_dummy%zzquad, address(16), ierr )
  call MPI_ADDRESS( node_dummy%xyquad, address(17), ierr )
  call MPI_ADDRESS( node_dummy%yzquad, address(18), ierr )
  call MPI_ADDRESS( node_dummy%zxquad, address(19), ierr )
  call MPI_ADDRESS( node_dummy%jx, address(20), ierr )
  call MPI_ADDRESS( node_dummy%jy, address(21), ierr )
  call MPI_ADDRESS( node_dummy%jz, address(22), ierr )
  call MPI_ADDRESS( node_dummy%magmx, address(23), ierr )
  call MPI_ADDRESS( node_dummy%magmy, address(24), ierr )
  call MPI_ADDRESS( node_dummy%magmz, address(25), ierr )

  displacements(1:nprops_multipole) = address(1:nprops_multipole) - send_base   !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_multipole, blocklengths, displacements, types, mpi_type_multipole, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_multipole, ierr)

  if (me==0 .and. db_level>0) write(*,*) '... done'
  max_prefetches = 0

end subroutine pepc_setup






