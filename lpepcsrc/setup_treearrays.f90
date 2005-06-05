subroutine pepc_setup(my_rank,n_cpu,npart_total,theta)
  use treevars
  implicit none
  include 'mpif.h'

  real, intent(in) :: theta     ! Multipole clumping parameter
  integer, intent(in) :: my_rank  ! MPI cpu rank
  integer, intent(in) :: n_cpu  ! MPI # CPUs
  integer, intent(in) :: npart_total  ! total (max) # simulation particles


  integer :: ibig, machinebits, maxleaf, maxtwig,k
  integer :: ierr,npsize
  integer :: mem_parts, mem_multipoles, mem_fields, mem_tree, mem_prefetch
  character(3) :: cme
  character(30) :: cfile

! copy call parameters to treevars module

  me = my_rank
  num_pe = n_cpu
  npart = npart_total


  npartm = npart 
  nppm = 1.5*max(npartm/num_pe,1000) ! allow 50% fluctuation
  nshortm = 2000    ! Max shortlist length: leave safety factor for nshort_list in FORCES

  ! Estimate of interaction list length - Hernquist expression
  if (theta >0 ) then
     nintmax = 2.*24*log(2.*npartm)/theta**2
  else
     nintmax = npartm
  endif
  max_list_length = 0 ! current max length of all interaction lists

  ! tree stuff

  nlev = 20                     ! max refinement level
  iplace = 2_8**(3*nlev)           ! place holder bit

  !  Space for # table and tree arrays
  !  TODO: need good estimate for max # branches
!   npsize=2.5*nppm
   npsize=nppm
!   size_tree = max(4*nintmax+npsize,1000000)+1
!   nbaddr = max(log(1.*size_tree)/log(2.) + 1,17.)
   nbaddr = 18   ! fixed address range
   maxaddress = 2**nbaddr
   size_tree=maxaddress+1
!   size_fetch = min(60*size_tree/num_pe,size_tree/2) 
   size_fetch=size_tree
   nbranch_max = size_tree/20
   if (num_pe==1) size_fetch=size_tree
!  maxaddress = 512
  hashconst = maxaddress-1
  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)
  if (me==0) then
    write(*,*) 'size_tree= ',size_tree
    write(*,*) 'size_fetch= ',size_fetch
  endif 


  ! array allocation

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), work(nppm), &
! Ex(nppm), Ey(nppm), Ez(nppm), pot(nppm), &
       Ax(nppm), Ay(nppm), Az(nppm), &
! Bx(nppm), By(nppm), Bz(nppm),  &
!       Axo(nppm), Ayo(nppm), Azo(nppm), &
       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE

  mem_parts = nppm*(22*8 + 2*4 + 8)


  allocate ( nterm(nshortm), intlist(nintmax,nshortm), nodelist(nintmax,nshortm) )! interaction key-, node-lists
  mem_tree = nshortm + nshortm*nintmax*(8+4)

  allocate ( htable(0:maxaddress), all_addr(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
       nbranches(num_pe+2), igap(num_pe+3), &
       treekey(size_tree), branch_key(nbranch_max), branch_owner(nbranch_max), &
       pebranch(nbranch_max), leaf_key(size_tree), twig_key(size_tree), &
       requested_keys(size_fetch, 0:num_pe-1), fetched_keys(size_fetch, 0:num_pe-1), &
       nreqs_total(0:num_pe-1), nfetch_total(0:num_pe-1) )

  mem_tree = mem_tree + maxaddress * (36 + 4 + 4 + 4) & ! # htable stuff
                      + num_pe * (4+4+4+4) + 2*size_fetch*num_pe*8 & ! request stuff
                      + size_tree * (3*8) & ! keys
                      + nbranch_max * (8 + 4 + 8)  ! branches

  all_addr = (/ (k,k=0,maxaddress) /)      ! List of all possible # table addresses
  htable%node = 0
  htable%key = 0
  htable%link = -1
  htable%leaves = 0
  htable%childcode = 0

  ! Allocate memory for tree node properties

  maxtwig = size_tree
  maxleaf = size_tree

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

  allocate ( pack_child(size_tree), get_child(size_tree) )    ! Multipole shipping buffers

  mem_multipoles = 2*size_tree * (8+2*4 + 23*8 + 8) 
 

! work balance arrays
  allocate  (work_loads(num_pe),npps(num_pe))  ! Work load & Particle distrib amoung PEs


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

  ! Create new contiguous datatype for shipping multipole properties (24 arrays)

  blocklengths(1:nprops_multipole) = 1   


  types(1) = MPI_INTEGER8
  types(2:3) = MPI_INTEGER
  types(4) = MPI_INTEGER8
  types(5:24) = MPI_REAL8

  call MPI_ADDRESS( node_dummy%key, send_base, ierr )  ! Base address for send buffer



  call MPI_ADDRESS( node_dummy%key, address(1), ierr )
  call MPI_ADDRESS( node_dummy%byte, address(2), ierr )
  call MPI_ADDRESS( node_dummy%leaves, address(3), ierr )
  call MPI_ADDRESS( node_dummy%next, address(4), ierr )
  call MPI_ADDRESS( node_dummy%q, address(5), ierr )
  call MPI_ADDRESS( node_dummy%absq, address(6), ierr )
  call MPI_ADDRESS( node_dummy%xcoc, address(7), ierr )
  call MPI_ADDRESS( node_dummy%ycoc, address(8), ierr )
  call MPI_ADDRESS( node_dummy%zcoc, address(9), ierr )
  call MPI_ADDRESS( node_dummy%xdip, address(10), ierr )
  call MPI_ADDRESS( node_dummy%ydip, address(11), ierr )
  call MPI_ADDRESS( node_dummy%zdip, address(12), ierr )
  call MPI_ADDRESS( node_dummy%xxquad, address(13), ierr )
  call MPI_ADDRESS( node_dummy%yyquad, address(14), ierr )
  call MPI_ADDRESS( node_dummy%zzquad, address(15), ierr )
  call MPI_ADDRESS( node_dummy%xyquad, address(16), ierr )
  call MPI_ADDRESS( node_dummy%yzquad, address(17), ierr )
  call MPI_ADDRESS( node_dummy%zxquad, address(18), ierr )
  call MPI_ADDRESS( node_dummy%jx, address(19), ierr )
  call MPI_ADDRESS( node_dummy%jy, address(20), ierr )
  call MPI_ADDRESS( node_dummy%jz, address(21), ierr )
  call MPI_ADDRESS( node_dummy%magmx, address(22), ierr )
  call MPI_ADDRESS( node_dummy%magmy, address(23), ierr )
  call MPI_ADDRESS( node_dummy%magmz, address(24), ierr )

  displacements(1:nprops_multipole) = address(1:nprops_multipole) - send_base   !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_multipole, blocklengths, displacements, types, mpi_type_multipole, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_multipole, ierr)

  mem_prefetch = size_fetch*(2*8*num_pe + 5*8) + num_pe*4*11 + size_fetch*(8+4)

  if (me==0) then
     write(*,'(//a/)') 'Initial memory allocation:'
     write(*,'(5(a15,f12.3,a3/)/)') 'Particles: ',mem_parts/1.e6,' MB', &
                               'Tree:',mem_tree/1.e6,' MB', &
                               'Prefetch:',mem_prefetch/1.e6,' MB', &
                               'Multipoles:',mem_multipoles/1.e6,' MB'
  endif

  cme = achar(me/100+48) // achar(mod(me/10,10)+48) // achar(mod(me,10)+48)  ! Convert 3-digit PE number into character string
  cfile="pe"//cme//"/dump."//cme
  ipefile = 20
  open(ipefile,file=cfile)


  sumprefetches = 0

end subroutine pepc_setup






