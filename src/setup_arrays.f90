subroutine setup_arrays
  use treevars
  use physvars
  use utils
  implicit none

  integer :: ibig, machinebits, maxleaf, maxtwig,k
  integer :: ierr

  !  npartm = npart + nt*np_beam  ! Max # particles permitted
  npartm = npart + np_beam  ! allow 50% fluctuation

  if (scheme==5 .or. target_dup) npartm=npartm*2  ! reserve extra space for electrons in ions-only mode
	                                        ! or double-target config

  nppm = 2.*max(npartm/num_pe,1000) ! allow 50% fluctuation
  nshortm = 2000    ! Max shortlist length: leave safety factor for nshort_list in FORCES

  ! Estimate of interaction list length - Hernquist expression
  if (theta >0 ) then
     nintmax = 2.5*24*log(2.*npartm)/theta**2
  else
     nintmax = npartm
  endif
  max_list_length = 0 ! current max length of all interaction lists

  ! tree stuff

  idim = 3               ! # dimensions (2 or 3)
  nlev = 20                     ! max refinement level
  iplace = 2_8**(idim*nlev)           ! place holder bit
  !  nbaddr = 15                  ! # bits for cell address in hash-table
  !  Space for # table and tree arrays
  !  TODO: need good estimate for max # branches

  size_tree = max(4*nintmax+8*nppm,2000)+1
  maxaddress = size_tree
  nbaddr = log(1.*maxaddress)/log(2.) + 2
  maxaddress = 2**nbaddr
!  maxaddress = 512
  hashconst = maxaddress-1

  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)




  ! array allocation

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), Ex(nppm), Ey(nppm), Ez(nppm), pot(nppm), &
       Ax(nppm), Ay(nppm), Az(nppm), Bx(nppm), By(nppm), Bz(nppm), work(nppm), &
       Axo(nppm), Ayo(nppm), Azo(nppm), &
       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE

  allocate ( xslice(nppm), yslice(nppm), zslice(nppm), uxslice(nppm), uyslice(nppm), uzslice(nppm), & 
       qslice(nppm), mslice(nppm) )    ! Reserve slice particle array space N/NPE

  allocate ( nterm(nshortm), intlist(nintmax,nshortm), nodelist(nintmax,nshortm) )      ! Space for interaction lists


  allocate ( htable(0:maxaddress), all_addr(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
       nbranches(num_pe+2), igap(num_pe+3), &
       treekey(size_tree), branch_key(size_tree), branch_owner(size_tree), &
       pebranch(size_tree), leaf_key(nppm), twig_key(nppm), &
       requested_keys(size_tree, 0:num_pe-1), fetched_keys(size_tree, 0:num_pe-1), &
       nreqs_total(0:num_pe-1), nfetch_total(0:num_pe-1) )

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
  allocate (rhoi(0:ngx+1,0:ngy+1,0:ngz+1),rhoe(0:ngx+1,0:ngy+1,0:ngz+1))  ! global field arrays

! local field arrays for cycle-averages
  allocate (rhoe_loc(0:ngx+1,0:ngy+1,0:ngz+1), rhoi_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      ex_loc(0:ngx+1,0:ngy+1,0:ngz+1), ey_loc(0:ngx+1,0:ngy+1,0:ngz+1),  ez_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      bx_loc(0:ngx+1,0:ngy+1,0:ngz+1), by_loc(0:ngx+1,0:ngy+1,0:ngz+1),  bz_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      jxe_loc(0:ngx+1,0:ngy+1,0:ngz+1), jye_loc(0:ngx+1,0:ngy+1,0:ngz+1), jze_loc(0:ngx+1,0:ngy+1,0:ngz+1) )   



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

end subroutine setup_arrays
