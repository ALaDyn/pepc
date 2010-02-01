subroutine pepc_setup(my_rank,n_cpu,npart_total,theta,db_level,t_np_mult,t_fetch_mult,init_mb,nppm_ori)
  use treevars
  use tree_utils
  implicit none
  include 'mpif.h'

  real, intent(in) :: theta,t_np_mult     ! Multipole clumping parameter
  integer, intent(in) :: my_rank,t_fetch_mult  ! MPI cpu rank
  integer, intent(in) :: n_cpu  ! MPI # CPUs
  integer, intent(in) :: npart_total  ! total (max) # simulation particles
  integer, intent(in) :: db_level
  integer, intent(out) :: init_mb,nppm_ori

  integer :: ibig, machinebits,k
  integer :: ierr,npsize,i
  integer :: mem_parts, mem_multipoles, mem_fields, mem_tree, mem_prefetch, mem_tot, mem_lists
  character(3) :: cme
  character(30) :: cfile
  real, parameter :: mb=2.**20

  type (particle) :: ship_props_a, get_props_a
  type (results) :: ship_props_b, get_props_b

  integer, parameter :: nprops_particle=15, &    ! # particle properties to ship
  			nprops_multipole=25, &      ! Number of multipole properties to ship
                        nprops_results=6       ! # results to ship
  integer, dimension(nprops_multipole) :: blocklengths, displacements, types

  ! address calculation, 8 byte 
  integer*8, dimension(nprops_multipole) :: address
  integer*8 :: send_base, receive_base

! copy call parameters to treevars module
  
  np_mult = t_np_mult
  fetch_mult = t_fetch_mult
  me = my_rank
  num_pe = n_cpu
  npart = npart_total
  ipefile = 20

  force_debug=.false.
  tree_debug=.false.
  build_debug=.false.
  domain_debug = .false.
  branch_debug=.false.
  prefetch_debug=.false.
  walk_debug=.false.
  walk_summary=.false.
  dump_tree=.false.


  if (db_level==1) then
!      domain_debug = .true.
	
  else if (db_level==2) then
      tree_debug=.true.  ! location information only

  else if (db_level==3) then
      tree_debug=.true.
      force_debug=.false.
      walk_summary=.true.
      prefetch_debug=.false. 
      domain_debug = .true.

  else if (db_level==4) then
     tree_debug=.true.
     domain_debug = .true.
     build_debug=.true.
     branch_debug=.true.
     props_debug=.true.
     prefetch_debug=.true.
     walk_debug=.true.
     walk_summary=.true.
     force_debug=.true.

  else if (db_level==5) then
     dump_tree=.true.

  else if (db_level==6) then
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
  if (num_pe.eq.1) then
    nppm=1.5*npart + 1000  ! allow for additional ghost particles for field plots
!  else if (np_mult<0) then 
!    nppm = abs(np_mult)*max(npartm/num_pe,1000) ! allow 50% fluctuation
  else
    nppm = 2*max(npartm/num_pe,1000) ! allow 50% fluctuation
  endif
  
  nppm_ori = nppm

  nshortm = 2000    ! Max shortlist length: leave safety factor for nshort_list in FORCES
  nlev = 20                     ! max refinement level
  iplace = 2_8**(3*nlev)           ! place holder bit
  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)
  work_local = 1 ! Initial value for local load

  ! array allocation

  init_mb = init_mb + nppm*(22*8 + 2*4 + 8) + 8*num_pe*4

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), work(nppm), &
       Ex(nppm), Ey(nppm), Ez(nppm), pot(nppm), &
       Ax(nppm), Ay(nppm), Az(nppm), &
       Bx(nppm), By(nppm), Bz(nppm),  &
       Axo(nppm), Ayo(nppm), Azo(nppm), &
       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE

  allocate (work_loads(num_pe),npps(num_pe),pivots(num_pe+1))  ! Work load & Particle distrib amoung PEs

  allocate (nbranches(num_pe+2), igap(num_pe+3), nreqs_total(0:num_pe-1), nfetch_total(0:num_pe-1) )

  nreqs_total(0:num_pe-1) = 0   ! Zero cumulative fetch/ship counters for non-local nodes
  nfetch_total(0:num_pe-1) = 0  

  ! Create new contiguous datatype for shipping particle properties (15 arrays)

  blocklengths(1:nprops_particle) = 1   


  types(1:12) = MPI_REAL8
  types(13) = MPI_INTEGER8
  types(14:15) = MPI_INTEGER

!  receive_base=LOC(get_props_a%x)
  call LOCADDRESS( get_props_a%x, receive_base, ierr )  ! Base address for receive buffer
  call LOCADDRESS( ship_props_a%x, send_base, ierr )  ! Base address for send buffer

!  if (me==0) write(*,'(a30,o21)') 'Particle address base:',receive_base
!  call MPI_GET_ADDRESS( get_props_a%x, receive_base, ierr )  ! Base address for receive buffer
!  call MPI_GET_ADDRESS( ship_props_a%x, send_base, ierr )  ! Base address for send buffer

  call LOCADDRESS( ship_props_a%x, address(1), ierr )
  call LOCADDRESS( ship_props_a%y, address(2), ierr )
  call LOCADDRESS( ship_props_a%z, address(3), ierr )
  call LOCADDRESS( ship_props_a%ux, address(4), ierr )
  call LOCADDRESS( ship_props_a%uy, address(5), ierr )
  call LOCADDRESS( ship_props_a%uz, address(6), ierr )
  call LOCADDRESS( ship_props_a%q, address(7), ierr )
  call LOCADDRESS( ship_props_a%m, address(8), ierr )
  call LOCADDRESS( ship_props_a%work, address(9), ierr )
  call LOCADDRESS( ship_props_a%ax, address(10), ierr )
  call LOCADDRESS( ship_props_a%ay, address(11), ierr )
  call LOCADDRESS( ship_props_a%az, address(12), ierr )
  call LOCADDRESS( ship_props_a%key, address(13), ierr )
  call LOCADDRESS( ship_props_a%label, address(14), ierr )
  call LOCADDRESS( ship_props_a%pid, address(15), ierr )

  displacements(1:nprops_particle) = address(1:nprops_particle) - send_base  !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, mpi_type_particle, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_particle, ierr)

  if (me==0 .and. db_level>1) then
! Check addresses for MPI particle structure
     write(*,'(a30/(o28,i8))') 'Particle addresses:',(address(i),displacements(i),i=1,15)
  endif 

  ! Create new contiguous datatype for shipping result props (6 arrays)

  blocklengths(1:nprops_results) = 1   

  types(1:5) = MPI_REAL8
  types(6) = MPI_INTEGER

  call LOCADDRESS( get_props_b%Ex, receive_base, ierr )  ! Base address for receive buffer
  call LOCADDRESS( ship_props_b%Ex, send_base, ierr )  ! Base address for send buffer

  call LOCADDRESS( ship_props_b%Ex, address(1), ierr )
  call LOCADDRESS( ship_props_b%Ey, address(2), ierr )
  call LOCADDRESS( ship_props_b%Ez, address(3), ierr )
  call LOCADDRESS( ship_props_b%pot, address(4), ierr )
  call LOCADDRESS( ship_props_b%work, address(5), ierr )
  call LOCADDRESS( ship_props_b%label, address(6), ierr )

  displacements(1:nprops_results) = address(1:nprops_results) - send_base  !  Addresses relative to start of results (receive) data

  call MPI_TYPE_STRUCT( nprops_results, blocklengths, displacements, types, mpi_type_results, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_results, ierr)

  if (me==0 .and. db_level>1) then
     ! Check addresses for MPI results structure
     write(*,'(a30/(o28,i8))') 'Results addresses:',(address(i),displacements(i),i=1,6)
  endif 


  ! Create new contiguous datatype for shipping multipole properties (25 arrays)

  blocklengths(1:nprops_multipole) = 1   


  types(1) = MPI_INTEGER8
  types(2:4) = MPI_INTEGER
  types(5) = MPI_INTEGER8
  types(6:25) = MPI_REAL8

  call LOCADDRESS( node_dummy%key, send_base, ierr )  ! Base address for send buffer



  call LOCADDRESS( node_dummy%key, address(1), ierr )
  call LOCADDRESS( node_dummy%byte, address(2), ierr )
  call LOCADDRESS( node_dummy%leaves, address(3), ierr )
  call LOCADDRESS( node_dummy%owner, address(4), ierr )
  call LOCADDRESS( node_dummy%next, address(5), ierr )
  call LOCADDRESS( node_dummy%q, address(6), ierr )
  call LOCADDRESS( node_dummy%absq, address(7), ierr )
  call LOCADDRESS( node_dummy%xcoc, address(8), ierr )
  call LOCADDRESS( node_dummy%ycoc, address(9), ierr )
  call LOCADDRESS( node_dummy%zcoc, address(10), ierr )
  call LOCADDRESS( node_dummy%xdip, address(11), ierr )
  call LOCADDRESS( node_dummy%ydip, address(12), ierr )
  call LOCADDRESS( node_dummy%zdip, address(13), ierr )
  call LOCADDRESS( node_dummy%xxquad, address(14), ierr )
  call LOCADDRESS( node_dummy%yyquad, address(15), ierr )
  call LOCADDRESS( node_dummy%zzquad, address(16), ierr )
  call LOCADDRESS( node_dummy%xyquad, address(17), ierr )
  call LOCADDRESS( node_dummy%yzquad, address(18), ierr )
  call LOCADDRESS( node_dummy%zxquad, address(19), ierr )
  call LOCADDRESS( node_dummy%jx, address(20), ierr )
  call LOCADDRESS( node_dummy%jy, address(21), ierr )
  call LOCADDRESS( node_dummy%jz, address(22), ierr )
  call LOCADDRESS( node_dummy%magmx, address(23), ierr )
  call LOCADDRESS( node_dummy%magmy, address(24), ierr )
  call LOCADDRESS( node_dummy%magmz, address(25), ierr )

  displacements(1:nprops_multipole) = address(1:nprops_multipole) - send_base   !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_multipole, blocklengths, displacements, types, mpi_type_multipole, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_multipole, ierr)

  max_prefetches = 0

end subroutine pepc_setup







