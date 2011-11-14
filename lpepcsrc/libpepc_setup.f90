!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Initializes debug level
!> Creates and registers user-defined MPI types
!> Calls initialization for periodic framework
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine libpepc_setup(my_rank,n_cpu,db_level)
  use treevars
  use module_fmm_framework
  use module_branching
  implicit none
  integer, intent(in) :: db_level, my_rank, n_cpu

  ! copy call parameters to treevars module
  me     = my_rank
  num_pe = n_cpu

  force_debug=.false.
  tree_debug=.false.
  build_debug=.false.
  domain_debug = .false.
  branch_debug=.false.
  walk_summary=.false.
  dump_tree=.false.
  periodic_debug=.false.

  load_file_debug=.false.
  timing_file_debug=.false.


  if (db_level==1) then
!      domain_debug = .true.
      tree_debug=.true.  ! location information only
	
  else if (db_level==2) then
      tree_debug=.true.  ! location information only
      walk_summary=.true.

  else if (db_level==3) then
      tree_debug=.true.
      force_debug=.false.
      walk_summary=.true.
      domain_debug = .true.
      periodic_debug=.true.

  else if (db_level==4) then
     tree_debug=.true.
     domain_debug = .true.
     build_debug=.true.
     branch_debug=.true.
     props_debug=.true.
     walk_summary=.true.
     force_debug=.true.
     periodic_debug=.true.

  else if (db_level==5) then
     dump_tree=.true.

  else if (db_level==6) then
     tree_debug=.true.
     domain_debug = .true.
     build_debug=.true.
     branch_debug=.true.
     walk_summary=.true.
     force_debug=.true.
     dump_tree=.true.
     periodic_debug=.true.
     load_file_debug=.true.
     timing_file_debug=.true.

  else
    ! all off by default
  endif

  ! create and register mpi types
  call libpepc_register_mpi_types(db_level)

  ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
  call fmm_framework_init(me, wellsep = 1)

  ! initialize data structures in module_branches
  call branches_initialize()

end subroutine libpepc_setup



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Creates and registers user-defined MPI types
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine libpepc_register_mpi_types(db_level)
  use treevars
  use treetypes
  implicit none
  include 'mpif.h'
  integer, intent(in) :: db_level

  integer :: ierr,i

  type (t_particle) :: ship_props_a, get_props_a
  type (t_particle_results) :: ship_props_b, get_props_b

  integer, dimension(nprops_multipole) :: blocklengths, displacements, types

  ! address calculation, 8 byte
  integer(KIND=MPI_ADDRESS_KIND), dimension(nprops_multipole) :: address
  integer(KIND=MPI_ADDRESS_KIND) :: send_base, receive_base

  ! Create new contiguous datatype for shipping particle properties (15 arrays)
  blocklengths(1:nprops_particle) = 1   

  types(1:8)   = MPI_REAL8
  types(9)     = MPI_INTEGER8
  types(10:11) = MPI_INTEGER

  call MPI_GET_ADDRESS( get_props_a%x, receive_base, ierr )  ! Base address for receive buffer
  call MPI_GET_ADDRESS( ship_props_a%x, send_base, ierr )  ! Base address for send buffer

!  if (me==0) write(*,'(a30,o21)') 'Particle address base:',receive_base

  call MPI_GET_ADDRESS( ship_props_a%x(1), address(1), ierr )
  call MPI_GET_ADDRESS( ship_props_a%x(2), address(2), ierr )
  call MPI_GET_ADDRESS( ship_props_a%x(3), address(3), ierr )
  call MPI_GET_ADDRESS( ship_props_a%u(1), address(4), ierr )
  call MPI_GET_ADDRESS( ship_props_a%u(2), address(5), ierr )
  call MPI_GET_ADDRESS( ship_props_a%u(3), address(6), ierr )
  call MPI_GET_ADDRESS( ship_props_a%q, address(7), ierr )
  call MPI_GET_ADDRESS( ship_props_a%work, address(8), ierr )
  call MPI_GET_ADDRESS( ship_props_a%key, address(9), ierr )
  call MPI_GET_ADDRESS( ship_props_a%label, address(10), ierr )
  call MPI_GET_ADDRESS( ship_props_a%pid, address(11), ierr )

  displacements(1:nprops_particle) = int(address(1:nprops_particle) - send_base)  !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, mpi_type_particle, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_particle, ierr)

  if (me==0 .and. db_level>2) then
  ! Check addresses for MPI particle structure
     write(*,'(a30/(o28,i8))') 'Particle addresses:',(address(i),displacements(i),i=1,15)
  endif 

  ! Create new contiguous datatype for shipping result props (6 arrays)

  blocklengths(1:nprops_particle_results) = 1   

  types(1:5) = MPI_REAL8
  types(6) = MPI_INTEGER

  call MPI_GET_ADDRESS( get_props_b%e(1), receive_base, ierr )  ! Base address for receive buffer
  call MPI_GET_ADDRESS( ship_props_b%e(1), send_base, ierr )  ! Base address for send buffer

  call MPI_GET_ADDRESS( ship_props_b%e(1), address(1), ierr )
  call MPI_GET_ADDRESS( ship_props_b%e(2), address(2), ierr )
  call MPI_GET_ADDRESS( ship_props_b%e(3), address(3), ierr )
  call MPI_GET_ADDRESS( ship_props_b%pot,  address(4), ierr )
  call MPI_GET_ADDRESS( ship_props_b%work, address(5), ierr )

  displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - send_base)  !  Addresses relative to start of results (receive) data

  call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_results, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_results, ierr)

  if (me==0 .and. db_level>2) then
     ! Check addresses for MPI results structure
     write(*,'(a30/(o28,i8))') 'Results addresses:',(address(i),displacements(i),i=1,6)
  endif 


  ! Create new contiguous datatype for shipping multipole properties (25 arrays)

  blocklengths(1:nprops_multipole) = 1   

  types(1)      = MPI_INTEGER8
  types(2:4)    = MPI_INTEGER
  types(5:9)    = MPI_REAL8
  types(10)     = MPI_INTEGER
  types(11:22)  = MPI_REAL8


  call MPI_GET_ADDRESS( node_dummy%key, send_base, ierr )  ! Base address for send buffer

  call MPI_GET_ADDRESS( node_dummy%key,        address( 1), ierr )
  call MPI_GET_ADDRESS( node_dummy%byte,       address( 2), ierr )
  call MPI_GET_ADDRESS( node_dummy%leaves,     address( 3), ierr )
  call MPI_GET_ADDRESS( node_dummy%owner,      address( 4), ierr )
  call MPI_GET_ADDRESS( node_dummy%charge,     address( 5), ierr )
  call MPI_GET_ADDRESS( node_dummy%abs_charge, address( 6), ierr )
  call MPI_GET_ADDRESS( node_dummy%coc(1),     address( 7), ierr )
  call MPI_GET_ADDRESS( node_dummy%coc(2),     address( 8), ierr )
  call MPI_GET_ADDRESS( node_dummy%coc(3),     address( 9), ierr )
  call MPI_GET_ADDRESS( node_dummy%level,      address(10), ierr )
  call MPI_GET_ADDRESS( node_dummy%xdip,       address(11), ierr )
  call MPI_GET_ADDRESS( node_dummy%ydip,       address(12), ierr )
  call MPI_GET_ADDRESS( node_dummy%zdip,       address(13), ierr )
  call MPI_GET_ADDRESS( node_dummy%xxquad,     address(14), ierr )
  call MPI_GET_ADDRESS( node_dummy%yyquad,     address(15), ierr )
  call MPI_GET_ADDRESS( node_dummy%zzquad,     address(16), ierr )
  call MPI_GET_ADDRESS( node_dummy%xyquad,     address(17), ierr )
  call MPI_GET_ADDRESS( node_dummy%yzquad,     address(18), ierr )
  call MPI_GET_ADDRESS( node_dummy%zxquad,     address(19), ierr )
  call MPI_GET_ADDRESS( node_dummy%xshift,     address(20), ierr )
  call MPI_GET_ADDRESS( node_dummy%yshift,     address(21), ierr )
  call MPI_GET_ADDRESS( node_dummy%zshift,     address(22), ierr )

  displacements(1:nprops_multipole) = int(address(1:nprops_multipole) - send_base)   !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_multipole, blocklengths, displacements, types, mpi_type_multipole, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_multipole, ierr)

end subroutine libpepc_register_mpi_types





