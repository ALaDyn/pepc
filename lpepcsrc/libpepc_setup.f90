!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Initializes debug level
!> Creates and registers user-defined MPI types
!> Calls initialization for periodic framework
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine libpepc_setup(my_rank,n_cpu,db_level)
  use treevars
  use treetypes
  use module_fmm_framework
  use module_branching
  implicit none
  integer, intent(in) :: db_level, my_rank, n_cpu

  call status('SETUP')

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
  call register_lpepc_mpi_types()

  ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
  call fmm_framework_init(me, wellsep = 1)

  ! initialize data structures in module_branches
  call branches_initialize()

end subroutine libpepc_setup
