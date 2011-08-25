!
!   Header module for treemp: global arrays and variables
!
!   $ Revision:   $ 
!   All key variables defined 64 bit (8 byte)
!
 
module treevars
  
  use treetypes

  implicit none

  ! Constants
  integer, dimension(0:7) :: bitarr = (/ 0,1,2,3,4,5,6,7 /)    ! Array of bit positions

  !  Associated MPI stuff
  integer :: me       ! Rank of current task
  integer :: num_pe   ! # cpus used by program

  type (multipole) :: node_dummy

  integer :: mpi_type_particle, mpi_type_multipole, mpi_type_results

  !  tree variables
  integer*8, allocatable :: &
                                treekey(:), &       ! keys of all twig and leaf nodes
                                branch_key(:), &    ! keys of branch nodes covering all domains
                                pebranch(:), &	    ! keys of branch nodes covering local domain
	                            twig_key(:)         ! local twig keys

  integer, allocatable :: &
                                nbranches(:), &       ! # branches in local domain
                                branch_owner(:), &    ! owners of branch nodes covering all domains
                                free_addr(:), &    ! List of free #table addresses (for HASHENTRY routine)
                                point_free(:)      ! Pointer to free address index

  real*8, allocatable  :: &                ! Tree node properties:
                               charge(:), &                          ! charge
                               abs_charge(:), &                      ! absolute charge
                               xcoc(:), ycoc(:), zcoc(:), &          ! centre of charge 
                               xshift(:), yshift(:), zshift(:), &    ! shift vector
                               xdip(:), ydip(:), zdip(:), &          ! dipole moment
                               xxquad(:), yyquad(:), zzquad(:), &    ! quadrupole moment
                               xyquad(:), yzquad(:), zxquad(:)       !

  integer, allocatable ::   node_level(:)       ! refinement level
 

  integer*8 ::  max_req_list_length, & ! maximum length of request queue
                 cum_req_list_length, & ! cumulative length of request queue
                 comm_loop_iterations(3)! number of comm loop iterations (total, sending, receiving)

  integer,   parameter :: nlev = 20 ! max refinement level
  integer*8, parameter :: iplace = 2_8**(3*nlev) ! value of place holder bit = 2^(idim*nlev)

  integer :: &
             nbaddr, &         ! # bits in hashing function
             nleaf, &          ! total # leaf nodes in local #table 
             ntwig, &          ! total # twig nodes in local #table
             nleaf_me, &       ! total # leaves in local domain
             ntwig_me, &       ! total # twigs in local domain
             nlist, &          ! # particles/PE + boundary bodies (1 or 2)
             nnodes, &         ! nleaf + ntwig
             nbranch, &        ! min # branch nodes covering local domain
             nbranch_sum, &    ! total # branch nodes covering all domains
             nintmax, &        ! max # terms allowed in interaction list
             maxaddress, &     ! max address allowed in #table
             maxleaf, &     ! max leaf allowed in #table
             maxtwig, &     ! max twig allowed in #table
             size_tree, &      ! array space needed for local tree
             maxships, &       ! max # multipole ships per traversal 
             sum_ships, &      ! total # multipole ships per iteration  
             sum_fetches, &    ! total # key fetches  per iteration  
             free_lo, &        ! min address allowed for resolving collisions
             tablehigh, &      ! highest current address in #table
             sum_unused, &     ! # free addresses
             npartm, &         ! absolute max # particles
             npart, &          ! actual # particles
             nppm, &           ! max # particles/PE
             npp, &            !  actual  # particles/PE
             iused          ! counter for collision resolution array free_addr()

  integer :: ipefile = 20 ! local O/P stream
  integer :: nkeys_total=1 ! total # keys in local tree
  integer :: proc_debug=0     ! Debug rank: set to -1 for all
  real*8 :: xmin, xmax    ! box limits
  real*8 :: ymin, ymax  
  real*8 :: zmin, zmax
  real*8 :: boxsize       ! length of box
  integer, parameter :: CHILDCODE_NODE_TOUCHED = 11 !< this bit is used inside the childcode to notify of nodes, that already contain valid multipole information and may not be set to zero in tree_global
  real*8 :: interactions_local = 0. !< number of interactions that have been processed locally
  real*8 :: mac_evaluations_local = 0.!< number of mac evaluations that have been processed locally
  real*8 :: thread_workload(-4:4) !< stores average particles and runtime per thread for diagnostic purposes, entry 0 contains number of worker threads
  ! Debugging switches (all off by default)
  logical :: tree_debug=.false.
  logical :: build_debug=.false.
  logical :: domain_debug = .false.
  logical :: branch_debug=.false.
  logical :: props_debug=.false.
  logical :: walk_debug=.false.
  logical :: walk_summary=.false.
  logical :: force_debug=.false.
  logical :: dump_tree=.false.
  logical :: timing_file_debug=.false.
  logical :: load_file_debug=.false.
  

  !  particle data - dimensions filled in when #PEs known

  real*8, allocatable :: work(:), &        ! interaction work load
                            x(:),   y(:),   z(:), &     ! position
                           ux(:),  uy(:),  uz(:), &     ! velocity
                            q(:)                ! charge

  integer*8, allocatable ::  pekey(:)  ! local particle keys

  integer, allocatable ::    pepid(:), & ! owner
                                pelabel(:)   ! particle label

! Memory control
  real :: np_mult=1.5

end module treevars



