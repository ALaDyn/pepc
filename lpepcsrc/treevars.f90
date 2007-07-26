!
!   Header module for treemp: global arrays and variables
!
!   $ Revision:   $ 
!   All key variables defined 64 bit (8 byte)
!
 
module treevars

!  implicit none


! fixed array sizes for debugging
!  integer, parameter :: size_tree = 10000, &
!                        maxaddress=32768, &
!                        nppm=2000, &
!                        nbranch_max=size_tree/10

 ! Constants

  integer, dimension(0:7) :: bitarr = (/ 0,1,2,3,4,5,6,7 /)    ! Array of bit positions

!  Associated MPI stuff

  integer :: me       ! Rank of current task
  integer :: num_pe   ! # cpus used by program

 ! Hash table datatype - 36 bytes per entry
 
  type hash
     integer   :: node          ! Address of particle/pseudoparticle data
     integer*8 :: key           ! Key
     integer   :: link          ! Pointer to next empty address in table in case of collision
     integer   :: leaves        ! # leaves contained within twig (=1 for leaf, npart for root)
     integer   :: childcode     ! Byte code indicating position of children (twig node); particle label (leaf node)
     integer*8 :: next          ! Pointer to next key to examine in tree-walk
     integer   :: owner         ! Node owner (for branches)
  end type hash

  !  hash table

  type (hash), allocatable :: htable(:)
  
  
 ! Data structure for shipping single particles
  
  type particle
     real*8 :: x    ! coords
     real*8 :: y
     real*8 :: z
     real*8 :: ux    ! momenta
     real*8 :: uy
     real*8 :: uz 
     real*8 :: q     ! charge
     real*8 :: m     ! mass
     real*8 :: work  ! work load from force sum
     real*8 :: ax   ! 'forces' 
     real*8 :: ay
     real*8 :: az
     integer*8 :: key           ! Key
     integer :: label    ! label
     integer :: pid      ! owner
  end type particle

  type (particle) :: ship_props, get_props



 ! Data structure for shipping multiple moments of child nodes
  
  type multipole
     integer*8 :: key     ! key
     integer   :: byte    ! byte code
     integer   :: leaves  ! # leaves contained
     integer :: owner    ! owner where multipole resides
     integer*8 :: next    ! next key on walk
     real*8 :: q    	! net charge sum
     real*8 :: absq  	!  absolute charge sum
     real*8 :: xcoc  	! center of charge
     real*8 :: ycoc    
     real*8 :: zcoc
     real*8 :: xdip  	! dipole moment
     real*8 :: ydip
     real*8 :: zdip
     real*8 :: xxquad  	! quadrupole moment
     real*8 :: yyquad  
     real*8 :: zzquad
     real*8 :: xyquad
     real*8 :: yzquad
     real*8 :: zxquad
     real*8 :: jx        ! current
     real*8 :: jy
     real*8 :: jz
     real*8 :: magmx      ! magnetic moment
     real*8 :: magmy
     real*8 :: magmz
  end type multipole

  type (multipole), allocatable :: pack_child(:), get_child(:)
  type (multipole) :: node_dummy

  integer :: mpi_type_particle, mpi_type_multipole

  !  tree variables

  integer*8, allocatable :: &
                                requested_keys(:), &  ! Remote nodes requested during tree walk/prefetch
                                fetched_keys(:), &  ! Remote nodes fetched during tree walk/prefetch
                                treekey(:), &       ! keys of all twig and leaf nodes
                                branch_key(:), &    ! keys of branch nodes covering all domains
                                pebranch(:), &	    ! keys of branch nodes covering local domain
				leaf_key(:), & 	    ! local leaf keys
	                        twig_key(:)         ! local twig keys

  integer, allocatable :: &
                                nbranches(:), &       ! # branches in local domain
                                igap(:), &    !  stride lengths of local branch arrays
                                branch_owner(:), &    ! owners of branch nodes covering all domains
                                all_addr(:), &  ! List of all possible #table addresses
                                free_addr(:), &    ! List of free #table addresses (for HASHENTRY routine)
                                point_free(:), &   ! Pointer to free address index
   				fetched_owner(:), & ! Owners of nonlocal keys fetched
   				requested_owner(:), & ! Recipients of multipole ships 
                                nreqs_total(:), &    ! total # nodes requested from local PE during tree walk
                                nfetch_total(:)   ! total # non-local nodes fetched during tree walk 

  real*8, allocatable  :: &                ! Tree node properties:
                               charge(:), &                          ! charge
                               abs_charge(:), &                      ! absolute charge
                               xcoc(:), ycoc(:), zcoc(:), &          ! center of charge 
                               xshift(:), yshift(:), zshift(:), &    ! shift vector
                               xdip(:), ydip(:), zdip(:), &          ! dipole moment
                               xxquad(:), yyquad(:), zzquad(:), &    ! quadrupole moment
                               xyquad(:), yzquad(:), zxquad(:), &    !
                               magmx(:), magmy(:), magmz(:), &       ! magnetic dipole moment
                               jx(:), jy(:), jz(:)                   ! current
  integer*8, allocatable ::    first_child(:)   ! key of first child

  integer, allocatable ::       n_children(:), &     ! # children
                                node_level(:)       ! refinement level
 

  integer*8 ::  hashconst, &   ! hashing constants
	        hashchild=7_8, &
                iplace         ! value of place holder bit = 2^(2*nlev)

  integer :: &
             nlev, &           ! max refinement level
             nbaddr, &         ! # bits in hashing function
             nleaf, &          ! total # leaf nodes in local #table 
             ntwig, &          ! total # twig nodes in local #table
             nleaf_me, &       ! total # leaves in local domain
             ntwig_me, &       ! total # twigs in local domain
             nlist, &          ! # particles/PE + boundary bodies (1 or 2)
             nnodes, &         ! nleaf + ntwig
             nnodes_me, &      ! nleaf_me + ntwig_me
             nnodes_pw, &      ! # nnodes prior to tree-walk
             nleaf_pw, &       ! # leaves prior to walk
             ntwig_pw, &       ! # twigs prior to walk
             nbranch, &        ! min # branch nodes covering local domain
             nbranch_sum, &    ! total # branch nodes covering all domains
             nintmax, &        ! max # terms allowed in interaction list
             max_list_length, & ! current max list length
             maxaddress, &     ! max address allowed in #table
             maxleaf, &     ! max leaf allowed in #table
             maxtwig, &     ! max twig allowed in #table
             size_tree, &      ! array space needed for local tree
             size_fetch, &      ! array space needed for fetch/request arrays 
             maxtraverse, &      ! max # traversals per iteration 
             maxships, &       ! max # multipole ships per traversal 
             sum_ships, &      ! total # multipole ships per iteration  
             sum_fetches, &    ! total # key fetches  per iteration  
             max_prefetches, &      ! total # multipole ships per prefetch 
             nbranch_max, &    ! array space needed for branches
             free_lo, &        ! min address allowed for resolving collisions
	     tablehigh, &      ! highest current address in #table 
             sum_unused, &     ! # free addresses
             ipefile, &             ! local O/P stream
             pe_debug, &       ! rank for debug O/P
             npartm, &         ! absolute max # particles
             npart, &          ! actual # particles
             nppm, &           ! max # particles/PE
             npp, &            !  actual  # particles/PE
             nshortm, &        ! shortlist length
             iused          ! counter for collision resolution array free_addr()
  real :: work_imbal=0.
  real :: part_imbal=0.
  real :: work_imbal_max, work_imbal_min  ! load stats
  integer ::  part_imbal_max, part_imbal_min
  real :: work_local ! Total local work load (=sum of interaction lists)
  real*8 :: xmin, xmax    ! box limits
  real*8 :: ymin, ymax  
  real*8 :: zmin, zmax
  real*8 :: boxsize       ! length of box

 ! Force control
  integer :: load_balance   ! Balances particles in || sort according to work load
  logical :: use_multipoles = .true.   ! Use of multipoles? 
  logical :: prefetch = .false.  ! Prefetch multipole info prior to walk

  ! Debugging switches (all off by default)
  logical :: tree_debug=.false.
  logical :: build_debug=.false.
  logical :: domain_debug = .false.
  logical :: branch_debug=.false.
  logical :: prefetch_debug=.false.
  logical :: walk_debug=.false.
  logical :: walk_summary=.false.
  logical :: force_debug=.false.
  logical :: dump_tree=.false.
  logical :: perf_anal=.false.  ! Performance analysis mode: turns off all diagnostic routines

  !  particle data - dimensions filled in when #PEs known

  real*8, allocatable :: x(:),  y(:),  z(:), &     ! position
                      ux(:), uy(:), uz(:), &     ! velocity
                              q(:),  m(:), &     ! charge and mass
			Ex(:), Ey(:), Ez(:), &   ! E-field
			pot(:), &	         ! scalar potential
                        Axo(:), Ayo(:), Azo(:), &   ! vector potential
                        Ax(:), Ay(:), Az(:), &   ! vector potential
                        Bx(:), By(:), Bz(:), &   ! B-field
			work(:)                  ! interaction work load 

  integer*8, allocatable ::   pekey(:), &  ! local particle keys                             
                              intlist(:,:) ! interaction key-list

  integer, allocatable ::     pepid(:), & ! owner
                              pelabel(:), &   ! particle label
                              nterm(:), &  ! # interactions 
                              nodelist(:,:) ! interaction key-list

  real, allocatable ::  work_loads(:)  ! Load balance array
  integer, allocatable :: npps(:)  ! Particle distrib amoung PEs
  integer*8, allocatable ::  pivots(:)  ! Pivot buffer for sort

end module treevars



