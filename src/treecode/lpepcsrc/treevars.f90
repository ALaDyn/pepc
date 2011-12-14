!
!   Header module for treemp: global arrays and variables
!
!   $ Revision:   $ 
!   All key variables defined 64 bit (8 byte)
!
 
module treevars
  
  use module_interaction_specific_types
  use module_pepc_types

  implicit none

  !  Associated MPI stuff
  integer :: me       !< Rank of current task
  integer :: num_pe   !< # cpus used by program

  !  tree variables
  integer*8, allocatable :: &
                                branch_key(:), &    !< keys of branch nodes covering all domains
                                pebranch(:)         !< keys of branch nodes covering local domain

  integer, allocatable :: &
                                nbranches(:), &       !< # branches in local domain
                                branch_owner(:)       !< owners of branch nodes covering all domains

  type(t_tree_node_interaction_data), target, allocatable  :: tree_nodes(:)                 !< Tree node properties TODO: move to module_tree

  integer,   parameter :: nlev = 20 !< max refinement level
  integer*8, parameter :: iplace = 2_8**(3*nlev) !< value of place holder bit = 2^(idim*nlev)

  integer :: &
             nleaf, &          ! total # leaf nodes in local #table 
             ntwig, &          ! total # twig nodes in local #table
             nleaf_me, &       ! total # leaves in local domain
             ntwig_me, &       ! total # twigs in local domain
             nlist, &          ! # particles/PE + boundary bodies (1 or 2)
             nbranch, &        ! min # branch nodes covering local domain
             nbranch_sum, &    ! total # branch nodes covering all domains
             nintmax, &        ! max # terms allowed in interaction list
             maxleaf, &        ! max leaf allowed in #table
             maxtwig, &        ! max twig allowed in #table
             size_tree, &      ! array space needed for local tree
             maxships, &       ! max # multipole ships per traversal 
             sum_ships, &      ! total # multipole ships per iteration  
             sum_fetches, &    ! total # key fetches  per iteration  
             npart, &          ! actual # particles (total)
             npp               !  actual  # particles/PE

  integer :: nkeys_total=1 ! total # keys in local tree
  real*8 :: xmin, xmax    ! box limits
  real*8 :: ymin, ymax  
  real*8 :: zmin, zmax
  real*8 :: boxsize       ! length of box
  real*8 :: interactions_local = 0. !< number of interactions that have been processed locally
  real*8 :: mac_evaluations_local = 0.!< number of mac evaluations that have been processed locally

! Memory control
  real :: np_mult=1.5


end module treevars



