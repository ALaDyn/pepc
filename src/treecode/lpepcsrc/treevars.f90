!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates all global variables for lpepc
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module treevars
  
  use module_interaction_specific_types
  use module_pepc_types

  implicit none

  !  Associated MPI stuff
  integer :: me       !< Rank of current task
  integer :: num_pe   !< # cpus used by program
  integer :: MPI_COMM_lpepc !> communicator that has been supplied to or created by pepc_initialize

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
             npp, &            ! actual  # particles/PE
             idim = 3          ! dimension of the system (currently only used in tree_domains)

  real*8 :: boxmin(3)  ! box min limits
  real*8 :: boxmax(3)  ! box max limits
  real*8 :: boxsize(3) ! box extension
  real*8 :: interactions_local = 0. !< number of interactions that have been processed locally
  real*8 :: mac_evaluations_local = 0.!< number of mac evaluations that have been processed locally

! Memory control
  real    :: np_mult=1.5
  integer :: defer_list_length_factor = 1 !< factor for increasing defer_list_length in case of respective warning (e.g. for very inhomogeneous or 2D cases set to 2..8)


end module treevars


