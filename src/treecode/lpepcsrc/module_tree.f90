! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!>
!> Defines a derived type `t_tree` that represents distributed hashed octrees
!> and associated procedures.
!>
module module_tree
    use module_htable, only: t_htable
    use module_box, only: t_box
    use module_comm_env, only: t_comm_env
    use module_domains, only: t_decomposition
    use module_atomic_ops, only: t_atomic_int, t_critical_section
    use module_pepc_types
    use, intrinsic :: iso_c_binding
    implicit none
    private

    integer(kind_key), public, parameter :: TREE_KEY_ROOT    =  1_kind_key

    integer, public, parameter :: TREE_COMM_ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
    integer, public, parameter :: TREE_COMM_REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue

    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STOPPED  = 1
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STARTING = 2
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_STARTED  = 3
    integer, public, parameter :: TREE_COMM_THREAD_STATUS_WAITING  = 4

    !> data type for tree communicator
    type, public :: t_tree_communicator
      ! counters and timers
      integer*8 :: comm_loop_iterations(3) !< number of comm loop iterations (total, sending, receiving)
      real*8 :: timings_comm(3) !< array for storing internal timing information
      integer(kind_node) :: sum_ships   !< total number of node ships
      integer(kind_node) :: sum_fetches !< total number of node fetches

      ! thread data
      type(c_ptr) :: comm_thread
      type(t_atomic_int), pointer :: thread_status
      type(t_critical_section), pointer :: cs_request
      integer :: processor_id
    end type t_tree_communicator

    !>
    !> A derived type representing a distributed hashed octree over a collection
    !> of particles.
    !>
    type, public :: t_tree
      integer(kind_particle) :: npart       !< number of particles across all ranks
      integer(kind_particle) :: npart_me    !< number of particles on local rank

      integer(kind_node) :: nleaf       !< number of leaves stored locally
      integer(kind_node) :: ntwig       !< number of twigs stored locally
      integer(kind_node) :: nleaf_me    !< number of leaves that originated on this rank
      integer(kind_node) :: ntwig_me    !< number of twigs that originated on this rank
      
      integer(kind_node)   :: nbranch     !< number of branch nodes in tree
      integer(kind_node)   :: nbranch_me  !< number of branch nodes that originated on this rank
      integer(kind_node)   :: nbranch_max_me !< upper limit estimate for number of local branch nodes

      integer(kind_particle) :: nintmax     !< maximum number of interactions
      
      real*8, allocatable :: boxlength2(:) !< precomputed square of maximum edge length of boxes for different levels - used for MAC evaluation
      
      type(t_box) :: bounding_box               !< bounding box enclosing all particles contained in the tree
      type(t_htable) :: node_storage            !< hash table in which tree nodes are stored for rapid retrieval
      type(t_tree_node), pointer :: nodes(:)    !< array of tree nodes, shorthand to node_sorage%values
      integer(kind_node) :: node_root            !< index of the root node in nodes-array
      type(t_comm_env) :: comm_env              !< communication environment over which the tree is distributed
      type(t_decomposition) :: decomposition    !< permutation of particles inserted into the tree
      type(t_tree_communicator) :: communicator !< associated communicator structure
    end type t_tree

    public tree_create
    public tree_allocated
    public tree_insert_node
    public tree_insert_or_update_node
    public tree_contains_key
    public tree_node_get_parent
    public tree_node_connect_children
    public tree_check
    public tree_dump
    public tree_stats
    public tree_destroy

    public tree_lookup_node
    public tree_lookup_node_i
    public tree_lookup_node_critical_i
    public tree_lookup_node_critical
    public tree_lookup_node_p
    public tree_lookup_node_critical_p

    interface tree_lookup_node
      module procedure tree_lookup_node_i, tree_lookup_node_p
    end interface

    interface tree_lookup_node_critical
      module procedure tree_lookup_node_critical_i, tree_lookup_node_critical_p
    end interface

    contains

    !>
    !> Create a tree (allocates memory but does not fill the tree)
    !> 
    !> Uses particle numbers (local and global) to estimate the memory needed
    !> for node storage.
    !> A communication environment over which the tree is distributed can be
    !> supplied as an MPI communicator `comm` or a `t_comm_env` in `comm_env`.
    !> If no environment is supplied, a duplicate of the PEPC global environment
    !> is used.
    !>
    subroutine tree_create(t, nl, n, comm, comm_env)
      use treevars, only: interaction_list_length_factor, MPI_COMM_lpepc, np_mult, nlev
      use module_htable, only: htable_create
      use module_interaction_specific, only: get_number_of_interactions_per_particle
      use module_comm_env, only: comm_env_dup, comm_env_mirror
      use module_timings
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t !< The tree
      integer(kind_particle), intent(in) :: nl !< Number of local particles to be inserted into the tree
      integer(kind_particle), intent(in) :: n !< Total number of particles across communication ranks
      integer, optional, intent(in) :: comm !< An MPI communicator
      type(t_comm_env), optional, intent(in) :: comm_env !< A communication environment

      integer(kind_node) :: maxaddress
      integer :: i

      call pepc_status('ALLOCATE TREE')
      DEBUG_ASSERT(.not. tree_allocated(t))

      ! initialize tree communication environment
      if (present(comm_env)) then
        call comm_env_mirror(comm_env, t%comm_env)
      else if (present(comm)) then
        call comm_env_mirror(comm, t%comm_env)
      else
        call comm_env_dup(MPI_COMM_lpepc, t%comm_env)
      end if

      t%npart = n
      t%npart_me = nl
      t%nleaf = 0
      t%ntwig = 0
      t%nleaf_me = 0
      t%ntwig_me = 0
      t%nbranch = 0
      t%nbranch_me = 0
      t%nbranch_max_me = 0
      t%nintmax = 0

      call timer_start(t_allocate)

      call get_number_of_interactions_per_particle(t%npart, t%nintmax)
      t%nintmax = interaction_list_length_factor * t%nintmax

      ! Space for hash table
      if (np_mult > 0) then
        maxaddress = max(30_kind_node * t%nintmax + 4_kind_node * t%npart_me, 10000_kind_node)
      else
        maxaddress = int(abs(np_mult) * 10000._8, kind_node)
      end if

      call htable_create(t%node_storage, maxaddress)
      t%nodes => t%node_storage%values

      if (maxaddress < t%npart_me + 2) then
        DEBUG_ERROR('("maxaddress = ", I0, " < npp + 2 = ", I0, ".", / , "You should increase np_mult.")', maxaddress, t%npart_me + 2)
      end if

      call tree_communicator_create(t%communicator)
      
      ! Preprocessed box sizes for each level
      allocate(t%boxlength2(0:nlev))
      t%boxlength2(0) = maxval(t%bounding_box%boxsize)**2
      do i = 1, nlev
        t%boxlength2(i) =  t%boxlength2(i-1)/4._8
      end do

      call timer_stop(t_allocate)
    end subroutine tree_create


    !>
    !> destroy a tree, freeing all memory used
    !>
    subroutine tree_destroy(t)
      use module_htable, only: htable_destroy
      use module_domains, only: decomposition_allocated, decomposition_destroy
      use module_comm_env, only: comm_env_destroy
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t !< the tree

      call pepc_status('DEALLOCATE TREE')
      DEBUG_ASSERT(tree_allocated(t))

      deallocate(t%boxlength2)
      call tree_communicator_destroy(t%communicator)
      t%nodes => null()
      call htable_destroy(t%node_storage)
      call comm_env_destroy(t%comm_env)
      if (decomposition_allocated(t%decomposition)) then
        call decomposition_destroy(t%decomposition)
      end if
    end subroutine tree_destroy


    !>
    !> returns `.true.` if resources have been allocated for the tree `t`
    !>
    function tree_allocated(t)
      use module_htable, only: htable_allocated
      implicit none

      logical :: tree_allocated
      type(t_tree), intent(in) :: t

      tree_allocated = htable_allocated(t%node_storage) .or. tree_communicator_allocated(t%communicator)
    end function tree_allocated


    !>
    !> allocates resources for the tree communicator `c` and initializes its fields
    !>
    subroutine tree_communicator_create(c)
      use, intrinsic :: iso_c_binding
      use pthreads_stuff, only: pthreads_alloc_thread
      use module_atomic_ops, only: atomic_allocate_int, atomic_store_int, critical_section_allocate
      use module_debug
      implicit none

      type(t_tree_communicator), intent(inout) :: c

      DEBUG_ASSERT(.not. tree_communicator_allocated(c))
      c%timings_comm = 0.
      
      call atomic_allocate_int(c%thread_status)
      call critical_section_allocate(c%cs_request)
      if (.not. (associated(c%thread_status) .or. associated(c%cs_request))) then
        DEBUG_ERROR(*, "atomic_allocate_int() failed!")
      end if
      call atomic_store_int(c%thread_status, TREE_COMM_THREAD_STATUS_STOPPED)

      c%sum_ships = 0
      c%sum_fetches = 0

      c%comm_thread = c_null_ptr
      c%comm_thread = pthreads_alloc_thread()
      if (.not. c_associated(c%comm_thread)) then
        DEBUG_ERROR(*, "pthreads_alloc_thread() failed!")
      end if
    end subroutine tree_communicator_create


    !>
    !> deallocates the resources of the tree communicator `c`
    !>
    subroutine tree_communicator_destroy(c)
      use, intrinsic :: iso_c_binding
      use module_atomic_ops, only: critical_section_deallocate, atomic_deallocate_int
      use pthreads_stuff, only: pthreads_free_thread
      use module_atomic_ops, only: atomic_load_int
      use module_debug
      implicit none

      type(t_tree_communicator), intent(inout) :: c

      DEBUG_ASSERT(tree_communicator_allocated(c))
      ! TODO: we could just stop the running thread, if only it was not in another module
      if (atomic_load_int(c%thread_status) == TREE_COMM_THREAD_STATUS_STARTED) then
        DEBUG_ERROR(*, "tree_communicator_destroy() called with comm thread still running!")
      end if

      call pthreads_free_thread(c%comm_thread)
      c%comm_thread = c_null_ptr
      call atomic_deallocate_int(c%thread_status)
      call critical_section_deallocate(c%cs_request)
    end subroutine tree_communicator_destroy


    !>
    !> returns `.true.` if resources have been allocated for tree communicator `c`
    !>
    function tree_communicator_allocated(c)
      use, intrinsic :: iso_c_binding
      implicit none

      logical :: tree_communicator_allocated
      type(t_tree_communicator), intent(in) :: c

      tree_communicator_allocated = c_associated(c%comm_thread)
    end function tree_communicator_allocated


    !>
    !> inserts the tree node `n` into the tree `t`.
    !>
    !> returns `.true.` if successfull, `.false.` if `n` exists in `t` allready.
    !>
    function tree_insert_node(t, n, preexisting_node)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_htable, only: htable_add
      use module_tree_node, only: tree_node_is_leaf
      use module_debug
      implicit none

      logical :: tree_insert_node
      type(t_tree), intent(inout) :: t !< Tree into which to insert the node
      type(t_tree_node), intent(in) :: n !< The tree node to insert
      integer(kind_node), optional, intent(out) :: preexisting_node !< points to preexisting node

      DEBUG_ASSERT(tree_allocated(t))
      tree_insert_node = htable_add(t%node_storage, n%key, n, preexisting_node)
      if (tree_insert_node) then
        ! everything is fine - keep count of leaves / twigs
        if (tree_node_is_leaf(n)) then
          t%nleaf =  t%nleaf + 1
          if (n%owner == t%comm_env%rank) t%nleaf_me = t%nleaf_me + 1
        else
          t%ntwig =  t%ntwig + 1
          if (n%owner == t%comm_env%rank) t%ntwig_me = t%ntwig_me + 1
        end if
      end if
    end function tree_insert_node
    
    
    !>
    !> inserts the node `n` into the tree `t` or, if a node with the same key
    !> exists allready, updates that node's entry
    !>
    !> this routine cannot be used to change a tree_node from leaf to twig or similar
    !>
    subroutine tree_insert_or_update_node(t, n, ptr)
        use module_pepc_types, only: t_tree_node, kind_node
        use module_tree_node, only: tree_node_is_leaf
        use module_debug
        implicit none

        type(t_tree), intent(inout) :: t !< Tree into which to insert the node
        type(t_tree_node), intent(in) :: n !< The tree node to insert
        integer(kind_node), optional, intent(out) :: ptr

        integer(kind_node) :: node_entry
        type(t_tree_node), pointer :: node

        DEBUG_ASSERT(tree_allocated(t))
        if (.not. tree_insert_node(t, n, node_entry)) then
          ! the node already exist --> update

          node => t%nodes(node_entry)
          ! if we change the owner from something else to 'me', we have to keep track of the leaf/twig counters
          if ((node%owner .ne. t%comm_env%rank) .and. &
            (n%owner .eq. t%comm_env%rank)) then
            if (tree_node_is_leaf(node)) then
              t%nleaf_me = t%nleaf_me + 1
            else
              t%ntwig_me = t%ntwig_me + 1
            end if
          end if

          ! copy all data
          node = n
        end if

        if (present(ptr)) then; ptr = node_entry; end if
    end subroutine


    !>
    !> checks whether a node of key `k` is contained in tree `t`
    !>
    function tree_contains_key(t, k)
      use module_htable, only: htable_contains
      use module_debug
      implicit none

      logical :: tree_contains_key
      type(t_tree), intent(in) :: t !< the tree
      integer(kind_key), intent(in) :: k !< key to look up

      DEBUG_ASSERT(tree_allocated(t))
      tree_contains_key = htable_contains(t%node_storage, k)
    end function tree_contains_key

    
    !>
    !> looks up a node for key `k` in tree `t`,
    !> returns `.true.` in case a node is found and makes `n` point to it,
    !> `.false.` is returned otherwise
    !>
    function tree_lookup_node_i(t, k, n, caller)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_htable, only: htable_lookup
      use module_debug
      implicit none

      logical :: tree_lookup_node_i

      type(t_tree), intent(in) :: t !< the tree
      integer(kind_key), intent(in) :: k !< key to look up
      integer(kind_node), intent(out) :: n !< node that is identified by `k`
      character(LEN = *), optional, intent(in) :: caller

      DEBUG_ASSERT(tree_allocated(t))
      tree_lookup_node_i = htable_lookup(t%node_storage, k, n, caller)
    end function tree_lookup_node_i


    !>
    !> looks up a node for key `k` in tree `t` and makes `n` point to it if one
    !> is found, otherwise debug information is dumped and program execution is aborted
    !>
    subroutine tree_lookup_node_critical_i(t, k, n, caller)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_htable, only: htable_lookup_critical
      use module_debug
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      integer(kind_key), intent(in) :: k !< key to look up
      integer(kind_node), intent(out) :: n !< node that is identified by `k`
      character(LEN = *), intent(in) :: caller

      DEBUG_ASSERT(tree_allocated(t))
      call htable_lookup_critical(t%node_storage, k, n, caller)
    end subroutine tree_lookup_node_critical_i


    !>
    !> looks up a node for key `k` in tree `t`,
    !> returns `.true.` in case a node is found and makes `n` point to it,
    !> `.false.` is returned otherwise
    !>
    function tree_lookup_node_p(t, k, n, caller)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_htable, only: htable_lookup
      use module_debug
      implicit none

      logical :: tree_lookup_node_p

      type(t_tree), intent(in) :: t !< the tree
      integer(kind_key), intent(in) :: k !< key to look up
      type(t_tree_node), pointer, intent(out) :: n !< node that is identified by `k`
      character(LEN = *), optional, intent(in) :: caller

      integer(kind_node) :: in

      tree_lookup_node_p = tree_lookup_node_i(t, k, in, caller)
      if (tree_lookup_node_p) then
        n => t%nodes(in)
      else
        n => null()
      end if
      
    end function tree_lookup_node_p


    !>
    !> looks up a node for key `k` in tree `t` and makes `n` point to it if one
    !> is found, otherwise debug information is dumped and program execution is aborted
    !>
    subroutine tree_lookup_node_critical_p(t, k, n, caller)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_htable, only: htable_lookup_critical
      use module_debug
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      integer(kind_key), intent(in) :: k !< key to look up
      type(t_tree_node), pointer, intent(out) :: n !< node that is identified by `k`
      character(LEN = *), intent(in) :: caller
      
      integer(kind_node) :: in

      call tree_lookup_node_critical_i(t, k, in, caller)
      n => t%nodes(in)
      
    end subroutine tree_lookup_node_critical_p


    !>
    !> Returns the parent `p` of node `n` in tree `t`.
    !>
    !> If a parent exists, `.true.` is returned and `p` points to the parent.
    !> If `n` is the root node, no parent exists and `.false.` is returned and
    !> `p` points to `null()`.
    !>
    !> @todo Currently implemented via hash table lookups which could be replaced by
    !> pointers.
    function tree_node_get_parent(t, n, p)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_tree_node, only: tree_node_is_root, NODE_INVALID
      use module_spacefilling, only: parent_key_from_key
      use module_debug
      implicit none

      logical tree_node_get_parent
      type(t_tree), intent(in) :: t
      type(t_tree_node), intent(in) :: n
      integer(kind_node), intent(out) :: p

      DEBUG_ASSERT(tree_allocated(t))
      tree_node_get_parent = .false.
      p = NODE_INVALID

      if (.not. tree_node_is_root(n)) then
        tree_node_get_parent = .true.
        call tree_lookup_node_critical(t, parent_key_from_key(n%key), p, &
          "tree_node_get_parent")
      end if
    end function tree_node_get_parent      


    !>
    !> connects `next_sibling` pointers in list of
    !> given nodes `c` and attaches the first of them to
    !> `n` via `first_child` within tree `t`
    !>
    subroutine tree_node_connect_children(t, n, c)
      use module_pepc_types, only: kind_node
      use module_debug
      use module_tree_node, only: NODE_INVALID
      implicit none

      type(t_tree), intent(inout) :: t
      integer(kind_node), intent(in) :: n
      integer(kind_node), intent(in) :: c(:)
      
      integer :: ic, nc

      nc = size(c, kind=kind(nc)); DEBUG_ASSERT(nc > 0)
      ! // DEBUG_ASSERT_MSG(all((c(2:nc)%p%key - c(1:nc - 1)%p%key) >= 1), *, "children are not arranged as expected.")
      ! // DEBUG_ASSERT_MSG(2**idim >= c(nc)%p%key - c(1)%p%key, '("= ", I3, ". Children do not all belong to the same parent.")', c(nc)%p%key - c(1)%p%key)

      t%nodes(n)%first_child = c(1)
      
      do ic = 1, nc - 1
        t%nodes(c(ic))%parent       = n
        t%nodes(c(ic))%next_sibling = c(ic + 1)
      end do
      
      t%nodes(c(nc))%parent       = n
      t%nodes(c(nc))%next_sibling = NODE_INVALID

    end subroutine tree_node_connect_children


    !>
    !> Do some quick checks on the tree structure
    !>
    function tree_check(t, callpoint)
      use module_debug
      use module_pepc_types, only: t_tree_node, kind_node
      use module_debug
      implicit none

      logical :: tree_check
      type(t_tree), intent(in) :: t !< the tree
      character(*), intent(in) :: callpoint !< caller

      integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check

      call pepc_status('CHECK TREE')

      tree_check = .true.
      nleaf_check = 0
      nleaf_me_check = 0
      ntwig_check = 0
      ntwig_me_check = 0

      DEBUG_ASSERT(tree_allocated(t))
      call tree_check_helper(t%node_root)

      if (t%nleaf /= nleaf_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# leaves in table = ',nleaf_check,' vs ',t%nleaf,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      if (t%ntwig /= ntwig_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# twigs in table = ',ntwig_check,' vs ',t%ntwig,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      if (t%nleaf_me /= nleaf_me_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# own leaves in table = ',nleaf_me_check,' vs ',t%nleaf_me,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      if (t%ntwig_me /= ntwig_me_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank, '# own twigs in table = ',ntwig_me_check,' vs ',t%ntwig_me,' accumulated', 'Fixing and continuing for now..')
        tree_check = .false.
      end if

      contains

      recursive subroutine tree_check_helper(nid)
        use module_tree_node, only: tree_node_is_leaf, tree_node_get_first_child, &
          tree_node_get_next_sibling, NODE_INVALID
        use module_pepc_types, only: t_tree_node
        implicit none

        integer(kind_node), intent(in) :: nid
        type(t_tree_node), pointer :: n

        integer(kind_node) :: s, ns

        s  = NODE_INVALID
        ns = NODE_INVALID

        n => t%nodes(nid)
        if (tree_node_is_leaf(n)) then
          nleaf_check = nleaf_check + 1
          if (n%owner == t%comm_env%rank) then
            nleaf_me_check = nleaf_me_check + 1
          end if
          return
        else
          ntwig_check = ntwig_check + 1
          if (n%owner == t%comm_env%rank) then
            ntwig_me_check = ntwig_me_check + 1
          end if
        end if

        if (tree_node_get_first_child(n, s)) then
          do
            call tree_check_helper(s)
            if (.not. tree_node_get_next_sibling(t%nodes(s), ns)) then
              exit
            end if
            s = ns
          end do
        end if
      end subroutine tree_check_helper
    end function tree_check


    !>
    !> Print tree structure from hash table to ipefile
    !>
    subroutine tree_dump(t, particles)
      use module_pepc_types
      use module_htable
      implicit none
      type(t_tree), intent(in) :: t
      type(t_particle), optional, intent(in) :: particles(:)
      
      call htable_dump(t%node_storage, particles)

    end subroutine

    !>
    !> gather statistics on the tree structure and dump them to a file
    !>
    subroutine tree_stats(t, u)
      use module_htable, only: htable_entries, htable_maxentries
      use treevars, only: np_mult
      use module_debug
      implicit none
      include 'mpif.h'

      type(t_tree), intent(in) :: t
      integer, intent(in) :: u

      integer :: i, s
      integer(kind_default) :: ierr
      integer(kind_particle), allocatable :: nparticles(:)
      integer(kind_node), allocatable :: fetches(:), ships(:)
      integer(kind_node), allocatable :: total_keys(:), tot_nleaf(:), tot_ntwig(:)
      integer(kind_particle) :: total_part
      integer(kind_node) :: max_nbranch, min_nbranch, nbranch, branch_max_global
      integer(kind_node) :: nhashentries, gmax_keys
      real, save :: part_imbal = 0.
      integer(kind_particle) ::  part_imbal_max, part_imbal_min
      integer(kind_node) :: nkeys_total

      call pepc_status('STATISTICS')
      DEBUG_ASSERT(tree_allocated(t))

      s = t%comm_env%size
      allocate(nparticles(s), fetches(s), ships(s), total_keys(s), tot_nleaf(s), tot_ntwig(s))

      ! particle distrib
      call MPI_GATHER(t%npart_me,    1, MPI_KIND_PARTICLE, nparticles, 1, MPI_KIND_PARTICLE, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%ntwig_me,    1, MPI_KIND_NODE,     tot_ntwig,  1, MPI_KIND_NODE,     0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%nleaf_me,    1, MPI_KIND_NODE,     tot_nleaf,  1, MPI_KIND_NODE,     0,  t%comm_env%comm, ierr )
      nkeys_total = t%nleaf + t%ntwig
      call MPI_GATHER(nkeys_total,                         1, MPI_KIND_NODE, total_keys, 1, MPI_KIND_NODE,  0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%communicator%sum_fetches,          1, MPI_KIND_NODE, fetches,    1, MPI_KIND_NODE,  0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%communicator%sum_ships,            1, MPI_KIND_NODE, ships,      1, MPI_KIND_NODE,  0,  t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, max_nbranch,           1, MPI_KIND_NODE, MPI_MAX,    0, t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, min_nbranch,           1, MPI_KIND_NODE, MPI_MIN,    0, t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, nbranch,               1, MPI_KIND_NODE, MPI_SUM,    0, t%comm_env%comm, ierr)
      call MPI_REDUCE(t%nbranch_max_me, branch_max_global, 1, MPI_KIND_NODE, MPI_MAX,    0, t%comm_env%comm, ierr)
      nhashentries = htable_entries(t%node_storage)
      call MPI_REDUCE(nhashentries, gmax_keys,             1, MPI_KIND_NODE, MPI_MAX,    0, t%comm_env%comm, ierr )

      part_imbal_max = MAXVAL(nparticles)
      part_imbal_min = MINVAL(nparticles)
      part_imbal = (part_imbal_max - part_imbal_min) / 1.0 / t%npart * s
      total_part = sum(nparticles)

      if (t%comm_env%first) then
        write (u,'(a20,i7,a22)') 'Tree stats for CPU ', t%comm_env%rank, ' and global statistics'
        write (u,*) '######## GENERAL DATA #####################################################################'
        write (u,'(a50,1i12)') '# procs', s
        write (u,'(a50,i12,f12.2,i12)') 'nintmax, np_mult, maxaddress: ',t%nintmax, np_mult, htable_maxentries(t%node_storage)
        write (u,'(a50,2i12)') 'npp, npart: ', t%npart_me, t%npart
        write (u,'(a50,2i12)') 'total # nparticles, N/P: ', total_part, int(t%npart/s)
        write (u,'(a50,f12.3,2i12)')   'Particle imbalance ave,min,max: ',part_imbal,part_imbal_min,part_imbal_max
        write (u,*) '######## TREE STRUCTURES ##################################################################'
        write (u,'(a50,3i12)') 'local # leaves, twigs, keys: ', t%nleaf_me, t%ntwig_me, t%nleaf_me + t%ntwig_me
        write (u,'(a50,3i12)') 'non-local # leaves, twigs, keys: ',t%nleaf - t%nleaf_me, t%ntwig - t%ntwig_me, t%nleaf + t%ntwig - t%nleaf_me - t%ntwig_me
        write (u,'(a50,3i12,f12.1,a6,i12)') 'final # leaves, twigs, keys, (max): ', t%nleaf, t%ntwig, t%nleaf + t%ntwig, &
                  real((t%nleaf + t%ntwig), kind(0._8)) / (.01 * real(htable_maxentries(t%node_storage), kind(0._8))), ' % of ', htable_maxentries(t%node_storage)
        write (u,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # keys: ',gmax_keys, real(gmax_keys, kind(0._8))/(.01 * real(htable_maxentries(t%node_storage), kind(0._8))), ' % of  ', htable_maxentries(t%node_storage)
        write (u,*) '######## BRANCHES #########################################################################'
        write (u,'(a50,3i12)') '#branches local, max_global, min_global: ', t%nbranch_me, max_nbranch, min_nbranch
        write (u,'(a50,2i12)') '#branches global sum estimated, sum actual: ', branch_max_global, nbranch
        write (u,'(a50,2i12)') 'max res.space for local branches, global br.: ', t%nbranch_max_me, branch_max_global
        write (u,*) '######## WALK-COMMUNICATION ###############################################################'
        write (u,'(a50,2i12)') 'Max # multipole fetches/ships per cpu: ',maxval(fetches), maxval(ships)
        write (u,'(a50,2i12)') 'Min # multipole fetches/ships per cpu: ',minval(fetches), minval(ships)
        write (u,'(a50,2i12)') 'Local #  multipole fetches & ships: ', t%communicator%sum_fetches, t%communicator%sum_ships
        write (u,'(a50,3i12)') '# of comm-loop iterations (tot,send,recv): ', t%communicator%comm_loop_iterations(:)
        write (u,*) '######## DETAILED DATA ####################################################################'
        write (u,'(2a/(4i10,F8.4,4i15))') '        PE     parts     nleaf     ntwig   ratio        nl_keys', &
                  '       tot_keys        fetches          ships', &
                  (i-1,nparticles(i),tot_nleaf(i),tot_ntwig(i),1.0*tot_nleaf(i)/(1.0*tot_ntwig(i)), &
                  total_keys(i)-(tot_nleaf(i)+tot_ntwig(i)),total_keys(i),fetches(i),ships(i),i=1,s)
      end if

      deallocate(nparticles, fetches, ships, total_keys, tot_nleaf, tot_ntwig)
    end subroutine tree_stats
end module module_tree
