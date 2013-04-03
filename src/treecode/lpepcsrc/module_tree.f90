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
    use module_atomic_ops, only: t_atomic_int
    use module_pepc_types, only: t_tree_node
    use, intrinsic :: iso_c_binding
    implicit none
    private

    !> data type for communicator request queue
    type, public :: t_request_queue_entry
      integer*8 :: key
      integer   :: owner
      type(t_tree_node), pointer :: node
    end type

    integer, public, parameter :: TREE_COMM_ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
    integer, public, parameter :: TREE_COMM_REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue

    !> data type for tree communicator
    type, public :: t_tree_communicator
      ! request queue
      ! TODO: this used to be volatile, but that is not allowed for derived type components
      type(t_request_queue_entry) :: req_queue(TREE_COMM_REQUEST_QUEUE_LENGTH)
      integer :: req_queue_top !< position of queue top in array; pushed towards bottom by communicator only when sending
      type(t_atomic_int), pointer :: req_queue_bottom !< position of queue bottom in array; pushed away from top by tree users

      ! counters and timers
      integer*8 :: comm_loop_iterations(3) !< number of comm loop iterations (total, sending, receiving)
      real*8 :: timings_comm(3) !< array for storing internal timing information
      integer :: request_balance !< total (#requests - #answers), should be zero after complete traversal
      integer*8 :: sum_ships   !< total number of node ships
      integer*8 :: sum_fetches !< total number of node fetches

      ! thread data
      type(c_ptr) :: comm_thread
      logical :: comm_thread_running
      logical :: comm_thread_stopping
      logical :: comm_thread_stop_requested
    end type t_tree_communicator

    !>
    !> A derived type representing a distributed hashed octree over a collection
    !> of particles.
    !>
    type, public :: t_tree
      integer*8 :: npart       !< number of particles across all ranks
      integer*8 :: npart_me    !< number of particles on local rank

      integer*8 :: nleaf       !< number of leaves stored locally
      integer*8 :: ntwig       !< number of twigs stored locally
      integer*8 :: nleaf_me    !< number of leaves that originated on this rank
      integer*8 :: ntwig_me    !< number of twigs that originated on this rank
      
      integer   :: nbranch     !< number of branch nodes in tree
      integer   :: nbranch_me  !< number of branch nodes that originated on this rank
      integer   :: nbranch_max_me !< upper limit estimate for number of local branch nodes

      !TODO: do we need to keep this?
      integer*8 :: nintmax     !< maximum number of interactions
      
      type(t_box) :: bounding_box               !< bounding box enclosing all particles contained in the tree
      type(t_htable) :: node_storage            !< hash table in which tree nodes are stored for rapid retrieval
      type(t_comm_env) :: comm_env              !< communication environment over which the tree is distributed
      type(t_decomposition) :: decomposition    !< permutation of particles inserted into the tree
      type(t_tree_communicator) :: communicator !< associated communicator structure
    end type t_tree

    public tree_create
    public tree_allocated
    public tree_insert_node
    public tree_insert_or_update_node
    public tree_contains_key
    public tree_lookup_root
    public tree_lookup_node
    public tree_lookup_node_critical
    public tree_node_get_first_child
    public tree_node_get_next_sibling
    public tree_node_get_parent
    public tree_check
    public tree_stats
    public tree_destroy

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
      use treevars, only: interaction_list_length_factor, MPI_COMM_lpepc, np_mult
      use module_timings
      use module_debug
      use module_htable, only: htable_create
      use module_interaction_specific, only: get_number_of_interactions_per_particle
      use module_comm_env, only: comm_env_dup, comm_env_mirror
      implicit none

      type(t_tree), intent(inout) :: t !< The tree
      integer, intent(in) :: nl !< Number of local particles to be inserted into the tree
      integer*8, intent(in) :: n !< Total number of particles across communication ranks
      integer, optional, intent(in) :: comm !< An MPI communicator
      type(t_comm_env), optional, intent(in) :: comm_env !< A communication environment

      integer*8 :: maxaddress

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
      call pepc_status('ALLOCATE TREE')

      call get_number_of_interactions_per_particle(t%npart, t%nintmax)
      t%nintmax = interaction_list_length_factor * t%nintmax

      ! Space for hash table
      if (np_mult > 0) then
        maxaddress = max(30_8 * t%nintmax + 4_8 * t%npart_me, 10000_8)
      else
        maxaddress = int(abs(np_mult) * 10000_8, kind = 8)
      end if

      call htable_create(t%node_storage, maxaddress)

      if (maxaddress < t%npart_me + 2) then
        DEBUG_ERROR('("maxaddress = ", I0, " < npp + 2 = ", I0, ".", / , "You should increase np_mult.")', maxaddress, t%npart_me + 2)
      end if

      call tree_communicator_create(t%communicator)

      call timer_stop(t_allocate)
    end subroutine tree_create


    !>
    !> destroy a tree, freeing all memory used
    !>
    subroutine tree_destroy(t)
      use module_htable, only: htable_destroy
      use module_domains, only: decomposition_allocated, decomposition_destroy
      use module_debug, only: pepc_status
      use module_comm_env, only: comm_env_destroy
      implicit none

      type(t_tree), intent(inout) :: t !< the tree

      call pepc_status('DEALLOCATE TREE')

      call tree_communicator_destroy(t%communicator)
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

      tree_allocated = htable_allocated(t%node_storage)
    end function tree_allocated


    !>
    !> allocates resouces for the tree communicator `c` and initializes its fields
    !>
    subroutine tree_communicator_create(c)
      use, intrinsic :: iso_c_binding
      use module_pepc_types, only: mpi_type_tree_node
      use pthreads_stuff, only: pthreads_alloc_thread
      use module_atomic_ops, only: atomic_allocate_int, atomic_store_int
      use module_debug
      implicit none

      type(t_tree_communicator), intent(inout) :: c

      c%timings_comm = 0.
      
      call atomic_allocate_int(c%req_queue_bottom)
      if (.not. associated(c%req_queue_bottom)) then
        DEBUG_ERROR(*, "atomic_allocate_int() failed!")
      end if
      call atomic_store_int(c%req_queue_bottom, 0)

      c%req_queue_top =  0
      c%request_balance =  0
      c%req_queue(:)%owner = -1 ! used in send_requests() to ensure that only completely stored entries are sent form the list
      c%sum_ships = 0
      c%sum_fetches = 0

      c%comm_thread = c_null_ptr
      c%comm_thread = pthreads_alloc_thread()
      if (.not. c_associated(c%comm_thread)) then
        DEBUG_ERROR(*, "pthreads_alloc_thread() failed!")
      end if
      c%comm_thread_running = .false.
      c%comm_thread_stopping = .false.
      c%comm_thread_stop_requested = .false.
    end subroutine tree_communicator_create


    !>
    !> deallocates the resources of the tree communicator `c`
    !>
    subroutine tree_communicator_destroy(c)
      use module_atomic_ops, only: atomic_deallocate_int 
      use pthreads_stuff, only: pthreads_free_thread
      use module_debug
      implicit none

      type(t_tree_communicator), intent(inout) :: c

      ! TODO: we could just stop the running thread, if only it was not in another module
      if (c%comm_thread_running) then
        DEBUG_ERROR(*, "tree_communicator_destroy() called with comm thread still running!")
      end if

      call pthreads_free_thread(c%comm_thread)
      call atomic_deallocate_int(c%req_queue_bottom)
    end subroutine tree_communicator_destroy


    !>
    !> inserts the tree node `n` into the tree `t`.
    !>
    !> returns `.true.` if successfull, `.false.` if `n` exists in `t` allready.
    !>
    function tree_insert_node(t, n, preexisting_node)
      use module_pepc_types, only: t_tree_node
      use module_htable, only: htable_add
      use module_tree_node, only: tree_node_is_leaf
      use module_debug
      implicit none

      logical :: tree_insert_node
      type(t_tree), intent(inout) :: t !< Tree into which to insert the node
      type(t_tree_node), intent(in) :: n !< The tree node to insert
      type(t_tree_node), optional, pointer, intent(out) :: preexisting_node !< points to preexisting node

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
    subroutine tree_insert_or_update_node(t, n)
        use module_pepc_types, only: t_tree_node
        use module_tree_node, only: tree_node_is_leaf
        implicit none
        include 'mpif.h'

        type(t_tree), intent(inout) :: t !< Tree into which to insert the node
        type(t_tree_node), intent(in) :: n !< The tree node to insert

        type(t_tree_node), pointer :: preexisting_node

        if (.not. tree_insert_node(t, n, preexisting_node)) then
          ! the node already exist --> update

          ! if we change the owner from someting else to 'me', we have to keep track of the leaf/twig counters
          if ((preexisting_node%owner .ne. t%comm_env%rank) .and. &
            (n%owner .eq. t%comm_env%rank)) then
            if (tree_node_is_leaf(preexisting_node)) then
              t%nleaf_me = t%nleaf_me + 1
            else
              t%ntwig_me = t%ntwig_me + 1
            end if
          end if

          preexisting_node%leaves           = n%leaves
          preexisting_node%flags            = n%flags     
          preexisting_node%owner            = n%owner
          preexisting_node%interaction_data = n%interaction_data
        end if
    end subroutine


    !>
    !> checks whether a node of key `k` is contained in tree `t`
    !>
    function tree_contains_key(t, k)
      use module_htable, only: htable_contains
      implicit none

      logical :: tree_contains_key
      type(t_tree), intent(in) :: t !< the tree
      integer*8, intent(in) :: k !< key to look up

      tree_contains_key = htable_contains(t%node_storage, k)
    end function tree_contains_key


    !>
    !> look up the root node `r` of tree `t`
    !>
    subroutine tree_lookup_root(t, r, caller)
      use module_pepc_types, only: t_tree_node
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      type(t_tree_node), pointer, intent(out) :: r !< root node
      character(len = *), optional, intent(in) :: caller !< identifies the caller in case an error message is printed

      if (present(caller)) then
        call tree_lookup_node_critical(t, 1_8, r, caller)
      else 
        call tree_lookup_node_critical(t, 1_8, r, 'tree_lookup_root')
      end if
    end subroutine tree_lookup_root

    
    !>
    !> looks up a node for key `k` in tree `t`,
    !> returns `.true.` in case a node is found and makes `n` point to it,
    !> `.false.` is returned otherwise
    !>
    function tree_lookup_node(t, k, n)
      use module_pepc_types, only: t_tree_node
      use module_htable, only: htable_lookup
      implicit none

      logical :: tree_lookup_node

      type(t_tree), intent(in) :: t !< the tree
      integer*8, intent(in) :: k !< key to look up
      type(t_tree_node), pointer, intent(out) :: n !< node that is identified by `k`

      tree_lookup_node = htable_lookup(t%node_storage, k, n)
    end function tree_lookup_node


    !>
    !> looks up a node for key `k` in tree `t` and makes `n` point to it if one
    !> is found, otherwise debug information is dumped and program execution is aborted
    !>
    subroutine tree_lookup_node_critical(t, k, n, caller)
      use module_pepc_types, only: t_tree_node
      use module_htable, only: htable_lookup_critical
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      integer*8, intent(in) :: k !< key to look up
      type(t_tree_node), pointer, intent(out) :: n !< node that is identified by `k`
      character(LEN = *), intent(in) :: caller

      call htable_lookup_critical(t%node_storage, k, n, caller)
    end subroutine tree_lookup_node_critical


    !>
    !> Returns the first child of node `p` in tree `t`.
    !>
    !> If `p` has children that are locally available, returns `.true.` and `fc`
    !> points to the child.
    !> Otherwise, `.false.` is returned and `fc` points to `null()`.
    !>
    !> @todo Currently implemented via hash table lookups which could be replaced by
    !> pointers.
    function tree_node_get_first_child(t, p, fc)
      use module_pepc_types, only: t_tree_node
      use module_tree_node, only: tree_node_is_leaf, &
        tree_node_children_available, tree_node_has_child
      use module_spacefilling, only: child_key_from_parent_key
      use treevars, only: idim
      implicit none

      logical :: tree_node_get_first_child
      type(t_tree), intent(in) :: t
      type(t_tree_node), intent(in) :: p
      type(t_tree_node), pointer, intent(out) :: fc

      integer :: ic

      tree_node_get_first_child = .false.
      fc => null()

      if (.not. tree_node_is_leaf(p) .and. tree_node_children_available(p)) then
        do ic = 0, 2**idim - 1
          if (tree_node_has_child(p, ic)) then
            tree_node_get_first_child = .true.
            call tree_lookup_node_critical(t, child_key_from_parent_key(p%key, ic), &
              fc, "tree_node_get_first_child")
            exit
          end if
        end do
      end if
    end function tree_node_get_first_child


    !>
    !> Returns the next sibling of node `p` in tree `t`.
    !>
    !> If there is a next sibling to `n` in `t` returns `.true.` and `s` points
    !> to the sibling node.
    !> Otherwise `.false.` is returned and `s` points to `null()`.
    !>
    !> In this context, "next" is defined by the ordering of the node keys.
    !>
    !> @todo Currently implemented via hash table lookups which could be replaced by
    !> pointers.
    function tree_node_get_next_sibling(t, n, s)
      use module_pepc_types, only: t_tree_node
      use module_tree_node, only: tree_node_has_child
      use module_spacefilling, only: parent_key_from_key, &
        child_key_from_parent_key, child_number_from_key
      use treevars, only: idim
      implicit none

      logical :: tree_node_get_next_sibling
      type(t_tree), intent(in) :: t
      type(t_tree_node), intent(in) :: n
      type(t_tree_node), pointer, intent(out) :: s

      type(t_tree_node), pointer :: p
      integer :: is

      tree_node_get_next_sibling = .false.
      s => null()

      if (tree_node_get_parent(t, n, p)) then
        do is = child_number_from_key(n%key) + 1, 2**idim - 1    
          if (tree_node_has_child(p, is)) then
            tree_node_get_next_sibling = .true.
            call tree_lookup_node_critical(t, child_key_from_parent_key(p%key, is), &
              s, "tree_node_get_next_sibling: sibling")
            exit
          end if
        end do
      end if
    end function tree_node_get_next_sibling


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
      use module_pepc_types, only: t_tree_node
      use module_tree_node, only: tree_node_is_root
      use module_spacefilling, only: parent_key_from_key
      implicit none

      logical tree_node_get_parent
      type(t_tree), intent(in) :: t
      type(t_tree_node), intent(in) :: n
      type(t_tree_node), pointer, intent(out) :: p

      tree_node_get_parent = .false.
      p => null()

      if (.not. tree_node_is_root(n)) then
        tree_node_get_parent = .true.
        call tree_lookup_node_critical(t, parent_key_from_key(n%key), p, &
          "tree_node_get_parent")
      end if
    end function tree_node_get_parent      


    !>
    !> Do some quick checks on the tree structure
    !>
    subroutine tree_check(t, callpoint)
      use treevars, only: me
      use module_debug
      use module_pepc_types, only: t_tree_node
      use module_htable, only: htable_dump
      implicit none

      type(t_tree), intent(inout) :: t !< the tree
      character(*), intent(in) :: callpoint !< caller

      type(t_tree_node), pointer :: r
      integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check
      logical :: error

      call pepc_status('CHECK TREE')

      error = .false.
      nleaf_check = 0
      nleaf_me_check = 0
      ntwig_check = 0
      ntwig_me_check = 0

      call tree_lookup_root(t, r)
      call tree_check_helper(r)

      if (t%nleaf /= nleaf_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank,
        '# leaves in table = ',nleaf_check,' vs ',t%nleaf,' accumulated',
        'Fixing and continuing for now..')
        t%nleaf = nleaf_check ! TODO: does this really help anyone?
        error = .true.
      end if

      if (t%ntwig /= ntwig_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank,
        '# twigs in table = ',ntwig_check,' vs ',t%ntwig,' accumulated',
        'Fixing and continuing for now..')
        t%ntwig = ntwig_check
        error = .true.
      end if

      if (t%nleaf_me /= nleaf_me_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank,
        '# own leaves in table = ',nleaf_me_check,' vs ',t%nleaf_me,' accumulated',
        'Fixing and continuing for now..')
        t%nleaf_me = nleaf_me_check
        error = .true.
      end if

      if (t%ntwig_me /= ntwig_me_check) then
        DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',t%comm_env%rank,
        '# own twigs in table = ',ntwig_me_check,' vs ',t%ntwig_me,' accumulated',
        'Fixing and continuing for now..')
        t%ntwig_me = ntwig_me_check
        error = .true.
      end if

      if (error) then
        call htable_dump(t%node_storage)
      end if

      contains

      recursive subroutine tree_check_helper(n)
        use module_tree_node, only: tree_node_is_leaf
        implicit none

        type(t_tree_node), intent(in) :: n

        type(t_tree_node), pointer :: s, ns

        s => null()
        ns => null()

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

        if (tree_node_get_first_child(t, n, s)) then
          do
            call tree_check_helper(s)
            if (.not. tree_node_get_next_sibling(t, s, ns)) then
              exit
            end if
            s => ns
          end do
        end if
      end subroutine tree_check_helper
    end subroutine tree_check


    !>
    !> gather statistics on the tree structure (and also communication and also
    !> the last traversal) and dump them to a file
    !>
    !> @todo split this up!
    subroutine tree_stats(t, timestamp)
      use treevars
      !use module_walk, only : tree_walk_statistics
      use module_debug, only : pepc_status
      use module_htable, only: htable_entries, htable_maxentries
      use module_utils, only : create_directory
      implicit none
      include 'mpif.h'

      type(t_tree), intent(in) :: t
      integer, intent(in) :: timestamp

      integer :: i, s, ierr
      integer, allocatable :: nparticles(:)
      integer*8, allocatable :: fetches(:), ships(:), total_keys(:), tot_nleaf(:), tot_ntwig(:)
      real*8, allocatable ::  num_interactions(:), num_mac_evaluations(:)  ! Load balance arrays
      character*40 :: cfile
      integer :: total_part, max_nbranch, min_nbranch, nbranch, branch_max_global
      integer*8 :: nhashentries, gmax_keys
      real*8 :: average_interactions, average_mac_evaluations, total_interactions, total_mac_evaluations, max_interactions, max_mac_evaluations
      real, save :: part_imbal = 0.
      real*8, save :: work_imbal = 0.
      real*8 :: work_imbal_max, work_imbal_min  ! load stats
      integer ::  part_imbal_max, part_imbal_min
      integer*8 :: nkeys_total
      logical, save :: firstcall = .true.

      call pepc_status('STATISTICS')

      s = t%comm_env%size
      allocate(nparticles(s), fetches(s), ships(s), total_keys(s), tot_nleaf(s), &
        tot_ntwig(s), num_interactions(s), num_mac_evaluations(s))

      ! particle distrib
      call MPI_GATHER(t%npart_me, 1, MPI_INTEGER, nparticles, 1, MPI_INTEGER, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%ntwig_me,    1, MPI_INTEGER8, tot_ntwig,  1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%nleaf_me,    1, MPI_INTEGER8, tot_nleaf,  1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      nkeys_total = t%nleaf + t%ntwig
      call MPI_GATHER(nkeys_total,   1, MPI_INTEGER8, total_keys, 1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%communicator%sum_fetches, 1, MPI_INTEGER8, fetches,    1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%communicator%sum_ships,   1, MPI_INTEGER8, ships,      1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(interactions_local,    1, MPI_REAL8, num_interactions,      1, MPI_REAL8,   0,  t%comm_env%comm, ierr )
      call MPI_GATHER(mac_evaluations_local, 1, MPI_REAL8, num_mac_evaluations,   1, MPI_REAL8,   0,  t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, max_nbranch,     1, MPI_INTEGER, MPI_MAX, 0, t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, min_nbranch,     1, MPI_INTEGER, MPI_MIN, 0, t%comm_env%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, nbranch,         1, MPI_INTEGER, MPI_SUM, 0, t%comm_env%comm, ierr)
      call MPI_REDUCE(t%nbranch_max_me, branch_max_global, 1, MPI_INTEGER, MPI_MAX, 0, t%comm_env%comm, ierr)
      nhashentries = htable_entries(t%node_storage)
      call MPI_REDUCE(nhashentries, gmax_keys, 1, MPI_INTEGER8, MPI_MAX, 0, t%comm_env%comm, ierr )

      part_imbal_max = MAXVAL(nparticles)
      part_imbal_min = MINVAL(nparticles)
      part_imbal = (part_imbal_max-part_imbal_min)/1.0/t%npart*s

      total_interactions       = SUM(num_interactions)
      total_mac_evaluations    = SUM(num_mac_evaluations)
      max_interactions         = MAXVAL(num_interactions)
      max_mac_evaluations      = MAXVAL(num_mac_evaluations)
      average_interactions     = total_interactions    / s
      average_mac_evaluations  = total_mac_evaluations / s
      work_imbal_max = max_interactions/average_interactions
      work_imbal_min = MINVAL(num_interactions)/average_interactions
      work_imbal = 0.
      do i = 1, s
        work_imbal = work_imbal + abs(num_interactions(i) - average_interactions)/average_interactions/s
      end do

      total_part = sum(nparticles)

      if (t%comm_env%first) then
        if (firstcall) then
          call create_directory("stats")
          firstcall = .false.
        end if

        write(cfile,'("stats/stats.",i6.6)') timestamp
      
        open (60,file=trim(cfile))

        write (60,'(a20,i7,a22)') 'Tree stats for CPU ', t%comm_env%rank, ' and global statistics'
        write (60,*) '######## GENERAL DATA #####################################################################'
        write (60,'(a50,1i12)') '# procs', s
        write (60,'(a50,i12,f12.2,i12)') 'nintmax, np_mult, maxaddress: ',t%nintmax, np_mult, htable_maxentries(t%node_storage)
        write (60,'(a50,2i12)') 'npp, npart: ', t%npart_me, t%npart
        write (60,'(a50,2i12)') 'total # nparticles, N/P: ', total_part, int(t%npart/s)
        write (60,*) '######## TREE STRUCTURES ##################################################################'
        write (60,'(a50,3i12)') 'local # leaves, twigs, keys: ', t%nleaf_me, t%ntwig_me, t%nleaf_me + t%ntwig_me
        write (60,'(a50,3i12)') 'non-local # leaves, twigs, keys: ',t%nleaf - t%nleaf_me, t%ntwig - t%ntwig_me, t%nleaf + t%ntwig - t%nleaf_me - t%ntwig_me
        write (60,'(a50,3i12,f12.1,a6,i12)') 'final # leaves, twigs, keys, (max): ', t%nleaf, t%ntwig, t%nleaf + t%ntwig, &
                  (t%nleaf + t%ntwig) / (.01 * htable_maxentries(t%node_storage)), ' % of ', htable_maxentries(t%node_storage)
        write (60,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # keys: ',gmax_keys, gmax_keys/(.01 * htable_maxentries(t%node_storage)), ' % of  ', htable_maxentries(t%node_storage)
        write (60,*) '######## BRANCHES #########################################################################'
        write (60,'(a50,3i12)') '#branches local, max_global, min_global: ', t%nbranch_me, max_nbranch, min_nbranch
        write (60,'(a50,2i12)') '#branches global sum estimated, sum actual: ', branch_max_global, nbranch
        write (60,'(a50,2i12)') 'max res.space for local branches, global br.: ', t%nbranch_max_me, branch_max_global
        write (60,*) '######## TREE TRAVERSAL MODULE ############################################################'
      end if

      ! TODO: cannot call this from here!
      !call tree_walk_statistics(60, t%comm_env%first)

      if (t%comm_env%first) then
        write (60,*) '######## WALK-COMMUNICATION ###############################################################'
        write (60,'(a50,2i12)') 'Max # multipole fetches/ships per cpu: ',maxval(fetches), maxval(ships)
        write (60,'(a50,2i12)') 'Min # multipole fetches/ships per cpu: ',minval(fetches), minval(ships)
        write (60,'(a50,2i12)') 'Local #  multipole fetches & ships: ', t%communicator%sum_fetches, t%communicator%sum_ships
        write (60,*) '######## WORKLOAD AND WALK ################################################################'
        write (60,'(a50,3e12.4)')       'total/ave/max_local # interactions(work): ', total_interactions, average_interactions, max_interactions
        write (60,'(a50,3e12.4)')       'total/ave/max_local # mac evaluations: ', total_mac_evaluations, average_mac_evaluations, max_mac_evaluations
        write (60,'(a50,3f12.3)')       'Load imbalance percent,min,max: ',work_imbal,work_imbal_min,work_imbal_max
        write (60,'(a50,f12.3,2i12)')   'Particle imbalance ave,min,max: ',part_imbal,part_imbal_min,part_imbal_max
        write (60,*) '###########################################################################################'
        write (60,*) '######## DETAILED DATA ####################################################################'
        write (60,'(2a/(4i10,F8.4,6i15,F8.4))') '         PE     parts    nleaf     ntwig   ratio    nl_keys', &
                  '   tot_keys   fetches    ships    #interactions(work)   #mac_evals   rel.work*ncpu', &
                  (i-1,nparticles(i),tot_nleaf(i),tot_ntwig(i),1.0*tot_nleaf(i)/(1.0*tot_ntwig(i)), &
                  total_keys(i)-(tot_nleaf(i)+tot_ntwig(i)),total_keys(i),fetches(i),ships(i),int(num_interactions(i)),int(num_mac_evaluations(i)),&
                  num_interactions(i)/average_interactions,i=1,s)
        close(60)

      end if

      deallocate(nparticles, fetches, ships, total_keys, tot_nleaf, tot_ntwig, num_interactions, num_mac_evaluations)
    end subroutine tree_stats
end module module_tree
