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
!> Contains all tree specific helper routines and data fields
!>
module module_tree
    use module_htable, only: t_htable
    use module_box, only: t_box
    use module_comm_data, only: t_comm_data
    use module_domains, only: t_decomposition
    implicit none
    private

    type, public :: t_tree
      integer*8 :: npart       !< number of particles across all PE
      integer*8 :: npart_local !< number of particles on local PE
      integer*8 :: nleaf       !< number of leaves stored locally
      integer*8 :: ntwig       !< number of twigs stored locally
      !TODO: is it _me or _local? make this consistent!
      integer*8 :: nleaf_me    !< number of leaves that originated on this PE
      integer*8 :: ntwig_me    !< number of twigs that originated on this PE
      integer   :: nbranch     !< number of branch nodes
      integer   :: nbranch_me  !< number of branch nodes that originated on this PE
      integer   :: nbranch_max_me !< upper limit estimate for number of local branch nodes
      !TODO: do we need to keep this?
      integer*8 :: nintmax     !< maximum number of interactions
      !TODO: factor these out into tree comm statistics?
      integer*8 :: sum_ships   !< total number of node ships
      integer*8 :: sum_fetches !< total number of node fetches


      type(t_box) :: bounding_box
      type(t_htable) :: node_storage
      type(t_comm_data) :: comm_data
      type(t_decomposition) :: decomposition
    end type t_tree

    public tree_create
    public tree_insert_node
    public tree_update_or_insert_node
    public tree_exchange
    public tree_exchange_branches
    public tree_get_root
    public tree_lookup_node
    public tree_lookup_node_critical
    public tree_build_upwards
    public tree_build_from_particles
    public tree_destroy
    public tree_stats

    contains

    !>
    !> create a tree (allocates memory but does not fill the tree)
    !>
    subroutine tree_create(t, nl, n, comm, comm_data)
      use treevars, only: interaction_list_length_factor, MPI_COMM_lpepc, np_mult
      use module_timings
      use module_debug
      use module_htable, only: htable_create
      use module_interaction_specific, only: get_number_of_interactions_per_particle
      use module_comm_data, only: comm_data_create
      implicit none

      type(t_tree), intent(inout) :: t
      integer, intent(in) :: nl
      integer*8, intent(in) :: n
      integer, optional, intent(in) :: comm
      type(t_comm_data), optional, intent(in) :: comm_data

      integer*8 :: maxaddress

      ! initialize tree communication infrastructure
      if (present(comm_data)) then
        t%comm_data = comm_data
      else if (present(comm)) then
        call comm_data_create(t%comm_data, comm)
      else
        call comm_data_create(t%comm_data, MPI_COMM_lpepc)
      end if

      t%npart = n
      t%npart_local = nl
      t%nleaf = 0
      t%ntwig = 0
      t%nleaf_me = 0
      t%ntwig_me = 0
      t%nbranch = 0
      t%nbranch_me = 0
      t%nbranch_max_me = 0
      t%nintmax = 0
      t%sum_ships = 0
      t%sum_fetches = 0

      call timer_start(t_allocate)
      call pepc_status('ALLOCATE TREE')

      call get_number_of_interactions_per_particle(t%npart, t%nintmax)
      t%nintmax = interaction_list_length_factor * t%nintmax

      ! Space for hash table and tree arrays
      if (np_mult > 0) then
        maxaddress = max(30_8 * t%nintmax + 4_8 * t%npart_local, 10000_8)
      else
        maxaddress = int(abs(np_mult) * 10000_8, kind = 8)
      end if

      call htable_create(t%node_storage, maxaddress)

      if (maxaddress < t%npart_local + 2) then
        DEBUG_ERROR('("maxaddress = ", I0, " < npp + 2 = ", I0 ".", / , "You should increase np_mult.")', maxaddress, t%npart_local + 2)
      end if

      call timer_stop(t_allocate)
    end subroutine tree_create


    !>
    !> destroy a tree, freeing its memory
    !>
    subroutine tree_destroy(t)
      use module_htable, only: htable_destroy
      use module_domains, only: decomposition_allocated, decomposition_destroy
      use module_debug, only: pepc_status
      use module_comm_data, only: comm_data_destroy
      implicit none

      type(t_tree), intent(inout) :: t

      call pepc_status('DEALLOCATE TREE')

      t%npart = 0
      t%npart_local = 0
      t%nleaf = 0
      t%ntwig = 0
      t%nleaf_me = 0
      t%ntwig_me = 0
      t%nbranch = 0
      t%nbranch_me = 0
      t%nbranch_max_me = 0
      t%nintmax = 0
      t%sum_ships = 0
      t%sum_fetches = 0

      call htable_destroy(t%node_storage)
      call comm_data_destroy(t%comm_data)
      if (decomposition_allocated(t%decomposition)) then
        call decomposition_destroy(t%decomposition)
      end if
    end subroutine tree_destroy


    !>
    !> if an entry with tree_node%key already exists in htable, then
    !> updates htable and tree_nodes with new data
    !> otherwise: creates new entries
    !>
    !> this routine cannot be used to change a tree_node from leaf to twig or similar
    !>
    subroutine tree_update_or_insert_node(t, tree_node)
        use module_pepc_types
        use module_htable
        use module_tree_node
        implicit none
        include 'mpif.h'

        type(t_tree), intent(inout) :: t
        type(t_tree_node), intent(in) :: tree_node

        type(t_tree_node), pointer :: preexisting_node

        if (htable_lookup(t%node_storage, tree_node%key, preexisting_node)) then
          ! the htable-entry and node already exist --> update

          ! if we change the owner from someting else to 'me', we have to keep track of the leaf/twig counters
          if ((preexisting_node%owner .ne. t%comm_data%rank) .and. &
            (tree_node%owner .eq. t%comm_data%rank)) then
            if (tree_node_is_leaf(preexisting_node)) then
              t%nleaf_me = t%nleaf_me + 1
            else
              t%ntwig_me = t%ntwig_me + 1
            end if
          end if

          preexisting_node%leaves           = tree_node%leaves
          preexisting_node%flags            = tree_node%flags     
          preexisting_node%owner            = tree_node%owner
          preexisting_node%interaction_data = tree_node%interaction_data
        else
          ! create new htable and nodelist entry
          call tree_insert_node(t, tree_node)
        end if
    end subroutine


    !>
    !> Inserts a given tree node into the next free position in the tree ( -(ntwig+1) or (nleaf+1) )
    !>
    subroutine tree_insert_node(t, tree_node)
      use module_pepc_types
      use module_htable
      use module_debug
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_tree_node), intent(in) :: tree_node
      type(t_tree_node), pointer :: preexisting_entry

      if (htable_add(t%node_storage, tree_node%key, tree_node, preexisting_entry)) then
        ! anything is fine - we will have to assign a node number now
        if ( tree_node%leaves == 1 ) then ! TODO: this should be replaced by tree_node_is_leaf()
          t%nleaf =  t%nleaf + 1
          if (tree_node%owner == t%comm_data%rank) t%nleaf_me = t%nleaf_me + 1
        else if ( tree_node%leaves > 1 ) then
          ! twig
          t%ntwig =  t%ntwig + 1
          if (tree_node%owner == t%comm_data%rank) t%ntwig_me = t%ntwig_me + 1
        else
          DEBUG_ERROR(*, "Found a tree node with less than 1 leaf.")
        endif

      else
        ! entry with the same key is already existing, so we just overwrite it
        preexisting_entry%interaction_data = tree_node%interaction_data

        ! TODO: cleanup
        DEBUG_WARNING_ALL(*, "PE", t%comm_data%rank, "has found an already inserted entry while calling make_hashentry(", tree_node%key, tree_node%leaves, tree_node%flags, tree_node%owner, tree_node%level, ") - overwriting it")
      endif
    end subroutine tree_insert_node


    !>
    !> Accumulates properties of child nodes (given by keys) to parent node
    !>
    subroutine shift_nodes_up_key(t, parent, childkeys, parent_owner)
      use module_pepc_types
      use module_htable
      use module_spacefilling
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_tree_node), intent(inout) :: parent
      integer*8, intent(in) :: childkeys(:)
      integer, intent(in) :: parent_owner

      integer :: nchild, i
      type(t_tree_node) :: child_nodes(1:8)
      type(t_tree_node), pointer :: p
      integer :: childnumber(1:8)

      nchild = size(childkeys)

      do i = 1, nchild
        childnumber(i) = child_number_from_key(childkeys(i))
        call htable_lookup_critical(t%node_storage, childkeys(i), p, 'shift_nodes_up_key')
        child_nodes(i) = p
      end do

      call shift_nodes_up(parent, child_nodes(1:nchild), childnumber(1:nchild), parent_owner)
    end subroutine


    !>
    !> Accumulates properties of child nodes to parent node
    !>
    subroutine shift_nodes_up(parent, children, childnumber, parent_owner)
      use module_pepc_types
      use module_tree_node
      use module_interaction_specific, only : shift_multipoles_up
      use module_spacefilling
      use module_debug
      implicit none

        type(t_tree_node), intent(inout) :: parent
        type(t_tree_node), intent(in) :: children(:)
        integer, intent(in) :: childnumber(:)
        integer, intent(in) :: parent_owner

        integer*8 :: parent_keys(1:8)
        integer :: nchild, i, flags

        nchild = size(children)

        ! check if all keys fit to the same parent
        parent_keys(1:nchild) = parent_key_from_key(children(1:nchild)%key)

        if ( any(parent_keys(2:nchild) .ne. parent_keys(1))) then
          DEBUG_ERROR(*,"Error in shift nodes up: not all supplied children contribute to the same parent node")
        endif

        flags = 0
        do i = 1, nchild
          ! set bits for available children
          flags = ibset(flags, childnumber(i))
          ! parents of nodes with local contributions also contain local contributions
          if (btest(children(i)%flags, TREE_NODE_FLAG_HAS_LOCAL_CONTRIBUTIONS)) flags = ibset(flags, TREE_NODE_FLAG_HAS_LOCAL_CONTRIBUTIONS)
          ! parents of nodes with remote contributions also contain remote contributions
          if (btest(children(i)%flags, TREE_NODE_FLAG_HAS_REMOTE_CONTRIBUTIONS)) flags = ibset(flags, TREE_NODE_FLAG_HAS_REMOTE_CONTRIBUTIONS)
          ! parents of branch and fill nodes will also be fill nodes
          if (btest(children(i)%flags, TREE_NODE_FLAG_IS_FILL_NODE) .or. btest(children(i)%flags, TREE_NODE_FLAG_IS_BRANCH_NODE)) flags = ibset(flags, TREE_NODE_FLAG_IS_FILL_NODE)
        end do

        ! Set children_HERE flag parent since we just built it from its children
        flags =  ibset(flags, TREE_NODE_FLAG_CHILDREN_AVAILABLE)

        parent%key        = parent_keys(1)
        parent%flags      = flags
        parent%leaves     = sum(children(1:nchild)%leaves)
        parent%owner      = parent_owner
        parent%level      = level_from_key( parent_keys(1) )

        call shift_multipoles_up(parent%interaction_data, children(1:nchild)%interaction_data)
    end subroutine


    !>
    !> Exchanges tree nodes that are given in local_branch_keys with remote PEs
    !> incoming tree nodes are inserted into tree_nodes array and htable, but the
    !> tree above these nodes is not corrected
    !> outputs keys of new(and own) htable/tree_node entries in branch_keys
    !>
    subroutine tree_exchange(t, local_branch_keys, branch_keys)
        use module_pepc_types
        use module_debug, only : pepc_status
        use module_timings
        use module_htable
        use module_tree_node
        implicit none
        include 'mpif.h'

        type(t_tree), intent(inout) :: t
        integer*8, intent(in) :: local_branch_keys(:)
        integer*8, intent(inout), allocatable :: branch_keys(:)

        integer :: ierr
        integer :: i, nbranch, nbranch_sum
        type(t_tree_node), pointer :: branch_node
        type(t_tree_node), allocatable :: pack_mult(:), get_mult(:)
        integer, allocatable :: nbranches(:), igap(:)

        call timer_start(t_exchange_branches_pack)

        call pepc_status('EXCHANGE BRANCHES')

        nbranch = size(local_branch_keys)
        ! Pack local branches for shipping
        allocate(pack_mult(nbranch))
        do i = 1, nbranch
          call tree_lookup_node_critical(t, local_branch_keys(i), branch_node, 'EXCHANGE: info')

          ! additionally, we mark all local branches as branches since this is only done for remote branches during unpack (is used for fill node identification)
          branch_node%flags = ibset(branch_node%flags, TREE_NODE_FLAG_IS_BRANCH_NODE)
          pack_mult(i) = branch_node
        end do

        call timer_stop(t_exchange_branches_pack)
        call timer_start(t_exchange_branches_admininstrative)

        ! work out stride lengths so that partial arrays placed sequentially in global array
        allocate (nbranches(t%comm_data%size), igap(t%comm_data%size + 1))
        call mpi_allgather(nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, t%comm_data%comm, ierr)

        igap(1) = 0
        do i = 2, t%comm_data%size + 1
          igap(i) = igap(i - 1) + nbranches(i - 1)
        end do

        nbranch_sum = igap(t%comm_data%size + 1)

        allocate(get_mult(1:nbranch_sum), branch_keys(1:nbranch_sum))

        call timer_stop(t_exchange_branches_admininstrative)
        call timer_start(t_exchange_branches_allgatherv)

        ! actually exchange the branch nodes
        call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_tree_node, get_mult, nbranches, igap, MPI_TYPE_tree_node, &
          t%comm_data%comm, ierr)

        deallocate(pack_mult)
        deallocate(nbranches, igap)

        call timer_stop(t_exchange_branches_allgatherv)
        call timer_start(t_exchange_branches_integrate)

        ! Integrate remote branches into local tree
        do i = 1, nbranch_sum

          ! insert all remote branches into local data structures (this does *not* prepare the internal tree connections, but only copies multipole properties and creates the htable-entries)
          if (get_mult(i)%owner /= t%comm_data%rank) then
            ! delete all custom flags from incoming nodes (e.g. TREE_NODE_FLAG_CHILDREN_AVAILABLE)
            get_mult(i)%flags = IAND(get_mult(i)%flags, TREE_NODE_CHILDBYTE)
            ! after clearing all bits we have to set the flag for branches again to propagate this property upwards during global buildup
            get_mult(i)%flags = ibset(get_mult(i)%flags, TREE_NODE_FLAG_IS_BRANCH_NODE)
            ! additionally, we mark all remote branches as remote nodes (this information is propagated upwards later)
            get_mult(i)%flags = ibset(get_mult(i)%flags, TREE_NODE_FLAG_HAS_REMOTE_CONTRIBUTIONS)

            call tree_insert_node(t, get_mult(i))
          end if
          ! store branch key for later (global tree buildup)
          branch_keys(i) = get_mult(i)%key
        end do

        deallocate(get_mult)

        call timer_stop(t_exchange_branches_integrate)
    end subroutine tree_exchange


    subroutine tree_exchange_branches(t, p, bp, bk)
      use module_pepc_types
      use module_timings
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_particle), intent(in) :: p(:)
      type(t_particle), intent(in) :: bp(2)
      integer*8, allocatable :: bk(:)

      integer :: num_local_branch_keys
      integer*8, allocatable :: local_branch_keys(:)

      ! identification of branch nodes
      call timer_start(t_branches_find)
      call find_local_branches(t, p, bp, num_local_branch_keys, local_branch_keys)
      t%nbranch_me = num_local_branch_keys
      call timer_stop(t_branches_find)

      call tree_exchange(t, local_branch_keys(1:num_local_branch_keys), bk)
      t%nbranch = size(bk)

      deallocate(local_branch_keys)
      
      contains

      subroutine find_local_branches(t, p, bp, nb, bk)
        use module_pepc_types, only: t_particle
        use module_spacefilling, only: shift_key_by_level
        use module_math_tools, only: bpi
        use module_htable, only: htable_contains
        use treevars, only: idim, nlev
        use module_debug, only: pepc_status
        implicit none

        type(t_tree), intent(inout) :: t
        type(t_particle), intent(in) :: p(:)
        type(t_particle), intent(in) :: bp(2)
        integer, intent(out) :: nb
        integer*8, allocatable, intent(out) :: bk(:)

        integer :: ilevel
        integer*8 :: vld_llim, vld_rlim, L, D1, D2, pos, j, possible_branch, &
          branch_level(0:nlev), branch_level_D1(0:nlev), branch_level_D2(0:nlev)

        call pepc_status('FIND BRANCHES')
        call find_vld_limits(t, p, bp, vld_llim, vld_rlim)

        ! First find highest power in the Virtual Domain to ensure a correct branch definition
        L = bpi(vld_llim, vld_rlim)
        ! divide in two sub-domains
        ! only the last tasks must get 1 particle more
        ! because it s right limit is possibly not presentable
        D1 = L - vld_llim
        D2 = vld_rlim - L
        if (t%comm_data%last) then
          D2 = D2 + 1
        end if

        ! get estimate for number of local branches per level and total
        do ilevel = 0, nlev
          pos = idim * (nlev - ilevel)
          branch_level_D1(ilevel) = ibits(D1, pos, idim)
          branch_level_D2(ilevel) = ibits(D2, pos, idim)
          branch_level(ilevel) = branch_level_D1(ilevel) + branch_level_D2(ilevel)
        end do
        t%nbranch_max_me = sum(branch_level(:))
            
        allocate(bk(1:t%nbranch_max_me))
        
        nb = 0
        ! for D1
        pos = L
        do ilevel = 0, nlev
          do j = 1, branch_level_D1(ilevel)
            pos = pos - 2_8**(idim * (nlev - ilevel))
            possible_branch = shift_key_by_level(pos, -(nlev - ilevel))
          
            ! After local build hashtable should contain branch key
            ! otherwise branch does not exists
            ! if entry exists it is counted as branch
            ! otherwise discarded
            if (htable_contains(t%node_storage, possible_branch)) then ! entry exists
              nb = nb + 1
              bk(nb) = possible_branch
            end if
          end do
        end do

        ! for D2
        pos = L - 1
        do ilevel = 0, nlev
          do j = 1, int(branch_level_D2(ilevel))
            pos = pos + 2_8**(idim * (nlev - ilevel))
            possible_branch = shift_key_by_level(pos, -(nlev - ilevel))
            
            ! After build hashtable should contain branch key
            ! otherwise branch does not exists
            ! if entry exists it is counted as branch
            ! otherwise discarded
            if (htable_contains(t%node_storage, possible_branch)) then ! entry exists
              nb = nb + 1
              bk(nb) = possible_branch
            end if
          end do
        end do
      end subroutine find_local_branches


      subroutine find_vld_limits(t, p, bp, l, r)
        use module_pepc_types, only: t_particle
        use module_math_tools, only: bpi
        use treevars, only: idim, nlev
        implicit none

        type(t_tree), intent(in) :: t
        type(t_particle), intent(in) :: p(:)
        type(t_particle), intent(in) :: bp(2)
        integer*8, intent(out) :: l
        integer*8, intent(out) :: r

        integer*8 :: lme, rme, lb, rb
          
        ! get local key limits
        lme = p(1)%key
        rme = p(ubound(p, 1))%key

        ! get key limits for neighbor tasks
        ! and build virtual limits, so that a minimum set a branch nodes comes arround
        ! boundary tasks can access their boundary space fully only need one virtual limit
        l = 2_8**(nlev * idim)
        r = 2_8**(nlev * idim + 1) - 1

        if (.not. t%comm_data%first) then
          lb = bp(1)%key
          l = bpi(lb, lme)
        end if

        if (.not. t%comm_data%last) then
          rb = bp(2)%key
          r = bpi(rme, rb)
        end if
      end subroutine find_vld_limits

    end subroutine tree_exchange_branches


    !>
    !> Builds up the tree from the given start keys towards root
    !>  - expects, that the nodes that correspond to start_keys already
    !>    have been inserted into the tree
    !>  - missing nodes on the way towards root are added automatically
    !>  - already existing nodes are updated
    !>
    subroutine tree_build_upwards(t, keys)
        use module_debug, only : pepc_status, DBG_TREE, dbg
        use module_timings
        use module_utils
        use module_htable
        use module_spacefilling
        use module_interaction_specific
        use module_pepc_types
        implicit none

        type(t_tree), intent(inout) :: t
        integer*8, intent(in) :: keys(:)

        integer, allocatable :: key_level(:)
        integer*8, allocatable :: sub_key(:), parent_key(:)
        type(t_tree_node) :: parent_node

        integer :: numkeys, ilevel, maxlevel, nsub, groupstart, groupend, i, nparent, nuniq
        integer*8 :: current_parent_key

        call pepc_status('BUILD TOWARDS ROOT')

        numkeys = size(keys)
        allocate(key_level(numkeys))
        allocate(sub_key(0:numkeys + 1), parent_key(0:numkeys + 1))

        if (dbg(DBG_TREE)) call htable_check(t%node_storage, 'after make_branches ')

        ! get levels of branch nodes
        key_level(:) = level_from_key(keys(:))
        maxlevel = maxval( key_level(:) ) ! Find maximum level

        nparent = 0
        ! iterate through key levels
        do ilevel = maxlevel,1,-1 ! Start at finest level
            ! Collect all keys at this level
            nsub = 0
            do i = 1, numkeys
                if (key_level(i) == ilevel) then
                    nsub          = nsub + 1
                    sub_key(nsub) = keys(i)
                end if
            end do

            ! Augment list with parent keys checked at previous level
            sub_key(nsub+1:nsub+nparent) = parent_key(1:nparent)
            nsub                         = nsub + nparent

            call sort(sub_key(1:nsub)) ! Sort keys

            sub_key(0) = 0 ! remove all duplicates from the list
            nuniq = 0
            do i = 1, nsub
                if (sub_key(i) .ne. sub_key(i-1)) then
                    nuniq          = nuniq + 1
                    sub_key(nuniq) = sub_key(i)
                end if
            end do

            nsub = nuniq
            sub_key(nsub+1) = 0

            ! now, sub_key(1:nsub) contains a list of all keys (unique) at ilevel that
            ! (1) just have been inserted into the tree
            ! (2) have been modified due to additional child data
            ! tree_nodes() and htable()-entries exist for both cases
            ! --> their parents need to be created and/or updated
            i       = 1
            nparent = 0

            do while (i <= nsub)
              ! group keys with the same parent
              current_parent_key = parent_key_from_key(sub_key(i))

              groupstart = i
              do while ((parent_key_from_key(sub_key(i+1)) .eq. current_parent_key) .and. (i+1 <= nsub))
                i = i + 1
              end do
              groupend   = i

              call shift_nodes_up_key(t, parent_node, sub_key(groupstart:groupend), t%comm_data%rank)
              call tree_update_or_insert_node(t, parent_node)

              nparent             = nparent + 1
              parent_key(nparent) = current_parent_key

              ! go on with next group
              i = i + 1
            end do

        end do

        deallocate(key_level, sub_key, parent_key)
    end subroutine tree_build_upwards


    !>
    !> clears the htable and inserts all given particles at the correct level
    !> by recursively subdividing the cells if necessary
    !> after function execution, htable- and tree_node-entries for all twig- and
    !> leaf-keys exist. entries for leaves are completely valid while those
    !> for twigs have to be updated via a call to tree_build_upwards()
    !>
    !> upon exit, the key_leaf property of any particle_list entry
    !> is set appropriately
    !>
    !> warning: in contrast to tree_build_upwards(), this function may *not*
    !> be called several times to add further particles etc.
    !>
    subroutine tree_build_from_particles(t, p, bp)
      use treevars, only : idim, nlev
      use module_pepc_types
      use module_spacefilling
      use module_htable
      use module_interaction_specific, only : multipole_from_particle
      use module_timings
      use module_debug
      use module_tree_node
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_particle), intent(inout) :: p(:)
      type(t_particle), intent(in) :: bp(2)

      type t_keyidx
        integer*8 :: key
        integer :: idx
      end type

      type(t_keyidx), allocatable :: kidx(:)
      integer :: i, j

      call pepc_status('INSERT PARTICLES')

      call timer_start(t_build_pure)

      allocate(kidx(ubound(p, 1) + 2)) ! might not need those two, but what the heck
      i = 0
      if (.not. t%comm_data%first) then
        i = i + 1
        kidx(i) = t_keyidx( bp(1)%key, 0 )
      end if

      do j = 1, ubound(p, 1)
        i = i + 1
        kidx(i) = t_keyidx( p(j)%key, j )
      end do

      if (.not. t%comm_data%last) then
        i = i + 1
        kidx(i) = t_keyidx( bp(2)%key, 0 )
      end if

      p(:)%key_leaf = 0_8
      t%nleaf = i

      call timer_start(t_props_leafs) ! TODO: this timer does not make sense like this
      call insert_helper(1_8, level_from_key(1_8), kidx(1:i))
      call timer_stop(t_props_leafs)

      deallocate(kidx)
      call timer_stop(t_build_pure)

      t%nleaf_me = t%nleaf !  Retain leaves and twigs belonging to local PE
      t%ntwig_me = t%ntwig

      ! check if we did not miss any particles
      if (any(p(:)%key_leaf == 0_8)) then
        DEBUG_WARNING_ALL(*, ' did not incorporate all particles into its tree')
      end if

      contains

      recursive subroutine insert_helper(k, l, ki)
        implicit none

        integer*8, intent(in) :: k
        integer, intent(in) :: l
        type(t_keyidx), intent(in) :: ki(:)

        type(t_tree_node) :: this_node
        integer*8 :: childkey
        integer :: ip, ichild, childlevel, pstart, pend

        this_node%flags      = ibset(0, TREE_NODE_FLAG_HAS_LOCAL_CONTRIBUTIONS)
        this_node%owner      = t%comm_data%rank
        this_node%key        = k
        this_node%level      = l

        if (size(ki) == 1) then ! we have arrived at a leaf
          if (ki(1)%idx /= 0) then ! boundary particles have idx == 0, do not really need to be inserted
            this_node%leaves     = 1
            call multipole_from_particle(p(ki(1)%idx)%x, p(ki(1)%idx)%data, this_node%interaction_data)
            p(ki(1)%idx)%key_leaf = k
            if (.not. htable_add(t%node_storage, k, this_node)) then
              DEBUG_ERROR(*, "Leaf allready inserted, aborting.") ! TODO: tell me more!
            end if
          end if

        else ! more particles left, split twig
          if (l >= nlev) then ! no more levels left, cannot split
            DEBUG_WARNING_ALL('("Problem with tree: No more levels. Remaining particles 1..",I0,"  [i, key, label, x, y, z]:")', size(ki))
            DEBUG_ERROR_NO_HEADER('(I6,x,O22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, p(ki(ip)%idx)%key, p(ki(ip)%idx)%label, p(ki(ip)%idx)%x(1:3), ip=1,size(ki) ) )
          end if

          pstart = lbound(ki, dim = 1)
          do ichild = 0, 2**idim - 1
            childkey   = child_key_from_parent_key(k, ichild)
            childlevel = l + 1

            pend = pstart - 1
            do while (pend < ubound(ki, dim = 1))
              if (.not. is_ancestor_of_particle(ki(pend + 1)%key, childkey, childlevel)) exit
              pend = pend + 1
            end do
            
            if (pend >= pstart) then
              call insert_helper(childkey, childlevel, ki(pstart:pend))
              ! TODO: could shift up childrens' properties here
            end if

            pstart = pend + 1
          end do

          if (.not. pend == ubound(ki, dim = 1)) then
            DEBUG_WARNING_ALL('("Problem with tree: Could not distribute particles among children of ", I0, ". Remaining particles ",I0,"..",I0,"  [i, key, label, x, y, z]:")', k, pend + 1, ubound(ki, dim = 1) )
            DEBUG_ERROR_NO_HEADER('(I6,x,I22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, p(ki(ip)%idx)%key, p(ki(ip)%idx)%label, p(ki(ip)%idx)%x(1:3), ip=pend + 1, ubound(ki, dim = 1) ) )
          end if

          if (.not. htable_add(t%node_storage, k, this_node)) then
            DEBUG_ERROR(*, "Twig allready inserted, aborting.") ! TODO: tell me more!
          end if
          t%ntwig = t%ntwig + 1
        end if

      end subroutine insert_helper
    end subroutine tree_build_from_particles


    subroutine tree_get_root(t, r, caller)
      use module_pepc_types
      implicit none

      type(t_tree), intent(in) :: t
      type(t_tree_node), pointer, intent(out) :: r
      character(len = *), optional, intent(in) :: caller

      if (present(caller)) then
        call tree_lookup_node_critical(t, 1_8, r, caller)
      else 
        call tree_lookup_node_critical(t, 1_8, r, 'tree_get_root')
      end if
    end subroutine tree_get_root


    function tree_lookup_node(t, k, n)
      use module_pepc_types
      use module_htable
      implicit none

      logical :: tree_lookup_node

      type(t_tree), intent(in) :: t
      integer*8, intent(in) :: k
      type(t_tree_node), pointer, intent(out) :: n

      tree_lookup_node = htable_lookup(t%node_storage, k, n)
    end function tree_lookup_node


    subroutine tree_lookup_node_critical(t, k, n, caller)
      use module_pepc_types
      use module_htable
      implicit none

      type(t_tree), intent(in) :: t
      integer*8, intent(in) :: k
      type(t_tree_node), pointer, intent(out) :: n
      character(LEN = *), intent(in) :: caller

      call htable_lookup_critical(t%node_storage, k, n, caller)
    end subroutine tree_lookup_node_critical


    !>
    !> gather statistics on the tree structure (and also communication and also
    !> the last traversal) and dump them to a file
    !>
    ! TODO: split this up!
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

      s = t%comm_data%size
      allocate(nparticles(s), fetches(s), ships(s), total_keys(s), tot_nleaf(s), &
        tot_ntwig(s), num_interactions(s), num_mac_evaluations(s))

      ! particle distrib
      call MPI_GATHER(t%npart_local, 1, MPI_INTEGER, nparticles, 1, MPI_INTEGER, 0,  t%comm_data%comm, ierr )
      call MPI_GATHER(t%ntwig_me,    1, MPI_INTEGER8, tot_ntwig,  1, MPI_INTEGER8, 0,  t%comm_data%comm, ierr )
      call MPI_GATHER(t%nleaf_me,    1, MPI_INTEGER8, tot_nleaf,  1, MPI_INTEGER8, 0,  t%comm_data%comm, ierr )
      nkeys_total = t%nleaf + t%ntwig
      call MPI_GATHER(nkeys_total,   1, MPI_INTEGER8, total_keys, 1, MPI_INTEGER8, 0,  t%comm_data%comm, ierr )
      call MPI_GATHER(t%sum_fetches, 1, MPI_INTEGER8, fetches,    1, MPI_INTEGER8, 0,  t%comm_data%comm, ierr )
      call MPI_GATHER(t%sum_ships,   1, MPI_INTEGER8, ships,      1, MPI_INTEGER8, 0,  t%comm_data%comm, ierr )
      call MPI_GATHER(interactions_local,    1, MPI_REAL8, num_interactions,      1, MPI_REAL8,   0,  t%comm_data%comm, ierr )
      call MPI_GATHER(mac_evaluations_local, 1, MPI_REAL8, num_mac_evaluations,   1, MPI_REAL8,   0,  t%comm_data%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, max_nbranch,     1, MPI_INTEGER, MPI_MAX, 0, t%comm_data%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, min_nbranch,     1, MPI_INTEGER, MPI_MIN, 0, t%comm_data%comm, ierr )
      call MPI_REDUCE(t%nbranch_me, nbranch,         1, MPI_INTEGER, MPI_SUM, 0, t%comm_data%comm, ierr)
      call MPI_REDUCE(t%nbranch_max_me, branch_max_global, 1, MPI_INTEGER, MPI_MAX, 0, t%comm_data%comm, ierr)
      nhashentries = htable_entries(t%node_storage)
      call MPI_REDUCE(nhashentries, gmax_keys, 1, MPI_INTEGER8, MPI_MAX, 0, t%comm_data%comm, ierr )

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

      if (t%comm_data%first) then
        if (firstcall) then
          call create_directory("stats")
          firstcall = .false.
        end if

        write(cfile,'("stats/stats.",i6.6)') timestamp
      
        open (60,file=trim(cfile))

        write (60,'(a20,i7,a22)') 'Tree stats for CPU ', t%comm_data%rank, ' and global statistics'
        write (60,*) '######## GENERAL DATA #####################################################################'
        write (60,'(a50,1i12)') '# procs', s
        write (60,'(a50,i12,f12.2,i12)') 'nintmax, np_mult, maxaddress: ',t%nintmax, np_mult, htable_maxentries(t%node_storage)
        write (60,'(a50,2i12)') 'npp, npart: ', t%npart_local, t%npart
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

      !call tree_walk_statistics(60, t%comm_data%first)

      if (t%comm_data%first) then
        write (60,*) '######## WALK-COMMUNICATION ###############################################################'
        write (60,'(a50,2i12)') 'Max # multipole fetches/ships per cpu: ',maxval(fetches), maxval(ships)
        write (60,'(a50,2i12)') 'Min # multipole fetches/ships per cpu: ',minval(fetches), minval(ships)
        write (60,'(a50,2i12)') 'Local #  multipole fetches & ships: ', t%sum_fetches, t%sum_ships
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
