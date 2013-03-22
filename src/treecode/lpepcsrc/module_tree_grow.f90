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
!> Contains procedures to construct a tree from a collection of particles.
!>
module module_tree_grow
  implicit none
  private

  public tree_grow

  contains

  !>
  !> Builds the tree from the given particles, redistributes particles
  !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
  !>
  subroutine tree_grow(t, n, p, npl)
    use module_pepc_types, only: t_particle, t_tree_node
    use module_htable, only: htable_dump
    use module_timings
    use module_tree, only: t_tree, tree_create, tree_lookup_root
    use module_domains, only: domain_decompose
    use module_debug
    use module_comm_env, only: t_comm_env, comm_env_dup
    use module_box, only: box_create
    use treevars, only: MPI_COMM_lpepc
    use module_spacefilling, only: compute_particle_keys
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    integer*8, intent(in) :: n !< total number of simulation particles (across all MPI ranks)
    type(t_particle), allocatable, intent(inout) :: p(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
    integer, optional, intent(in) :: npl !< number of valid entries in p (local particles)

    type(t_particle) :: bp(2) 
    type(t_tree_node), pointer :: root_node
    type(t_comm_env) :: tree_comm_env
    integer*8, allocatable :: branch_keys(:)

    call pepc_status('GROW TREE')
    !call MPI_BARRIER( MPI_COMM_lpepc, ierr)  ! Wait for everyone to catch up
    call timer_start(t_all)
    call timer_start(t_fields_tree)

    call comm_env_dup(MPI_COMM_lpepc, tree_comm_env)

    ! determine the bounding box of the particle configuration
    call box_create(t%bounding_box, p, tree_comm_env)
    call timer_start(t_domains_keys)
    ! assign SFC coordinate to each particle
    call compute_particle_keys(t%bounding_box, p)
    call timer_stop(t_domains_keys)

    ! Domain decomposition: allocate particle keys to PEs
    call domain_decompose(t%decomposition, t%bounding_box, n, p, bp, tree_comm_env, npl = npl)

    ! allocate the tree structure
    call tree_create(t, size(p), n, comm_env = tree_comm_env)

    ! build local part of tree
    call timer_start(t_local)
    call tree_build_from_particles(t, p, bp)

    call tree_lookup_root(t, root_node, 'libpepc_grow_tree:root node')
    if (root_node%leaves .ne. size(p)) then
      call htable_dump(t%node_storage, p)
      DEBUG_ERROR(*, 'did not find all its particles inside the tree after local tree buildup: root_node%leaves =', root_node%leaves, ' but size(particles) =', size(p))
    end if
    call timer_stop(t_local)

    ! Should now have multipole information up to root list level(s) (only up to branch level, the information is correct)
    ! By definition, this is complete: each branch node is self-contained.
    ! This information has to be broadcast to the other PEs so that the top levels can be filled in.

    call timer_start(t_exchange_branches)
    call tree_exchange_branches(t, p, bp, branch_keys)
    call timer_stop(t_exchange_branches)

    ! build global part of tree
    call timer_start(t_global)
    call tree_build_upwards(t, branch_keys)
    call timer_stop(t_global)
    deallocate(branch_keys)

    if (root_node%leaves .ne. t%npart) then
      call htable_dump(t%node_storage, p)
      DEBUG_ERROR(*, 'did not find all particles inside the htable after global tree buildup: root_node%leaves =', root_node%leaves, ' but npart_total =', t%npart)
    endif

    call timer_stop(t_fields_tree)
    call pepc_status('TREE GROWN')
  end subroutine tree_grow


  !>
  !> Exchanges tree nodes that are given in `local_branch_keys` with remote ranks.
  !>
  !> Incoming tree nodes are inserted into tree `t`, but the tree above these nodes is not corrected.
  !> Outputs keys of new (and own) tree nodes in `branch_keys`.
  !>
  subroutine tree_exchange(t, local_branch_keys, branch_keys)
    use module_tree, only: t_tree, tree_lookup_node_critical, &
      tree_insert_or_update_node
    use module_pepc_types, only: t_tree_node, MPI_TYPE_tree_node
    use module_debug, only : pepc_status
    use module_timings
    use module_tree_node
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t !< the tree
    integer*8, intent(in) :: local_branch_keys(:) !< keys of all local branch nodes
    integer*8, intent(inout), allocatable :: branch_keys(:) !< keys of all branch nodes in the tree

    integer :: ierr
    integer :: i, nbranch, nbranch_sum
    type(t_tree_node), pointer :: branch_node
    type(t_tree_node), allocatable :: pack_mult(:), get_mult(:)
    integer, allocatable :: nbranches(:), igap(:)

    call timer_start(t_exchange_branches_pack)

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
    allocate (nbranches(t%comm_env%size), igap(t%comm_env%size + 1))
    call mpi_allgather(nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, t%comm_env%comm, ierr)

    igap(1) = 0
    do i = 2, t%comm_env%size + 1
      igap(i) = igap(i - 1) + nbranches(i - 1)
    end do

    nbranch_sum = igap(t%comm_env%size + 1)

    allocate(get_mult(1:nbranch_sum), branch_keys(1:nbranch_sum))

    call timer_stop(t_exchange_branches_admininstrative)
    call timer_start(t_exchange_branches_allgatherv)

    ! actually exchange the branch nodes
    call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_tree_node, get_mult, nbranches, igap, MPI_TYPE_tree_node, &
      t%comm_env%comm, ierr)

    deallocate(pack_mult)
    deallocate(nbranches, igap)

    call timer_stop(t_exchange_branches_allgatherv)
    call timer_start(t_exchange_branches_integrate)

    ! Integrate remote branches into local tree
    do i = 1, nbranch_sum

      ! insert all remote branches into local data structures (this does *not* prepare the internal tree connections, but only copies multipole properties and creates the htable-entries)
      if (get_mult(i)%owner /= t%comm_env%rank) then
        ! delete all custom flags from incoming nodes (e.g. TREE_NODE_FLAG_CHILDREN_AVAILABLE)
        get_mult(i)%flags = IAND(get_mult(i)%flags, TREE_NODE_CHILDBYTE)
        ! after clearing all bits we have to set the flag for branches again to propagate this property upwards during global buildup
        get_mult(i)%flags = ibset(get_mult(i)%flags, TREE_NODE_FLAG_IS_BRANCH_NODE)
        ! additionally, we mark all remote branches as remote nodes (this information is propagated upwards later)
        get_mult(i)%flags = ibset(get_mult(i)%flags, TREE_NODE_FLAG_HAS_REMOTE_CONTRIBUTIONS)

        call tree_insert_or_update_node(t, get_mult(i))
      end if
      ! store branch key for later (global tree buildup)
      branch_keys(i) = get_mult(i)%key
    end do

    deallocate(get_mult)

    call timer_stop(t_exchange_branches_integrate)
  end subroutine tree_exchange


  !>
  !> identifies branch nodes in tree `t` and makes them available on all communication ranks.
  !>
  !> the keys of all boundary nodes are returned in `bk`.
  !>
  subroutine tree_exchange_branches(t, p, bp, bk)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use module_timings
    use module_debug, only: pepc_status
    implicit none

    type(t_tree), intent(inout) :: t !< the tree
    type(t_particle), intent(in) :: p(:) !< list of local particles
    type(t_particle), intent(in) :: bp(2) !< boundary particles
    integer*8, allocatable, intent(inout) :: bk(:) !< keys of all branch nodes

    integer :: num_local_branch_keys
    integer*8, allocatable :: local_branch_keys(:)

    ! identification of branch nodes
    call timer_start(t_branches_find)
    call find_local_branches(t, p, bp, num_local_branch_keys, local_branch_keys)
    t%nbranch_me = num_local_branch_keys
    call timer_stop(t_branches_find)

    call pepc_status('EXCHANGE BRANCHES')

    call tree_exchange(t, local_branch_keys(1:num_local_branch_keys), bk)
    t%nbranch = size(bk)

    deallocate(local_branch_keys)
    
    contains

    !>
    !> identifies all `nb` local branch nodes in tree `t` and returns their
    !> keys in `bk`.
    !>
    subroutine find_local_branches(t, p, bp, nb, bk)
      use module_spacefilling, only: shift_key_by_level
      use module_math_tools, only: bpi
      use treevars, only: idim, nlev
      use module_debug, only: pepc_status
      use module_tree, only: tree_contains_key
      implicit none

      type(t_tree), intent(inout) :: t !< the tree
      type(t_particle), intent(in) :: p(:) !< list of local particles
      type(t_particle), intent(in) :: bp(2) !< boundary particles
      integer, intent(out) :: nb !< number of local branch nodes
      integer*8, allocatable, intent(out) :: bk(:) !< list of keys of local branch nodes

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
      if (t%comm_env%last) then
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
          if (tree_contains_key(t, possible_branch)) then ! entry exists
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
          if (tree_contains_key(t, possible_branch)) then ! entry exists
            nb = nb + 1
            bk(nb) = possible_branch
          end if
        end do
      end do
    end subroutine find_local_branches


    !>
    !> find the virtual limits of the local key domain
    !>
    subroutine find_vld_limits(t, p, bp, l, r)
      use module_math_tools, only: bpi
      use treevars, only: idim, nlev
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      type(t_particle), intent(in) :: p(:) !< list of local particles
      type(t_particle), intent(in) :: bp(2) !< boundary particles
      integer*8, intent(out) :: l !< left virtual key domain limit
      integer*8, intent(out) :: r !< right virtual key domain limit

      integer*8 :: lme, rme, lb, rb
        
      ! get local key limits
      lme = p(1)%key
      rme = p(ubound(p, 1))%key

      ! get key limits for neighbor tasks
      ! and build virtual limits, so that a minimum set of branch nodes comes arround
      ! boundary tasks can access their boundary space fully only need one virtual limit
      l = 2_8**(nlev * idim)
      r = 2_8**(nlev * idim + 1) - 1

      if (.not. t%comm_env%first) then
        lb = bp(1)%key
        l = bpi(lb, lme)
      end if

      if (.not. t%comm_env%last) then
        rb = bp(2)%key
        r = bpi(rme, rb)
      end if
    end subroutine find_vld_limits
  end subroutine tree_exchange_branches


  !>
  !> Builds up the tree `t` from the given start keys `keys` towards root
  !>  - expects, that the nodes that correspond to `keys` already
  !>    have been inserted into the tree
  !>  - missing nodes on the way towards root are added automatically
  !>  - already existing nodes are updated
  !>
  subroutine tree_build_upwards(t, keys)
    use module_tree, only: t_tree, tree_insert_or_update_node
    use module_debug, only: pepc_status, DBG_TREE, dbg
    use module_timings
    use module_utils, only: sort
    use module_htable, only: htable_check
    use module_spacefilling, only: level_from_key, parent_key_from_key
    use module_pepc_types, only: t_tree_node
    implicit none

    type(t_tree), intent(inout) :: t !< the tree
    integer*8, intent(in) :: keys(:) !< keys of nodes at which to start the build upwards

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

        call shift_nodes_up_key(t, parent_node, sub_key(groupstart:groupend), t%comm_env%rank)
        call tree_insert_or_update_node(t, parent_node)

        nparent             = nparent + 1
        parent_key(nparent) = current_parent_key

        ! go on with next group
        i = i + 1
        end do

    end do

    deallocate(key_level, sub_key, parent_key)
  end subroutine tree_build_upwards


  !>
  !> Inserts all given local particles `p` at the correct leaf level in tree `t`
  !> by recursively subdividing the nodes if necessary.
  !>
  !> After function execution, leaf nodes exist for all particles in `p` and
  !> all twig nodes above those leaves until the root node.
  !> Data for leaves are completely valid while those
  !> for twigs have to be updated via a call to `tree_build_upwards()`.
  !>
  !> Upon exit, the `key_leaf` property of any particle is set appropriately.
  !>
  !> @warning In contrast to `tree_build_upwards()`, this function may *not*
  !> be called several times to add further particles etc.
  !>
  subroutine tree_build_from_particles(t, p, bp)
    use module_tree, only: t_tree
    use treevars, only : idim, nlev
    use module_pepc_types, only: t_particle
    use module_spacefilling, only: level_from_key
    use module_timings
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t !< the tree
    type(t_particle), intent(inout) :: p(:) !< list of particles
    type(t_particle), intent(in) :: bp(2) !< boundary particles

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
    if (.not. t%comm_env%first) then
      i = i + 1
      kidx(i) = t_keyidx( bp(1)%key, 0 )
    end if

    do j = 1, ubound(p, 1)
      i = i + 1
      kidx(i) = t_keyidx( p(j)%key, j )
    end do

    if (.not. t%comm_env%last) then
      i = i + 1
      kidx(i) = t_keyidx( bp(2)%key, 0 )
    end if

    p(:)%key_leaf = 0_8

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

    !>
    !> Inserts the keys in `ki` below the key `k`.
    !>
    recursive subroutine insert_helper(k, l, ki, pi)
      use module_tree, only: tree_insert_node
      use module_pepc_types, only: t_tree_node
      use module_tree_node
      use module_spacefilling, only: child_key_from_parent_key, is_ancestor_of_particle
      use module_interaction_specific, only: multipole_from_particle
      implicit none

      integer*8, intent(in) :: k !< key blow which the `ki` are inserted
      integer, intent(in) :: l !< precomputed `level_from_key(k)`
      type(t_keyidx), intent(in) :: ki(:) !< keys (with particles) to be inserted
      type(t_tree_node), optional, pointer, intent(out) :: pi

      type(t_tree_node) :: this_node
      type(t_tree_node), pointer :: child_node
      type(t_tree_node) :: child_nodes(8)
      integer*8 :: childkey
      integer :: ip, ichild, childlevel, pstart, pend, nchild, child_number(8)

      if (present(pi)) then
        pi => null()
      end if

      if (size(ki) < 1) then ! no particles below this key
        ! do nothing
      else if (size(ki) == 1) then ! we have arrived at a leaf
        if (ki(1)%idx /= 0) then ! boundary particles have idx == 0, do not really need to be inserted
          this_node%flags      = ibset(0, TREE_NODE_FLAG_HAS_LOCAL_CONTRIBUTIONS)
          this_node%owner      = t%comm_env%rank
          this_node%key        = k
          this_node%level      = l
          this_node%leaves     = 1
          call multipole_from_particle(p(ki(1)%idx)%x, p(ki(1)%idx)%data, this_node%interaction_data)
          p(ki(1)%idx)%key_leaf = k
          if (.not. tree_insert_node(t, this_node, pi)) then
            DEBUG_ERROR(*, "Leaf allready inserted, aborting.") ! TODO: tell me more!
          end if
        end if
      else ! more particles left, split twig
        if (l >= nlev) then ! no more levels left, cannot split
          DEBUG_WARNING_ALL('("Problem with tree: No more levels. Remaining particles 1..",I0,"  [i, key, label, x, y, z]:")', size(ki))
          DEBUG_ERROR_NO_HEADER('(I6,x,O22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, p(ki(ip)%idx)%key, p(ki(ip)%idx)%label, p(ki(ip)%idx)%x(1:3), ip=1,size(ki) ) )
        end if

        nchild = 0
        pstart = lbound(ki, dim = 1)
        do ichild = 0, 2**idim - 1
          childkey   = child_key_from_parent_key(k, ichild)
          childlevel = l + 1

          pend = pstart - 1
          do while (pend < ubound(ki, dim = 1))
            if (.not. is_ancestor_of_particle(ki(pend + 1)%key, childkey, childlevel)) exit
            pend = pend + 1
          end do
          
          call insert_helper(childkey, childlevel, ki(pstart:pend), child_node)
          if (associated(child_node)) then
            nchild = nchild + 1
            child_nodes(nchild) = child_node
            child_number(nchild) = ichild
          end if

          pstart = pend + 1
        end do

        if (.not. pend == ubound(ki, dim = 1)) then
          DEBUG_WARNING_ALL('("Problem with tree: Could not distribute particles among children of ", I0, ". Remaining particles ",I0,"..",I0,"  [i, key, label, x, y, z]:")', k, pend + 1, ubound(ki, dim = 1) )
          DEBUG_ERROR_NO_HEADER('(I6,x,I22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, p(ki(ip)%idx)%key, p(ki(ip)%idx)%label, p(ki(ip)%idx)%x(1:3), ip=pend + 1, ubound(ki, dim = 1) ) )
        end if

        call shift_nodes_up(this_node, child_nodes(1:nchild), child_number(1:nchild), t%comm_env%rank)

        if (.not. tree_insert_node(t, this_node, pi)) then
          DEBUG_ERROR(*, "Twig allready inserted, aborting.") ! TODO: tell me more!
        end if
      end if
    end subroutine insert_helper
  end subroutine tree_build_from_particles


  !>
  !> accumulates properties of child nodes (given by keys) to parent node
  !>
  subroutine shift_nodes_up_key(t, parent, childkeys, parent_owner)
    use module_tree, only: t_tree, tree_lookup_node_critical
    use module_pepc_types, only: t_tree_node
    use module_spacefilling, only: child_number_from_key
    implicit none

    type(t_tree), intent(inout) :: t !< tree in which to find the nodes
    type(t_tree_node), intent(inout) :: parent !< parent node
    integer*8, intent(in) :: childkeys(:) !< keys of child nodes
    integer, intent(in) :: parent_owner !< communication rank that owns the parent

    integer :: nchild, i
    type(t_tree_node) :: child_nodes(1:8)
    type(t_tree_node), pointer :: p
    integer :: childnumber(1:8)

    nchild = size(childkeys)

    do i = 1, nchild
      childnumber(i) = child_number_from_key(childkeys(i))
      call tree_lookup_node_critical(t, childkeys(i), p, 'shift_nodes_up_key')
      child_nodes(i) = p
    end do

    call shift_nodes_up(parent, child_nodes(1:nchild), childnumber(1:nchild), parent_owner)
  end subroutine shift_nodes_up_key


  !>
  !> accumulates properties of child nodes to parent node
  !>
  subroutine shift_nodes_up(parent, children, childnumber, parent_owner)
    use module_pepc_types, only: t_tree_node
    use module_tree_node
    use module_interaction_specific, only : shift_multipoles_up
    use module_spacefilling, only: parent_key_from_key, level_from_key
    use module_debug
    implicit none

      type(t_tree_node), intent(inout) :: parent !< parent node
      type(t_tree_node), intent(in) :: children(:) !< child nodes
      integer, intent(in) :: childnumber(:) !< child number of each child node
      integer, intent(in) :: parent_owner !< communication rank that owns the parent

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
  end subroutine shift_nodes_up
end module module_tree_grow
