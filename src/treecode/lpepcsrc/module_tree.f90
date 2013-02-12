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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all tree specific helper routines and data fields
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_tree
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public tree_insert_node
    public tree_exchange
    public tree_build_upwards
    public tree_build_from_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> type for storing key and nideindex together in tree_build_from_particles
    type t_keyidx
      integer*8 :: key
      integer :: idx
    end type

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> if an entry with tree_node%key already exists in htable, then
    !> updates htable and tree_nodes with new data
    !> otherwise: creates new entries
    !>
    !> this routine cannot be used to change a tree_node from leaf to twig or similar
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_update_or_insert_node(tree_node)
        use treevars
        use module_pepc_types
        use module_htable
        use module_tree_node
        implicit none
        include 'mpif.h'

        type(t_tree_node), intent(in) :: tree_node

        type(t_tree_node), pointer :: preexisting_node

        if (htable_lookup(global_htable, tree_node%key, preexisting_node)) then
          ! the htable-entry and node already exist --> update

          ! if we change the owner from someting else to 'me', we have to keep track of the leaf/twig counters
          if ((preexisting_node%owner .ne. me) .and. (tree_node%owner .eq. me)) then
            if (tree_node_is_leaf(preexisting_node)) then
              nleaf_me = nleaf_me + 1
            else
              ntwig_me = ntwig_me + 1
            endif
          endif

          preexisting_node%leaves           = tree_node%leaves
          preexisting_node%info_field       = tree_node%info_field
          preexisting_node%owner            = tree_node%owner
          preexisting_node%interaction_data = tree_node%interaction_data
        else
          ! create new htable and nodelist entry
          call tree_insert_node(tree_node)
        endif

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Inserts a given tree node into the next free position in the tree ( -(ntwig+1) or (nleaf+1) )
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_insert_node(tree_node)
      use treevars
      use module_pepc_types
      use module_htable
      use module_debug
      implicit none

      type(t_tree_node), intent(in) :: tree_node
      type(t_tree_node), pointer :: preexisting_entry

      if (htable_add(global_htable, tree_node%key, tree_node, preexisting_entry)) then
        ! anything is fine - we will have to assign a node number now
        if ( tree_node%leaves == 1 ) then ! TODO: this should be replaced by tree_node_is_leaf()
          nleaf =  nleaf + 1
          if (tree_node%owner == me) nleaf_me = nleaf_me+1
        else if ( tree_node%leaves > 1 ) then
          ! twig
          ntwig =  ntwig + 1
          if (tree_node%owner == me) ntwig_me = ntwig_me+1
        else
          DEBUG_ERROR(*, "Found a tree node with less than 1 leaf.")
        endif

      else
        ! entry with the same key is already existing, so we just overwrite it
        preexisting_entry%interaction_data = tree_node%interaction_data

        ! TODO: cleanup
        DEBUG_WARNING_ALL(*, "PE", me, "has found an already inserted entry while calling make_hashentry(", tree_node%key, tree_node%leaves, tree_node%info_field, tree_node%owner, tree_node%level, ") - overwriting it")
      endif

    end subroutine tree_insert_node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Accumulates properties of child nodes (given by keys) to parent node
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine shift_nodes_up_key(parent, childkeys, parent_owner)
      use module_pepc_types
      use module_htable
      use module_spacefilling
      implicit none
      type(t_tree_node), intent(inout) :: parent
      integer*8, intent(in) :: childkeys(:)
      integer, intent(in) :: parent_owner

      integer :: nchild, i
      type(t_tree_node) :: child_nodes(1:8)
      type(t_tree_node), pointer :: p
      integer :: childnumber(1:8)

      nchild = size(childkeys)

      do i=1,nchild
        childnumber(i) = child_number_from_key(childkeys(i))
        call htable_lookup_critical(global_htable, childkeys(i), p, 'shift_nodes_up_key')
        child_nodes(i) = p
      end do

      call shift_nodes_up(parent, child_nodes(1:nchild), childnumber(1:nchild), parent_owner)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Accumulates properties of child nodes to parent node
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine shift_nodes_up(parent, children, childnumber, parent_owner)
      use module_pepc_types
      use module_tree_node
      use module_interaction_specific, only : shift_multipoles_up
      use module_spacefilling
      use module_debug
      use module_htable
      implicit none
        type(t_tree_node), intent(inout) :: parent
        type(t_tree_node), intent(in) :: children(:)
        integer, intent(in) :: childnumber(:)
        integer, intent(in) :: parent_owner
        integer*8 :: parent_keys(1:8)

        integer :: nchild, i, info_field

        nchild = size(children)

        ! check if all keys fit to the same parent
        parent_keys(1:nchild) = parent_key_from_key(children(1:nchild)%key)

        if ( any(parent_keys(2:nchild) .ne. parent_keys(1))) then
          DEBUG_ERROR(*,"Error in shift nodes up: not all supplied children contribute to the same parent node")
        endif

        info_field = 0
        do i = 1,nchild
          ! set bits for available children
          info_field = IBSET(info_field, childnumber(i))
          ! parents of nodes with local contributions also contain local contributions
          if (btest(children(i)%info_field, CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS)) info_field = ibset(info_field, CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS)
          ! parents of nodes with remote contributions also contain remote contributions
          if (btest(children(i)%info_field, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)) info_field = ibset(info_field, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)
          ! parents of branch and fill nodes will also be fill nodes
          if (btest(children(i)%info_field, CHILDCODE_BIT_IS_FILL_NODE) .or. btest(children(i)%info_field, CHILDCODE_BIT_IS_BRANCH_NODE)) info_field = ibset(info_field, CHILDCODE_BIT_IS_FILL_NODE)
        end do


        ! Set children_HERE flag parent since we just built it from its children
        info_field =  IBSET( info_field, CHILDCODE_BIT_CHILDREN_AVAILABLE )

        parent%key        = parent_keys(1)
        parent%info_field = info_field
        parent%leaves     = sum(children(1:nchild)%leaves)
        parent%owner      = parent_owner
        parent%level      = level_from_key( parent_keys(1) )

        call shift_multipoles_up(parent%interaction_data, children(1:nchild)%interaction_data)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Exchanges tree nodes that are given in local_branch_keys with remote PEs
    !> incoming tree nodes are inserted into tree_nodes array and htable, but the
    !> tree above these nodes is not corrected
    !> outputs keys of new(and own) htable/tree_node entries in branch_keys
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_exchange(local_branch_keys, nbranch, branch_keys, nbranch_sum)

        use treevars, only : me, num_pe, nbranches, MPI_COMM_lpepc
        use module_debug, only : pepc_status
        use module_pepc_types
        use module_timings
        use module_htable
        use module_tree_node
        implicit none
        include 'mpif.h'

        integer*8, intent(in) :: local_branch_keys(1:nbranch)
        integer, intent(in) :: nbranch
        integer*8, intent(inout), allocatable :: branch_keys(:)
        integer, intent(out) :: nbranch_sum

        integer :: i,ierr
        type(t_tree_node), pointer :: branch_node
        type(t_tree_node), allocatable :: pack_mult(:), get_mult(:)
        integer, allocatable :: igap(:)    !  stride lengths of local branch arrays

        if (allocated(branch_keys)) deallocate(branch_keys)

        call timer_start(t_exchange_branches_pack)

        call pepc_status('EXCHANGE BRANCHES')

        ! Pack local branches for shipping
        allocate(pack_mult(nbranch))
        do i=1,nbranch
            call htable_lookup_critical(global_htable, local_branch_keys(i), branch_node, 'EXCHANGE: info')

            ! additionally, we mark all local branches as branches since this is only done for remote branches during unpack (is used for fill node identification)
            branch_node%info_field = ibset(branch_node%info_field, CHILDCODE_BIT_IS_BRANCH_NODE)
            pack_mult(i) = branch_node
        end do

        call timer_stop(t_exchange_branches_pack)
        call timer_start(t_exchange_branches_admininstrative)

        call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_lpepc, ierr )

        ! work out stride lengths so that partial arrays placed sequentially in global array
        allocate (igap(num_pe+3))

        igap(1) = 0
        do i=2,num_pe+1
            igap(i) = igap(i-1) + nbranches(i-1)
        end do

        nbranch_sum = igap(num_pe+1)

        allocate(get_mult(1:nbranch_sum), branch_keys(1:nbranch_sum))

        call timer_stop(t_exchange_branches_admininstrative)
        call timer_start(t_exchange_branches_allgatherv)

        ! actually exchange the branch nodes
        call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_tree_node, get_mult, nbranches, igap, MPI_TYPE_tree_node, &
          MPI_COMM_lpepc, ierr)

        deallocate(pack_mult)
        deallocate (igap)

        call timer_stop(t_exchange_branches_allgatherv)
        call timer_start(t_exchange_branches_integrate)

        ! Integrate remote branches into local tree
        do i = 1,nbranch_sum

            ! insert all remote branches into local data structures (this does *not* prepare the internal tree connections, but only copies multipole properties and creates the htable-entries)
            if (get_mult(i)%owner /= me) then
              ! delete all custom flags from incoming nodes (e.g. CHILDCODE_BIT_CHILDREN_AVAILABLE)
              get_mult(i)%info_field = IAND(get_mult(i)%info_field, CHILDCODE_CHILDBYTE)
              ! after clearing all bits we have to set the flag for branches again to propagate this property upwards during global buildup
              get_mult(i)%info_field = ibset(get_mult(i)%info_field, CHILDCODE_BIT_IS_BRANCH_NODE)
              ! additionally, we mark all remote branches as remote nodes (this information is propagated upwards later)
              get_mult(i)%info_field = ibset(get_mult(i)%info_field, CHILDCODE_BIT_HAS_REMOTE_CONTRIBUTIONS)

              call tree_insert_node(get_mult(i))
            endif
            ! store branch key for later (global tree buildup)
            branch_keys(i) = get_mult(i)%key
        end do

        deallocate(get_mult)

        call timer_stop(t_exchange_branches_integrate)

    end subroutine tree_exchange


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds up the tree from the given start keys towards root
    !>  - expects, that the nodes that correspond to start_keys already
    !>    have been inserted into htable and tree_nodes array
    !>  - missing nodes on the way towards root are added automatically
    !>  - already existing nodes are updated
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_build_upwards(start_keys, numkeys)

        use treevars, only : me
        use module_debug, only : pepc_status, DBG_TREE, dbg
        use module_timings
        use module_utils
        use module_htable
        use module_spacefilling
        use module_interaction_specific
        use module_pepc_types
        implicit none

        integer*8, intent(in) :: start_keys(1:numkeys)
        integer, intent(in) :: numkeys
        integer, dimension(1:numkeys) :: branch_level
        integer*8, dimension(0:numkeys+1) :: sub_key, parent_key
        type(t_tree_node) :: parent_node

        integer :: ilevel, maxlevel, nsub, groupstart, groupend, i, nparent, nuniq
        integer*8 :: current_parent_key

        call pepc_status('BUILD TOWARDS ROOT')

        if (dbg(DBG_TREE)) call htable_check(global_htable, 'after make_branches ')

        ! get levels of branch nodes
        branch_level(1:numkeys) = level_from_key(start_keys(1:numkeys))
        maxlevel = maxval( branch_level(1:numkeys) )        ! Find maximum level

        nparent = 0

        ! iterate through branch levels
        do ilevel = maxlevel,1,-1                                            ! Start at finest branch level
            ! Collect all branches at this level
            nsub = 0
            do i=1,numkeys
                if (branch_level(i) == ilevel) then
                    nsub          = nsub + 1
                    sub_key(nsub) = start_keys(i)
                endif
            end do

            ! Augment list with parent keys checked at previous level
            sub_key(nsub+1:nsub+nparent) = parent_key(1:nparent)
            nsub                         = nsub + nparent

            call sort(sub_key(1:nsub))                                        ! Sort keys

            sub_key(0)   = 0                                                  ! remove all duplicates from the list
            nuniq = 0
            do i=1,nsub
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

              call shift_nodes_up_key(parent_node, sub_key(groupstart:groupend), me)
              call tree_update_or_insert_node(parent_node)

              nparent             = nparent + 1
              parent_key(nparent) = current_parent_key

              ! go on with next group
              i = i + 1
            end do

        end do

    end subroutine tree_build_upwards


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_build_from_particles(particle_list, nparticles, neighbour_particles, leaf_keys)
      use treevars, only : nleaf, ntwig, nlev, me, nleaf_me, ntwig_me, idim, num_pe
      use module_pepc_types
      use module_spacefilling
      use module_htable
      use module_interaction_specific, only : multipole_from_particle
      use module_timings
      use module_debug
      use module_tree_node
      implicit none
      type(t_particle), intent(inout) :: particle_list(1:nparticles + neighbour_particles)
      integer, intent(in) :: nparticles
      integer, intent(in) :: neighbour_particles
      integer*8, intent(out) :: leaf_keys(1:nparticles + neighbour_particles)

      type(t_keyidx), allocatable :: keyidx(:)
      integer :: i, j, numkeys

      call pepc_status('INSERT PARTICLES')

      call timer_start(t_build_pure)

      numkeys = nparticles + neighbour_particles
      allocate(keyidx(numkeys))
      i = 1

      if (me /= 0) then
        keyidx(i) = t_keyidx( particle_list(nparticles + neighbour_particles)%key, nparticles + neighbour_particles )
        i = i + 1
      end if

      do j = 1, nparticles
        keyidx(i) = t_keyidx( particle_list(j)%key, j )
        i = i + 1
      end do

      if (me /= num_pe - 1) then
        keyidx(i) = t_keyidx( particle_list(nparticles + 1)%key, nparticles + 1 )
      end if

      particle_list(1:numkeys)%key_leaf = 0_8
      nleaf = numkeys

      call timer_start(t_props_leafs) ! TODO: this timer does not make sense like this
      call insert_helper(1_8, level_from_key(1_8), keyidx(1:numkeys))
      call timer_stop(t_props_leafs)

      deallocate(keyidx)
      call timer_stop(t_build_pure)

      nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
      ntwig_me = ntwig

      ! copy leaf keys
      leaf_keys(1:numkeys) = particle_list(1:numkeys)%key_leaf

      ! check if we did not miss any particles
      if (any(leaf_keys(1:numkeys) == 0_8)) then
        DEBUG_WARNING_ALL(*, ' did not incorporate all particles into its leaf_keys array')
      endif

      contains


      recursive subroutine insert_helper(k, l, ki)
        implicit none

        integer*8, intent(in) :: k
        integer, intent(in) :: l
        type(t_keyidx), intent(inout) :: ki(:)

        type(t_tree_node) :: this_node
        integer*8 :: childkey
        integer :: ip, ichild, childlevel, pstart, pend

        this_node%info_field = ibset(0, CHILDCODE_BIT_HAS_LOCAL_CONTRIBUTIONS)
        this_node%owner      = me
        this_node%key        = k
        this_node%level      = l

        if (size(ki) == 1) then ! we have arrived at a leaf
          this_node%leaves     = 1
          call multipole_from_particle(particle_list(ki(1)%idx)%x, particle_list(ki(1)%idx)%data, this_node%interaction_data)
          particle_list(ki(1)%idx)%key_leaf = k
          if (.not. htable_add(global_htable, k, this_node)) then
            DEBUG_ERROR(*, "Leaf allready inserted, aborting.") ! TODO: tell me more!
          end if

        else ! more particles left, split twig
          if (l >= nlev) then ! no more levels left, cannot split
            DEBUG_WARNING_ALL('("Problem with tree: No more levels. Remaining particles 1..",I0,"  [i, key, label, x, y, z]:")', size(ki))
            DEBUG_ERROR_NO_HEADER('(I6,x,O22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, particle_list(ki(ip)%idx)%key, particle_list(ki(ip)%idx)%label, particle_list(ki(ip)%idx)%x(1:3), ip=1,size(ki) ) )
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
            DEBUG_ERROR_NO_HEADER('(I6,x,I22.22,x,I16,g20.12,x,g20.12,x,g20.12,x)', ( ip, particle_list(ki(ip)%idx)%key, particle_list(ki(ip)%idx)%label, particle_list(ki(ip)%idx)%x(1:3), ip=pend + 1, ubound(ki, dim = 1) ) )
          end if

          if (.not. htable_add(global_htable, k, this_node)) then
            DEBUG_ERROR(*, "Twig allready inserted, aborting.") ! TODO: tell me more!
          end if
          ntwig = ntwig + 1
        end if

      end subroutine insert_helper

    end subroutine




end module module_tree
