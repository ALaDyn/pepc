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
    implicit none
    private

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
      
      !TODO: factor these out into tree comm statistics?
      integer*8 :: sum_ships   !< total number of node ships
      integer*8 :: sum_fetches !< total number of node fetches

      type(t_box) :: bounding_box            !< bounding box enclosing all particles contained in the tree
      type(t_htable) :: node_storage         !< hash table in which tree nodes are stored for rapid retrieval
      type(t_comm_env) :: comm_env         !< communication environment over which the tree is distributed
      type(t_decomposition) :: decomposition !< permutation of particles inserted into the tree
    end type t_tree

    public tree_create
    public tree_insert_node
    public tree_insert_or_update_node
    public tree_contains_key
    public tree_lookup_root
    public tree_lookup_node
    public tree_lookup_node_critical
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
      t%sum_ships = 0
      t%sum_fetches = 0

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
        DEBUG_ERROR('("maxaddress = ", I0, " < npp + 2 = ", I0 ".", / , "You should increase np_mult.")', maxaddress, t%npart_me + 2)
      end if

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

      type(t_tree), intent(inout) :: t

      call pepc_status('DEALLOCATE TREE')

      t%npart = 0
      t%npart_me = 0
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
      call comm_env_destroy(t%comm_env)
      if (decomposition_allocated(t%decomposition)) then
        call decomposition_destroy(t%decomposition)
      end if
    end subroutine tree_destroy


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
      use module_pepc_types
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
      use module_pepc_types
      use module_htable
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
      use module_pepc_types
      use module_htable
      implicit none

      type(t_tree), intent(in) :: t !< the tree
      integer*8, intent(in) :: k !< key to look up
      type(t_tree_node), pointer, intent(out) :: n !< node that is identified by `k`
      character(LEN = *), intent(in) :: caller

      ! TODO: move mpi_abort here!
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

      s = t%comm_env%size
      allocate(nparticles(s), fetches(s), ships(s), total_keys(s), tot_nleaf(s), &
        tot_ntwig(s), num_interactions(s), num_mac_evaluations(s))

      ! particle distrib
      call MPI_GATHER(t%npart_me, 1, MPI_INTEGER, nparticles, 1, MPI_INTEGER, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%ntwig_me,    1, MPI_INTEGER8, tot_ntwig,  1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%nleaf_me,    1, MPI_INTEGER8, tot_nleaf,  1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      nkeys_total = t%nleaf + t%ntwig
      call MPI_GATHER(nkeys_total,   1, MPI_INTEGER8, total_keys, 1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%sum_fetches, 1, MPI_INTEGER8, fetches,    1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
      call MPI_GATHER(t%sum_ships,   1, MPI_INTEGER8, ships,      1, MPI_INTEGER8, 0,  t%comm_env%comm, ierr )
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
