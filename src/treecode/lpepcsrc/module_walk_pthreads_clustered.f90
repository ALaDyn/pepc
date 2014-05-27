! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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
!>  Perform tree walk for all local particles
!>  in a hybrid parallelization scheme using
!>  linux pthreads
!>
!>  Algorithm follows the implementation of
!>  Warren & Salmon`s 'latency-hiding' concept,
!>  retaining list-based tree-walk from vectorised
!>  code by Pfalzner & Gibbon.
!>
!>
!>  Structure:
!>    * the main thread upon entering `tree_walk` spawns a number of
!>      worker threads and waits for them to complete.
!>    * each worker thread (`walk_worker_thread`)
!>      grabs a number of particles and performs their
!>      individual walks.
!>    * if a potential interaction partner of a particle is not available
!>      locally, it is requested via `tree_node_fetch_children()`
!>      and the walk for that particle is deferred until later
!>    * as soon as one particle has
!>      finished walking, the worker thread takes an
!>      additional particle as long as there are still
!>      unprocessed particles available
!>    * when all walks of all particles are finished, the
!>      worker threads are terminated and the main thread
!>      continues execution
!>
!>
!>  Structure of individual walk_work_threads:
!>  ------------------------------------------
!>      do while (particles_active .or. particles_available)
!>
!>        particles_active = .false.
!>
!>        do i=1,max_clusters_per_thread
!>
!>          if ( (my_particles(i) == -1) .and. (particles_available)) then         ! i.e. the place for a particle is unassigned
!>            my_clusters(i) = get_first_unassigned_cluster()
!>          end if
!>
!>          call walk_single_cluster(my_cluster(i))
!>
!>        end do
!>      end do
!>
!>
!>  Structure of walk_single_cluster(cluster):
!>  ------------------------------------------
!>
!>     if (.not.finished(cluster)) then
!>       num_unfinished = num_unfinished + 1
!>
!>       check defer_list entries:
!>         if (requested children available)
!>            put them onto todo_list
!>
!>       do while (can take entry form todo_list)
!>           if (MAC OK)
!>                immedeately interact with node
!>           else
!>                if (node locally available)
!>                    resolve node
!>                    put all children to front of todo_list
!>                else
!>                    post_request(parentkey, owner)
!>                    put node on defer_list
!>                end if
!>           end if
!>       end do
!>     end if
!>
module module_walk
  use, intrinsic :: iso_c_binding
  use module_tree, only: t_tree
  use module_pepc_types
  use module_atomic_ops, only: t_atomic_int
  use pthreads_stuff, only: t_pthread_with_type
  implicit none
  private

  integer, parameter :: NUM_THREAD_COUNTERS                = 4
  integer, parameter :: THREAD_COUNTER_PROCESSED_CLUSTERS  = 1
  integer, parameter :: THREAD_COUNTER_INTERACTIONS        = 2
  integer, parameter :: THREAD_COUNTER_MAC_EVALUATIONS     = 3
  integer, parameter :: THREAD_COUNTER_POST_REQUEST        = 4

  !> type for input and return values of walk_threads
  type :: t_threaddata
    integer :: id !< just a running number to distinguish the threads, currently unused
    logical :: is_on_shared_core !< thread output value: is set to true if the thread detects that it shares its processor with the communicator thread
    integer :: coreid !< thread output value: id of thread's processor
    logical :: finished !< will be set to .true. when the thread has finished
    integer(kind_node) :: counters(NUM_THREAD_COUNTERS)
  end type t_threaddata

  type(t_pthread_with_type), target, allocatable :: thread_handles(:)
  type(t_threaddata), allocatable, target :: threaddata(:)
  integer :: num_walk_threads = -1 !< number of worker threads, value is set in tree_walk_init()
  real :: work_on_communicator_cluster_number_factor = 0.1 !< factor for reducing max_clusters_per_thread for thread which share their processor with the communicator
  ! variables for adjusting the thread's workload
  integer, public :: max_clusters_per_thread = 2000 !< maximum number of particles that will in parallel be processed by one workthread

  real*8 :: vbox(3)
  integer :: todo_list_length, defer_list_length, num_clusters
  type(t_particle), pointer, dimension(:) :: particle_data
  type(t_tree), pointer :: walk_tree

  type t_cluster_data
    integer(kind_particle) :: first !< index of the cluster`s first particle in the particle_list
    integer(kind_particle) :: last  !< index of the cluster`s last  particle in the particle_list
  end type

  type t_thread_cluster_data
    type(t_cluster_data) :: d
    real*8 :: cluster_centre(3) !< cluster center, already shifted by vbox
  end type

  type(t_cluster_data), allocatable, dimension(:) :: particle_clusters
  integer(kind_particle) :: max_particles_per_cluster = 32
  type(t_atomic_int), pointer :: next_unassigned_cluster

  namelist /walk_para_pthreads_clustered/ max_clusters_per_thread, max_particles_per_cluster

  public tree_walk_run
  public tree_walk_init
  public tree_walk_uninit
  public tree_walk_finalize
  public tree_walk_prepare
  public tree_walk_statistics
  public tree_walk_read_parameters
  public tree_walk_write_parameters

  contains

  !>
  !> writes walk-specific data to I/O unit `u`
  !>
  subroutine tree_walk_statistics(u)
    use treevars, only: me, num_pe, MPI_COMM_lpepc
    implicit none
    include 'mpif.h'

    integer, intent(in) :: u

    integer :: i
    integer(kind_default) :: ierr
    real*8 :: average_interactions, average_mac_evaluations, total_interactions, total_mac_evaluations, max_interactions, &
      max_mac_evaluations
    real*8 :: work_imbal = 0.
    real*8 :: work_imbal_max, work_imbal_min  ! load stats
    integer(kind_node) :: local_counters(NUM_THREAD_COUNTERS)
    integer(kind_node), allocatable :: global_counters(:,:)
    integer :: c

    allocate(global_counters(NUM_THREAD_COUNTERS, num_pe))

    do c=1,NUM_THREAD_COUNTERS
      local_counters(c) = sum(threaddata(:)%counters(c))
    end do

    call MPI_GATHER(local_counters, NUM_THREAD_COUNTERS, MPI_KIND_NODE, &
                   global_counters, NUM_THREAD_COUNTERS, MPI_KIND_NODE, 0, MPI_COMM_lpepc, ierr)

    total_interactions       = SUM(global_counters(THREAD_COUNTER_INTERACTIONS,:))
    total_mac_evaluations    = SUM(global_counters(THREAD_COUNTER_MAC_EVALUATIONS,:))
    max_interactions         = MAXVAL(global_counters(THREAD_COUNTER_INTERACTIONS,:))
    max_mac_evaluations      = MAXVAL(global_counters(THREAD_COUNTER_MAC_EVALUATIONS,:))
    average_interactions     = total_interactions    / num_pe
    average_mac_evaluations  = total_mac_evaluations / num_pe
    work_imbal_max = max_interactions / average_interactions
    work_imbal_min = MINVAL(global_counters(THREAD_COUNTER_INTERACTIONS,:)) / average_interactions
    work_imbal = 0.
    do i = 1, num_pe
      work_imbal = work_imbal + abs(global_counters(THREAD_COUNTER_INTERACTIONS,i) - average_interactions) / average_interactions / num_pe
    end do

    if (0 == me) then
      write (u,*) '######## WORKLOAD AND WALK ################################################################'
      write (u,'(a50,3e12.4)')       'total/ave/max_local # interactions(work): ', total_interactions, average_interactions, max_interactions
      write (u,'(a50,3e12.4)')       'total/ave/max_local # mac evaluations: ', total_mac_evaluations, average_mac_evaluations, max_mac_evaluations
      write (u,'(a50,3f12.3)')       'Load imbalance percent,min,max: ',work_imbal,work_imbal_min,work_imbal_max
      write (u,*) '######## TREE TRAVERSAL MODULE ############################################################'
      write (u,'(a50,2i12)') 'walk_threads, max_nclusters_per_thread: ', num_walk_threads, max_clusters_per_thread
      write (u,*) '######## DETAILED DATA ####################################################################'
      write (u,'(a)') '        PE  #interactions     #mac_evals    #posted_req  rel.work'
      do i = 1, num_pe
        write (u,'(i10,3i15,F10.4)') i-1, global_counters(THREAD_COUNTER_INTERACTIONS,i), global_counters(THREAD_COUNTER_MAC_EVALUATIONS,i), &
          global_counters(THREAD_COUNTER_POST_REQUEST,i), 1._8*global_counters(THREAD_COUNTER_INTERACTIONS,i) / average_interactions
      end do
    end if
  end subroutine tree_walk_statistics

  !>
  !> reads walk specific parameters from file
  !>
  subroutine tree_walk_read_parameters(filehandle)
    use module_debug
    implicit none
    integer, intent(in) :: filehandle

    call pepc_status("READ PARAMETERS, section walk_para_pthreads")
    read(filehandle, NML=walk_para_pthreads_clustered)
  end subroutine


  !>
  !> writes walk specific parameters to file
  !>
  subroutine tree_walk_write_parameters(filehandle)
    use module_debug, only: pepc_status
    implicit none
    integer, intent(in) :: filehandle

    write(filehandle, NML=walk_para_pthreads_clustered)
  end subroutine


  !>
  !> computes derived parameters for tree walk
  !>
  subroutine tree_walk_prepare()
    implicit none
  end subroutine


  !>
  !> finalizes walk, currently this is not needed by this walk-type,
  !> but needs to be implemented in the module_walk
  !>
  subroutine tree_walk_finalize()
    implicit none
    deallocate(threaddata) ! this cannot be deallocated in tree_walk_uninit since tree_walk_statistics might want to access the data
  end subroutine tree_walk_finalize


  !>
  !> calculates forces due to the sources contained in tree `t` on the particles
  !> `p` by performing a B-H tree traversal
  !>
  subroutine tree_walk_run(vbox_)
    use, intrinsic :: iso_c_binding
    use pthreads_stuff, only: pthreads_createthread, pthreads_jointhread, THREAD_TYPE_WORKER
    use module_atomic_ops, only: atomic_store_int
    use module_pepc_types
    use module_timings
    use module_debug
    implicit none
    include 'mpif.h'

    real*8, intent(in) :: vbox_(3) !< real space shift vector of box to be processed

    integer :: ith
    integer(kind_particle) :: num_processed_clusters

    call pepc_status('WALK HYBRID')
    ! box shift vector
    vbox = vbox_

    DEBUG_ASSERT(size(threaddata)==num_walk_threads)

    threaddata(1:num_walk_threads)%finished = .false. ! we do not do this within the following loop because all (!) entries have to be .false. before the first (!) thread starts

    call atomic_store_int(next_unassigned_cluster, 1)

    ! start the worker threads...
    do ith = 1, num_walk_threads
      threaddata(ith)%id = ith
      ERROR_ON_FAIL_MSG(pthreads_createthread(thread_handles(ith), c_funloc(walk_worker_thread), c_loc(threaddata(ith)), thread_type = THREAD_TYPE_WORKER, counter = ith), "Consider setting environment variable BG_APPTHREADDEPTH=2 if you are using BG/P.")
    end do

    ! ... and wait for work thread completion
    do ith = 1, num_walk_threads
      ERROR_ON_FAIL(pthreads_jointhread(thread_handles(ith)))

      if (dbg(DBG_WALKSUMMARY)) then
        DEBUG_INFO(*, "Hybrid walk finished for thread", ith, ". Returned data = ", threaddata(ith))
      end if
    end do

    ! check wether all particles really have been processed
    num_processed_clusters = sum(threaddata(:)%counters(THREAD_COUNTER_PROCESSED_CLUSTERS))
    if (num_processed_clusters .ne. num_clusters) then
      DEBUG_ERROR(*, "Serious issue on PE", walk_tree%comm_env%rank, ": all walk threads have terminated, but obviously not all clusters are finished with walking: num_processed_clusters =", num_processed_clusters, " num_clusters =", num_clusters)
    end if
  end subroutine tree_walk_run


  subroutine tree_walk_init(t, p, num_threads)
    use, intrinsic :: iso_c_binding
    use module_atomic_ops, only: atomic_allocate_int, atomic_store_int
    use module_debug
    implicit none
    type(t_tree), target, intent(inout) :: t !< a B-H tree
    type(t_particle), target, intent(in) :: p(:) !< a list of particles
    integer, intent(in) :: num_threads !< number of traversal threads to be used

    call pepc_status('WALK INIT')

    ! we have to have at least one walk thread
    num_walk_threads = max(num_threads, 1)

    ! to keep this alive for tree_walk_statistics(), this is only deallocated in tree_walk_finalize() -> we do not know if it is allocated here or not
    if (allocated(threaddata)) deallocate(threaddata)
    allocate(threaddata(num_walk_threads))

    particle_data => p
    walk_tree     => t

    ! defer-list per particle (estimations) - will set total_defer_list_length = defer_list_length*my_max_particles_per_thread later
    defer_list_length = max(int(t%nintmax), 10)
    ! in worst case, each entry in the defer list can spawn 8 children in the todo_list
    todo_list_length  = 8 * defer_list_length

    ! find particle clusters
    call identify_particle_clusters()

    ! initialize atomic variables
    call atomic_allocate_int(next_unassigned_cluster)

    if (.not. associated(next_unassigned_cluster)) then
      DEBUG_ERROR(*, "atomic_allocate_int() failed!")
    end if

    ! evenly balance particles to threads if there are less than the maximum
    max_clusters_per_thread = max(min(num_clusters/num_walk_threads, max_clusters_per_thread),1)
    ! allocate storage for thread handles
    allocate(thread_handles(num_walk_threads))
  end subroutine tree_walk_init


  subroutine tree_walk_uninit(t, p)
    use module_atomic_ops, only: atomic_deallocate_int
    use module_debug
    implicit none
    type(t_tree), target, intent(inout) :: t !< a B-H tree
    type(t_particle), target, intent(inout) :: p(:) !< a list of particles

    call pepc_status('WALK UNINIT')

    call atomic_deallocate_int(next_unassigned_cluster)
    deallocate(thread_handles)
    deallocate(particle_clusters)
  end subroutine tree_walk_uninit


  subroutine identify_particle_clusters()
    use module_pepc_types, only: kind_particle, kind_key
    use module_spacefilling, only : level_from_key
    use module_vtk
    use module_vtk_helpers
    use module_math_tools, only : bpi
    implicit none
    include 'mpif.h'

    integer(kind_particle) :: i, k, np
    integer(kind_key) :: cut_key

    ! FIXME: we have to assure, that all particle keys are available and already have been sorted
    !        we can do this now if we want to, but should not forget to restore original order afterwards!
    np = size(particle_data, kind=kind(np))

    do k=1,2
      ! 1st loop run: count the number of clusters to be identified
      ! 2nd loop run: put them onto the particle_clusters list
      num_clusters   = 0
      i = 1
      do while (i<=np)
        num_clusters = num_clusters + 1
        if (k==2) particle_clusters(num_clusters)%first = i
        if (i<np) then
          cut_key = bpi(particle_data(i)%key, particle_data(min(i+max_particles_per_cluster,np))%key)
          do while ((iand(particle_data(i)%key, cut_key) .ne. cut_key) .and. i<np)
            i = i + 1
          end do
        else
          i = i+1
        end if
        if (k==2) particle_clusters(num_clusters)%last = i - 1
      end do
      ! allocate enough space for storing them
      if (k==1) allocate(particle_clusters(num_clusters))
    end do
  end subroutine


  function walk_worker_thread(arg) bind(c)
    use, intrinsic :: iso_c_binding
    use pthreads_stuff
    use module_interaction_specific
    use module_debug
    use module_atomic_ops
    use module_pepc_types
    use treevars, only: main_thread_processor_id
    implicit none
    include 'mpif.h'

    type(c_ptr) :: walk_worker_thread
    type(c_ptr), value :: arg

    integer, dimension(:), allocatable :: thread_cluster_indices
    integer(kind_node), dimension(:), allocatable :: partner_leaves ! list for storing number of interaction partner leaves
    integer(kind_node), dimension(:), pointer :: defer_list_old, defer_list_new, ptr_defer_list_old, ptr_defer_list_new
    integer, dimension(:), allocatable :: defer_list_start_pos
    integer :: defer_list_entries_new, defer_list_entries_old, total_defer_list_length
    integer :: defer_list_new_tail
    integer(kind_node), dimension(:), allocatable :: todo_list
    integer :: i
    logical :: clusters_available
    logical :: clusters_active
    type(t_threaddata), pointer :: my_threaddata
    logical :: shared_core
    integer :: my_max_clusters_per_thread
    integer :: my_processor_id
    logical :: cluster_has_finished

    integer(kind_node), dimension(1), target :: defer_list_root_only ! start at root node (addr, and key)
    defer_list_root_only(1) = walk_tree%node_root

    my_processor_id = get_my_core()
    shared_core = (my_processor_id == walk_tree%communicator%processor_id) .or. &
                  (my_processor_id == main_thread_processor_id)

    if ((shared_core) .and. (num_walk_threads > 1)) then
          my_max_clusters_per_thread = max(int(work_on_communicator_cluster_number_factor * max_clusters_per_thread), 1)
    else
          my_max_clusters_per_thread = max_clusters_per_thread
    end if

    call c_f_pointer(arg, my_threaddata)
    my_threaddata%is_on_shared_core = shared_core
    my_threaddata%coreid = my_processor_id
    my_threaddata%finished = .false.
    my_threaddata%counters = 0

    if (my_max_clusters_per_thread > 0) then
      total_defer_list_length = defer_list_length*my_max_clusters_per_thread

      allocate(thread_cluster_indices(my_max_clusters_per_thread), &
                      defer_list_start_pos(my_max_clusters_per_thread+1), &
                          partner_leaves(my_max_clusters_per_thread))
      allocate(defer_list_old(1:total_defer_list_length), &
                defer_list_new(1:total_defer_list_length) )
      allocate(todo_list(0:todo_list_length - 1))

      thread_cluster_indices(:) = -1     ! no clusters assigned to this thread
      clusters_available        = .true. ! but there might be clusters to be picked by the thread
      clusters_active           = .false.

      do while (clusters_active .or. clusters_available)

        call swap_defer_lists() ! swap _old and _new - lists
                                ! we will always read entries from _old and write/copy entries to _new and swap again later

        clusters_active = .false.

        ! after processing a number of clusters: handle control to other (possibly comm) thread
        if (shared_core) then
          ERROR_ON_FAIL(pthreads_sched_yield())
        end if

        do i=1,my_max_clusters_per_thread

          if (contains_cluster(i)) then
            call setup_defer_list(i)
          else
            call get_new_cluster_and_setup_defer_list(i)
          end if

          if (contains_cluster(i)) then

            ptr_defer_list_new      => defer_list_new(defer_list_new_tail:total_defer_list_length)
            defer_list_start_pos(i) =  defer_list_new_tail

            cluster_has_finished  = walk_single_cluster(particle_clusters(thread_cluster_indices(i)), &
                                      ptr_defer_list_old, defer_list_entries_old, &
                                      ptr_defer_list_new, defer_list_entries_new, &
                                      todo_list, partner_leaves(i), my_threaddata)

            if (cluster_has_finished) then
              ! walk for cluster i has finished
              ! check whether the cluster really interacted with all other clusters
              if (partner_leaves(i) .ne. walk_tree%npart) then
                write(*,'("Algorithmic problem on PE", I7, ": cluster ", I10)') walk_tree%comm_env%rank, thread_cluster_indices(i)
                write(*,'("should have been interacting (directly or indirectly) with", I16," leaves (clusters), but did with", I16)') walk_tree%npart, partner_leaves(i)
                write(*,*) "Its force and potential will be wrong due to some algorithmic error during tree traversal. Continuing anyway"
                call debug_mpi_abort()
              end if

              ! mark cluster entry i as free
              thread_cluster_indices(i)                = -1
              ! count total processed clusters for this thread
              my_threaddata%counters(THREAD_COUNTER_PROCESSED_CLUSTERS) = my_threaddata%counters(THREAD_COUNTER_PROCESSED_CLUSTERS) + 1
            else
              ! walk for cluster i has not been finished
              defer_list_new_tail = defer_list_new_tail + defer_list_entries_new
              clusters_active    = .true.
            end if

            if (defer_list_new_tail > total_defer_list_length) then
              DEBUG_ERROR('("defer_list is full for cluster ", I20, " defer_list_length =", I6, ", total =", I0," is too small (you should increase interaction_list_length_factor)")', i, defer_list_length, total_defer_list_length)
            end if
          else
            ! there is no cluster to process at position i, set the corresponding defer list to size 0
            defer_list_start_pos(i) = defer_list_new_tail
          end if
        end do ! i=1,my_max_clusters_per_thread

        defer_list_start_pos(my_max_clusters_per_thread+1) = defer_list_new_tail ! this entry is needed to store the length of the (max_clusters_per_thread)th clusters defer_list
      end do

      deallocate(thread_cluster_indices, defer_list_start_pos, partner_leaves)
      deallocate(defer_list_old, defer_list_new)
      deallocate(todo_list)
    end if

    my_threaddata%finished = .true.

    walk_worker_thread = c_null_ptr
    ERROR_ON_FAIL(pthreads_exitthread())

    contains

    subroutine swap_defer_lists()
      use module_pepc_types, only: kind_node
      implicit none
      integer(kind_node), dimension(:), pointer :: tmp_list

      tmp_list       => defer_list_old
      defer_list_old => defer_list_new
      defer_list_new => tmp_list

      defer_list_new_tail = 1 ! position of first free entry in defer_list_new (i.e. it is considered as empty now)
    end subroutine


    logical function contains_cluster(idx)
      implicit none
      integer, intent(in) :: idx

      contains_cluster = ( thread_cluster_indices(idx) .ge. 0 )
    end function contains_cluster


    function get_first_unassigned_cluster()
      use module_atomic_ops, only: atomic_fetch_and_increment_int, atomic_store_int
      implicit none
      integer :: get_first_unassigned_cluster

      integer :: next_unassigned_cluster_local

      next_unassigned_cluster_local = atomic_fetch_and_increment_int(next_unassigned_cluster)

      if (next_unassigned_cluster_local < num_clusters + 1) then
        get_first_unassigned_cluster = next_unassigned_cluster_local
      else
        get_first_unassigned_cluster = -1
      end if
    end function get_first_unassigned_cluster


    subroutine get_new_cluster_and_setup_defer_list(idx)
      implicit none
      integer, intent(in) :: idx

      if (clusters_available) then
        thread_cluster_indices(idx) = get_first_unassigned_cluster()

        if (contains_cluster(idx)) then
          ! for clusters that we just inserted into our list, we start with only one defer_list_entry: the root node
          ptr_defer_list_old      => defer_list_root_only
          defer_list_entries_old  =  1
          partner_leaves(idx)     =  0 ! no interactions yet
        else
          clusters_available     = .false.
        end if ! contains_cluster(idx)
      end if ! clusters_available
    end subroutine get_new_cluster_and_setup_defer_list


    subroutine setup_defer_list(idx)
      implicit none
      integer, intent(in) :: idx

      ptr_defer_list_old     => defer_list_old(defer_list_start_pos(idx):defer_list_start_pos(idx+1)-1)
      defer_list_entries_old =  defer_list_start_pos(idx+1) - defer_list_start_pos(idx)
    end subroutine setup_defer_list
  end function walk_worker_thread


  function walk_single_cluster(cluster, defer_list_old, defer_list_entries_old, &
                                        defer_list_new, defer_list_entries_new, &
                                        todo_list, partner_leaves, my_threaddata)
    use module_tree_node
    use module_tree_communicator, only: tree_node_fetch_children
    use module_interaction_specific
    use module_debug
    #ifndef NO_SPATIAL_INTERACTION_CUTOFF
    use module_mirror_boxes, only : spatial_interaction_cutoff
    #endif
    use module_atomic_ops
    use module_pepc_types
    implicit none
    include 'mpif.h'

    type(t_cluster_data), intent(in) :: cluster
    integer(kind_node), dimension(:), pointer, intent(in) :: defer_list_old
    integer, intent(in) :: defer_list_entries_old
    integer(kind_node), dimension(:), pointer, intent(out) :: defer_list_new
    integer, intent(out) :: defer_list_entries_new
    integer(kind_node), intent(inout) :: todo_list(0:todo_list_length-1)
    integer(kind_node), intent(inout) :: partner_leaves
    type(t_threaddata), intent(inout) :: my_threaddata
    logical :: walk_single_cluster !< function will return .true. if this cluster has finished its walk

    interface t_calc_force
      subroutine t_calc_force(particle, node, node_idx, delta, dist2, vbox)
        use module_pepc_types
        use treevars
        implicit none
        type(t_tree_node_interaction_data), intent(in) :: node
        integer(kind_node), intent(in) :: node_idx
        type(t_particle), intent(inout) :: particle
        real*8, intent(in) :: vbox(3), delta(3), dist2
      end subroutine
    end interface

    integer :: todo_list_entries
    type(t_tree_node), pointer :: walk_node
    integer(kind_node) :: walk_node_idx
    real*8 :: dist2, delta(3)
    logical :: is_leaf
    integer(kind_node) :: num_interactions, num_mac_evaluations, num_post_request
    integer(kind_particle) :: p

    todo_list_entries      = 0
    num_interactions       = 0
    num_mac_evaluations    = 0
    num_post_request       = 0
    walk_node_idx          = NODE_INVALID

    ! for each entry on the defer list, we check, whether children are already available and put them onto the todo_list
    ! another mac-check for each entry is not necessary here, since due to having requested the children, we already know,
    ! that the node has to be resolved
    ! if the defer_list is empty, the call returns without doing anything
    call defer_list_parse_and_compact()

    ! read all todo_list-entries and start further traversals there
    todo_list_loop: do while (todo_list_pop(walk_node_idx))
      walk_node => walk_tree%nodes(walk_node_idx)

      ! we may not interact with the cluster itself or its ancestors
      ! if we are in the central box
      ! interaction with ancestor nodes should be prevented by the MAC
      ! but this does not always work (i.e. if theta > 0.7 or if keys and/or coordinates have
      ! been modified due to 'duplicate keys'-error)
      is_leaf = tree_node_is_leaf(walk_node)


      if (is_leaf) then
        partner_leaves = partner_leaves + 1
        call calc_force_for_cluster(calc_force_per_interaction_with_leaf, calc_force_per_interaction_with_self)
      else ! not a leaf, evaluate MAC
        do p=cluster%first,cluster%last
          delta = particle_data(p)%x - vbox - walk_node%interaction_data%coc ! Separation vector
          dist2 = DOT_PRODUCT(delta, delta)

          num_mac_evaluations = num_mac_evaluations + 1
          if (.not. mac(walk_node%interaction_data, dist2, walk_tree%boxlength2(walk_node%level))) then ! MAC negative, resolve
            call resolve(particle_data(p))
            cycle todo_list_loop ! we do not have to evaluate the MAC for any other particles in this cluster
          endif

        end do

        ! MAC was positive for all particles in the cluster, interact
        partner_leaves = partner_leaves + walk_node%leaves
        call calc_force_for_cluster(calc_force_per_interaction_with_twig, calc_force_per_interaction_with_twig)
      end if
    end do todo_list_loop ! (while (todo_list_pop(walk_key)))

    ! if todo_list and defer_list are now empty, the walk has finished
    walk_single_cluster = (todo_list_entries == 0) .and. (defer_list_entries_new == 0)

    my_threaddata%counters(THREAD_COUNTER_INTERACTIONS)    = my_threaddata%counters(THREAD_COUNTER_INTERACTIONS)    + num_interactions
    my_threaddata%counters(THREAD_COUNTER_MAC_EVALUATIONS) = my_threaddata%counters(THREAD_COUNTER_MAC_EVALUATIONS) + num_mac_evaluations
    my_threaddata%counters(THREAD_COUNTER_POST_REQUEST)    = my_threaddata%counters(THREAD_COUNTER_POST_REQUEST)    + num_post_request

    contains

    ! FIXME: self with twig does not make too much sense but makes the code much more beautiful
    ! TODO: this should be vectorized or whatever. put this into interaction-specific bla
    subroutine calc_force_for_cluster(calc_force, calc_force_self)
      implicit none
      procedure(t_calc_force) :: calc_force, calc_force_self
      real*8 :: pdist2, pdelta(3)
      integer(kind_particle) :: ipart

      do ipart=cluster%first,cluster%last
        pdelta = particle_data(ipart)%x - vbox - walk_node%interaction_data%coc ! Separation vector
        pdist2 = DOT_PRODUCT(pdelta, pdelta)

        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (any(abs(pdelta) >= spatial_interaction_cutoff)) cycle
        #endif

        if (pdist2 > 0.0_8) then ! not self, interact
          call calc_force(particle_data(ipart), walk_node%interaction_data, walk_node_idx, pdelta, pdist2, vbox)
          else ! self, count as interaction partner, otherwise ignore
          call calc_force_self(particle_data(ipart), walk_node%interaction_data, walk_node_idx, pdelta, pdist2, vbox)
        endif
        num_interactions = num_interactions + 1
        particle_data(ipart)%work = particle_data(ipart)%work + 1._8
      end do
    end subroutine


    subroutine resolve(particle)
      implicit none
      type(t_particle), intent(in) :: particle

      integer(kind_node) :: n

      ! resolve
      n = tree_node_get_first_child(walk_node)
      if (n /= NODE_INVALID) then
        ! children for twig are present
        ! --> resolve cell & put all children in front of todo_list
        if (.not. todo_list_push_siblings(n)) then
          ! the todo_list is full --> put parent back onto defer_list
          call defer_list_push(walk_node_idx)
        end if
      else
        ! children for twig are _absent_
        ! --> put node on REQUEST list and put walk_key on bottom of todo_list
        ! eager requests
        call tree_node_fetch_children(walk_tree, walk_node, walk_node_idx, particle, particle%x - vbox) ! fetch children from remote
        ! simple requests
        ! call tree_node_fetch_children(walk_tree, walk_node, walk_node_idx)
        num_post_request = num_post_request + 1
        ! if posting the request failed, this is not a problem, since we defer the particle anyway
        ! since it will not be available then, the request will simply be repeated
        call defer_list_push(walk_node_idx) ! Deferred list of nodes to search, pending request
                                            ! for data from nonlocal PEs
      end if
    end subroutine resolve


    ! helper routines for todo_list manipulation
    function todo_list_pop(node)
      implicit none

      logical :: todo_list_pop
      integer(kind_node), intent(out) :: node

      todo_list_pop = (todo_list_entries > 0)

      if (todo_list_pop) then
        todo_list_entries = todo_list_entries - 1
        node = todo_list(todo_list_entries)
      end if
    end function


    function todo_list_push_siblings(node) result(res)
      implicit none

      integer(kind_node), intent(in) :: node

      integer(kind_node) :: n
      logical :: res

      ! check for enough space on todo_list
      res = (todo_list_entries + 8 <= todo_list_length)
      if (res) then
        n = node
        do
          todo_list(todo_list_entries) = n
          todo_list_entries = todo_list_entries + 1
          n = tree_node_get_next_sibling(walk_tree%nodes(n))
          if (n == NODE_INVALID) exit
        end do
      else
        DEBUG_WARNING_ALL('("todo_list is full for one cluster, todo_list_length =", I6, " is too small (you should increase interaction_list_length_factor). Putting particles back onto defer_list. Programme will continue without errors.")', todo_list_length)
      end if
    end function


    subroutine defer_list_push(node)
      implicit none

      integer(kind_node), intent(in) :: node

      defer_list_entries_new = defer_list_entries_new + 1
      defer_list_new(defer_list_entries_new) = node
    end subroutine


    subroutine defer_list_parse_and_compact()
      implicit none

      integer(kind_node) :: n
      integer :: iold

      defer_list_entries_new = 0
      iold = 1
      do
        if (iold > defer_list_entries_old) return
        n = tree_node_get_first_child(walk_tree%nodes(defer_list_old(iold)))
        if (n /= NODE_INVALID) then
          ! children for deferred node have arrived --> put children onto todo_list
          if (.not. todo_list_push_siblings(n)) exit
        else
          ! children for deferred node are still unavailable - put onto defer_list_new (do not use defer_list_push for performance reasons)
          defer_list_entries_new                 = defer_list_entries_new + 1
          defer_list_new(defer_list_entries_new) = defer_list_old(iold)
        end if
        iold = iold + 1
      end do

      do
        if (iold > defer_list_entries_old) return
        defer_list_entries_new                 = defer_list_entries_new + 1
        defer_list_new(defer_list_entries_new) = defer_list_old(iold)
        iold = iold + 1
      end do
    end subroutine
  end function walk_single_cluster
end module module_walk
