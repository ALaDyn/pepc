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
!>        do i=1,max_particles_per_thread
!>
!>          if ( (my_particles(i) == -1) .and. (particles_available)) then         ! i.e. the place for a particle is unassigned
!>            my_particles(i) = get_first_unassigned_particle()
!>          end if
!>
!>          call walk_single_particle(my_particles(i))
!>
!>        end do
!>      end do
!>
!>
!>  Structure of walk_single_particle(particle):
!>  ------------------------------------------
!>
!>     if (.not.finished(particle)) then
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
!>
module module_walk_pthreads_commutils
  use module_atomic_ops, only: t_atomic_int
  implicit none
  private

  !> debug flags - cannot be modified at runtime due to performance reasons
  logical, parameter, public  :: walk_debug     = .false.

  ! variables for adjusting the thread's workload
  real, public :: work_on_communicator_particle_number_factor = 0.1 !< factor for reducing max_particles_per_thread for thread which share their processor with the communicator
  integer, public :: particles_per_yield = 500 !< number of particles to process in a work_thread before it shall call sched_yield to hand the processor over to some other thread
  ! atomic variables for regulating concurrent access
  type(t_atomic_int), pointer, public :: next_unassigned_particle
  type(t_atomic_int), pointer, public :: threads_finished

  ! local walktime (i.e. from comm_loop start until send_walk_finished() )
  real*8, public, pointer :: twalk_loc

  integer, public, parameter :: NUM_THREAD_TIMERS                 = 4
  integer, public, parameter :: THREAD_TIMER_TOTAL                = 1
  integer, public, parameter :: THREAD_TIMER_POST_REQUEST         = 2
  integer, public, parameter :: THREAD_TIMER_GET_NEW_PARTICLE     = 3
  integer, public, parameter :: THREAD_TIMER_WALK_SINGLE_PARTICLE = 4

  integer, public, parameter :: NUM_THREAD_COUNTERS                = 4
  integer, public, parameter :: THREAD_COUNTER_PROCESSED_PARTICLES = 1
  integer, public, parameter :: THREAD_COUNTER_INTERACTIONS        = 2
  integer, public, parameter :: THREAD_COUNTER_MAC_EVALUATIONS     = 3
  integer, public, parameter :: THREAD_COUNTER_POST_REQUEST        = 4

  !> type for input and return values of walk_threads
  type, public :: t_threaddata
    integer :: id                         !< just a running number to distinguish the threads, currently unused
    logical :: is_on_shared_core !< thread output value: is set to true if the thread detects that it shares its processor with the communicator thread
    integer :: coreid !< thread output value: id of thread`s processor
    logical :: finished !< will be set to .true. when the thread has finished
    real*8 :: timers(NUM_THREAD_TIMERS)
    integer*8 :: counters(NUM_THREAD_COUNTERS)
  end type t_threaddata

  type(t_threaddata), public, allocatable, target :: threaddata(:)

  real*8, public :: thread_timers_nonshared_avg(NUM_THREAD_TIMERS)
  real*8, public :: thread_timers_nonshared_dev(NUM_THREAD_TIMERS)
  real*8, public :: thread_timers_shared_avg(NUM_THREAD_TIMERS)
  real*8, public :: thread_timers_shared_dev(NUM_THREAD_TIMERS)
  real*8, public :: thread_counters_nonshared_avg(NUM_THREAD_COUNTERS)
  real*8, public :: thread_counters_nonshared_dev(NUM_THREAD_COUNTERS)
  real*8, public :: thread_counters_shared_avg(NUM_THREAD_COUNTERS)
  real*8, public :: thread_counters_shared_dev(NUM_THREAD_COUNTERS)

  integer, public :: num_nonshared_threads, num_shared_threads

  public comm_sched_yield
  public retval
  public init_commutils
  public uninit_commutils

  contains
  
  subroutine init_commutils()
    use module_atomic_ops, only: atomic_allocate_int, atomic_store_int
    use module_debug
    implicit none

    ! initialize atomic variables
    call atomic_allocate_int(next_unassigned_particle)
    call atomic_allocate_int(threads_finished)
    if (.not. (associated(next_unassigned_particle) .and. associated(threads_finished))) then
      DEBUG_ERROR(*, "atomic_allocate_int() failed!")
    end if

    call atomic_store_int(next_unassigned_particle, 1)
    call atomic_store_int(threads_finished,         0)
  end subroutine init_commutils


  subroutine uninit_commutils
    use module_atomic_ops, only: atomic_deallocate_int
    implicit none

    call atomic_deallocate_int(threads_finished)
    call atomic_deallocate_int(next_unassigned_particle)
  end subroutine uninit_commutils


  subroutine retval(iret, msg)
    use module_debug
    use, intrinsic :: iso_c_binding
    implicit none
    integer( kind= c_int) :: iret
    character(*), intent(in) :: msg

    if (iret .ne. 0) then
      DEBUG_ERROR('("[",a,"] iret = ", I0)',msg, iret)
    end if
  end subroutine retval


  subroutine comm_sched_yield()
    use pthreads_stuff
    implicit none

    call retval(pthreads_sched_yield(), "pthreads_sched_yield()")
  end subroutine comm_sched_yield
end module module_walk_pthreads_commutils



module module_walk
  use, intrinsic :: iso_c_binding
  use module_tree, only: t_tree
  use module_interaction_specific
  use pthreads_stuff
  use module_pepc_types
  implicit none

  private
  
  integer, public :: max_particles_per_thread = 2000 !< maximum number of particles that will in parallel be processed by one workthread

  integer, public :: num_walk_threads = -1 !< number of worker threads, default value is set to treevars%num_threads in tree_walk_read_parameters()
  integer, private :: primary_processor_id = 0

  real*8, dimension(:), allocatable :: boxlength2
  real*8 :: vbox(3)
  logical :: in_central_box
  integer*8 :: num_interaction_leaves
  integer :: todo_list_length, defer_list_length

  integer :: num_particles
  type(t_particle), pointer, dimension(:) :: particle_data
  type(t_tree), pointer :: walk_tree

  type(c_ptr), allocatable :: thread_handles(:)

  type t_defer_list_entry
    type(t_tree_node), pointer :: node
    integer*8 :: key
  end type t_defer_list_entry
  
  namelist /walk_para_pthreads/ num_walk_threads, max_particles_per_thread

  public tree_walk
  public tree_walk_finalize
  public tree_walk_prepare
  public tree_walk_statistics
  public tree_walk_read_parameters
  public tree_walk_write_parameters

  contains

  !>
  !> writes walk-specific data to file steam ifile
  !>
  subroutine tree_walk_statistics(t, ifile, perform_output)
    use module_walk_pthreads_commutils
    use module_tree, only: t_tree
    implicit none
    include 'mpif.h'

    type(t_tree), intent(in) :: t
    integer, intent(in) :: ifile !< file stream to write to
    logical, intent(in) :: perform_output !< if set to false, output is disabled (e.g. for MPI ranks that shall not print anything)

    integer :: ierr
    real*8 :: global_thread_timers_nonshared_avg(NUM_THREAD_TIMERS)
    real*8 :: global_thread_timers_nonshared_dev(NUM_THREAD_TIMERS)
    real*8 :: global_thread_timers_shared_avg(NUM_THREAD_TIMERS)
    real*8 :: global_thread_timers_shared_dev(NUM_THREAD_TIMERS)
    real*8 :: global_thread_counters_nonshared_avg(NUM_THREAD_COUNTERS)
    real*8 :: global_thread_counters_nonshared_dev(NUM_THREAD_COUNTERS)
    real*8 :: global_thread_counters_shared_avg(NUM_THREAD_COUNTERS)
    real*8 :: global_thread_counters_shared_dev(NUM_THREAD_COUNTERS)

    call MPI_REDUCE(thread_timers_nonshared_avg(:), global_thread_timers_nonshared_avg(:), NUM_THREAD_TIMERS, MPI_REAL8, MPI_SUM, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_timers_shared_avg(:), global_thread_timers_shared_avg(:), NUM_THREAD_TIMERS, MPI_REAL8, MPI_SUM, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_counters_nonshared_avg(:), global_thread_counters_nonshared_avg(:), NUM_THREAD_COUNTERS, MPI_REAL8, MPI_SUM, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_counters_shared_avg(:), global_thread_counters_shared_avg(:), NUM_THREAD_COUNTERS, MPI_REAL8, MPI_SUM, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_timers_nonshared_dev(:), global_thread_timers_nonshared_dev(:), NUM_THREAD_TIMERS, MPI_REAL8, MPI_MAX, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_timers_shared_dev(:), global_thread_timers_shared_dev(:), NUM_THREAD_TIMERS, MPI_REAL8, MPI_MAX, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_counters_nonshared_dev(:), global_thread_counters_nonshared_dev(:), NUM_THREAD_COUNTERS, MPI_REAL8, MPI_MAX, 0, t%comm_env%comm, ierr)
    call MPI_REDUCE(thread_counters_shared_dev(:), global_thread_counters_shared_dev(:), NUM_THREAD_COUNTERS, MPI_REAL8, MPI_MAX, 0, t%comm_env%comm, ierr)

    global_thread_timers_nonshared_avg   = global_thread_timers_nonshared_avg / t%comm_env%size
    global_thread_timers_shared_avg      = global_thread_timers_shared_avg / t%comm_env%size
    global_thread_counters_nonshared_avg = global_thread_counters_nonshared_avg / t%comm_env%size
    global_thread_counters_shared_avg    = global_thread_counters_shared_avg / t%comm_env%size

    if (perform_output) then
      write (ifile,'(a50,2i12)') 'walk_threads, max_nparticles_per_thread: ', num_walk_threads, max_particles_per_thread
      ! TODO: this goes somewhere else
      !write (ifile,'(a50,3i12)') '# of comm-loop iterations (tot,send,recv): ', comm_loop_iterations(:)
      write (ifile,*) '######## WALK-WORKER-THREAD WORKLOAD ######################################################'
      write (ifile,'(a50)')              'average # processed nparticles per thread    '
      write (ifile,'(a50,3f12.3)')       '  threads on exclusive cores, shared cores: ', &
                                          thread_counters_nonshared_avg(THREAD_COUNTER_PROCESSED_PARTICLES), &
                                          thread_counters_shared_avg(THREAD_COUNTER_PROCESSED_PARTICLES)
      write (ifile,'(a50,3f12.3)')       '  maximum relative deviation: ', &
                                          thread_counters_nonshared_dev(THREAD_COUNTER_PROCESSED_PARTICLES), &
                                          thread_counters_shared_dev(THREAD_COUNTER_PROCESSED_PARTICLES)
      write (ifile,'(a50)')              'average wallclocktime per thread    '
      write (ifile,'(a50,3f12.3)')       '  threads on exclusive cores, shared cores: ', &
                                          thread_timers_nonshared_avg(THREAD_TIMER_TOTAL), &
                                          thread_timers_shared_avg(THREAD_TIMER_TOTAL)
      write (ifile,'(a50,3f12.3)')       '  maximum relative deviation: ', &
                                          thread_timers_nonshared_dev(THREAD_TIMER_TOTAL), &
                                          thread_timers_shared_dev(THREAD_TIMER_TOTAL)
    end if

  end subroutine


  !>
  !> reads walk specific parameters from file
  !>
  subroutine tree_walk_read_parameters(filehandle)
    use module_debug
    use treevars, only: num_threads
    implicit none
    integer, intent(in) :: filehandle

    call pepc_status("READ PARAMETERS, section walk_para_pthreads")
    read(filehandle, NML=walk_para_pthreads)

    if (num_walk_threads > 0) then ! it has been set through the namelist
      DEBUG_WARNING(*,  'Setting num_walk_threads through the walk_para_pthreads namelist directly is deprecated and will be removed soon. Please switch to parameter num_threads in namelist libpepc in your parameter files.')
      num_threads = num_walk_threads
    else          
      num_walk_threads = num_threads
    end if
  end subroutine


  !>
  !> writes walk specific parameters to file
  !>
  subroutine tree_walk_write_parameters(filehandle)
    use module_debug, only: pepc_status
    implicit none
    integer, intent(in) :: filehandle

    write(filehandle, NML=walk_para_pthreads)
  end subroutine


  !>
  !> computes derived parameters for tree walk
  !>
  subroutine tree_walk_prepare()
    use module_walk_pthreads_commutils
    implicit none
    ! nothing to do here

    !if (me == 0) then
    !  write(*,'("MPI-PThreads walk: Using ", I0," worker-threads in treewalk on each processor (i.e. per MPI rank)")') num_walk_threads
    !  write(*,'("Maximum number of particles per work_thread = ", I0)') max_particles_per_thread
    !end if
  end subroutine


  !>
  !> finilizes walk, currently this is not needed by this walk-type,
  !> but needs to be implemented in the module_walk
  !>
  subroutine tree_walk_finalize()
    implicit none
  end subroutine tree_walk_finalize


  subroutine tree_walk(t, p, twalk, twalk_loc_, vbox_)
    use, intrinsic :: iso_c_binding
    use module_pepc_types
    use module_timings
    use module_walk_pthreads_commutils
    use module_debug, only : pepc_status
    implicit none
    include 'mpif.h'

    type(t_tree), target, intent(inout) :: t
    type(t_particle), target, intent(in) :: p(:)
    real*8, intent(in) :: vbox_(3) !< real space shift vector of box to be processed
    real*8, target, intent(inout) :: twalk, twalk_loc_

    call pepc_status('WALK HYBRID')

    num_particles = size(p)
    particle_data => p
    walk_tree => t
    ! box shift vector
    vbox = vbox_
    ! defer-list per particle (estimations) - will set total_defer_list_length = defer_list_length*my_max_particles_per_thread later
    defer_list_length = max(int(t%nintmax), 10)
    ! in worst case, each entry in the defer list can spawn 8 children in the todo_list
    todo_list_length  = 8 * defer_list_length

    ! pure local walk time (i.e. from start of communicator till send_walk_finished)
    twalk_loc => twalk_loc_

    twalk  = MPI_WTIME()

    call init_walk_data(t)

    call init_commutils()

    call walk_hybrid()

    call uninit_commutils()

    call uninit_walk_data()

    twalk = MPI_WTIME() - twalk
  end subroutine tree_walk


  subroutine walk_hybrid()
    use module_debug
    use module_walk_pthreads_commutils
    use module_atomic_ops
    use treevars, only: interactions_local, mac_evaluations_local
    use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'

    integer :: ith
    integer*8 :: num_processed_particles

    allocate(threaddata(num_walk_threads))

    threaddata(1:num_walk_threads)%finished = .false. ! we do not do this within the following loop because all (!) entries have to be .false. before the first (!) thread starts

    twalk_loc = MPI_WTIME()

    ! start the worker threads...
    do ith = 1, num_walk_threads
      threaddata(ith)%id = ith
      call retval(pthreads_createthread(thread_handles(ith), c_funloc(walk_worker_thread), c_loc(threaddata(ith))), &
        "walk_schedule_thread_inner:pthread_create. Consider setting environment variable BG_APPTHREADDEPTH=2 if you are using BG/P.")
    end do

    ! ... and wait for work thread completion
    do ith = 1, num_walk_threads
      call retval(pthreads_jointhread(thread_handles(ith)), "walk_schedule_thread_inner:pthread_join")

      if (dbg(DBG_WALKSUMMARY)) then
        DEBUG_INFO(*, "Hybrid walk finished for thread", ith, ". Returned data = ", threaddata(ith))
      end if
    end do

    twalk_loc = MPI_WTIME() - twalk_loc

    ! store workload data
    call collect_thread_counters_timers()

    ! check wether all particles really have been processed
    num_processed_particles = sum(threaddata(:)%counters(THREAD_COUNTER_PROCESSED_PARTICLES))
    if (num_processed_particles .ne. num_particles) then
      DEBUG_ERROR(*, "Serious issue on PE", walk_tree%comm_env%rank, ": all walk threads have terminated, but obviously not all particles are finished with walking: num_processed_particles =",
                          num_processed_particles, " num_particles =", num_particles)
    end if

    deallocate(threaddata)

    contains

    subroutine collect_thread_counters_timers()
      implicit none

      integer :: icounter, itimer

      num_shared_threads = count(threaddata(:)%is_on_shared_core)
      num_nonshared_threads = num_walk_threads - num_shared_threads

      do itimer = 1,NUM_THREAD_TIMERS
        thread_timers_shared_avg(itimer)    =  sum(threaddata(:)%timers(itimer), mask =       threaddata(:)%is_on_shared_core) / num_shared_threads
        thread_timers_nonshared_avg(itimer) =  sum(threaddata(:)%timers(itimer), mask = .not. threaddata(:)%is_on_shared_core) / num_nonshared_threads

        thread_timers_shared_dev(itimer)    = (maxval(threaddata(:)%timers(itimer), mask =       threaddata(:)%is_on_shared_core) - thread_timers_shared_avg(itimer))    / thread_timers_shared_avg(itimer)
        thread_timers_nonshared_dev(itimer) = (maxval(threaddata(:)%timers(itimer), mask = .not. threaddata(:)%is_on_shared_core) - thread_timers_nonshared_avg(itimer)) / thread_timers_nonshared_avg(itimer)
      end do

      do icounter = 1,NUM_THREAD_COUNTERS
        thread_counters_shared_avg(icounter) = sum(threaddata(:)%counters(icounter), mask = threaddata(:)%is_on_shared_core) / real(num_shared_threads, kind = 8)
        thread_counters_nonshared_avg(icounter) = sum(threaddata(:)%counters(icounter), mask = .not. threaddata(:)%is_on_shared_core) / real(num_nonshared_threads, kind = 8)

        thread_counters_shared_dev(icounter) = (maxval(threaddata(:)%counters(icounter), mask = threaddata(:)%is_on_shared_core) - thread_counters_shared_avg(icounter)) / thread_counters_shared_avg(icounter)
        thread_counters_nonshared_dev(icounter) = (maxval(threaddata(:)%counters(icounter), mask = .not. threaddata(:)%is_on_shared_core) - thread_counters_nonshared_avg(icounter)) / thread_counters_nonshared_avg(icounter)
      end do

      interactions_local    = sum(threaddata(:)%counters(THREAD_COUNTER_INTERACTIONS))
      mac_evaluations_local = sum(threaddata(:)%counters(THREAD_COUNTER_MAC_EVALUATIONS))
    end subroutine collect_thread_counters_timers
  end subroutine walk_hybrid


  subroutine init_walk_data(t)
    use, intrinsic :: iso_c_binding
    use module_walk_pthreads_commutils
    use module_tree, only: t_tree
    use module_debug
    use treevars, only: nlev
    use pthreads_stuff, only: pthreads_alloc_thread
    implicit none

    type(t_tree), intent(in) :: t

    integer :: i

    ! we have to have at least one walk thread
    num_walk_threads = max(num_walk_threads, 1)
    ! evenly balance particles to threads if there are less than the maximum
    max_particles_per_thread = max(min(num_particles/num_walk_threads, max_particles_per_thread),1)
    ! allocate storage for thread handles
    allocate(thread_handles(num_walk_threads))
    thread_handles = c_null_ptr
    do i = 1, num_walk_threads
      thread_handles(i) = pthreads_alloc_thread()
      if (.not. c_associated(thread_handles(i))) then
        DEBUG_ERROR(*, "pthreads_alloc_thread() failed!")
      end if
    end do

    ! we will only want to reject the root node and the particle itself if we are in the central box
    in_central_box = (dot_product(vbox,vbox) == 0)
    ! every particle has directly or indirectly interact with each other, and outside the central box even with itself
    num_interaction_leaves = walk_tree%npart
    ! Preprocessed box sizes for each level
    allocate(boxlength2(0:nlev))
    boxlength2(0)=maxval(t%bounding_box%boxsize)**2
    do i = 1, nlev
      boxlength2(i) =  boxlength2(i-1)/4.
    end do

    ! store ID of primary (comm-thread) processor
    primary_processor_id = get_my_core()
  end subroutine init_walk_data


  subroutine uninit_walk_data()
    use module_walk_pthreads_commutils
    use pthreads_stuff, only: pthreads_free_thread
    implicit none

    integer :: i
    deallocate(boxlength2)

    do i = 1, num_walk_threads
      call pthreads_free_thread(thread_handles(i))
    end do
    deallocate(thread_handles)
  end subroutine uninit_walk_data


  function get_first_unassigned_particle()
    use module_walk_pthreads_commutils
    use module_atomic_ops
    implicit none
    integer :: get_first_unassigned_particle

    integer :: next_unassigned_particle_local

    next_unassigned_particle_local = atomic_fetch_and_increment_int(next_unassigned_particle)

    if (next_unassigned_particle_local < num_particles + 1) then
      get_first_unassigned_particle = next_unassigned_particle_local
    else
      call atomic_store_int(next_unassigned_particle, num_particles + 1)
      get_first_unassigned_particle = -1
    end if
  end function get_first_unassigned_particle


  function walk_worker_thread(arg) bind(c)
    use, intrinsic :: iso_c_binding
    use module_walk_pthreads_commutils
    use pthreads_stuff
    use module_interaction_specific
    use module_debug
    use module_atomic_ops
    use module_tree, only: tree_lookup_root
    implicit none
    include 'mpif.h'

    type(c_ptr) :: walk_worker_thread
    type(c_ptr), value :: arg

    integer, dimension(:), allocatable :: thread_particle_indices
    type(t_particle), dimension(:), allocatable :: thread_particle_data
    integer*8, dimension(:), allocatable :: partner_leaves ! list for storing number of interaction partner leaves
    type(t_defer_list_entry), dimension(:), pointer :: defer_list_old,           defer_list_new, ptr_defer_list_old, ptr_defer_list_new
    integer, dimension(:), allocatable :: defer_list_start_pos
    integer :: defer_list_entries_new, defer_list_entries_old, total_defer_list_length
    integer :: defer_list_new_tail
    integer*8, dimension(:), allocatable :: todo_list
    integer :: i
    logical :: particles_available
    logical :: particles_active
    type(t_threaddata), pointer :: my_threaddata
    integer :: particles_since_last_yield
    logical :: same_core_as_communicator
    integer :: my_max_particles_per_thread
    integer :: my_processor_id, num_finished
    logical :: particle_has_finished
    real*8  :: t_get_new_particle, t_walk_single_particle

    type(t_defer_list_entry), dimension(1), target :: defer_list_root_only ! start at root node (addr, and key)
    defer_list_root_only(1)%key = 1_8
    call tree_lookup_root(walk_tree, defer_list_root_only(1)%node, 'walk_worker_thread:root node')

    t_get_new_particle = 0._8
    t_walk_single_particle = 0._8

    my_processor_id = get_my_core()
    same_core_as_communicator = (my_processor_id == primary_processor_id)

    if ((same_core_as_communicator) .and. (num_walk_threads > 1)) then
          my_max_particles_per_thread = max(int(work_on_communicator_particle_number_factor * max_particles_per_thread), 1)
    else
          my_max_particles_per_thread = max_particles_per_thread
    end if

    call c_f_pointer(arg, my_threaddata)
    my_threaddata%is_on_shared_core = same_core_as_communicator
    my_threaddata%coreid = my_processor_id
    my_threaddata%finished = .false.
    my_threaddata%timers(THREAD_TIMER_TOTAL) = - MPI_WTIME()
    my_threaddata%timers(THREAD_TIMER_POST_REQUEST) = 0
    my_threaddata%counters = 0

    if (my_max_particles_per_thread > 0) then
      total_defer_list_length = defer_list_length*my_max_particles_per_thread

      allocate(thread_particle_indices(my_max_particles_per_thread), &
                    thread_particle_data(my_max_particles_per_thread), &
                      defer_list_start_pos(my_max_particles_per_thread+1), &
                          partner_leaves(my_max_particles_per_thread))
      allocate(defer_list_old(1:total_defer_list_length), &
                defer_list_new(1:total_defer_list_length) )
      allocate(todo_list(0:todo_list_length - 1))

      thread_particle_indices(:) = -1     ! no particles assigned to this thread
      particles_available        = .true. ! but there might be particles to be picked by the thread
      particles_active           = .false.
      particles_since_last_yield =  0

      do while (particles_active .or. particles_available)

        call swap_defer_lists() ! swap _old and _new - lists
                                ! we will always read entries from _old and write/copy entries to _new and swap again later

        particles_active = .false.

        do i=1,my_max_particles_per_thread

          if (contains_particle(i)) then
            call setup_defer_list(i)
          else
            t_get_new_particle = t_get_new_particle - MPI_WTIME()
            call get_new_particle_and_setup_defer_list(i)
            t_get_new_particle = t_get_new_particle + MPI_WTIME()
          end if

          if (contains_particle(i)) then

            call do_sched_yield_if_necessary()

            ptr_defer_list_new      => defer_list_new(defer_list_new_tail:total_defer_list_length)
            defer_list_start_pos(i) =  defer_list_new_tail

            t_walk_single_particle = t_walk_single_particle - MPI_WTIME()
            particle_has_finished  = walk_single_particle(thread_particle_data(i), &
                                      ptr_defer_list_old, defer_list_entries_old, &
                                      ptr_defer_list_new, defer_list_entries_new, &
                                      todo_list, partner_leaves(i), my_threaddata)
            t_walk_single_particle = t_walk_single_particle + MPI_WTIME()

            if (particle_has_finished) then
              ! walk for particle i has finished
              if (walk_debug) then
                  DEBUG_INFO('("PE", I6, " particle ", I12, " obviously finished walking around :-)")', walk_tree%comm_env%rank, i)
              end if

              ! check whether the particle really interacted with all other particles
              if (partner_leaves(i) .ne. num_interaction_leaves) then
                write(*,'("Algorithmic problem on PE", I7, ": Particle ", I10, " label ", I16)') walk_tree%comm_env%rank, thread_particle_indices(i), thread_particle_data(i)%label
                write(*,'("should have been interacting (directly or indirectly) with", I16," leaves (particles), but did with", I16)') num_interaction_leaves, partner_leaves(i)
                write(*,*) "Its force and potential will be wrong due to some algorithmic error during tree traversal. Continuing anyway"
                call debug_mpi_abort()
              end if

              ! copy forces and potentials back to thread-global array
              particle_data(thread_particle_indices(i)) = thread_particle_data(i)
              ! mark particle entry i as free
              thread_particle_indices(i)                = -1
              ! count total processed particles for this thread
              my_threaddata%counters(THREAD_COUNTER_PROCESSED_PARTICLES) = my_threaddata%counters(THREAD_COUNTER_PROCESSED_PARTICLES) + 1
            else
              ! walk for particle i has not been finished
              defer_list_new_tail = defer_list_new_tail + defer_list_entries_new
              particles_active    = .true.
            end if

            if (defer_list_new_tail > total_defer_list_length) then
              DEBUG_ERROR('("defer_list is full for particle ", I20, " defer_list_length =", I6, ", total =", I0," is too small (you should increase interaction_list_length_factor)")', i, defer_list_length, total_defer_list_length)
            end if
          else
            ! there is no particle to process at position i, set the corresponding defer list to size 0
            defer_list_start_pos(i) = defer_list_new_tail
          end if
        end do ! i=1,my_max_particles_per_thread

        defer_list_start_pos(my_max_particles_per_thread+1) = defer_list_new_tail ! this entry is needed to store the length of the (max_particles_per_thread)th particles defer_list
      end do

      deallocate(thread_particle_indices, thread_particle_data, defer_list_start_pos, partner_leaves)
      deallocate(defer_list_old, defer_list_new)
      deallocate(todo_list)
    end if

    my_threaddata%timers(THREAD_TIMER_TOTAL) = my_threaddata%timers(THREAD_TIMER_TOTAL) + MPI_WTIME()
    my_threaddata%timers(THREAD_TIMER_GET_NEW_PARTICLE) = t_get_new_particle
    my_threaddata%timers(THREAD_TIMER_WALK_SINGLE_PARTICLE) = t_walk_single_particle

    my_threaddata%finished = .true.
    num_finished = atomic_fetch_and_increment_int(threads_finished) + 1

    ! tell rank 0 that we are finished with our walk
    if (num_finished == num_walk_threads .and. walk_debug) then
      DEBUG_INFO(*, "PE", walk_tree%comm_env%rank, "has finished walking")
    end if

    walk_worker_thread = c_null_ptr
    call retval(pthreads_exitthread(), "walk_worker_thread:pthread_exit")

    contains
  
    subroutine swap_defer_lists()
      implicit none
      type(t_defer_list_entry), dimension(:), pointer :: tmp_list

      tmp_list       => defer_list_old
      defer_list_old => defer_list_new
      defer_list_new => tmp_list

      defer_list_new_tail = 1 ! position of first free entry in defer_list_new (i.e. it is considered as empty now)
    end subroutine


    logical function contains_particle(idx)
      implicit none
      integer, intent(in) :: idx

      contains_particle = ( thread_particle_indices(idx) .ge. 0 )
    end function contains_particle


    subroutine get_new_particle_and_setup_defer_list(idx)
      implicit none
      integer, intent(in) :: idx

      if (particles_available) then
        thread_particle_indices(idx) = get_first_unassigned_particle()

        if (contains_particle(idx)) then
          ! we make a copy of all particle data to avoid thread-concurrent access to particle_data array
          thread_particle_data(idx) = particle_data(thread_particle_indices(idx))
          ! for particles that we just inserted into our list, we start with only one defer_list_entry: the root node
          ptr_defer_list_old      => defer_list_root_only
          defer_list_entries_old  =  1
          partner_leaves(idx)     =  0 ! no interactions yet
        else
          particles_available     = .false.
        end if ! contains_particle(idx)
      end if ! particles_available
    end subroutine get_new_particle_and_setup_defer_list


    subroutine setup_defer_list(idx)
      implicit none
      integer, intent(in) :: idx

      ptr_defer_list_old     => defer_list_old(defer_list_start_pos(idx):defer_list_start_pos(idx+1)-1)
      defer_list_entries_old =  defer_list_start_pos(idx+1) - defer_list_start_pos(idx)
    end subroutine setup_defer_list

    subroutine do_sched_yield_if_necessary()
      implicit none
      ! after processing a number of particles: handle control to other (possibly comm) thread
      if (same_core_as_communicator) then
        if (particles_since_last_yield >= particles_per_yield) then
          call comm_sched_yield()
          particles_since_last_yield = 0
        else
          particles_since_last_yield = particles_since_last_yield + 1
        end if
      end if
    end subroutine
  end function walk_worker_thread


  function walk_single_particle(particle, defer_list_old, defer_list_entries_old, &
                                          defer_list_new, defer_list_entries_new, &
                                          todo_list, partner_leaves, my_threaddata)
    use module_walk_pthreads_commutils
    use module_tree_node
    use module_tree, only: tree_lookup_node_critical
    use module_tree_communicator, only: tree_node_fetch_children
    use module_interaction_specific
    use module_spacefilling, only : is_ancestor_of_particle
    use module_debug
    use module_mirror_boxes, only : spatial_interaction_cutoff
    use module_atomic_ops
    implicit none
    include 'mpif.h'

    type(t_particle), intent(inout) :: particle
    type(t_defer_list_entry), dimension(:), pointer, intent(in) :: defer_list_old
    integer, intent(in) :: defer_list_entries_old
    type(t_defer_list_entry), dimension(:), pointer, intent(out) :: defer_list_new
    integer, intent(out) :: defer_list_entries_new
    integer*8, intent(inout) :: todo_list(0:todo_list_length-1) 
    integer*8, intent(inout) :: partner_leaves
    type(t_threaddata), intent(inout) :: my_threaddata
    logical :: walk_single_particle !< function will return .true. if this particle has finished its walk

    integer :: todo_list_entries
    integer*8 :: walk_key, childlist(8)
    integer :: childnum, walk_level
    type(t_tree_node), pointer :: walk_node
    real*8 :: dist2
    real*8 :: delta(3), shifted_particle_position(3)
    logical :: same_particle, same_particle_or_parent_node, mac_ok, ignore_node
    integer*8 :: num_interactions, num_mac_evaluations, num_post_request
    real*8 :: t_post_request

    todo_list_entries      = 0
    num_interactions       = 0
    num_mac_evaluations    = 0
    t_post_request         = 0._8
    num_post_request       = 0
    walk_node              => null()
    shifted_particle_position = particle%x - vbox ! precompute shifted particle position to avoid subtracting vbox in every loop iteration below

    ! for each entry on the defer list, we check, whether children are already available and put them onto the todo_list
    ! another mac-check for each entry is not necessary here, since due to having requested the children, we already know,
    ! that the node has to be resolved
    ! if the defer_list is empty, the call reurns without doing anything
    call defer_list_parse_and_compact()
    call atomic_read_barrier()

    ! read all todo_list-entries and start further traversals there
    do while (todo_list_pop(walk_key))
      call tree_lookup_node_critical(walk_tree, walk_key, walk_node, 'WALK:walk_single_particle')
      walk_level = walk_node%level

      delta = shifted_particle_position - walk_node%interaction_data%coc  ! Separation vector
      dist2 = DOT_PRODUCT(delta, delta)

      if (tree_node_is_leaf(walk_node)) then
        mac_ok = .true.
      else
        mac_ok = mac(particle, walk_node%interaction_data, dist2, boxlength2(walk_level))
        num_mac_evaluations = num_mac_evaluations + 1
      end if

      ! we may not interact with the particle itself or its ancestors
      ! if we are in the central box
      ! interaction with ancestor nodes should be prevented by the MAC
      ! but this does not always work (i.e. if theta > 0.7 or if keys and/or coordinates have
      ! been modified due to 'duplicate keys'-error)
      same_particle_or_parent_node  = (in_central_box) .and. ( is_ancestor_of_particle(particle%key, walk_key, walk_level))
      ! set ignore flag if leaf node corresponds to particle itself
      same_particle = same_particle_or_parent_node .and. tree_node_is_leaf(walk_node)

      ! ignore interactions with the particle itself (this is the place for possible other exclusion options)
      ignore_node = same_particle

      if (.not. ignore_node) then
        !  always accept leaf-nodes since they cannot be refined any further
        !  further resolve ancestor nodes if we are in the central box
        mac_ok = tree_node_is_leaf(walk_node) .or. ( mac_ok .and. (.not. same_particle_or_parent_node))

        ! ========= Possible courses of action:
        if (mac_ok) then
          ! 1) leaf node or MAC test OK ===========
          !    --> interact with cell if it does not lie outside the cutoff box
          if (all(abs(delta) < spatial_interaction_cutoff)) then
            call calc_force_per_interaction(particle, walk_node%interaction_data, walk_key, delta, dist2, vbox, tree_node_is_leaf(walk_node))

            num_interactions = num_interactions + 1
          end if

          partner_leaves = partner_leaves + walk_node%leaves
        else
          ! 2) MAC fails for twig node ============
          if ( tree_node_children_available(walk_node) ) then
            ! 2a) children for twig are present --------
            ! --> resolve cell & put all children in front of todo_list
            call tree_node_get_childkeys(walk_node, childnum, childlist)
            if (.not. todo_list_push(childnum, childlist)) then
              ! the todo_list is full --> put parent back onto defer_list
              call defer_list_push(walk_key, walk_node)
            end if
          else
            ! 2b) children for twig are _absent_ --------
            ! --> put node on REQUEST list and put walk_key on bottom of todo_list
            t_post_request = t_post_request - MPI_WTIME()
            call tree_node_fetch_children(walk_tree, walk_node) ! fetch children from remote
            t_post_request = t_post_request + MPI_WTIME()
            num_post_request = num_post_request + 1
            ! if posting the request failed, this is not a problem, since we defer the particle anyway
            ! since it will not be available then, the request will simply be repeated
            call defer_list_push(walk_key, walk_node) ! Deferred list of nodes to search, pending request
                                                      ! for data from nonlocal PEs
            if (walk_debug) then
              DEBUG_INFO('("PE ", I6, " adding nonlocal key to defer_list, defer_list_entries=", I6)',  walk_tree%comm_env%rank, defer_list_entries_new)
            end if
          end if
        end if
      else !(ignore_node)
        partner_leaves = partner_leaves + walk_node%leaves
      end if !(.not. ignore_node)
    end do ! (while (todo_list_pop(walk_key)))

    ! if todo_list and defer_list are now empty, the walk has finished
    walk_single_particle = (todo_list_entries == 0) .and. (defer_list_entries_new == 0)

    my_threaddata%counters(THREAD_COUNTER_INTERACTIONS) = my_threaddata%counters(THREAD_COUNTER_INTERACTIONS) + num_interactions
    my_threaddata%counters(THREAD_COUNTER_MAC_EVALUATIONS) = my_threaddata%counters(THREAD_COUNTER_MAC_EVALUATIONS) + num_mac_evaluations
    my_threaddata%timers(THREAD_TIMER_POST_REQUEST) = my_threaddata%timers(THREAD_TIMER_POST_REQUEST) + t_post_request
    my_threaddata%counters(THREAD_COUNTER_POST_REQUEST) = my_threaddata%counters(THREAD_COUNTER_POST_REQUEST) + num_post_request

    contains

    ! helper routines for todo_list manipulation
    function todo_list_pop(key)
      implicit none
      logical :: todo_list_pop
      integer*8, intent(out) :: key

      todo_list_pop = (todo_list_entries > 0)

      if (todo_list_pop) then
        todo_list_entries = todo_list_entries - 1
        key               = todo_list(todo_list_entries)
      end if
    end function


    function todo_list_push(numkeys, keys)
      use module_debug
      implicit none
      integer, intent(in) :: numkeys
      integer*8, dimension(numkeys), intent(in) :: keys
      logical :: todo_list_push
      integer :: i

      if (todo_list_entries + numkeys > todo_list_length) then
        DEBUG_WARNING_ALL('("todo_list is full for particle with label ", I20, " todo_list_length =", I6, " is too small (you should increase interaction_list_length_factor). Putting particles back onto defer_list. Programme will continue without errors.")', particle%label, todo_list_length)
        todo_list_push = .false.
      else
        do i=1,numkeys
          todo_list(todo_list_entries) = keys(i)
          todo_list_entries            = todo_list_entries + 1
        end do

        todo_list_push = .true.
      end if
    end function


    ! helper routines for defer_list manipulation
    subroutine defer_list_push(key_, node_)
      use module_debug
      implicit none
      integer*8, intent(in) :: key_
      type(t_tree_node), pointer, intent(in) :: node_

      defer_list_entries_new                 = defer_list_entries_new + 1
      defer_list_new(defer_list_entries_new) = t_defer_list_entry( node_, key_ )
    end subroutine


    subroutine defer_list_parse_and_compact()
      use module_htable
      implicit none
      integer :: iold, cnum
      integer*8 :: clist(8)

      defer_list_entries_new = 0
      do iold = 1,defer_list_entries_old
        if ( tree_node_children_available(defer_list_old(iold)%node) ) then
          ! children for deferred node have arrived --> put children onto todo_list
          call tree_node_get_childkeys(defer_list_old(iold)%node, cnum, clist)
          if (.not. todo_list_push(cnum, clist)) then
            ! the todo_list is full --> put parent back onto defer_list
            defer_list_entries_new                 = defer_list_entries_new + 1
            defer_list_new(defer_list_entries_new) = defer_list_old(iold)
          end if
        else
          ! children for deferred node are still unavailable - put onto defer_list_new (do not use defer_list_push for performance reasons)
          defer_list_entries_new                 = defer_list_entries_new + 1
          defer_list_new(defer_list_entries_new) = defer_list_old(iold)
        end if
      end do
    end subroutine
  end function walk_single_particle
end module module_walk

