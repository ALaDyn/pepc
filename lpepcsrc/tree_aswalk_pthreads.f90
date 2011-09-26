! ===========================================
!
!           TREE_WALK
!
!  Perform tree walk for all local particles
!  in a hybrid parallelization scheme using
!  linux pthreads
!
!  Algorithm follows the implementation of
!  Warren & Salmon`s 'latency-hiding' concept,
!  retaining list-based tree-walk from vectorised
!  code by Pfalzner & Gibbon.
!
!
!  Structure:
!    * primary thread (aka 'do_communication_loop')
!      performs any MPI-communication with other
!      compute nodes and inserts all received multipole
!      information into the local tree structure
!    * secondary threads ('walk_worker_thread')
!      grab a number of particles and perform their
!      individual walks. as soon as one particle has
!      finished walking, the worker_thread takes an
!      additional particle as long as there are still
!      unprocessed particles available
!    * communication requests are inhibited from the
!      schedule- and work-threads via
!           - post_request (req. child data for some parent key)
!           - notify_walk_finished (notify MPI-rank 0, that
!                this PE has finished)
!    * while the worker threads simply exit
!      after their work is finished, the communicator stays
!      active until the walks on all MPI ranks have been
!      completed. the information about this status
!      is distributed by MPI rank 0 to all other PEs
!      using a special MPI message
!    * this results in having effectively only one pass of tree_walk.
!      the communicator can stay active until anything is
!      finished on all PEs and the workload distribution is much easier.
!    * additionally, this allows to initialize all walk data only
!      once per run instead of once per chunk. we just use
!      nintmax and max_particles_per_thread as limits for maximum
!      number of interactions and maximum chunk length
!
!
!  Structure of individual walk_work_threads:
!  ------------------------------------------
!      do while (particles_active .or. particles_available)
!
!        particles_active = .false.
!
!        do i=1,max_particles_per_thread
!
!          if ( (my_particles(i) == -1) .and. (particles_available)) then         ! i.e. the place for a particle is unassigned
!            my_particles(i) = get_first_unassigned_particle()
!          end if
!
!          call walk_single_particle(my_particles(i))
!
!        end do
!      end do
!
!
!  Structure of walk_single_particle(particle):
!  ------------------------------------------
!
!     if (.not.finished(particle)) then
!       num_unfinished = num_unfinished + 1
!
!       check defer_list entries:
!         if (requested children available)
!            put them onto todo_list
!
!       do while (can take entry form todo_list)
!           if (MAC OK)
!                immedeately interact with node
!           else
!                if (node locally available)
!                    resolve node
!                    put all children to front of todo_list
!                else
!                    post_request(parentkey, owner)
!                    put node on defer_list
!                end if
!           end if
!       end do
!     end if
!
!
!
!  Structure of communicator:
!  --------------------------
!    do while (not all PEs finished)
!
!        if (requests have been posted)
!            send all requests
!
!        if (all my walks finished)
!            send_walk_finished(to rank 0)
!
!        while (received MPI-message)
!          case (message tag) of
!              TAG_REQUEST_KEY:    send child_data(for parent key we just received to sender)
!              TAG_REQUESTED_DATA: unpack received child data, insert into local tree structure and mark as locally available
!              TAG_FINISHED_PE:    mark sender as having its walk finished (this msg is only received by rank 0)
!              TAG_FINISHED_ALL:   exit communicator loop
!        end while
!
!    end while
!
!
! ===========================================
module tree_walk_pthreads_commutils
  use tree_walk_communicator
  use pthreads_stuff
  implicit none
  private

    !> debug flags - cannot be modified at runtime due to performance reasons
    logical, parameter :: rwlock_debug    = .false.
    logical, parameter, public  :: walk_debug     = .false.

    ! variables for adjusting the thread`s workload
    real, public :: work_on_communicator_particle_number_factor = 0.1 !< factor for reducing max_particles_per_thread for thread which share their processor with the communicator
    integer, public :: particles_per_yield = 500 !< number of particles to process in a work_thread before it shall call sched_yield to hand the processor over to some other thread
    integer, parameter :: MAX_MESSAGES_PER_ITERATION = 20
    integer, parameter :: MIN_MESSAGES_PER_ITERATION = 5

    integer, parameter :: ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
    integer, parameter :: REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue

    type(t_request_queue_entry), volatile, private, target :: req_queue(REQUEST_QUEUE_LENGTH)
    integer, private, volatile :: req_queue_top, req_queue_bottom ! we will insert data at bottom and take from top
    integer, private :: request_balance !< total (#requests - #answers), should be zero after complete traversal

    ! rwlocks for regulating concurrent access
    integer, private, parameter :: NUM_RWLOCKS = 2
    integer, public, parameter :: RWLOCK_REQUEST_QUEUE      = 1
    integer, public, parameter :: RWLOCK_NEXT_FREE_PARTICLE = 2

    ! internal initialization status
    logical, private :: initialized = .false.

    ! local walktime (i.e. from comm_loop start until send_walk_finished() )
    real*8, public, pointer :: twalk_loc

    !> type for input and return values of walk_threads
    type, public :: t_threaddata
      integer :: id                         !< just a running number to distinguish the threads, currently unused
      integer :: num_processed_particles  !< thread output value: number of particles that it has processed
      real*8  :: num_interactions !< thread output value: number of interactions that were performed
      real*8  :: num_mac_evaluations !< thread output value: number of mac evaluations that have been performed
      logical :: is_on_shared_core !< thread output value: is set to true if the thread detects that it shares its processor with the communicator thread
      integer :: coreid !< thread output value: id of thread`s processor
      real*8 :: runtime_seconds !< thread wallclock-runtime in seconds, measured with MPI_WTIME()
      logical :: finished !< will be set to .true. when the thread has finished
    end type t_threaddata

    type(t_threaddata), public, allocatable, target :: threaddata(:)


    public run_communication_loop
    public send_requests
    public post_request
    public rwlock_rdlock
    public rwlock_wrlock
    public rwlock_unlock
    public comm_sched_yield
    public retval
    public init_commutils
    public uninit_commutils

  contains
  
    subroutine init_commutils()
        use treevars
        implicit none

        if (.not. initialized) then

          call init_comm_data(REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH)

          ! initialize rwlock objects
          call retval(rwlocks_init(NUM_RWLOCKS), "rwlocks_init")

          req_queue_top    = 0
          req_queue_bottom = 0
          request_balance  = 0

          initialized = .true.

        endif

      end subroutine init_commutils



      subroutine uninit_commutils
        use tree_walk_communicator
        implicit none

        initialized = .false.

        call uninit_comm_data

        ! free the conditional variable objects
        call retval(pthreads_conds_uninit(), "pthreads_conds_uninit")
        ! free the rwlock objects
        call retval(rwlocks_uninit(), "rwlocks_uninit")

      end subroutine uninit_commutils



    subroutine check_comm_finished()
      use treevars
      implicit none
      ! if request_balance contains positive values, not all requested data has arrived
      ! this means, that there were algorithmically unnecessary requests, which would be a bug
      ! negative values denote more received datasets than requested, which implies
      ! a bug in bookkeeping on the sender`s side
        if (req_queue_top .ne. req_queue_bottom) then
           write(*,*) "PE", me, " has finished its walk, but the request list is not empty"
           write(*,*) "obviously, there is an error in the todo_list bookkeeping"
           write(*,*) "Trying to recover from that"
           flush(6)
        elseif (request_balance .ne. 0) then
           write(*,*) "PE", me, " finished walking but is are still waiting for requested data: request_balance =", request_balance
           write(*,*) "This should never happen - obviously, the request_queue or some todo_list is corrupt"
           write(*,*) "Trying to recover from this situation anyway: Waiting until everything arrived although we will not interact with it."
           flush(6)
        else
           walk_status = WALK_ALL_MSG_DONE
        end if
    end subroutine



      subroutine retval(iret, msg)
        use treevars
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'
        integer( kind= c_int) :: iret
        character(*) :: msg
        integer :: ierr

        if (iret .ne. 0) then
          write(*,*)       "PE:", me, "[", msg, "] iret == ", iret
          write(ipefile,*) "PE:", me, "[", msg, "] iret == ", iret
          flush(6)
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

      end subroutine retval



      ! this routine is thread-safe to prevent concurrent write access to the queue
      subroutine post_request(request_key, request_addr)
        use treevars
        use module_htable
        use tree_walk_communicator
        implicit none
        include 'mpif.h'
        integer*8, intent(in) :: request_key
        integer, intent(in) :: request_addr
        integer :: local_queue_bottom
        logical, save :: warned = .false.

        ! check wether the node has already been requested
        ! this if-construct has to be secured against synchronous invocation (together with the modification while receiving data)
        ! otherwise it will be possible that two walk threads can synchronously post a prticle to the request queue
        if (BTEST( htable(request_addr)%childcode, CHILDCODE_BIT_REQUEST_POSTED ) ) then
          return
        endif

        ! we first flag the particle as having been already requested to prevent other threads from doing it while
        ! we are inside this function
        !call rwlock_wrlock(RWLOCK_CHILDBYTE, "walk_single_particle")
        htable(request_addr)%childcode   =  IBSET( htable(request_addr)%childcode, CHILDCODE_BIT_REQUEST_POSTED ) ! Set requested flag
        !call rwlock_unlock(RWLOCK_CHILDBYTE, "walk_single_particle")

        call rwlock_wrlock(RWLOCK_REQUEST_QUEUE, "post_request")

        ! use a thread safe list to put the requests onto
        local_queue_bottom = mod(req_queue_bottom, REQUEST_QUEUE_LENGTH) + 1

        if (local_queue_bottom == req_queue_top) then
          if (.not. warned) then
            write(*,*) "Issue with request sending queue: REQUEST_QUEUE_LENGTH is too small: ", REQUEST_QUEUE_LENGTH, ". Will try again later."
            warned = .true.
          end if

          call rwlock_unlock(RWLOCK_REQUEST_QUEUE, "post_request")
          ! since posting the request failed due to a too short queue, we have to flag the particle as to be procssed again
          !call rwlock_wrlock(RWLOCK_CHILDBYTE, "walk_single_particle")
          htable(request_addr)%childcode   =  IBCLR( htable(request_addr)%childcode, CHILDCODE_BIT_REQUEST_POSTED )
          !call rwlock_unlock(RWLOCK_CHILDBYTE, "walk_single_particle")
          return
        end if

        req_queue(local_queue_bottom)   = t_request_queue_entry( request_key, request_addr, htable( request_addr )%owner )

        if (walk_comm_debug) then
          write(ipefile,'("PE", I6, " posting request. local_queue_bottom=", I5, ", request_key=", O22, ", request_owner=", I6, " request_addr=", I12)') &
                         me, local_queue_bottom, request_key, htable( request_addr )%owner, request_addr
        end if

        ! now, we can tell the communicator that there is new data available
        req_queue_bottom = local_queue_bottom

        call rwlock_unlock(RWLOCK_REQUEST_QUEUE, "post_request")

    end subroutine post_request



    subroutine send_requests()
      use tree_walk_communicator
      use treevars
      implicit none
      include 'mpif.h'

      integer :: local_queue_bottom ! buffer for avoiding interference with threads that post data to the queue
      real*8 :: tsend
      integer*8 :: req_queue_length

      ! send all requests from our thread-safe list

      local_queue_bottom = req_queue_bottom

      if (req_queue_top .ne. local_queue_bottom) then

        req_queue_length    = modulo(local_queue_bottom-req_queue_top, REQUEST_QUEUE_LENGTH)
        max_req_list_length = max(max_req_list_length,  req_queue_length)
        cum_req_list_length =     cum_req_list_length + req_queue_length
        comm_loop_iterations(2) = comm_loop_iterations(2) + 1

        tsend = MPI_WTIME()

        do while (req_queue_top .ne. local_queue_bottom)

          req_queue_top = mod(req_queue_top, REQUEST_QUEUE_LENGTH) + 1

          if (walk_comm_debug) then
            write(ipefile,'("PE", I6, " sending request.      req_queue_top=", I5, ", request_key=", O22, ", request_owner=", I6)') &
                         me, req_queue_top, req_queue(req_queue_top)%key, req_queue(req_queue_top)%owner
          end if

          if (send_request(req_queue(req_queue_top))) then
            request_balance = request_balance + 1
          endif

        end do

       timings_comm(TIMING_SENDREQS) = timings_comm(TIMING_SENDREQS) + ( MPI_WTIME() - tsend )

     end if

    end subroutine send_requests



      subroutine run_communication_loop(max_particles_per_thread)
        use pthreads_stuff
        use treevars
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'
        integer, intent(in) :: max_particles_per_thread
        integer, dimension(mintag:maxtag) :: nummessages
        integer :: messages_per_iteration !< tracks current number of received and transmitted messages per commloop iteration for adjusting particles_per_yield
        logical :: walk_finished(num_pe) ! will hold information on PE 0 about which processor
                                          ! is still working and which ones are finished
                                          ! to emulate a non-blocking barrier
        integer :: ierr
        nummessages = 0

        if (me==0) walk_finished = .false.

        twalk_loc = MPI_WTIME()

        ! check whether initialization has correctly been performed
        if (.not. initialized) then
           write(*,*) "Serious issue in PE", me, ": walk_communicator has not been initialized. Call init_comm_data() before run_communication_loop(..)"
           flush(6)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

        if (walk_comm_debug) write(ipefile,'("PE", I6, " run_communication_loop start. walk_status = ", I6)') me, walk_status

        timings_comm(TIMING_COMMLOOP) = MPI_WTIME()

        do while (walk_status < WALK_ALL_FINISHED)

          comm_loop_iterations(1) = comm_loop_iterations(1) + 1

          ! send our requested keys
          call send_requests()

          ! check whether we are still waiting for data or some other communication
          if (walk_status == WALK_IAM_FINISHED) call check_comm_finished()

          ! process any incoming answers
          call run_communication_loop_inner(walk_finished, nummessages)

          messages_per_iteration = messages_per_iteration + sum(nummessages)
          request_balance = request_balance - nummessages(TAG_REQUESTED_DATA)
          nummessages(TAG_REQUESTED_DATA) = 0

          ! adjust the sched_yield()-timeout for the thread that shares its processor with the communicator
          if (messages_per_iteration > MAX_MESSAGES_PER_ITERATION) then
            particles_per_yield = int(max(0.75 * particles_per_yield, 0.01*max_particles_per_thread))
            if (walk_debug) write(ipefile,'("messages_per_iteration = ", I6, " > ", I6, " --> Decreased particles_per_yield to", I10)') messages_per_iteration, MAX_MESSAGES_PER_ITERATION, particles_per_yield
          elseif ((particles_per_yield < max_particles_per_thread) .and. (messages_per_iteration < MIN_MESSAGES_PER_ITERATION)) then
            particles_per_yield = int(min(1.5 * particles_per_yield, 1.*max_particles_per_thread))
            if (walk_debug) write(ipefile,'("messages_per_iteration = ", I6, " < ", I6, " --> Increased particles_per_yield to", I10)') messages_per_iteration, MIN_MESSAGES_PER_ITERATION, particles_per_yield
          endif

          ! currently, there is no further communication request --> other threads may do something interesting
          call comm_sched_yield()

        end do ! while (walk_status .ne. WALK_ALL_FINISHED)

        if (walk_comm_debug) write(ipefile,'("PE", I6, " run_communication_loop end.   walk_status = ", I6)') me, walk_status

        timings_comm(TIMING_COMMLOOP) = MPI_WTIME() - timings_comm(TIMING_COMMLOOP)

    end subroutine run_communication_loop



      subroutine rwlock_rdlock(idx, reason)
        implicit none
        integer, intent(in) :: idx
        character(*), intent(in) :: reason

        call retval(rwlocks_rdlock(idx), "pthread_rwlock_rdlock:"//reason)

        if (rwlock_debug) then
          write(*,*) "pthread_rwlock_rdlock:"//reason, ", idx =", idx
          flush(6)
        end if

      end subroutine rwlock_rdlock



      subroutine rwlock_wrlock(idx, reason)
        implicit none
        integer, intent(in) :: idx
        character(*), intent(in) :: reason

        call retval(rwlocks_wrlock(idx), "pthread_rwlock_wrlock:"//reason)

        if (rwlock_debug) then
          write(*,*) "pthread_rwlock_wrlock:"//reason, ", idx =", idx
          flush(6)
        end if

      end subroutine rwlock_wrlock



      subroutine rwlock_unlock(idx, reason)
        implicit none
        integer, intent(in) :: idx
        character(*), intent(in) :: reason

        call retval(rwlocks_unlock(idx), "pthread_rwlock_unlock:"//reason)

        if (rwlock_debug) then
          write(*,*) "pthread_rwlock_unlock:"//reason, ", idx =", idx
          flush(6)
        end if

      end subroutine rwlock_unlock


      subroutine comm_sched_yield()
        use pthreads_stuff
        implicit none

        call retval(pthreads_sched_yield(), "pthreads_sched_yield()")
      end subroutine comm_sched_yield


end module tree_walk_pthreads_commutils






module tree_walk_pthreads
  use treetypes
  use treevars
  use pthreads_stuff
  implicit none

  private

    integer, public :: max_particles_per_thread = 2000 !< maximum number of particles that will in parallel be processed by one workthread

    integer, public :: num_walk_threads = 3 !< number of worker threads
    integer, private :: primary_processor_id = 0
    real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8
    real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8

    real*8, dimension(:), allocatable :: boxlength2
    integer :: mac
    real*8 :: theta2
    real*8 :: vbox(3)
    type(calc_force_params) :: cf_par
    logical :: in_central_box
    integer*8 :: num_interaction_leaves
    integer :: todo_list_length, defer_list_length


    integer :: np_local
    real*8, dimension(:), pointer :: work_per_particle => NULL()

    integer :: next_unassigned_particle !< index of next particle that has not been assigned to a work thread

    type t_defer_list_entry
      integer  :: addr
      integer*8 :: key
    end type t_defer_list_entry


    public tree_walk

  contains


    subroutine tree_walk(np_local_,theta_,cf_par_,itime,mac_,twalk,twalk_loc_,vbox_,work_per_particle_, tcomm)
      use, intrinsic :: iso_c_binding
      use treetypes
      use tree_utils
      use timings
      use tree_walk_pthreads_commutils
      use tree_walk_communicator
      implicit none
      include 'mpif.h'

      real, intent(in) :: theta_  ! MAC angle
      type(calc_force_params), intent(in) :: cf_par_
      integer, intent(in) :: np_local_
      integer, intent(in) :: itime
      integer, intent(in) :: mac_
      real*8, intent(in) :: vbox_(3) !< real space shift vector of box to be processed
      real*8, dimension(nppm), intent(inout), target :: work_per_particle_
      real*8, target, intent(inout) :: twalk, twalk_loc_
      real*8, target, intent(out), dimension(3) :: tcomm

      if (me.eq.0 .and. walk_summary) write(*,'(2(a,i6))') 'LPEPC | TREE WALK (HYBRID) for timestep ',itime

      np_local = np_local_
      ! allow global access to workload stuff
      work_per_particle => work_per_particle_
      ! box shift vector
      vbox = vbox_
      ! force calculation parameters
      cf_par = cf_par_
      ! mac selection parameter
      mac = mac_
      ! Clumping parameter**2 for MAC
      theta2 = theta_**2
      ! length of todo- and defer-list per particle (estimations)
      todo_list_length  = nintmax
      defer_list_length = max(todo_list_length / 8, 50)
      ! pure local walk time (i.e. from start of communicator till sned_walk_finished)
      twalk_loc => twalk_loc_

      twalk  = MPI_WTIME()

      call init_walk_data()

      call init_commutils()

      call walk_hybrid()

      call uninit_commutils()

      call uninit_walk_data()

      twalk = MPI_WTIME() - twalk
      tcomm = timings_comm

    end subroutine tree_walk



    subroutine walk_hybrid()
      use tree_walk_pthreads_commutils
      use, intrinsic :: iso_c_binding
      implicit none
      include 'mpif.h'

      integer :: ith, displ, ierr

      allocate(threaddata(num_walk_threads))

      ! start the worker threads...
      do ith = 1,num_walk_threads
        threaddata(ith)%id = ith
        threaddata(ith)%finished = .false.
        call retval(pthreads_createthread(ith, c_funloc(walk_worker_thread), c_loc(threaddata(ith))), "walk_schedule_thread_inner:pthread_create")
      end do

      call run_communication_loop(num_walk_threads)

      ! ... and wait for work thread completion
      thread_workload = 0.

      do ith = 1,num_walk_threads
        call retval( pthreads_jointhread( ith ), "walk_schedule_thread_inner:pthread_join" )

        if (threaddata(ith)%is_on_shared_core) then
          displ = 2
        else
          displ = 0
          thread_workload(0) = thread_workload(0) + 1
        endif

        thread_workload( displ+1) =     thread_workload( displ+1) + 1._8*threaddata(ith)%num_processed_particles
        thread_workload( displ+2) = max(thread_workload( displ+2),  1._8*threaddata(ith)%num_processed_particles)
        thread_workload(-displ-1) =     thread_workload(-displ-1) + 1._8*threaddata(ith)%runtime_seconds
        thread_workload(-displ-2) = max(thread_workload(-displ-2),  1._8*threaddata(ith)%runtime_seconds)


        if (walk_summary) then
          write(ipefile,*) "Hybrid walk finished for thread", ith, ". Returned data = ", threaddata(ith)
        end if
      end do

      ! compute relative deviation
      thread_workload( 2) = abs(thread_workload( 2) - thread_workload( 1)) / thread_workload( 1)
      thread_workload( 4) = abs(thread_workload( 4) - thread_workload( 1)) / thread_workload( 1)
      thread_workload( 1) =     thread_workload( 1)                        / thread_workload( 0)
      thread_workload( 3) =     thread_workload( 3)                        / (num_walk_threads - thread_workload(0))
      thread_workload(-2) = abs(thread_workload(-2) - thread_workload(-1)) / thread_workload(-1)
      thread_workload(-4) = abs(thread_workload(-4) - thread_workload(-1)) / thread_workload(-1)
      thread_workload(-1) =     thread_workload(-1)                        / thread_workload( 0)
      thread_workload(-3) =     thread_workload(-3)                        / (num_walk_threads - thread_workload(0))
      ! store workload data
      interactions_local    = sum(threaddata(:)%num_interactions)
      mac_evaluations_local = sum(threaddata(:)%num_mac_evaluations)

      ! check wether all particles really have been processed
      if (next_unassigned_particle .ne. np_local + 1) then
        write(*,*) "Serious issue on PE", me, ": all walk threads have terminated, but obviously not all particles are finished with walking: next_unassigned_particle =", &
                            next_unassigned_particle, " np_local =", np_local
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if

      deallocate(threaddata)

    end subroutine walk_hybrid



    subroutine init_walk_data()
      use, intrinsic :: iso_c_binding
      use tree_walk_pthreads_commutils
      implicit none
      integer :: i

      ! we have to have at least one walk thread
      num_walk_threads = max(num_walk_threads, 1)
      ! evenly balance particles to threads if there are less than the maximum
      max_particles_per_thread = max(min(np_local/num_walk_threads, max_particles_per_thread),1)
      ! allocate storage for thread handles, the 0th entry is the walk scheduler thread, the other ones are the walk worker threads
      call retval(pthreads_init(num_walk_threads + 1), "init_walk_data:pthreads_init")
      ! we will only want to reject the root node and the particle itself if we are in the central box
      in_central_box = (dot_product(vbox,vbox) == 0)
      ! every particle has directly or indirectly interact with each other, and outside the central box even with itself
      if (in_central_box) then
        num_interaction_leaves = npart - 1
      else
        num_interaction_leaves = npart
      endif
      ! Preprocessed box sizes for each level
      allocate(boxlength2(0:nlev))
      boxlength2(0)=boxsize**2
      do i=1,nlev
         boxlength2(i) =  boxlength2(i-1)/4.
      end do

      next_unassigned_particle = 1

      ! store ID of primary (comm-thread) processor
      primary_processor_id = get_my_core()

    end subroutine init_walk_data



    subroutine uninit_walk_data()
      use tree_walk_pthreads_commutils
      implicit none
      deallocate(boxlength2)
      call retval(pthreads_uninit(), "uninit_walk_data:pthreads_uninit")
    end subroutine uninit_walk_data



    function get_first_unassigned_particle(success)
      use tree_walk_pthreads_commutils
      implicit none
      integer :: get_first_unassigned_particle
      logical, intent(out) :: success

      call rwlock_wrlock(RWLOCK_NEXT_FREE_PARTICLE,"get_first_unassigned_particle")
      if (next_unassigned_particle < np_local + 1) then
        get_first_unassigned_particle = next_unassigned_particle
        next_unassigned_particle      = next_unassigned_particle + 1
        success = .true.
      else
        get_first_unassigned_particle = -1
        success = .false.
      end if
      call rwlock_unlock(RWLOCK_NEXT_FREE_PARTICLE,"get_first_unassigned_particle")

    end function get_first_unassigned_particle



    function walk_worker_thread(arg) bind(c)
      use, intrinsic :: iso_c_binding
      use tree_walk_pthreads_commutils
      use tree_walk_communicator
      use pthreads_stuff
      implicit none
      include 'mpif.h'
      type(c_ptr) :: walk_worker_thread
      type(c_ptr), value :: arg

      integer, dimension(:), allocatable :: my_particles
      integer*8, dimension(:), allocatable :: partner_leaves ! list for storing number of interaction partner leaves
      type(t_defer_list_entry), dimension(:,:), allocatable :: defer_list
      integer, dimension(:), allocatable :: defer_list_entries
      integer :: i
      logical :: particles_available
      logical :: particles_active
      logical :: process_particle
      type(t_threaddata), pointer :: my_threaddata
      integer :: particles_since_last_yield
      logical :: same_core_as_communicator
      integer :: my_max_particles_per_thread
      integer :: my_processor_id

      my_processor_id = get_my_core()
      same_core_as_communicator = (my_processor_id == primary_processor_id)

      if (same_core_as_communicator) then
            my_max_particles_per_thread = int(work_on_communicator_particle_number_factor * max_particles_per_thread)
      else
            my_max_particles_per_thread = max_particles_per_thread
      endif

      call c_f_pointer(arg, my_threaddata)
      my_threaddata = t_threaddata(my_threaddata%id, 0, 0._8, 0._8, same_core_as_communicator, my_processor_id, MPI_WTIME(), .false.)

      if (my_max_particles_per_thread > 0) then

          allocate(my_particles(my_max_particles_per_thread),                        &
                                defer_list_entries(my_max_particles_per_thread),      &
                                partner_leaves(my_max_particles_per_thread))
          allocate(defer_list(0:defer_list_length-1,my_max_particles_per_thread))

          ! every particle will start at the root node (one entry per todo_list, no particle is finished)
          my_particles(:)     = -1 ! no particles assigned to this thread
          defer_list_entries  =  1 ! one entry in defer_list:
          defer_list(0,:)     =  t_defer_list_entry(1, 1_8) !     start at root node (addr, and key)
          partner_leaves      =  0 ! no interactions yet
          particles_available = .true.
          particles_active    = .true.
          particles_since_last_yield = 0


          do while (particles_active .or. particles_available)

            particles_active = .false.

            do i=1,my_max_particles_per_thread

              process_particle = (my_particles(i) .ne. -1)

              if ((.not. process_particle) .and. (particles_available)) then ! i.e. the place for a particle is unassigned
                my_particles(i) = get_first_unassigned_particle(process_particle)
    !            write(ipefile,*) "PE", me, getfullid(), "my_particles(",i, ") :=", my_particles(i)
              end if

              if (process_particle) then

                ! after processing a number of particles: handle control to other (possibly comm) thread
                if (same_core_as_communicator) then
                  if (particles_since_last_yield >= particles_per_yield) then
                    call comm_sched_yield()
                    particles_since_last_yield = 0
                  else
                    particles_since_last_yield = particles_since_last_yield + 1
                  endif
                endif

                if (walk_single_particle(i, my_particles(i), defer_list, defer_list_entries(i), &
                                       my_max_particles_per_thread, partner_leaves(i), my_threaddata)) then
                  ! this particle`s walk has finished
    !              write(ipefile,*) "PE", me, getfullid(), " walk for particle", i, " nodeidx =", my_particles(i), " has finished"

                  ! walk for particle i has finished
                  ! check whether it really interacted with all other particles (in central box: excluding itself)
                  if (partner_leaves(i) .ne. num_interaction_leaves) then
                    write(*,'("Algorithmic problem on PE", I7, ": Particle ", I10, " label ", I16)') me, my_particles(i), particles(my_particles(i))%label
                    write(*,'("should have been interacting (directly or indirectly) with", I16," leaves (particles), but did with", I16)') num_interaction_leaves, partner_leaves(i)
                    write(*,*) "Its force and potential will be wrong due to some algorithmic error during tree traversal. Continuing anyway"
                  endif

                  !remove entries from defer_list
                  my_particles(i)            = -1
                  defer_list_entries(i)      =  1
                  defer_list(0,i)            =  t_defer_list_entry(1, 1_8)
                  partner_leaves(i)          =  0
                  my_threaddata%num_processed_particles = my_threaddata%num_processed_particles + 1
                else
                  particles_active = .true.
                end if
              else
                particles_available = .false.
              end if

            end do

          end do

          deallocate(my_particles, defer_list_entries, partner_leaves)
          deallocate(defer_list)

      endif

      my_threaddata%finished = .true.

      ! tell rank 0 that we are finished with our walk
      if (all(threaddata(:)%finished)) then
        call notify_walk_finished()

        twalk_loc = MPI_WTIME() - twalk_loc

        if (walk_debug) then
          write(ipefile,*) "PE", me, "has finished walking"
        end if

      endif

      walk_worker_thread = c_null_ptr

      my_threaddata%runtime_seconds = MPI_WTIME() - my_threaddata%runtime_seconds

      call retval(pthreads_exitthread(), "walk_worker_thread:pthread_exit")

    end function walk_worker_thread




   function walk_single_particle(myidx, nodeidx, defer_list, defer_list_entries, listlengths, partner_leaves, my_threaddata)
      use tree_walk_pthreads_commutils
      use module_htable
      use module_calc_force
      use module_spacefilling, only : level_from_key
      implicit none
      integer, intent(in) :: nodeidx, myidx, listlengths
      integer*8 :: todo_list(0:todo_list_length-1)
      integer*8, intent(inout) :: partner_leaves
      type(t_defer_list_entry), intent(inout) :: defer_list(0:defer_list_length-1,1:listlengths)
      integer, intent(inout) :: defer_list_entries
      integer :: todo_list_entries
      type(t_threaddata), intent(inout) :: my_threaddata
      logical :: walk_single_particle !< function will return .true. if this particle has finished its walk

      integer*8 :: walk_key, childlist(8)
      integer :: walk_addr, walk_node, childnum

      real*8 :: dist2
      real*8 :: delta(3)
      logical :: same_particle, mac_ok
      integer :: ierr

      todo_list_entries = 0

      ! for each entry on the defer list, we check, whether children are already available and put them onto the todo_list
      ! another mac-check for each entry is not necessary here, since due to having requested the children, we already know,
      ! that the node has to be resolved
      ! if the defer_list is empty, the call reurns without doing anything
      call defer_list_parse_and_compact()

      ! read all todo_list-entries and start further traversals there
      do while (todo_list_pop(walk_key))

          walk_addr = key2addr( walk_key, 'WALK:walk_single_particle' )  ! get htable address
          walk_node = htable( walk_addr )%node            ! Walk node index - points to multipole moments

          delta     = particles(nodeidx)%x - ([tree_nodes(walk_node)%xcoc, &
                                               tree_nodes(walk_node)%ycoc, &
                                               tree_nodes(walk_node)%zcoc] + vbox)  ! Separations

          dist2 = DOT_PRODUCT(delta, delta)

          if (mac == 0) then
              ! Barnes-Hut-MAC
              mac_ok = (theta2 * dist2 > boxlength2(tree_nodes(walk_node)%level))
          else
              ! TODO:  BH MAC is also implemented in mac_choose routine,
              !        which needs tuning and inlining to avoid excessive parameter-passing overhead
              !        Other MACS need additional preprocessed info on tree nodes (box lengths etc)
              !        before they can be used for performance tests.
              !        mac=1,2 only useful for accuracy tests at the moment (interaction list comparisons)
              !             call mac_choose(pshort(p),ex_nps(p),ey_nps(p),ez_nps(p),&
              !                  walk_node,walk_key(idx),abs_charge(walk_node),boxlength2(tree_nodes(walk_node)%level), &
              !                  theta2,mac,mac_ok, vbox)
              mac_ok = .false.
          endif

          !  always accept leaf-nodes since they cannot be refined any further
          !  reject root node if we are in the central box
          mac_ok = (walk_node > 0) .or. ( mac_ok .and. ((.not. in_central_box) .or. walk_key.ne.1 ) )

          work_per_particle(nodeidx)        = work_per_particle(nodeidx)        + WORKLOAD_PENALTY_MAC
          my_threaddata%num_mac_evaluations = my_threaddata%num_mac_evaluations + 1._8

          ! set ignore flag if leaf node corresponds to particle itself (number in pshort)
          ! NB: this uses local leaf #, not global particle label
          ! only ignore particle if we are processing the central box
          ! when walking through neighbour boxes, it has to be included
          same_particle = (in_central_box) .and. ( nodeidx == walk_node )

          if (.not. same_particle) then

          ! ========= Possible courses of action:

              if (mac_ok) then
                  ! 1) leaf node or MAC test OK ===========
                  !    --> interact with cell
                  call calc_force_per_interaction(nodeidx, walk_node, delta, dist2, vbox, cf_par)
                  work_per_particle(nodeidx)     = work_per_particle(nodeidx)     + WORKLOAD_PENALTY_INTERACTION
                  partner_leaves                 = partner_leaves                 + htable(walk_addr)%leaves
                  my_threaddata%num_interactions = my_threaddata%num_interactions + 1._8

              else
                  ! 2) MAC fails for twig node ============
                  if ( children_available(walk_addr) ) then
                      ! 2a) children for twig are present --------
                      ! --> resolve cell & put all children in front of todo_list
                      call get_childkeys(walk_addr, childnum, childlist)
                      call todo_list_push(childnum, childlist)
                  else
                      ! 2b) children for twig are _absent_ --------
                      ! --> put node on REQUEST list and put walk_key on bottom of todo_list
                      call post_request(walk_key, walk_addr)        ! tell the communicator about our needs
                      ! if posting the request failed, this is not a problem, since we defer the particle anyway
                      ! since it will not be available then, the request will simply be repeated
                      call defer_list_push(walk_key, walk_addr) ! Deferred list of nodes to search, pending request
                                                                     ! for data from nonlocal PEs
                      if (walk_debug) then
                          write(ipefile,'("PE ", I6, " adding nonlocal key to tail for particle ", I20, " defer_list_entries=", I6)') me, nodeidx, defer_list_entries
                      end if
                  end if
              endif
          endif !(.not. same_particle)
      end do ! (while (todo_list_pop_front(walk_key)))

      ! if todo_list and defer_list are now empty, the walk has finished
      walk_single_particle = (todo_list_entries == 0) .and. (defer_list_entries == 0)

      if (walk_single_particle .and. walk_debug) then
          write(ipefile,'("PE", I6, " particle ", I12, " has an empty todo_list: particle obviously finished walking around :-)")') me, nodeidx
      end if


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
       endif
     end function

     subroutine todo_list_push(numkeys, keys)
       implicit none
       integer, intent(in) :: numkeys
       integer*8, dimension(numkeys), intent(in) :: keys
       integer :: i

       if (todo_list_entries + numkeys >= todo_list_length) call todo_list_full()

       do i=1,numkeys
         todo_list(todo_list_entries) = keys(i)
         todo_list_entries            = todo_list_entries + 1
       end do
     end subroutine

     subroutine todo_list_full()
       implicit none
       include 'mpif.h'
        write(*,'("Rather serious issue on PE ", I6, ": todo_list is full for particle ", I20, " todo_list_length =", I6, " is too small (you should increase nintmax)")') me, nodeidx, todo_list_length
        write(*,*) "We could skip the current particle until some child_data has arrive, but this can result in a deadlock... --> aborting."
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end subroutine

     ! helper routines for defer_list manipulation
     subroutine defer_list_push(key_, addr_)
       implicit none
       integer, intent(in) :: addr_
       integer*8, intent(in) :: key_

       if (defer_list_entries == defer_list_length) call defer_list_full()

       defer_list(defer_list_entries, myidx) = t_defer_list_entry(addr_, key_)
       defer_list_entries                    = defer_list_entries + 1
     end subroutine

     subroutine defer_list_parse_and_compact()
       use module_htable
       implicit none
       integer :: iold, inew, cnum
       integer*8 :: clist(8)

       inew = 0
       do iold = 0,defer_list_entries-1
         if ( children_available(defer_list(iold, myidx)%addr) ) then
           ! children for deferred node have arrived --> put children onto todo_list
           call get_childkeys(defer_list(iold, myidx)%addr, cnum, clist)
           call todo_list_push(cnum, clist)
         else
           ! children for deferred node are still unavailable - put back onto defer_list
           defer_list(inew, myidx) = defer_list(iold, myidx)
           inew                    = inew + 1
         end if
       end do

       defer_list_entries = inew
     end subroutine

     subroutine defer_list_full()
       implicit none
       include 'mpif.h'
        write(*,'("Rather serious issue on PE ", I6, ": defer_list is full for particle ", I20, " defer_list_length =", I6, " is too small (you should increase nintmax)")') me, nodeidx, defer_list_length
        write(*,*) "We could skip the current particle until some child_data has arrive, but this can result in a deadlock... --> aborting."
        flush(6)
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end subroutine

    end function walk_single_particle

end module tree_walk_pthreads

