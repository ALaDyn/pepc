! ===========================================
!
!           TREE_WALK
!
!   $Revision: 1461 $
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
!      additional one as long as there are still
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
!       if (MAC OK)
!            immedeately interact with node
!       else
!            if (node locally available)
!                resolve node
!                put first child to front of todo_list
!            else
!                post_request(parentkey, owner)
!                put node on end of todo_list
!            end if
!       end if
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

module tree_walk_communicator
  use treevars
  use pthreads_stuff
  implicit none

  private

    integer, public :: walk_status


    logical, parameter :: walk_comm_debug = .false.
    logical, parameter :: rwlock_debug = .false.

    ! tags to be used in communication
    integer, parameter :: TAG_REQUEST_KEY    = 1257 !< message tag for walk communication: message for requesting child data for a certain key
    integer, parameter :: TAG_REQUESTED_DATA = 1258 !< message tag for walk communication: message that contains requested child data
    integer, parameter :: TAG_FINISHED_PE    = 1259 !< message tag for walk communication: message to rank 0, that the sender has finished with walking, i.e. its walk status == WALK_IAM_FINISHED
    integer, parameter :: TAG_FINISHED_ALL   = 1260 !< message tag for walk communication: message from rank 0, that all PEs are finished with walking, i.e. we can set walk_status := WALK_ALL_FINISHED

    ! internal status
    integer, public, parameter :: WALK_STILL_RUNNING = 0 !< value of walk_status: there are still particles left to be processed
    integer, public, parameter :: WALK_IAM_FINISHED  = 1 !< value of walk_status: all particles on this PE have finished their walk
    integer, public, parameter :: WALK_ALL_MSG_DONE  = 2 !< value of walk_status: there are no more pending requests or answers in the message queue
    integer, public, parameter :: WALK_I_NOTIFIED_0  = 3 !< value of walk_status: rank 0 has been informed that this PE is finished
    integer, public, parameter :: WALK_ALL_FINISHED  = 4 !< value of walk_status: rank 0 told this pe, that every PE is finished with walking, i.e. no more communication is necessary

    ! bits in childcode to be set when children are requested, the request has been sent, and they have arrived
    integer, public, parameter :: CHILDCODE_BIT_REQUEST_POSTED     =  8
    integer, public, parameter :: CHILDCODE_BIT_REQUEST_SENT       = 10
    integer, public, parameter :: CHILDCODE_BIT_CHILDREN_AVAILABLE =  9

    ! variables for adjusting the thread`s workload
    integer, public :: max_particles_per_thread = 2000 !< maximum number of particles that will in parallel be processed by one workthread
    real, public :: work_on_communicator_particle_number_factor = 0.1 !< factor for reducing max_particles_per_thread for thread which share their processor with the communicator
    integer, public :: particles_per_yield = 500 !< number of particles to process in a work_thread before it shall call sched_yield to hand the processor over to some other thread
    integer, private :: messages_per_iteration !< tracks current number of received and transmitted messages per commloop iteration for adjusting particles_per_yield
    integer, parameter :: MAX_MESSAGES_PER_ITERATION = 20
    integer, parameter :: MIN_MESSAGES_PER_ITERATION = 5
    logical, public :: comm_on_shared_processor = .false.

    ! internal communication variables - not to be touched from outside the module
    integer, parameter :: ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
    integer, parameter :: REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue
    integer*8, private, volatile :: req_queue_keys(REQUEST_QUEUE_LENGTH)
    integer, private, volatile :: req_queue_addrs(REQUEST_QUEUE_LENGTH)
    integer, private, volatile :: req_queue_owners(REQUEST_QUEUE_LENGTH)
    integer, private, volatile :: req_queue_top, req_queue_bottom ! we will insert data at bottom and take from top
    integer*1, private, allocatable, target :: bsend_buffer(:) !< buffer for bsend-alls
    integer, private :: comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)
    integer, private, allocatable :: request_balance(:) !< (#requests - #answers) per PE
    real*8, public :: timings_comm(3) !< array for storing internal timing information
    integer, parameter :: cond_comm_timeout = 1000 !< timeout for pthread_cond_wait() in commloop in microseconds TODO: tune this parameter automatically during runtime

    ! rwlocks for regulating concurrent access
    integer, private, parameter :: NUM_RWLOCKS = 4
    integer, public, parameter :: RWLOCK_REQUEST_QUEUE      = 1
    integer, public, parameter :: RWLOCK_NEXT_FREE_PARTICLE = 2
    integer, public, parameter :: RWLOCK_FINISHED_THREADS   = 3
    integer, public, parameter :: RWLOCK_CHILDBYTE          = 4

    ! pthread conditional variables
    integer, private, parameter :: NUM_CONDS = 1
    integer, public, parameter :: PTHREAD_COND_COMMBLOCK = 1

    ! IDs for internal timing measurement
    integer, public, parameter :: TIMING_COMMLOOP = 1
    integer, public, parameter :: TIMING_RECEIVE  = 2
    integer, public, parameter :: TIMING_SENDREQS = 3


    ! internal initialization status
    logical, private :: initialized = .false.

    public init_comm_data
    public run_communication_loop
    public run_communication_loop_inner
    public uninit_comm_data
    public notify_walk_finished
    public post_request
    public retval
    public rwlock_rdlock
    public rwlock_wrlock
    public rwlock_unlock
    public comm_sched_yield

  contains

      subroutine retval(iret, msg)
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

      ! initializes bsend buffer and rwlock objects
      ! returns size of bsend buffer in buffsize in bytes
      subroutine init_comm_data()
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'
        integer :: msg_size_request, msg_size_data
        integer :: buffsize !< size of bsend buffer in bytes
        integer :: ierr

        if (.not. initialized) then

          ! compute upper bounds for request and data message size
          call MPI_PACK_SIZE(1, MPI_INTEGER8, MPI_COMM_WORLD, msg_size_request, ierr)
          msg_size_request = msg_size_request + MPI_BSEND_OVERHEAD

          call MPI_PACK_SIZE(8, MPI_type_multipole, MPI_COMM_WORLD, msg_size_data, ierr)
          msg_size_data = msg_size_data + MPI_BSEND_OVERHEAD

          buffsize = (REQUEST_QUEUE_LENGTH * msg_size_request + ANSWER_BUFF_LENGTH * msg_size_data)
          ! reserve memory for buffered mpi communication
          allocate(bsend_buffer(buffsize))
          ! and tell mpi , where it can be found
          call MPI_BUFFER_ATTACH(bsend_buffer, buffsize, ierr)

          ! field for request balance per pe
          allocate(request_balance(num_pe))

          ! initialize rwlock objects
          call retval(rwlocks_init(NUM_RWLOCKS), "rwlocks_init")

          ! initialize conditional variable objects
          call retval(pthreads_conds_init(NUM_CONDS), "pthreads_conds_init")

        end if

        req_queue_top    = 0
        req_queue_bottom = 0

        walk_status = WALK_STILL_RUNNING
        request_balance = 0

        timings_comm = 0.

        initialized = .true.

      end subroutine init_comm_data





      subroutine run_communication_loop()
        use pthreads_stuff
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'
        integer( kind= c_int) :: iret
        logical :: walk_finished(num_pe) ! will hold information on PE 0 about which processor
                                          ! is still working and which ones are finished
                                          ! to emulate a non-blocking barrier
        integer :: ierr
        integer :: timeout = cond_comm_timeout

        if (me==0) walk_finished = .false.

        ! check whether initialization has correctly been performed
        if (.not. initialized) then
           write(*,*) "Serious issue in PE", me, ": walk_communicator has not been initialized. Call init_comm_data() before run_communication_loop(..)"
           flush(6)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

        if (walk_comm_debug) write(ipefile,'("PE", I6, " run_communication_loop start. walk_status = ", I6)') me, walk_status

        timings_comm(TIMING_COMMLOOP) = MPI_WTIME()

        do while (walk_status < WALK_ALL_FINISHED)

          messages_per_iteration = 0;

          comm_loop_iterations(1) = comm_loop_iterations(1) + 1
          call run_communication_loop_inner(walk_finished)

          ! adjust the sched_yield()-timeout for the thread that shares its processor with the communicator
          if (messages_per_iteration > MAX_MESSAGES_PER_ITERATION) then
            particles_per_yield = int(max(0.75 * particles_per_yield, 0.01*max_particles_per_thread))
            if (walk_debug) write(ipefile,'("messages_per_iteration = ", I6, " > ", I6, " --> Decreased particles_per_yield to", I10)') messages_per_iteration, MAX_MESSAGES_PER_ITERATION, particles_per_yield
          elseif ((particles_per_yield < max_particles_per_thread) .and. (messages_per_iteration < MIN_MESSAGES_PER_ITERATION)) then
            particles_per_yield = int(min(1.5 * particles_per_yield, 1.*max_particles_per_thread))
            if (walk_debug) write(ipefile,'("messages_per_iteration = ", I6, " < ", I6, " --> Increased particles_per_yield to", I10)') messages_per_iteration, MIN_MESSAGES_PER_ITERATION, particles_per_yield
          endif

          ! currently, there is no further communication request --> other threads may do something interesting
          if (comm_on_shared_processor) then
             timeout = cond_comm_timeout
          else
             timeout = cond_comm_timeout / 100
          endif

          iret = pthreads_conds_timedwait(PTHREAD_COND_COMMBLOCK, cond_comm_timeout)
          if (iret == -1) then
            ! Timeout occured, hence no pending send-requests,
            ! but we have to check whether there are incoming messages
          else
            call retval(iret, "pthreads_conds_timedwait(PTHREAD_COND_COMMBLOCK, cond_comm_timeout)")
          endif

        end do ! while (walk_status .ne. WALK_ALL_FINISHED)

        if (walk_comm_debug) write(ipefile,'("PE", I6, " run_communication_loop end.   walk_status = ", I6)') me, walk_status

        timings_comm(TIMING_COMMLOOP) = MPI_WTIME() - timings_comm(TIMING_COMMLOOP)

    end subroutine run_communication_loop



    subroutine run_communication_loop_inner(walk_finished)
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'
        logical :: msg_avail
        integer :: ierr
        integer :: stat(MPI_STATUS_SIZE)
        integer :: reqhandle
        integer*8 :: requested_key
        type (multipole) :: child_data(8) ! child data to be received - maximum up to eight children per particle
        integer :: num_children
        integer :: ipe_sender
        logical, intent(inout) :: walk_finished(num_pe)
        integer :: dummy, i
        real*8 :: tcomm

!         if (walk_comm_debug) then
!           write(ipefile,'("PE", I6, " run_communication_loop_inner. walk_status = ", I6)') &
!                          me, walk_status
!         end if

          ! send our requested keys
          call send_requests()

          ! check whether we are still waiting for data or some other communication
          if (walk_status == WALK_IAM_FINISHED) call check_comm_finished()

          ! notify rank 0 if we completed our traversal
          if (walk_status == WALK_ALL_MSG_DONE) call send_walk_finished()

          ! probe for incoming messages
          ! TODO: could be done with a blocking probe, but
          ! since my Open-MPI Version is *not thread-safe*,
          ! we will have to guarantee, that there is only one thread
          ! doing all the communication during walk...
          ! otherwise, we cannot send any messages while this thread is idling in a blocking probe
          ! and hence cannot abort this block
          ! if a blocking probe is used,
          ! the calls to send_requests() and send_walk_finished() have
          ! to be performed asynchonously (i.e. from the walk threads)
          call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_avail, stat, ierr)
          if (msg_avail) then

            comm_loop_iterations(3) = comm_loop_iterations(3) + 1
            tcomm = MPI_WTIME()

          do while (msg_avail)

          messages_per_iteration = messages_per_iteration + 1

          ipe_sender = stat(MPI_SOURCE)

          select case (stat(MPI_TAG))
             ! another PE requested child data for a certain key
             case (TAG_REQUEST_KEY)
                ! actually receive this request... ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
                call MPI_RECV( requested_key, 1, MPI_INTEGER8, ipe_sender, TAG_REQUEST_KEY, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                ! ... and answer it
                call send_data(requested_key, ipe_sender)

             ! some PE answered our request and sends
             case (TAG_REQUESTED_DATA)
                ! actually receive the data... ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
                call MPI_GET_COUNT(stat, MPI_type_multipole, num_children, ierr)
                call MPI_RECV( child_data, num_children, MPI_type_multipole, ipe_sender, TAG_REQUESTED_DATA, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                ! ... and put it into the tree and all other data structures
                call unpack_data(child_data, num_children, ipe_sender)

             ! rank 0 does bookkeeping about which PE is already finished with its walk
             ! no one else will ever receive this message tag
             case (TAG_FINISHED_PE)
                ! actually receive the data (however, we are not interested in it here) ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
                call MPI_RECV( dummy, 1, MPI_INTEGER, ipe_sender, TAG_FINISHED_PE, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

                if (walk_comm_debug) then
                  write(ipefile,'("PE", I6, " has been told that PE", I6, "has finished walking")') &
                          me, ipe_sender
                end if

                walk_finished(ipe_sender+1) = .true.

                if ( all(walk_finished) ) then
                  ! all PEs have to be informed
                  ! TODO: need better idea here...
                  if (walk_comm_debug) then
                    write(ipefile,'("PE", I6, " has found out that all PEs have finished walking - telling them to exit now")') &
                            me
                  end if

                  do i=0,num_pe-1
                    call MPI_IBSEND(comm_dummy, 1, MPI_INTEGER, i, TAG_FINISHED_ALL, &
                       MPI_COMM_WORLD, reqhandle, ierr )
                    call MPI_REQUEST_FREE( reqhandle, ierr)
                  end do
                end if


             ! all PEs have finished their walk
             case (TAG_FINISHED_ALL)
                call MPI_RECV( dummy, 1, MPI_INTEGER, ipe_sender, TAG_FINISHED_ALL, &
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication

                if (walk_comm_debug) then
                  write(ipefile,'("PE", I6, " has been told to terminate by PE", I6, " since all walks on all PEs are finished")') &
                          me, ipe_sender
                end if

                walk_status = WALK_ALL_FINISHED

          end select

          call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_avail, stat, ierr)

        end do ! while (msg_avail)

          timings_comm(TIMING_RECEIVE) = timings_comm(TIMING_RECEIVE) +  (MPI_WTIME() - tcomm)

        endif

    end subroutine





    subroutine check_comm_finished()
      implicit none
      integer :: i
      ! if request_balance contains positive values, not all requested data has arrived
      ! this means, that there were algorithmically unnecessary requests, which would be a bug
      ! negative values denote more received datasets than requested, which implies
      ! a bug in bookkeeping on the sender`s side
        if (req_queue_top .ne. req_queue_bottom) then
           write(*,*) "PE", me, " has finished its walk, but the request list is not empty"
           write(*,*) "obviously, there is an error in the todo_list bookkeeping"
           write(*,*) "Trying to recover from that"
           flush(6)
        elseif (any(request_balance .ne. 0)) then
           write(*,*) "PE", me, " finished walking but is are still waiting for requested data"
           write(*,*) "This should never happen - obviously, the request_queue or some todo_list is corrupt"

           do i=0,num_pe-1
             if (request_balance(i+1) .ne. 0) write(*,*) "        partner PE:", i, "pending requests: ", request_balance(i+1)
           end do

           write(*,*) "Trying to recover from this situation anyway."
           flush(6)
        else
           walk_status = WALK_ALL_MSG_DONE
        end if
    end subroutine




      subroutine uninit_comm_data()
        use, intrinsic :: iso_c_binding
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer :: buffsize
        integer :: dummy

        initialized = .false.

        ! free our buffer that was reserved for buffered communication
        call MPI_BUFFER_DETACH(dummy, buffsize, ierr) ! FIXME: what is the dummy thought for?
        deallocate(bsend_buffer)
        ! free the conditional variable objects
        call retval(pthreads_conds_uninit(), "pthreads_conds_uninit")
        ! free the rwlock objects
        call retval(rwlocks_uninit(), "rwlocks_uninit")

        deallocate(request_balance)

      end subroutine uninit_comm_data



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


      subroutine notify_walk_finished()
        use, intrinsic :: iso_c_binding
        use treevars
        implicit none

        walk_status = max(walk_status, WALK_IAM_FINISHED)

        if (walk_debug) then
          write(ipefile,*) "PE", me, "has finished walking"
        end if

      end subroutine notify_walk_finished





      subroutine send_walk_finished()
        use treevars
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer :: reqhandle

        ! notify rank 0 that we are finished with our walk
        call MPI_IBSEND(comm_dummy, 1, MPI_INTEGER, 0, TAG_FINISHED_PE, &
                        MPI_COMM_WORLD, reqhandle, ierr )
        call MPI_REQUEST_FREE( reqhandle, ierr)

        walk_status = WALK_I_NOTIFIED_0

      end subroutine send_walk_finished





      ! this routine is thread-safe to prevent concurrent write access to the queue
      subroutine post_request(request_key, request_addr)
        use, intrinsic :: iso_c_binding
        use treevars
        implicit none
        include 'mpif.h'
        integer*8, intent(in) :: request_key
        integer, intent(in) :: request_addr
        integer :: local_queue_bottom
        logical, save :: warned = .false.

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

        req_queue_keys(local_queue_bottom)   = request_key
        req_queue_owners(local_queue_bottom) = htable( request_addr )%owner
        req_queue_addrs(local_queue_bottom)  = request_addr

        if (walk_comm_debug) then
          write(ipefile,'("PE", I6, " posting request. local_queue_bottom=", I5, ", request_key=", O22, ", request_owner=", I6, " request_addr=", I12)') &
                         me, local_queue_bottom, request_key, htable( request_addr )%owner, request_addr
        end if

        ! now, we can tell the communicator that there is new data available
        req_queue_bottom = local_queue_bottom

        call rwlock_unlock(RWLOCK_REQUEST_QUEUE, "post_request")

        ! since there is at least one active communication request, we can wakeup the comm thread
        !call retval(pthreads_conds_signal(PTHREAD_COND_COMMBLOCK), "post_request() - pthreads_cond_signal()")

    end subroutine post_request





    subroutine send_requests()
      use, intrinsic :: iso_c_binding
      implicit none
      include 'mpif.h'

      integer :: local_queue_bottom ! buffer for avoiding interference with threads that post data to the queue
      integer :: reqhandle, ierr
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

        call rwlock_rdlock(RWLOCK_REQUEST_QUEUE, "send_requests")

        do while (req_queue_top .ne. local_queue_bottom)

          messages_per_iteration = messages_per_iteration + 1

          req_queue_top = mod(req_queue_top, REQUEST_QUEUE_LENGTH) + 1

          if (walk_comm_debug) then
            write(ipefile,'("PE", I6, " sending request.      req_queue_top=", I5, ", request_key=", O22, ", request_owner=", I6)') &
                           me, req_queue_top, req_queue_keys(req_queue_top), req_queue_owners(req_queue_top)
          end if

          if (.not. BTEST( htable(req_queue_addrs(req_queue_top))%childcode, CHILDCODE_BIT_REQUEST_SENT ) ) then
            ! send a request to PE req_queue_owners(req_queue_top)
            ! telling, that we need child data for particle request_key(req_queue_top)
            call MPI_IBSEND(req_queue_keys(req_queue_top), 1, MPI_INTEGER8, &
                            req_queue_owners(req_queue_top), TAG_REQUEST_KEY, &
                                MPI_COMM_WORLD, reqhandle, ierr )
            call MPI_REQUEST_FREE( reqhandle, ierr)

            htable(req_queue_addrs(req_queue_top))%childcode = ibset(htable(req_queue_addrs(req_queue_top))%childcode, CHILDCODE_BIT_REQUEST_SENT )

            request_balance(req_queue_owners(req_queue_top)+1) = request_balance(req_queue_owners(req_queue_top)+1) + 1

          end if
        end do

       call rwlock_unlock(RWLOCK_REQUEST_QUEUE, "send_requests")

       timings_comm(TIMING_SENDREQS) = timings_comm(TIMING_SENDREQS) + ( MPI_WTIME() - tsend )

     end if

    end subroutine send_requests







    subroutine send_data(requested_key, ipe_sender)
      implicit none
      include 'mpif.h'
      integer*8, intent(in) :: requested_key
      integer, intent(in) :: ipe_sender
      integer :: process_addr, childbyte
      type(multipole) :: children_to_send(8)
      integer :: reqhandle
      integer*8, dimension(8) :: sub_key, key_child, next_child
      integer, dimension(8) :: addr_child, node_child, byte_child, leaves_child
      integer :: j, ic, ierr, nchild

      integer :: key2addr        ! Mapping function to get hash table address from key

      if (walk_comm_debug) then
        write(ipefile,'("PE", I6, " answering request.                         request_key=", O22, ",        sender=", I6)') &
                       me, requested_key, ipe_sender
      end if

      j = 0
      process_addr = key2addr( requested_key,'WALK:send_data:parentkey')       ! get htable addresses
      childbyte    = htable( process_addr )%childcode              ! children byte-code

      nchild = int(SUM( (/ (ibits(childbyte,j,1),j=0,7) /) ))        ! Get # children

      sub_key(1:nchild)      = pack( bitarr, mask=(/ (btest(childbyte,j),j=0,7) /) )      ! Extract sub key from byte code
      key_child(1:nchild)    = IOR( ishft( requested_key,3 ), sub_key(1:nchild) )    ! Construct keys of children
      addr_child(1:nchild)   = (/( key2addr( key_child(j),'WALK:send_data:childkey' ),j=1,nchild)/)  ! Table address of children
      node_child(1:nchild)   = htable( addr_child(1:nchild) )%node                        ! Child node index
      byte_child(1:nchild)   = IAND( htable( addr_child(1:nchild) )%childcode,255 )       ! Catch lowest 8 bits of childbyte - filter off requested and here flags
      leaves_child(1:nchild) = htable( addr_child(1:nchild) )%leaves                      ! # contained leaves
      next_child(1:nchild-1) = htable( addr_child(1:nchild-1) )%next                      ! # next-node pointer
      next_child(nchild)     = -1  ! Last child gets pointed back to _parent_ for non-local nodes
                                   ! This is used to distinguish particles` walks during 'defer' phase
      ! Package children properties into user-defined multipole array for shipping
      do ic = 1,nchild
         children_to_send(ic) = multipole ( key_child(ic), &
                                      byte_child(ic), &
                                      leaves_child(ic), &
                                      me, &
                                      next_child(ic), &
                                      charge( node_child(ic) ), &
                                      abs_charge( node_child(ic) ), &
                                      xcoc( node_child(ic)), &
                                      ycoc( node_child(ic)), &
                                      zcoc( node_child(ic)), &
                                      xdip( node_child(ic)), &
                                      ydip( node_child(ic)), &
                                      zdip( node_child(ic)), &
                                      xxquad( node_child(ic)), &
                                      yyquad( node_child(ic)), &
                                      zzquad( node_child(ic)), &
                                      xyquad( node_child(ic)), &
                                      yzquad( node_child(ic)), &
                                      zxquad( node_child(ic)), &
                                      xshift( node_child(ic)), &
                                      yshift( node_child(ic)), &
                                      zshift( node_child(ic)) )
      end do

      ! Ship child data back to PE that requested it
      call MPI_IBSEND( children_to_send(1:nchild), nchild, MPI_type_multipole, ipe_sender, TAG_REQUESTED_DATA, &
              MPI_COMM_WORLD, reqhandle, ierr )
      call MPI_REQUEST_FREE(reqhandle, ierr)

      ! statistics on number of sent children-packages
      sum_ships = sum_ships + 1

    end subroutine send_data





    subroutine unpack_data(child_data, num_children, ipe_sender)
      implicit none
      include 'mpif.h'
      type (multipole) :: child_data(8) !< child data that has been received
      integer :: num_children !< actual number of valid children in dataset
      integer, intent(in) :: ipe_sender
      integer*8, dimension(8) :: sub_key
      integer*8 :: kchild, kparent(8), nxchild
      integer :: node_addr, hashaddr, lchild, nchild, nodchild, bchild, ownerchild
      integer :: j, ic, ierr

      integer :: key2addr        ! Mapping function to get hash table address from key
      integer*8 :: next_node   ! Function to get next node key for local tree walk

      request_balance(ipe_sender+1) = request_balance(ipe_sender+1) - 1

      do ic = 1, num_children
        kchild      = child_data(ic)%key
        kparent(ic) = ishft( kchild,-3 )
        bchild      = child_data(ic)%byte
        lchild      = child_data(ic)%leaves
        nxchild     = child_data(ic)%next
        ownerchild  = child_data(ic)%owner

        if (walk_comm_debug) then
          write(ipefile,'("PE", I6, " received answer.                            parent_key=", O22, ",        sender=", I6, ",        owner=", I6, ", kchild=", O22)') &
                         me, kparent(ic), ipe_sender, ownerchild, kchild
        end if

        if (lchild == 1 ) then
           nleaf = nleaf + 1
           nodchild = nleaf
           first_child(nodchild) = kchild
           ! Array bound checks
           if (nleaf>=maxleaf) then
             write (6,*) 'LPEPC | WARNING: tree arrays full on CPU ',me,' leaves',nleaf,' / ',maxleaf
             flush(6)
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           end if

        else if (lchild > 1) then
           ntwig = ntwig + 1
           nodchild = -ntwig
           nchild = int(SUM( (/ (ibits(bchild,j,1),j=0,7) /) ))   ! Get # children
           sub_key(1:nchild) = pack( bitarr(0:7), mask=(/ (btest(bchild,j),j=0,7) /) )  ! Extract child sub-keys from byte code
           first_child( nodchild ) = IOR( ishft( kchild,3), sub_key(1) )                ! Construct key of 1st (grand)child
           ! Array bound checks
           if (ntwig>=maxtwig) then
             write (6,*) 'LPEPC | WARNING: tree arrays full on CPU ',me,' twigs ',ntwig,' / ',maxtwig
             flush(6)
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           end if

        else
           write(ipefile,'(a,o15,a,i7)') '# leaves <= 0 for received child node ',kchild,' from PE ',ipe_sender
        endif

        ! Insert new node into local #-table
        call make_hashentry( kchild, nodchild, lchild, bchild, ownerchild, hashaddr, ierr )

        select case (ierr)
          case (0)
           ! anything is fine
          case (1)
           ! entry with the same key is already existing, so we just overwrite it
           write(*,*) "PE", me, "has found an already inserted entry unpack_data while calling make_hashentry(", kchild, nodchild, lchild, bchild, ipe_sender, hashaddr, ierr, ") - overwriting it"
           flush(6)
          case (2)
           ! some serious issue happened
           write(*,*) "PE", me, "has encountered problems in unpack_data while calling make_hashentry(", kchild, nodchild, lchild, bchild, ipe_sender, hashaddr, ierr, ")"
           flush(6)
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end select

        !  Determine 'next_node' pointers for 'last child' list & update hash table (tree).
        !  This  ensures that traversals in next pass treat already-fetched nodes as local,
        !  avoiding deferral list completely.
        if (nxchild == -1) then
          htable( hashaddr )%next = next_node(kchild)  !   Get next sibling, uncle, great-uncle in local tree
        else
          htable(hashaddr)%next = nxchild           ! Fill in special next-node pointer for non-local children
        endif

        node_level( nodchild ) = int(log(1.*kchild)/log(8.))  ! get level from keys and prestore as node property

        ! Physical properties
        charge( nodchild ) = child_data(ic)%q
        abs_charge( nodchild ) = child_data(ic)%absq
        xcoc( nodchild )   = child_data(ic)%xcoc
        ycoc( nodchild )   = child_data(ic)%ycoc
        zcoc( nodchild )   = child_data(ic)%zcoc
        xdip( nodchild )   = child_data(ic)%xdip
        ydip( nodchild )   = child_data(ic)%ydip
        zdip( nodchild )   = child_data(ic)%zdip
        xxquad( nodchild ) = child_data(ic)%xxquad
        yyquad( nodchild ) = child_data(ic)%yyquad
        zzquad( nodchild ) = child_data(ic)%zzquad
        xyquad( nodchild ) = child_data(ic)%xyquad
        yzquad( nodchild ) = child_data(ic)%yzquad
        zxquad( nodchild ) = child_data(ic)%zxquad
        xshift( nodchild ) = child_data(ic)%xshift
        yshift( nodchild ) = child_data(ic)%yshift
        zshift( nodchild ) = child_data(ic)%zshift

        !  Add child key to list of fetched nodes
        sum_fetches=sum_fetches+1
     end do

     ! first check, if all received particles share the same parent node
     do ic=2,num_children
       if (kparent(ic) .ne. kparent(1)) then
         write(*,*) "PE", me, "received child data for different parent nodes within one message - this is not allowed"
         write(*,*) "child_data =", child_data
         flush(6)
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
     end do
     ! mark the parent node inside the hashtable: it`s children are now accessible
     node_addr = key2addr( kparent(1),'WALK:unpack_data' )

     !call rwlock_wrlock(RWLOCK_CHILDBYTE, "unpack_data")
     htable( node_addr )%childcode = IBSET(  htable( node_addr )%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node
     !call rwlock_unlock(RWLOCK_CHILDBYTE, "unpack_data")

    end subroutine unpack_data

end module tree_walk_communicator









module tree_walk_utils
  use treetypes
  use treevars
  use pthreads_stuff
  implicit none

  private
    integer, public :: num_walk_threads = 3
    integer, private :: primary_processor_id = 0
    real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8
    real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8

    real*8, dimension(:), allocatable :: boxlength2
    integer :: mac
    real*8 :: theta2
    real*8 :: vbox(3)
    type(calc_force_params) :: cf_par
    logical :: in_central_box


    integer :: np_local
    real*8, dimension(:), pointer :: work_per_particle => NULL()

    integer, allocatable :: pshort_thread(:,:) !< start- and endpoint of shortlist-chunks to be processed by the distinct walk threads

    integer :: next_unassigned_particle !< index of next particle that has not been assigned to a work thread
    integer :: finished_threads !< number of worker threads that have terminated. if equal to num_walk_threads, all particles have been completely processed on this PE

    !> type for input and return values of walk_threads
    type t_threaddata
      integer :: id                         !< just a running number to distinguish the threads, currently unused
      integer :: num_processed_particles  !< thread output value: number of particles that it has processed
      real*8  :: num_interactions !< thread output value: number of interactions that were performed
      real*8  :: num_mac_evaluations !< thread output value: number of mac evaluations that have been performed
      logical :: is_on_shared_core !< thread output value: is set to true if the thread detects that it shares its processor with the communicator thread
      integer :: coreid !< thread output value: id of thread`s processor
      real*8 :: runtime_seconds !< thread wallclock-runtime in seconds, measured with MPI_WTIME()
    end type t_threaddata


  public tree_walk

  contains


    subroutine tree_walk(np_local_,theta_,cf_par_,itime,mac_,twalk,tfetch,vbox_,work_per_particle_, tcomm)
      use, intrinsic :: iso_c_binding
      use treetypes
      use tree_utils
      use timings
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
      real*8, target, intent(inout) :: twalk, tfetch
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

      tfetch = 0.
      twalk  = MPI_WTIME()

      call init_walk_data()

      call init_comm_data()

      call walk_hybrid()

      call uninit_comm_data()

      call uninit_walk_data()

      twalk = MPI_WTIME() - twalk
      tcomm = timings_comm

    end subroutine tree_walk



    subroutine walk_hybrid()
      use tree_walk_communicator
      use, intrinsic :: iso_c_binding
      implicit none
      include 'mpif.h'

      integer :: ith, displ, ierr
      type(t_threaddata), target :: threaddata(num_walk_threads)

      ! we count how many threads already have finished
      finished_threads = 0

      ! start the worker threads...
      do ith = 1,num_walk_threads
        threaddata(ith)%id = ith
        call retval(pthreads_createthread(ith, c_funloc(walk_worker_thread), c_loc(threaddata(ith))), "walk_schedule_thread_inner:pthread_create")
      end do

      call run_communication_loop()

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

    end subroutine walk_hybrid



    subroutine init_walk_data()
      use, intrinsic :: iso_c_binding
      use tree_walk_communicator
      implicit none
      integer :: i

      ! we have to have at least one walk thread
      num_walk_threads = max(num_walk_threads, 1)
      ! evenly balance particles to threads if there are less than the maximum
      max_particles_per_thread = max(min(np_local/num_walk_threads, max_particles_per_thread),1)
      ! allocate storage for thread handles, the 0th entry is the walk scheduler thread, the other ones are the walk worker threads
      call retval(pthreads_init(num_walk_threads + 1), "init_walk_data:pthreads_init")
      allocate(pshort_thread(2, num_walk_threads))
      ! we will only want to reject the root node and the particle itself if we are in the central box
      in_central_box = (dot_product(vbox,vbox) == 0)
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
      use, intrinsic :: iso_c_binding
      use tree_walk_communicator
      implicit none
      deallocate(boxlength2)
      deallocate(pshort_thread)
      call retval(pthreads_uninit(), "uninit_walk_data:pthreads_uninit")
    end subroutine uninit_walk_data



    function get_first_unassigned_particle(success)
      use tree_walk_communicator
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
      use tree_walk_communicator
      use pthreads_stuff
      implicit none
      include 'mpif.h'
      type(c_ptr) :: walk_worker_thread
      type(c_ptr), value :: arg

      integer, dimension(:), allocatable :: my_particles
      integer*8, dimension(:,:), allocatable :: todo_list
      integer, dimension(:), allocatable :: todo_list_top, todo_list_bottom
      integer, dimension(:), allocatable :: todo_list_minlevel_next
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
            comm_on_shared_processor = .true.
      else
            my_max_particles_per_thread = max_particles_per_thread
      endif


      allocate(my_particles(my_max_particles_per_thread),                        &
                            todo_list_top(my_max_particles_per_thread),           &
                            todo_list_bottom(my_max_particles_per_thread),        &
                            todo_list_minlevel_next(my_max_particles_per_thread));
      allocate(todo_list(nintmax,my_max_particles_per_thread));

      ! every particle will start at the root node (one entry per todo_list, no particle is finished)
      todo_list_top       = 0
      todo_list_bottom    = 1
      todo_list           = 0
      todo_list(1,:)      = 1
      todo_list_minlevel_next = 0
      my_particles(:)     = -1
      particles_available = .true.
      particles_active    = .true.
      call c_f_pointer(arg, my_threaddata)
      my_threaddata = t_threaddata(my_threaddata%id, 0, 0._8, 0._8, same_core_as_communicator, my_processor_id, MPI_WTIME())
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

            if (walk_single_particle(my_particles(i), nintmax, todo_list(1:nintmax,i), todo_list_top(i), todo_list_bottom(i), todo_list_minlevel_next(i), my_threaddata)) then
              ! this particle`s walk has finished
!              write(ipefile,*) "PE", me, getfullid(), " walk for particle", i, " nodeidx =", my_particles(i), " has finished"

              todo_list_top(i)           =  0
              todo_list_bottom(i)        =  1
              todo_list(:,i)             =  0
              todo_list(1,i)             =  1
              todo_list_minlevel_next(i) =  0
              my_particles(i)            = -1
              my_threaddata%num_processed_particles = my_threaddata%num_processed_particles + 1
            else
              particles_active = .true.
            end if
          else
            particles_available = .false.
          end if

        end do

      end do

      deallocate(my_particles, todo_list_top, todo_list_bottom, todo_list_minlevel_next);
      deallocate(todo_list);

      ! bookekeeping about finished worker threads
      call rwlock_wrlock(RWLOCK_FINISHED_THREADS, "walk_worker_thread: increment finished_threads")

      finished_threads = finished_threads + 1

      if (finished_threads == num_walk_threads) then
        ! tell rank 0 that we are finished with our walk
        call notify_walk_finished()
      end if

      call rwlock_unlock(RWLOCK_FINISHED_THREADS, "walk_worker_thread: increment finished_threads")

      if (same_core_as_communicator) comm_on_shared_processor = .false.

      walk_worker_thread = c_null_ptr

      my_threaddata%runtime_seconds = MPI_WTIME() - my_threaddata%runtime_seconds

      call retval(pthreads_exitthread(), "walk_worker_thread:pthread_exit")

    end function walk_worker_thread




   function walk_single_particle(nodeidx, nintmax, todo_list, todo_list_top, todo_list_bottom, todo_list_minlevel_next, my_threaddata)
      use tree_walk_communicator
      use module_calc_force
      implicit none
      include 'mpif.h'
      integer, intent(in) :: nodeidx
      integer, intent(in) :: nintmax
      integer*8, intent(inout) :: todo_list(nintmax)
      integer, intent(inout) :: todo_list_top, todo_list_bottom, todo_list_minlevel_next
      type(t_threaddata), intent(inout) :: my_threaddata
      logical :: walk_single_particle !< function will return .true. if this particle has finished its walk
      logical :: result

      integer*8 :: walk_key, next_key
      integer :: walk_addr, walk_node

      integer :: key2addr        ! Mapping function to get hash table address from key
      integer :: newtop
      real*8 :: dist2
      real*8 :: delta(3)
      logical :: ignore, mac_ok
      integer :: ierr

      newtop    = todo_list_top ! todo_list_top( will only be updated if really necessary

      result = (newtop == todo_list_bottom)
      walk_single_particle = result

      if (result) then
        if (walk_debug) then
          write(ipefile,'("PE", I6, " particle ", I12, " has an empty todo_list: particle obviously finished walking around :-)")') me, nodeidx
        end if

        return ! nothing to walk here any more
      end if

      ! read next todo_list-entry
      newtop = mod(newtop, nintmax) + 1

      walk_key  = todo_list(newtop)
      walk_addr = key2addr( walk_key, 'WALK:walk_single_particle' )  ! get htable address
      walk_node = htable( walk_addr )%node            ! Walk node index - points to multipole moments
      next_key  = htable( walk_addr )%next            ! Next node pointer

      if (int(log(1.*next_key)/log(8.)) < todo_list_minlevel_next) then ! TODO: use pretabulated values for node level
        next_key = 1                                    ! special case for nodes from the deferral list to prevent restarting the walk with their siblings again and again
      end if


! TODO:  BH MAC is also implemented in mac_choose routine,
!        which needs tuning and inlining to avoid excessive parameter-passing overhead
!        Other MACS need additional preprocessed info on tree nodes (box lengths etc)
!        before they can be used for performance tests.
!        mac=1,2 only useful for accuracy tests at the moment (interaction list comparisons)

       if (mac==0) then
          ! BH-MAC
          delta(1) = x(nodeidx) - (xcoc( walk_node ) + vbox(1) )     ! Separations
          delta(2) = y(nodeidx) - (ycoc( walk_node ) + vbox(2) )
          delta(3) = z(nodeidx) - (zcoc( walk_node ) + vbox(3) )

          dist2 = DOT_PRODUCT(delta, delta)

          mac_ok = (theta2 * dist2 > boxlength2(node_level(walk_node)))
       else
!             call mac_choose(pshort(p),ex_nps(p),ey_nps(p),ez_nps(p),&
!                  walk_node,walk_key(idx),abs_charge(walk_node),boxlength2(node_level(walk_node)), &
!                  theta2,mac,mac_ok, vbox)
          mac_ok = .false.
       endif

       mac_ok = ( mac_ok .and. ((.not. in_central_box) .or. walk_key>1 ) )  !  reject root node if we are in the central box

       my_threaddata%num_mac_evaluations = my_threaddata%num_mac_evaluations + 1._8
       work_per_particle(nodeidx) = work_per_particle(nodeidx) + WORKLOAD_PENALTY_MAC

       ! set ignore flag if leaf node corresponds to particle itself (number in pshort)
       ! NB: this uses local leaf #, not global particle label
       ! only ignore particle if we are processing t.not. all_my_particles_finishedhe central box
       ! when walking through neighbour boxes, it has to be included
       ignore = (in_central_box) .and. ( nodeidx == walk_node )

       ! Possible courses of action:

       ! 1) MAC test OK, or leaf node, so put cell on interaction list and find next node for tree walk
       !    - reject self

       if ( (mac_ok .or. walk_node > 0) .and. .not.ignore) then

          call calc_force_per_interaction(nodeidx, walk_node, vbox, cf_par)
          work_per_particle(nodeidx) = work_per_particle(nodeidx) + WORKLOAD_PENALTY_INTERACTION
          my_threaddata%num_interactions = my_threaddata%num_interactions + 1._8

          ! if walk_key is already the rootnode, next_key is zero --> walk is finished
          ! if next_key == 1, then we reached root --> walk also finished
          todo_list(newtop) = next_key


       else  if ( .not.mac_ok .and. walk_node < 0 ) then

          if ( btest(htable( walk_addr )%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE) ) then
              ! 2) MAC fails at node for which children present, so resolve cell & put 1st child in front of todo_list
              ! if local put 1st child node on walk_list
              todo_list(newtop) = first_child( walk_node )
          else
              ! 3) MAC fails at node for which children _absent_, so put node on REQUEST list (flag with add=2) and
              !    put walk_key on bottom of todo_list
              todo_list_bottom  = mod(todo_list_bottom, nintmax) + 1

              if (walk_debug) then
                write(ipefile,'("PE ", I6, " adding nonlocal key to tail for particle ", I20, " todo_list_bottom=", I6)') me, nodeidx, todo_list_bottom
              end if

              if (todo_list_bottom == todo_list_top) then ! the "overfull todo_list" error can be avoided by simply deferring
                                                          ! the current walk_key and going on with the next one
                                                          ! TODO: abort this process if there is no hope anymore
                write(*,'("Rather serious issue on PE ", I6, ": todo_list is full for particle ", I20, " nintmax =", I6, " is too small")') me, nodeidx, nintmax
                write(*,*) "We could skip the current particle until some child_data has arrive, but this can result in a deadlock... --> aborting."
                flush(6)
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
              else
                todo_list( todo_list_bottom ) = walk_key ! Deferred list of nodes to search, pending request
                                                         ! for data from nonlocal PEs
                ! check wether the node has already been requested
                ! this if-construct has to be secured against synchronous invocation (together with the modification while receiving data)
                ! otherwise it will be possible that two walk threads can synchronously post a prticle to the request queue

                if (.not. BTEST( htable(walk_addr)%childcode, CHILDCODE_BIT_REQUEST_POSTED ) ) then  ! Check if node already requested
                    call post_request(walk_key, walk_addr)        ! tell the communicator about our needs
                    ! if posting the request failed, this is not a problem, since we defer the particle anyway
                    ! since it will not be available then, the request will simply be repeated
                end if

                todo_list(newtop) = next_key  ! Continue with walk for now
              end if

          endif

       else
          ! 4) particle and leaf node identical, so skip
          todo_list(newtop) = next_key  ! Continue with walk
       endif

       ! Possible special cases now:
       !  1. ) we reached last of nonlocal children (walk_key == -1)
       !  2a.) we reached root node (walk_key == 1) while still having pending nodes in our todo_list
       !  2b.) we just took a deferred node - its next_node has already been processed before, the next_node value has been set to 1 (root node) manually
       !  3. ) next_key was zero because we directly interacted with box root
       ! then we simply point to the next entry in the todo_list
       if (( todo_list(newtop) == -1 ) .or. ( todo_list(newtop) == 1 ) .or. ( todo_list(newtop) == 0 )) then
         ! skip the entry that is currently on top of the todo_list and continue with deferred keys
         todo_list_top           = newtop
         todo_list_minlevel_next = int(log(1.*todo_list(mod(newtop,nintmax)+1))/log(8.))+1
       end if


    end function walk_single_particle

end module tree_walk_utils















