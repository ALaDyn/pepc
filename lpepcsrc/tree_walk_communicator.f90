!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates helper functions that simplify communication during tree traversal
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module tree_walk_communicator
  use treevars
  use pthreads_stuff
  implicit none

  private

    integer, public :: walk_status

    !> debug flags - cannot be modified at runtime due to performance reasons
    logical, parameter :: walk_comm_debug = .false.
    logical, parameter :: rwlock_debug    = .false.
    logical, parameter, public  :: walk_debug     = .false.

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

    ! variables for adjusting the thread`s workload
    integer, public :: max_particles_per_thread = 2000 !< maximum number of particles that will in parallel be processed by one workthread
    real, public :: work_on_communicator_particle_number_factor = 0.1 !< factor for reducing max_particles_per_thread for thread which share their processor with the communicator
    integer, public :: particles_per_yield = 500 !< number of particles to process in a work_thread before it shall call sched_yield to hand the processor over to some other thread
    integer, private :: messages_per_iteration !< tracks current number of received and transmitted messages per commloop iteration for adjusting particles_per_yield
    integer, parameter :: MAX_MESSAGES_PER_ITERATION = 20
    integer, parameter :: MIN_MESSAGES_PER_ITERATION = 5
    logical, public :: comm_on_shared_processor = .false.

    !> data type for internal request queue
    type t_request_queue_entry
      integer*8 :: key
      integer   :: addr
      integer   :: owner
    end type

    ! internal communication variables - not to be touched from outside the module
    integer, parameter :: ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
    integer, parameter :: REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue
    type(t_request_queue_entry), volatile, private, target :: req_queue(REQUEST_QUEUE_LENGTH)
    integer, private, volatile :: req_queue_top, req_queue_bottom ! we will insert data at bottom and take from top
    integer*1, private, allocatable, target :: bsend_buffer(:) !< buffer for bsend-alls
    integer, private :: comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)
    integer, private, allocatable :: request_balance(:) !< (#requests - #answers) per PE
    real*8, public :: timings_comm(3) !< array for storing internal timing information
    integer, parameter :: cond_comm_timeout = 1000 !< timeout for pthread_cond_wait() in commloop in microseconds TODO: tune this parameter automatically during runtime

    ! rwlocks for regulating concurrent access
    integer, private, parameter :: NUM_RWLOCKS = 2
    integer, public, parameter :: RWLOCK_REQUEST_QUEUE      = 1
    integer, public, parameter :: RWLOCK_NEXT_FREE_PARTICLE = 2

    ! pthread conditional variables
    integer, private, parameter :: NUM_CONDS = 1
    integer, public, parameter :: PTHREAD_COND_COMMBLOCK = 1

    ! IDs for internal timing measurement
    integer, public, parameter :: TIMING_COMMLOOP = 1
    integer, public, parameter :: TIMING_RECEIVE  = 2
    integer, public, parameter :: TIMING_SENDREQS = 3


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
        implicit none
        include 'mpif.h'
        logical :: msg_avail
        integer :: ierr
        integer :: stat(MPI_STATUS_SIZE)
        integer :: reqhandle
        integer*8 :: requested_key
        type (multipole), allocatable :: child_data(:) ! child data to be received - maximum up to eight children per particle
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
                allocate(child_data(num_children))
                call MPI_RECV( child_data, num_children, MPI_type_multipole, ipe_sender, TAG_REQUESTED_DATA, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                ! ... and put it into the tree and all other data structures
                call unpack_data(child_data, num_children, ipe_sender)
                deallocate(child_data)

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
        use treevars
        implicit none

        if (all(threaddata(:)%finished)) then

          walk_status = max(walk_status, WALK_IAM_FINISHED)

          if (walk_debug) then
            write(ipefile,*) "PE", me, "has finished walking"
          end if

        endif

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

        twalk_loc = MPI_WTIME() - twalk_loc

      end subroutine send_walk_finished





      ! this routine is thread-safe to prevent concurrent write access to the queue
      subroutine post_request(request_key, request_addr)
        use treevars
        use module_htable
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

        ! since there is at least one active communication request, we can wakeup the comm thread
        !call retval(pthreads_conds_signal(PTHREAD_COND_COMMBLOCK), "post_request() - pthreads_cond_signal()")

    end subroutine post_request





    subroutine send_requests()
      use module_htable
      implicit none
      include 'mpif.h'

      integer :: local_queue_bottom ! buffer for avoiding interference with threads that post data to the queue
      integer :: reqhandle, ierr
      real*8 :: tsend
      integer*8 :: req_queue_length
      type(t_request_queue_entry), pointer :: req

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

          req=>req_queue(req_queue_top)

            if (walk_comm_debug) then
              write(ipefile,'("PE", I6, " sending request.      req_queue_top=", I5, ", request_key=", O22, ", request_owner=", I6)') &
                           me, req_queue_top, req%key, req%owner
            end if

            if (.not. BTEST( htable(req%addr)%childcode, CHILDCODE_BIT_REQUEST_SENT ) ) then
              ! send a request to PE req_queue_owners(req_queue_top)
              ! telling, that we need child data for particle request_key(req_queue_top)
              call MPI_IBSEND(req%key, 1, MPI_INTEGER8, req%owner, TAG_REQUEST_KEY, MPI_COMM_WORLD, reqhandle, ierr )
              call MPI_REQUEST_FREE( reqhandle, ierr)

              htable(req%addr)%childcode = ibset(htable(req%addr)%childcode, CHILDCODE_BIT_REQUEST_SENT )

              request_balance(req%owner+1) = request_balance(req%owner+1) + 1

            end if

        end do

       call rwlock_unlock(RWLOCK_REQUEST_QUEUE, "send_requests")

       timings_comm(TIMING_SENDREQS) = timings_comm(TIMING_SENDREQS) + ( MPI_WTIME() - tsend )

     end if

    end subroutine send_requests







    subroutine send_data(requested_key, ipe_sender)
      use module_htable
      implicit none
      include 'mpif.h'
      integer*8, intent(in) :: requested_key
      integer, intent(in) :: ipe_sender
      integer :: process_addr
      type(multipole) :: children_to_send(8)
      integer :: reqhandle
      integer*8, dimension(8) :: key_child
      integer, dimension(8) :: addr_child, node_child, byte_child, leaves_child, owner_child
      integer :: j, ic, ierr, nchild

      if (walk_comm_debug) then
        write(ipefile,'("PE", I6, " answering request.                         request_key=", O22, ",        sender=", I6)') &
                       me, requested_key, ipe_sender
      end if

      j = 0
      process_addr = key2addr( requested_key,'WALK:send_data:parentkey')       ! get htable addresses

      call get_childkeys(process_addr, nchild, key_child)

      addr_child(1:nchild)   = (/( key2addr( key_child(j),'WALK:send_data:childkey' ),j=1,nchild)/)  ! Table address of children
      node_child(1:nchild)   = htable( addr_child(1:nchild) )%node                        ! Child node index
      byte_child(1:nchild)   = int(IAND( htable( addr_child(1:nchild) )%childcode, CHILDCODE_CHILDBYTE ))! Catch lowest 8 bits of childbyte - filter off requested and here flags
      leaves_child(1:nchild) = htable( addr_child(1:nchild) )%leaves                      ! # contained leaves
      owner_child(1:nchild)  = htable( addr_child(1:nchild) )%owner                       ! real owner of child (does not necessarily have to be identical to me, at least after futural modifications)
      ! Package children properties into user-defined multipole array for shipping
      do ic = 1,nchild
         associate(c=>children_to_send(ic))
           c        = tree_nodes(node_child(ic))
           c%key    = key_child(ic)   ! TODO: this data is maybe not consistently stored in tree_nodes array
           c%byte   = byte_child(ic)  ! therefore, we have to take it directly form the htable --> repair this
           c%leaves = leaves_child(ic)
           c%owner  = owner_child(ic)
          end associate
      end do

      ! Ship child data back to PE that requested it
      call MPI_IBSEND( children_to_send(1:nchild), nchild, MPI_type_multipole, ipe_sender, TAG_REQUESTED_DATA, &
              MPI_COMM_WORLD, reqhandle, ierr )
      call MPI_REQUEST_FREE(reqhandle, ierr)

      ! statistics on number of sent children-packages
      sum_ships = sum_ships + 1

    end subroutine send_data





    subroutine unpack_data(child_data, num_children, ipe_sender)
      use module_htable
      use module_spacefilling
      implicit none
      include 'mpif.h'
      type (multipole) :: child_data(num_children) !< child data that has been received
      integer :: num_children !< actual number of valid children in dataset
      integer, intent(in) :: ipe_sender
      integer*8 :: kchild, kparent
      integer :: hashaddr, lchild,  nodchild, bchild, ownerchild, parent_addr(num_children)
      integer :: ic, ierr

      request_balance(ipe_sender+1) = request_balance(ipe_sender+1) - 1

      do ic = 1, num_children
        kchild      = child_data(ic)%key
        kparent     = ishft( kchild,-3 )
        bchild      = child_data(ic)%byte
        lchild      = child_data(ic)%leaves
        ownerchild  = child_data(ic)%owner
        ! save parent address - after (!) inserting all (!) children we can flag it: it`s children are then accessible
        parent_addr(ic) = key2addr( kparent, 'WALK:unpack_data() - get parent address' )

        if (walk_comm_debug) then
          write(ipefile,'("PE", I6, " received answer.                            parent_key=", O22, ",        sender=", I6, ",        owner=", I6, ", kchild=", O22)') &
                         me, kparent, ipe_sender, ownerchild, kchild
        end if

        if (lchild == 1 ) then
           nleaf = nleaf + 1
           nodchild = nleaf
           ! Array bound checks
           if (nleaf>=maxleaf) then
             write (6,*) 'LPEPC | WARNING: tree arrays full on CPU ',me,' leaves',nleaf,' / ',maxleaf
             flush(6)
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           end if

        else if (lchild > 1) then
           ntwig = ntwig + 1
           nodchild = -ntwig
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

        ! Physical properties
        tree_nodes( nodchild ) = child_data( ic )

        !  Add child key to list of fetched nodes
        sum_fetches=sum_fetches+1
     end do

     ! set 'children-here'-flag for all parent addresses
     ! may only be done *after inserting all* children, hence not(!) during the loop above
     do ic=1,num_children
         htable( parent_addr(ic) )%childcode = IBSET(  htable( parent_addr(ic) )%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node
     end do


    end subroutine unpack_data

end module tree_walk_communicator
