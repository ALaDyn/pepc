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
!> Defines communication routines for the distributed hashed k-ary tree.
!>
!  Structure of communicator:
!  --------------------------
!    do while (not all PEs finished)
!
!        if (requests have been posted)
!            send all requests
!
!        if (shutdown requested)
!            send_comm_finished(to rank 0)
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
module module_tree_communicator
  use module_tree, only: t_request_queue_entry, t_tree_communicator, &
    TREE_COMM_REQUEST_QUEUE_LENGTH, TREE_COMM_ANSWER_BUFF_LENGTH
  implicit none
  private

  !> debug flags - cannot be modified at runtime due to performance reasons
  logical, public, parameter :: tree_comm_debug = .false.

  integer, private :: tree_comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)

  ! tags to be used in communication
  integer, public, parameter :: TREE_COMM_TAG_REQUEST_KEY    = 1257 !< message requesting child data for a certain key
  integer, public, parameter :: TREE_COMM_TAG_REQUESTED_DATA = 1258 !< message that contains requested child data
  integer, public, parameter :: TREE_COMM_TAG_FINISHED_PE    = 1259 !< message to rank 0, to announce requested stop of communicator
  integer, public, parameter :: TREE_COMM_TAG_FINISHED_ALL   = 1260 !< message from rank 0, that communication has finished
  integer, public, parameter :: mintag = TREE_COMM_TAG_REQUEST_KEY
  integer, public, parameter :: maxtag = TREE_COMM_TAG_FINISHED_ALL
  
  ! IDs for internal timing measurement
  integer, public, parameter :: TREE_COMM_TIMING_COMMLOOP = 1
  integer, public, parameter :: TREE_COMM_TIMING_RECEIVE  = 2
  integer, public, parameter :: TREE_COMM_TIMING_SENDREQS = 3

  integer, parameter :: TREE_COMM_MAX_MESSAGES_PER_ITERATION = 20
  integer, parameter :: TREE_COMM_MIN_MESSAGES_PER_ITERATION = 5

  ! MPI buffer
  integer*1, allocatable, target :: tree_comm_bsend_buffer(:) !< buffer for bsend-calls
  
  public :: tree_communicator_start
  public :: tree_communicator_stop
  public :: tree_node_fetch_children
  public :: tree_communicator_prepare
  public :: tree_communicator_finalize

  contains

  !>
  !> One time initialization of MPI communication buffers to be called from
  !> pepc_prepare.
  !>
  !> @todo need more space for multiple communicators?
  !>
  subroutine tree_communicator_prepare()
    use treevars, only: mpi_comm_lpepc
    use module_pepc_types, only: mpi_type_tree_node_package
    implicit none
    include 'mpif.h'

    integer :: msg_size_request, msg_size_data, buffsize, ierr

    if (.not. allocated(tree_comm_bsend_buffer)) then
      ! compute upper bounds for request and data message size
      call MPI_PACK_SIZE(1, MPI_INTEGER8, mpi_comm_lpepc, msg_size_request, ierr)
      msg_size_request = msg_size_request + MPI_BSEND_OVERHEAD

      call MPI_PACK_SIZE(8, MPI_TYPE_tree_node_package, mpi_comm_lpepc, msg_size_data, ierr)
      msg_size_data = msg_size_data + MPI_BSEND_OVERHEAD

      buffsize = (TREE_COMM_REQUEST_QUEUE_LENGTH * msg_size_request + TREE_COMM_ANSWER_BUFF_LENGTH * msg_size_data)
      ! reserve memory for buffered mpi communication
      allocate(tree_comm_bsend_buffer(buffsize))
      ! and tell mpi , where it can be found
      call MPI_BUFFER_ATTACH(tree_comm_bsend_buffer, buffsize, ierr)
    end if
  end subroutine tree_communicator_prepare
    

  !>
  !> Deallocate MPI communication buffers in pepc_finalize.
  !>
  subroutine tree_communicator_finalize()
    implicit none
    include 'mpif.h'

    integer :: dummy, buffsize, ierr

    if (allocated(tree_comm_bsend_buffer)) then
      ! free our buffer that was reserved for buffered communication
      call MPI_BUFFER_DETACH(dummy, buffsize, ierr)
      deallocate(tree_comm_bsend_buffer)
    end if
  end subroutine tree_communicator_finalize


  !>
  !> Starts the communication thread associated with tree `t`.
  !>
  !> The thread answers incoming requests for tree nodes on the communicator in
  !> `t%comm_env` and sends out requests accumulated in the associated request
  !> queue.
  !>
  !> @note Smooth sailing is only expected as long as `t%comm_env` is not used
  !> concurrently by any other thread due to varying levels of thread safeness
  !> of MPI implementations.
  !>
  !> @todo Add an MPI barrier here for symmetry with `tree_communicator_stop()`?
  !>
  subroutine tree_communicator_start(t)
    use, intrinsic :: iso_c_binding
    use module_tree, only: t_tree
    use pthreads_stuff, only: pthreads_createthread
    use module_debug
    use module_timings
    implicit none

    type(t_tree), target, intent(inout) :: t
    type(c_ptr) :: tp
    integer(c_int) :: res

    DEBUG_ASSERT(.not. t%communicator%comm_thread_running)
    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " run_communication_loop start.")', t%comm_env%rank)
    end if
    ! TODO: in future, need to handle multiple communicators.
    call timer_reset(t_comm_total)
    call timer_reset(t_comm_recv)
    call timer_reset(t_comm_sendreqs)
    
    t%communicator%comm_thread_stopping = .false.
    t%communicator%comm_thread_stop_requested = .false.
    tp = c_loc(t)
    res = pthreads_createthread(t%communicator%comm_thread, c_funloc(run_communication_loop), tp)
    if (0 /= res) then
      DEBUG_ERROR(*, "pthreads_createthread() failed with return value: ", res)
    end if
    t%communicator%comm_thread_running = .true.
  end subroutine tree_communicator_start


  !>
  !> Stops the communication thread associated with tree `t`.
  !>
  !> Stopping tree communication is a global action across all ranks owning
  !> parts of the tree and can only complete once everyone has agreed that
  !> communications should come to an end. This routine blocks until consensus
  !> achieved an the communicator thread can actually be stopped.
  !>
  subroutine tree_communicator_stop(t)
    use module_tree, only: t_tree
    use pthreads_stuff, only: pthreads_jointhread
    use module_debug
    use module_timings
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t

    DEBUG_ASSERT(t%communicator%comm_thread_running)
    t%communicator%comm_thread_stop_requested = .true.
    if (0 /= pthreads_jointhread(t%communicator%comm_thread)) then
      DEBUG_ERROR(*, "pthreads_jointhread() failed!")
    end if
    t%communicator%comm_thread_stopping = .false.
    t%communicator%comm_thread_running = .false.

    call timer_add(t_comm_total,    t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP))
    call timer_add(t_comm_recv,     t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE))
    call timer_add(t_comm_sendreqs, t%communicator%timings_comm(TREE_COMM_TIMING_SENDREQS))

    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " run_communication_loop end.")', t%comm_env%rank)
    end if
  end subroutine tree_communicator_stop


  !>
  !> Request children of node `n` in tree `t` from the responsible remote rank.
  !> The local node that actually needs the remote node is `n_targ`.
  !>
  !> This routine returns immediately as the actual communication with the
  !> remote rank is handled by the communicator thread. The caller is then free
  !> to continue working on something different. Later on, `tree_node_children_available`
  !> can be used to check whether the requested data has arrived in the
  !> meantime.
  !> Information from `n_targ` will be used for the eager send algorithm.
  !> If `n_targ` is omitted, no eager sending will be requested.
  !>
  subroutine tree_node_fetch_children(t, n, pos_targ)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_tree_node
    use module_atomic_ops, only: atomic_mod_increment_and_fetch_int, &
      atomic_write_barrier
    use module_tree_node
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    type(t_tree_node), target, intent(inout) :: n
    real*8, optional, intent(in) :: pos_targ(3)
    
    integer :: local_queue_bottom

    ! check wether the node has already been requested
    ! this if-construct has to be secured against synchronous invocation (together with the modification while receiving data)
    ! otherwise it will be possible that two walk threads can synchronously post a particle to the request queue
    if (btest(n%flags, TREE_NODE_FLAG_REQUEST_POSTED)) then
      return
    end if

    ! we first flag the particle as having been already requested to prevent other threads from doing it while
    ! we are inside this function
    n%flags = ibset(n%flags, TREE_NODE_FLAG_REQUEST_POSTED) ! Set requested flag

    ! thread-safe way of reserving storage for our request
    local_queue_bottom = atomic_mod_increment_and_fetch_int(t%communicator%req_queue_bottom, TREE_COMM_REQUEST_QUEUE_LENGTH)

    if (local_queue_bottom == t%communicator%req_queue_top) then
      DEBUG_ERROR(*, "Issue with request sending queue: TREE_COMM_REQUEST_QUEUE_LENGTH is too small: ", TREE_COMM_REQUEST_QUEUE_LENGTH)
    end if

    ! the communicator will check validity of the request and will only proceed as soon as the entry is valid -- this actually serializes the requests
    t%communicator%req_queue(local_queue_bottom)%request%key             =  n%key
    if (present(pos_targ)) then
      t%communicator%req_queue(local_queue_bottom)%request%pos_target    = pos_targ
      t%communicator%req_queue(local_queue_bottom)%request%eager_request = .true.
    else
      t%communicator%req_queue(local_queue_bottom)%request%eager_request = .false.
    endif
    t%communicator%req_queue(local_queue_bottom)%node => n
    call atomic_write_barrier() ! make sure the above information is actually written before flagging the entry valid by writing the owner
    t%communicator%req_queue(local_queue_bottom)%owner = n%owner

    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " posting request. local_queue_bottom=", I5, ", request_key=", O22, ", request_owner=", I6)', t%comm_env%rank, local_queue_bottom, n%key, n%owner)
    end if
  end subroutine tree_node_fetch_children


  !>
  !> Notify all ranks about a communication shutdown.
  !>
  subroutine broadcast_comm_finished(t)
    use module_tree, only: t_tree
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(in) :: t
    integer :: i, ierr

    ! all PEs have to be informed
    ! TODO: need better idea here...
    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " has found out that all PEs have finished walking - telling them to exit now")', t%comm_env%rank)
    end if

    do i = 0, t%comm_env%size - 1
      call MPI_BSEND(tree_comm_dummy, 1, MPI_INTEGER, i, TREE_COMM_TAG_FINISHED_ALL, &
        t%comm_env%comm, ierr)
    end do
  end subroutine


  !>
  !> Answer a request for node data.
  !>
  subroutine answer_request(t, request, ipe_sender)
    use module_tree, only: t_tree, tree_lookup_node_critical
    use module_tree_node, only: tree_node_pack, tree_node_get_first_child, tree_node_get_next_sibling
    use module_pepc_types, only: t_tree_node, t_tree_node_package, MPI_TYPE_tree_node_package, t_request
    use module_debug
    use module_tree_node
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    type(t_Request), intent(in) :: request
    integer, intent(in) :: ipe_sender

    type(t_tree_node_package) :: children_to_send(8)
    type(t_tree_node), pointer :: n, nn
    integer :: nchild, ierr

    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " answering request.                         request%key=", O22, ",        sender=", I6)', t%comm_env%rank, request%key, ipe_sender )
    end if

    call tree_lookup_node_critical(t, request%key, n, 'WALK:answer_request_simple:parentkey')
    if (tree_node_get_first_child(n, nn)) then
      nchild = 0

      do
        n => nn
        nchild = nchild + 1

        call tree_node_pack(n, children_to_send(nchild))
        children_to_send(nchild)%flags = int(iand( children_to_send(nchild)%flags, TREE_NODE_CHILDBYTE ))! Catch lowest 8 bits of childbyte - filter off requested and here flags

        if (.not. tree_node_get_next_sibling(n, nn)) then; exit; end if
      end do
      
      ! Ship child data back to PE that requested it
      call MPI_BSEND(children_to_send(1:nchild), nchild, MPI_TYPE_tree_node_package, &
        ipe_sender, TREE_COMM_TAG_REQUESTED_DATA, t%comm_env%comm, ierr)

      ! statistics on number of sent children-packages
      t%communicator%sum_ships = t%communicator%sum_ships + 1
    end if
  end subroutine answer_request


  !>
  !> Send a request for data
  !>
  function send_request(t, req)
    use module_tree, only: t_tree
    use module_tree_node
    use module_pepc_types, only : t_request, MPI_TYPE_request
    implicit none
    include 'mpif.h'
    
    logical :: send_request
    type(t_tree), intent(inout) :: t
    type(t_request_queue_entry), intent(in) :: req

    integer :: ierr

    if (.not. btest( req%node%flags, TREE_NODE_FLAG_REQUEST_SENT ) ) then
      ! send a request to PE req_queue_owners(req_queue_top)
      ! telling, that we need child data for particle request_key(req_queue_top)
      call MPI_BSEND(req%request, 1, MPI_TYPE_REQUEST, req%owner, TREE_COMM_TAG_REQUEST_KEY, &
        t%comm_env%comm, ierr)

      req%node%flags = ibset(req%node%flags, TREE_NODE_FLAG_REQUEST_SENT)

      send_request = .true.
    else
      send_request = .false.
    end if
  end function


  !>
  !> Insert incoming data into the tree.
  !>
  subroutine unpack_data(t, child_data, num_children, ipe_sender)
    use module_tree, only: t_tree, tree_insert_node, tree_lookup_node_critical
    use module_pepc_types, only: t_tree_node, t_tree_node_ptr, t_tree_node_package
    use module_tree_node
    use module_spacefilling, only: parent_key_from_key
    use module_atomic_ops, only: atomic_write_barrier
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    type(t_tree_node_package) :: child_data(num_children) !< child data that has been received
    integer :: num_children !< actual number of valid children in dataset
    integer, intent(in) :: ipe_sender

    type(t_tree_node), pointer :: parent_node
    type(t_tree_node), pointer :: prev_sibling
    type(t_tree_node) :: unpack_node
    integer :: ic

    DEBUG_ASSERT(num_children > 0)
    
    prev_sibling => null()

    call tree_lookup_node_critical(t, parent_key_from_key(child_data(1)%key), parent_node, &
        'TREE_COMMUNICATOR:unpack_data(): - get parent node')

    do ic = num_children,1,-1
      call tree_node_unpack(child_data(ic), unpack_node)
      unpack_node%first_child  => null()
      ! FIXME: originallz, we should have been using tree_node_connect_children() but
      ! for the sake of nonsense, bg_xlf produces a relieably segfaulting code
      ! if an array of t_tree_node_ptr is onvolved here :-(
      unpack_node%next_sibling => prev_sibling
      
      ! tree nodes coming from remote PEs are flagged for easier identification
      unpack_node%flags = ibset(unpack_node%flags, TREE_NODE_FLAG_HAS_REMOTE_CONTRIBUTIONS)

      if (tree_comm_debug) then
        DEBUG_INFO('("PE", I6, " received answer. parent_key=", O22, ",  sender=", I6, ",  owner=", I6, ",  kchild=", O22)', t%comm_env%rank, parent_node%key, ipe_sender, unpack_node%owner, unpack_node%key)
      end if

      if (.not. tree_insert_node(t, unpack_node, prev_sibling)) then
        DEBUG_WARNING_ALL(*, "Received a node that is already present.")
      end if

      ! count number of fetched nodes
      t%communicator%sum_fetches = t%communicator%sum_fetches+1
      
    end do

    parent_node%first_child => prev_sibling

    call atomic_write_barrier() ! make sure children are actually inserted before indicating their presence
    ! set 'children-here'-flag for all parent addresses
    ! may only be done *after inserting all* children, hence not(!) during the loop above
    parent_node%flags = ibset(parent_node%flags, TREE_NODE_FLAG_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node
  end subroutine unpack_data


  !>
  !> send all requests from our thread-safe list until we find an invalid one
  !>
  subroutine send_requests(t)
    use module_tree, only: t_tree
    use module_atomic_ops, only: atomic_load_int, atomic_read_barrier
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t

    real*8 :: tsend
    integer :: tmp_top

    tsend = MPI_WTIME()

    do while (t%communicator%req_queue_top .ne. atomic_load_int(t%communicator%req_queue_bottom))
      tmp_top = mod(t%communicator%req_queue_top, TREE_COMM_REQUEST_QUEUE_LENGTH) + 1

      ! first check whether the entry is actually valid	  
      if (t%communicator%req_queue(tmp_top)%owner >= 0) then
        t%communicator%req_queue_top = tmp_top
        call atomic_read_barrier() ! make sure that reads of parts of the queue entry occurr in the correct order

        if (tree_comm_debug) then
          DEBUG_INFO('("PE", I6, " sending request.      req_queue_top=", I5, ", request_key=", O22, ", request_owner=", I6)', t%comm_env%rank, t%communicator%req_queue_top, t%communicator%req_queue(t%communicator%req_queue_top)%request% key, t%communicator%req_queue(t%communicator%req_queue_top)%owner)
        end if

        if (send_request(t, t%communicator%req_queue(t%communicator%req_queue_top))) then
          t%communicator%request_balance = t%communicator%request_balance + 1
        end if

        ! we have to invalidate this request queue entry. this shows that we actually processed it and prevents it from accidentially being resent after the req_queue wrapped around
        t%communicator%req_queue(t%communicator%req_queue_top)%owner = -1
      else
        ! the next entry is not valid (obviously it has not been stored completely until now -> we abort here and try again later
        exit
      end if
    end do

    t%communicator%timings_comm(TREE_COMM_TIMING_SENDREQS) = t%communicator%timings_comm(TREE_COMM_TIMING_SENDREQS) &
      + ( MPI_WTIME() - tsend )
  end subroutine send_requests


  !>
  !> main routine of the communicator thread.
  !>
  !> Repeatedly checks the request queue for new requests to send out
  !> and calls MPI_IPROBE to check for requests to answer, then relinquishes the
  !> CPU.
  !>
  !> @todo Factor out thread scheduling code below and reactivate it.
  !>
  function run_communication_loop(arg) bind(c)
    use, intrinsic :: iso_c_binding
    use module_tree, only: t_tree
    use pthreads_stuff, only: pthreads_sched_yield, get_my_core
    use module_debug
    implicit none
    include 'mpif.h'

    type(c_ptr) :: run_communication_loop
    type(c_ptr), value :: arg

    type(t_tree), pointer :: t
    !integer, intent(in) :: max_particles_per_thread
    integer, dimension(mintag:maxtag) :: nummessages
    integer :: messages_per_iteration !< tracks current number of received and transmitted messages per commloop iteration for adjusting particles_per_yield
    integer :: ierr
    logical, allocatable :: comm_finished(:) ! will hold information on PE 0 about which processor
                                             ! is still communicating and which ones are finished
                                             ! to emulate a non-blocking barrier

    t => null()
    call c_f_pointer(arg, t)
    DEBUG_ASSERT(associated(t))
    
    ! store ID of comm-thread processor
    t%communicator%processor_id = get_my_core()

    nummessages            = 0
    messages_per_iteration = 0

    allocate(comm_finished(t%comm_env%size))
    comm_finished = .false.

    t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP) = MPI_WTIME()

    do while (.not. t%communicator%comm_thread_stopping)
      t%communicator%comm_loop_iterations(1) = t%communicator%comm_loop_iterations(1) + 1

      ! send our requested keys
      call send_requests(t)

      if (t%communicator%comm_thread_stop_requested) then
        t%communicator%comm_thread_stop_requested = .false.
        ! notify rank 0 that we are finished with our walk
        call MPI_BSEND(tree_comm_dummy, 1, MPI_INTEGER, 0, TREE_COMM_TAG_FINISHED_PE, &
          t%comm_env%comm, ierr)
      end if

      ! check whether we are still waiting for data or some other communication
      !if (walk_status == WALK_IAM_FINISHED) call check_comm_finished(t)

      ! process any incoming answers
      call run_communication_loop_inner(t, comm_finished, nummessages)

      messages_per_iteration = messages_per_iteration + sum(nummessages)
      t%communicator%request_balance = t%communicator%request_balance - nummessages(TREE_COMM_TAG_REQUESTED_DATA)
      nummessages(TREE_COMM_TAG_REQUESTED_DATA) = 0

      ! adjust the sched_yield()-timeout for the thread that shares its processor with the communicator
      !if (messages_per_iteration > MAX_MESSAGES_PER_ITERATION) then
      !  particles_per_yield = int(max(0.75 * particles_per_yield, 0.01*max_particles_per_thread))
      !  if (walk_debug) then
      !    DEBUG_INFO('("messages_per_iteration = ", I6, " > ", I6, " --> Decreased particles_per_yield to", I10)', messages_per_iteration, MAX_MESSAGES_PER_ITERATION, particles_per_yield)
      !  end if
      !else if ((particles_per_yield < max_particles_per_thread) .and. (messages_per_iteration < MIN_MESSAGES_PER_ITERATION)) then
      !  particles_per_yield = int(min(1.5 * particles_per_yield, 1. * max_particles_per_thread))
      !  if (walk_debug) then
      !    DEBUG_INFO('("messages_per_iteration = ", I6, " < ", I6, " --> Increased particles_per_yield to", I10)', messages_per_iteration, MIN_MESSAGES_PER_ITERATION, particles_per_yield)
      !  end if
      !end if

      ! currently, there is no further communication request --> other threads may do something interesting
      ierr = pthreads_sched_yield()
      if (.not. ierr == 0) then
        DEBUG_ERROR(*, "pthreads_sched_yield failed with status: ", ierr)
      end if

    end do ! while (.not. t%communicator%comm_thread_stopping)

    deallocate(comm_finished)

    t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP) = MPI_WTIME() - t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP)

    run_communication_loop = c_null_ptr
  end function run_communication_loop


  !>
  !> Checks for incoming MPI messages and acts on them.
  !>
  subroutine run_communication_loop_inner(t, comm_finished, nummessages)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_tree_node_package, MPI_TYPE_tree_node_package, t_request, MPI_TYPE_REQUEST
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    logical, intent(inout) :: comm_finished(:)
    integer, intent(inout), dimension(mintag:maxtag) :: nummessages

    logical :: msg_avail
    integer :: ierr
    integer :: stat(MPI_STATUS_SIZE)
    type(t_request) :: request
    type (t_tree_node_package) :: child_data(8) ! child data to be received
    integer :: num_children
    integer :: ipe_sender, msg_tag
    integer :: dummy
    real*8 :: tcomm

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
    call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, t%comm_env%comm, msg_avail, stat, ierr)

    if (msg_avail) then
      t%communicator%comm_loop_iterations(3) = t%communicator%comm_loop_iterations(3) + 1
      tcomm = MPI_WTIME()

      do while (msg_avail)
        ipe_sender = stat(MPI_SOURCE)
        msg_tag    = stat(MPI_TAG)

        ! the functions returns the number of received messages of any tag
        nummessages(msg_tag) = nummessages(msg_tag) + 1

        select case (msg_tag)
          ! another PE requested child data for a certain key
          case (TREE_COMM_TAG_REQUEST_KEY)
            ! actually receive this request...
            ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
            call MPI_RECV(request, 1, MPI_TYPE_request, ipe_sender, TREE_COMM_TAG_REQUEST_KEY, &
                    t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
            ! ... and answer it
            call answer_request(t, request, ipe_sender)

          ! some PE answered our request and sends
          case (TREE_COMM_TAG_REQUESTED_DATA)
            ! actually receive the data...
            ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
            call MPI_GET_COUNT(stat, MPI_TYPE_tree_node_package, num_children, ierr)
            call MPI_RECV(child_data, num_children, MPI_TYPE_tree_node_package, ipe_sender, TREE_COMM_TAG_REQUESTED_DATA, &
                    t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
            ! ... and put it into the tree and all other data structures
            call unpack_data(t, child_data, num_children, ipe_sender)

          ! rank 0 does bookkeeping about which PE is already finished with its walk
          ! no one else will ever receive this message tag
          case (TREE_COMM_TAG_FINISHED_PE)
            ! actually receive the data (however, we are not interested in it here)
            ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
            call MPI_RECV(dummy, 1, MPI_INTEGER, ipe_sender, TREE_COMM_TAG_FINISHED_PE, &
                    t%comm_env%comm, MPI_STATUS_IGNORE, ierr)

            DEBUG_ASSERT_MSG(t%comm_env%rank == 0, *, "this kind of message is only expected at rank 0!")
            if (tree_comm_debug) then
              DEBUG_INFO('("PE", I6, " has been told that PE", I6, " has finished walking")', t%comm_env%rank, ipe_sender)
              DEBUG_INFO(*, 'comm_finished = ', comm_finished)
              DEBUG_INFO('("nummessages(TAG_FINISHED_PE) = ", I6, ", count(comm_finished) = ", I6)', nummessages(TREE_COMM_TAG_FINISHED_PE), count(comm_finished))
            end if

            if (.not. comm_finished(ipe_sender + 1)) then
              comm_finished(ipe_sender + 1) = .true.

              if (all(comm_finished)) then
                call broadcast_comm_finished(t)
                if (tree_comm_debug) then
                  DEBUG_INFO(*, 'BCWF: comm_finished = ', comm_finished)
                  DEBUG_INFO('("BCWF: nummessages(TREE_COMM_TAG_FINISHED_PE) = ", I6, ", count(comm_finished) = ", I6)', nummessages(TREE_COMM_TAG_FINISHED_PE), count(comm_finished))
                end if
              end if
            else
              DEBUG_WARNING_ALL('("PE", I6, " has been told that PE", I6, " has finished walking, but already knew that. Obviously received duplicate TAG_FINISHED_PE, ignoring.")', t%comm_env%rank, ipe_sender)
            end if

          ! all PEs have finished their walk
          case (TREE_COMM_TAG_FINISHED_ALL)
            call MPI_RECV( dummy, 1, MPI_INTEGER, ipe_sender, TREE_COMM_TAG_FINISHED_ALL, &
                        t%comm_env%comm, MPI_STATUS_IGNORE, ierr) ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication

            if (tree_comm_debug) then
              DEBUG_INFO('("PE", I6, " has been told to terminate by PE", I6, " since all walks on all PEs are finished")', t%comm_env%rank, ipe_sender)
            end if

            t%communicator%comm_thread_stopping = .true.
        end select

        call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, t%comm_env%comm, msg_avail, stat, ierr)
      end do ! while (msg_avail)

      t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE) = t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE) +  (MPI_WTIME() - tcomm)
    end if ! msg_avail
  end subroutine run_communication_loop_inner
end module module_tree_communicator
