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
  use module_tree, only: t_tree_communicator, TREE_COMM_REQUEST_QUEUE_LENGTH, TREE_COMM_ANSWER_BUFF_LENGTH, &
    TREE_COMM_THREAD_STATUS_STOPPED, TREE_COMM_THREAD_STATUS_STARTING, TREE_COMM_THREAD_STATUS_STARTED, &
    TREE_COMM_THREAD_STATUS_WAITING
  use module_pepc_types
  implicit none
  private

  !> debug flags - cannot be modified at runtime due to performance reasons
  logical, public, parameter :: tree_comm_debug = .false.

  integer, private :: tree_comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)
  integer, private, save :: tree_comm_thread_counter = 0

  ! tags to be used in communication
  integer, private, parameter :: TREE_COMM_TAG_REQUEST_KEY         = 1257 !< message requesting child data for a certain key
  integer, private, parameter :: TREE_COMM_TAG_REQUEST_KEY_EAGER   = 1258 !< message requesting child and grnadchild data for a certain key using the eager send algorithm
  integer, private, parameter :: TREE_COMM_TAG_REQUESTED_DATA      = 1259 !< message that contains requested child data
  integer, private, parameter :: TREE_COMM_TAG_FINISHED_PE         = 1260 !< message to rank 0, to announce requested stop of communicator
  integer, private, parameter :: TREE_COMM_TAG_FINISHED_ALL        = 1261 !< message from rank 0, that communication has finished
  integer, private, parameter :: TREE_COMM_TAG_DUMP_TREE_AND_ABORT = 1262 !< message from some rank telling us, that we shall dump the tree and abort
  integer, private, parameter :: mintag = TREE_COMM_TAG_REQUEST_KEY
  integer, private, parameter :: maxtag = TREE_COMM_TAG_DUMP_TREE_AND_ABORT

  ! IDs for internal timing measurement
  integer, public, parameter :: TREE_COMM_TIMING_COMMLOOP = 1
  integer, public, parameter :: TREE_COMM_TIMING_RECEIVE  = 2

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

    integer(kind_default) :: msg_size_request, msg_size_data, buffsize, ierr

    if (.not. allocated(tree_comm_bsend_buffer)) then
      ! compute upper bounds for request and data message size
      call MPI_PACK_SIZE(1, MPI_KIND_KEY, mpi_comm_lpepc, msg_size_request, ierr)
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

    integer(kind_default) :: dummy, buffsize, ierr

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
    use pthreads_stuff, only: pthreads_createthread, pthreads_sched_yield, THREAD_TYPE_COMMUNICATOR
    use module_atomic_ops, only: atomic_load_int, atomic_store_int
    use module_debug
    use module_timings
    implicit none

    type(t_tree), target, intent(inout) :: t
    type(c_ptr) :: tp

    DEBUG_ASSERT(atomic_load_int(t%communicator%thread_status) == TREE_COMM_THREAD_STATUS_STOPPED)
    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " run_communication_loop start.")', t%comm_env%rank)
    end if
    ! TODO: in future, need to handle multiple communicators.
    call timer_reset(t_comm_total)
    call timer_reset(t_comm_recv)
    
    tp = c_loc(t)
    call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_STARTING)
    tree_comm_thread_counter = tree_comm_thread_counter + 1
    DEBUG_ERROR_ON_FAIL(pthreads_createthread(t%communicator%comm_thread, c_funloc(run_communication_loop), tp, thread_type = THREAD_TYPE_COMMUNICATOR, counter = tree_comm_thread_counter))
    
    ! we have to wait here until the communicator has really started to find out its processor id
    do while (atomic_load_int(t%communicator%thread_status) /= TREE_COMM_THREAD_STATUS_STARTED)
      DEBUG_ERROR_ON_FAIL(pthreads_sched_yield())
    end do
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
    use module_atomic_ops, only: atomic_load_int, atomic_store_int
    use module_debug
    use module_timings
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    integer(kind_Default) :: ierr

    DEBUG_ASSERT(atomic_load_int(t%communicator%thread_status) == TREE_COMM_THREAD_STATUS_STARTED \
            .or. atomic_load_int(t%communicator%thread_status) == TREE_COMM_THREAD_STATUS_STARTING)
            
    ! notify rank 0 that we are finished with our walk
    call MPI_BSEND(tree_comm_dummy, 1, MPI_INTEGER, 0, TREE_COMM_TAG_FINISHED_PE, &
      t%comm_env%comm, ierr)
    call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_WAITING)

    DEBUG_ERROR_ON_FAIL(pthreads_jointhread(t%communicator%comm_thread))
    tree_comm_thread_counter = tree_comm_thread_counter - 1

    call timer_add(t_comm_total,    t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP))
    call timer_add(t_comm_recv,     t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE))

    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " run_communication_loop end.")', t%comm_env%rank)
    end if
  end subroutine tree_communicator_stop


  !>
  !> Request children of node `n` in tree `t` from the responsible remote rank.
  !> The local node that actually needs the remote node is `n_targ`.
  !>
  !> This routine returns after sending the request. The caller is then free
  !> to continue working on something different. Later on, `tree_node_children_available`
  !> can be used to check whether the requested data has arrived in the
  !> meantime.
  !> Information from `particle` will be used for the eager send algorithm.
  !> If `particle` is omitted, no eager sending will be requested.
  !> `pos` can be used to explicitly override the particle position that
  !> is reported to the eager sending algorithm (eg for periodic boxes)
  !>
  subroutine tree_node_fetch_children(t, n, particle, pos)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_tree_node, t_particle
    use module_atomic_ops, only: atomic_write_barrier, &
      critical_section_enter, critical_section_leave
    use module_tree_node
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    type(t_tree_node), target, intent(inout) :: n
    type(t_particle), intent(in), optional :: particle
    real*8, intent(in), optional :: pos(3)
    integer(kind_default) :: ierr
    type(t_request_eager) :: req
    logical(kind_byte) :: dontsend
    
    
    ! check wether the node has already been requested -- then there is no need for waiting for the critical section
    if (.not. n%request_sent) then
      ! the following part has to be secured against synchronous invocation
      ! otherwise it will be possible that two walk threads can synchronously send a node request for identical nodes
      call critical_section_enter(t%communicator%cs_request)
        dontsend = n%request_sent
        n%request_sent = .true.
        call atomic_write_barrier()
      call critical_section_leave(t%communicator%cs_request)

      if (.not. dontsend) then
        ! send a request to PE n%owner
        ! telling, that we need child data for particle request_key(req_queue_top)
        if (present(particle)) then
          req%key = n%key
          req%particle = particle
          if (present(pos)) then 
            req%particle%x = pos
          end if

          call MPI_BSEND(req, 1, MPI_TYPE_request_eager, n%owner, TREE_COMM_TAG_REQUEST_KEY_EAGER, &
            t%comm_env%comm, ierr)
        else
          call MPI_BSEND(n%key, 1, MPI_KIND_KEY, n%owner, TREE_COMM_TAG_REQUEST_KEY, &
            t%comm_env%comm, ierr)
        endif
      end if
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
    integer(kind_pe) :: i
    integer(kind_default) :: ierr

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
  !> Notify all ranks about a serious error:
  !> they shall dump their tree and abort
  !>
  subroutine broadcast_dump_tree_and_abort(t)
    use module_tree, only: t_tree
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(in) :: t
    integer(kind_pe) :: i
    integer(kind_default) :: ierr

    ! all PEs have to be informed
    ! TODO: need better idea here...
    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " is telling all ranks to dump their tree and exit now")', t%comm_env%rank)
    end if

    do i = 0, t%comm_env%size - 1
      call MPI_BSEND(tree_comm_dummy, 1, MPI_INTEGER, i, TREE_COMM_TAG_DUMP_TREE_AND_ABORT, &
        t%comm_env%comm, ierr)
    end do
  end subroutine


  !>
  !> Notify all ranks about a serious error:
  !> they shall dump their tree and abort
  !>
  subroutine dump_tree_and_abort(t, sender)
    use module_tree, only: t_tree, tree_check, tree_dump
    use module_debug
    implicit none

    type(t_tree), intent(in) :: t
    integer(kind_pe), intent(in) :: sender
    logical :: ret

    DEBUG_INFO('("PE", I6, " has told us to dump our tree and exit now")', sender)
    ret = tree_check(t, 'dump_tree_and_abort()')
    call tree_dump(t)
    call debug_barrier() ! make sure that every rank completed its dump
    call debug_mpi_abort() ! usually, continuing would not make any sense at this stage

  end subroutine


  !>
  !> send node data
  !>
  subroutine send_data(t, nodes, numnodes, adressee) 
    use module_tree, only: t_tree 
    use module_pepc_types, only: t_tree_node, t_tree_node_package, MPI_TYPE_tree_node_package 
    implicit none 
    include 'mpif.h' 

    type(t_tree), intent(inout) :: t 
    type(t_tree_node_package), contiguous, intent(in) :: nodes(:) 
    integer, intent(in) :: numnodes 
    integer(kind_pe), intent(in) :: adressee 
    integer(kind_default) :: ierr 

    ! Ship child data back to PE that requested it 
    call MPI_BSEND(nodes, numnodes, MPI_TYPE_tree_node_package, & 
      adressee, TREE_COMM_TAG_REQUESTED_DATA, t%comm_env%comm, ierr) 

    ! statistics on number of sent children-packages 
    t%communicator%sum_ships = t%communicator%sum_ships + 1 
 
  end subroutine send_data 

  
  !>
  !> Simply collect all child nodes for node with `key`
  !>
  subroutine answer_request_simple(t, key, ipe_sender)
    use module_tree, only: t_tree, tree_lookup_node
    use module_pepc_types, only : t_tree_node, t_tree_node_package, kind_node
    use module_tree_node
    use module_debug
    implicit none
    type(t_tree), intent(inout) :: t
    integer(kind_key), intent(in) :: key
    integer(kind_pe), intent(in) :: ipe_sender
    
    type(t_tree_node_package) :: children_to_send(8) ! for an octtree, there will never be more than 8 direct children - no need for counting in advance
    integer :: nchild
    
    integer(kind_node) :: n
    
    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " answering request.                         key=", O22, ",        sender=", I6)', t%comm_env%rank, key, ipe_sender )
    end if

    if (.not. tree_lookup_node(t, key, n, 'WALK:answer_request_simple:parentkey')) then
      call broadcast_dump_tree_and_abort(t)
      return
    endif

    nchild = 0
    if (tree_node_get_first_child(t%nodes(n), n)) then
      do
        nchild = nchild + 1

        call tree_node_pack(t%nodes(n), children_to_send(nchild))

        if (.not. tree_node_get_next_sibling(t%nodes(n), n)) then; exit; end if
      end do
    
      call send_data(t, children_to_send, nchild, ipe_sender)
    else
      DEBUG_WARNING_ALL('("PE", I6, " got a request from PE", I6, " for key ", O8," but did not find any children_to_send. - Will not answer.")', t%comm_env%rank, ipe_sender, key)
    end if
    
  end subroutine answer_request_simple

  
  !>
  !> Collect all child and grandchild nodes for `parent` that
  !> are needed to fulfill the mac() from given position `pos`
  !>
  subroutine answer_request_eager(t, request, ipe_sender)
    use module_tree, only: t_tree, tree_lookup_node
    use module_pepc_types, only : t_tree_node, &
      t_tree_node_package, t_request_eager, kind_node
    use module_tree_node
    use module_debug
    implicit none
    type(t_tree), intent(inout) :: t
    type(t_request_eager), intent(in) :: request
    integer(kind_pe), intent(in) :: ipe_sender
    
    integer(kind_node) :: parent_idx, firstchild
    type(t_tree_node), pointer :: parent
    type(t_tree_node_package), allocatable :: children_to_send(:)
    integer :: nchild, i
   
    integer(kind_node), allocatable :: child_nodes(:)
    
    if (tree_comm_debug) then
      DEBUG_INFO('("PE", I6, " answering eager request.                   key=", O22, ",        sender=", I6)', t%comm_env%rank, request%key, ipe_sender )
    end if

    if (.not. tree_lookup_node(t, request%key, parent_idx, 'WALK:answer_request_eager:parentkey')) then
      call broadcast_dump_tree_and_abort(t)
      return
    endif
    parent => t%nodes(parent_idx)

    nchild = 0

    if (tree_node_get_first_child(parent, firstchild)) then
      ! first, we only collect pointers to all children that have to be sent and count them
    
      allocate(child_nodes(parent%descendants)) ! enough space to keep all descendants in this array
      
      call eager_collect_traverse(firstchild)

      ! collect real data from pointers in child_nodes
      if (nchild > 0) then
        allocate(children_to_send(nchild))
        
        do i=1,nchild
          call tree_node_pack(t%nodes(child_nodes(i)), children_to_send(i))
        end do
      
        call send_data(t, children_to_send, nchild, ipe_sender)
        
        deallocate(children_to_send)
      endif
    else
      DEBUG_WARNING_ALL('("PE", I6, " got a request from PE", I6, " for key ", O8," but did not find any children_to_send. - Will not answer.")', t%comm_env%rank, ipe_sender, request%key)
    endif

    deallocate(child_nodes)
    
    contains
    
      recursive subroutine eager_collect_traverse(node)
        use module_interaction_specific, only : mac
        use module_pepc_types, only: kind_node, t_tree_node
        implicit none
        integer(kind_node), intent(in) :: node
        
        integer(kind_node) :: n, s
        type(t_tree_node), pointer :: n_ptr
        real*8 :: dist2
        real*8 :: delta(3)
        
        n = node
        
        do
          nchild = nchild + 1
          child_nodes(nchild) = n
          
          n_ptr => t%nodes(n)
        
          delta = request%particle%x - n_ptr%interaction_data%coc  ! Separation vector
          dist2 = dot_product(delta, delta)

          ! check MAC
          if (.not. (mac(request%particle, n_ptr%interaction_data, dist2, t%boxlength2(n_ptr%level)) .or. tree_node_is_leaf(n_ptr))) then
            ! resolve the node
            if (tree_node_get_first_child(n_ptr, s)) then
              call eager_collect_traverse(s)
            end if
          else
            ! according to the MAC, the node does not have to be resolved any further (or it is a leaf) - there is nothing to do here
          end if
          
          ! traverse to next sibling
          if (.not. tree_node_get_next_sibling(n_ptr, n)) then; exit; end if
        end do          
      end subroutine
    
  end subroutine answer_request_eager
  

  !>
  !> Insert incoming data into the tree.
  !>
  subroutine unpack_data(t, child_data, num_children, ipe_sender)
    use module_tree, only: t_tree, tree_insert_node, tree_lookup_node_critical
    use module_pepc_types, only: t_tree_node, t_tree_node_package, kind_node
    use module_tree_node
    use module_spacefilling, only: parent_key_from_key
    use module_atomic_ops, only: atomic_write_barrier
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    type(t_tree_node_package) :: child_data(num_children) !< child data that has been received
    integer :: num_children !< actual number of valid children in dataset
    integer(kind_pe), intent(in) :: ipe_sender

    integer(kind_node) :: parent_node
    integer :: ic
    
    DEBUG_ASSERT(num_children > 0)
    
    call tree_lookup_node_critical(t, parent_key_from_key(child_data(1)%key), parent_node, &
        'TREE_COMMUNICATOR:unpack_data(): - get parent node')

    if (tree_node_children_available(t%nodes(parent_node))) then
      DEBUG_WARNING_ALL(*, 'Received some node data but parent flags indicate that respective children are already present. Ignoring these nodes.')
    else
      ic           = 1
      call unpack_children(parent_node, .true.)
    
      ! count number of fetched nodes
      t%communicator%sum_fetches = t%communicator%sum_fetches + num_children

      DEBUG_ASSERT(num_children == ic-1) ! otherwise, the received list of children was not in appropriate traversal order
    end if

    contains
    
    recursive subroutine unpack_children(parent_idx, toplvl)
      use module_pepc_types, only: kind_node
      use module_tree_node, only: NODE_INVALID
      use module_tree, only: tree_node_connect_children
      implicit none
      integer(kind_node), intent(in) :: parent_idx
      logical, intent(in) :: toplvl
      
      integer :: lvl
      type(t_tree_node), pointer :: parent
      integer(kind_node) :: newnode
      type(t_tree_node) :: unpack_node
      integer(kind_node) :: child_nodes(1:8) 
      integer :: nchild    
      
      if (ic > num_children) then; return; endif

      lvl     =  child_data(ic)%level
      nchild  =  0
      newnode =  NODE_INVALID ! in case of an algorithmic error, this invalid pointer should be catched
      parent  => t%nodes(parent_idx)
      
      do
        if (child_data(ic)%level < lvl) then
          ! all nodes on this level (and below) have been inserted
          exit
        else if (child_data(ic)%level > lvl) then
          call unpack_children(newnode, .false.)
        else
          call tree_node_unpack(child_data(ic), unpack_node)
          unpack_node%first_child = NODE_INVALID
          ! tree nodes coming from remote PEs are flagged for easier identification
          unpack_node%flags_local = ibset(unpack_node%flags_local, TREE_NODE_FLAG_LOCAL_HAS_REMOTE_CONTRIBUTIONS)
          
          if (.not. tree_insert_node(t, unpack_node, newnode)) then
            DEBUG_WARNING_ALL(*, 'Received a node that is already present.')
          end if
          
          ic     = ic     + 1
          nchild = nchild + 1
          child_nodes(nchild) = newnode

          if (tree_comm_debug) then
            DEBUG_INFO('("PE", I6, " received answer. parent_key=", O22, ",  sender=", I6, ",  owner=", I6, ",  kchild=", O22)', t%comm_env%rank, parent%key, ipe_sender, unpack_node%owner, unpack_node%key)
          end if
        endif
        
        if (ic > num_children) then; exit; endif
      end do
      
      call tree_node_connect_children(t, parent_idx, child_nodes(1:nchild))

      if (toplvl) then
        ! make sure children are actually inserted before indicating their presence
        ! in fact, this is only relevant for the topmost parent node as the others cannot be traversed 
        ! before its direct children are present
        call atomic_write_barrier()
      endif
      ! set 'children-here'-flag for all parent addresses
      ! may only be done *after inserting all* children, hence not(!) during the loop above
      parent%flags_local = ibset(parent%flags_local, TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node
      
    end subroutine
    
  end subroutine unpack_data


  !>
  !> main routine of the communicator thread.
  !>
  !> Repeatedly calls MPI_PROBE to check for requests to answer, then relinquishes the
  !> CPU.
  !>
  !> @todo Factor out thread scheduling code below and reactivate it.
  !>
  function run_communication_loop(arg) bind(c)
    use, intrinsic :: iso_c_binding
    use module_tree, only: t_tree
    use pthreads_stuff, only: pthreads_sched_yield, get_my_core, pthreads_exitthread
    use module_atomic_ops, only: atomic_load_int, atomic_store_int, atomic_write_barrier
    use module_debug
    implicit none
    include 'mpif.h'

    type(c_ptr) :: run_communication_loop
    type(c_ptr), value :: arg

    type(t_tree), pointer :: t
    !integer, intent(in) :: max_particles_per_thread
    logical, allocatable :: comm_finished(:) ! will hold information on PE 0 about which processor
                                             ! is still communicating and which ones are finished
                                             ! to emulate a non-blocking barrier

    t => null()
    call c_f_pointer(arg, t)
    DEBUG_ASSERT(associated(t))
    
    ! store ID of comm-thread processor
    t%communicator%processor_id = get_my_core()
    call atomic_write_barrier()

    ! signal successfull start
    call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_STARTED)

    allocate(comm_finished(t%comm_env%size))
    comm_finished = .false.

    t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP) = MPI_WTIME()

    do while (atomic_load_int(t%communicator%thread_status) /= TREE_COMM_THREAD_STATUS_STOPPED)
      t%communicator%comm_loop_iterations(1) = t%communicator%comm_loop_iterations(1) + 1

      ! process any incoming messages
      call run_communication_loop_inner(t, comm_finished)

      ! currently, there is no further communication request --> other threads may do something interesting
      DEBUG_ERROR_ON_FAIL(pthreads_sched_yield())

    end do ! while (.not. t%communicator%comm_thread_stopping)

    deallocate(comm_finished)

    t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP) = MPI_WTIME() - t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP)

    run_communication_loop = c_null_ptr
    DEBUG_ERROR_ON_FAIL(pthreads_exitthread())
  end function run_communication_loop


  !>
  !> Checks for incoming MPI messages and acts on them.
  !>
  subroutine run_communication_loop_inner(t, comm_finished)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_tree_node_package, MPI_TYPE_tree_node_package, t_request_eager, MPI_TYPE_request_eager
    use module_atomic_ops, only: atomic_store_int
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), intent(inout) :: t
    logical, intent(inout) :: comm_finished(:)

    integer(kind_default) :: ierr
    integer :: stat(MPI_STATUS_SIZE)
    type(t_request_eager) :: request
    integer(kind_key) :: req_key
    type (t_tree_node_package), allocatable :: child_data(:) ! child data to be received
    integer :: num_children
    integer(kind_pe) :: ipe_sender
    integer ::msg_tag
    real*8 :: tcomm

    ! probe for incoming messages
    call MPI_PROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, t%comm_env%comm, stat, ierr)

    t%communicator%comm_loop_iterations(3) = t%communicator%comm_loop_iterations(3) + 1
    tcomm = MPI_WTIME()

    ipe_sender = stat(MPI_SOURCE)
    msg_tag    = stat(MPI_TAG)

    select case (msg_tag)
      ! another PE requested child data for a certain key
      case (TREE_COMM_TAG_REQUEST_KEY)
        ! actually receive this request...
        call MPI_RECV(req_key, 1, MPI_KIND_KEY, ipe_sender, TREE_COMM_TAG_REQUEST_KEY, &
                t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
        ! ... and answer it
        call answer_request_simple(t, req_key, ipe_sender)
      ! another PE requested child data for a certain key and wants us to do some additional work
      case (TREE_COMM_TAG_REQUEST_KEY_EAGER)
        ! actually receive this request...
        call MPI_RECV(request, 1, MPI_TYPE_request_eager, ipe_sender, TREE_COMM_TAG_REQUEST_KEY_EAGER, &
                t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
        ! ... and answer it
        call answer_request_eager(t, request, ipe_sender)

      ! some PE answered our request and sends
      case (TREE_COMM_TAG_REQUESTED_DATA)
        ! actually receive the data...
        call MPI_GET_COUNT(stat, MPI_TYPE_tree_node_package, num_children, ierr)
        allocate(child_data(num_children))
        call MPI_RECV(child_data, num_children, MPI_TYPE_tree_node_package, ipe_sender, TREE_COMM_TAG_REQUESTED_DATA, &
                t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
        ! ... and put it into the tree and all other data structures
        call unpack_data(t, child_data, num_children, ipe_sender)
        deallocate(child_data)
            
      ! rank 0 does bookkeeping about which PE is already finished with its walk
      ! no one else will ever receive this message tag
      case (TREE_COMM_TAG_FINISHED_PE)
        ! actually receive the data (however, we are not interested in it here)
        ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
        call MPI_RECV(tree_comm_dummy, 1, MPI_INTEGER, ipe_sender, TREE_COMM_TAG_FINISHED_PE, &
                t%comm_env%comm, MPI_STATUS_IGNORE, ierr)

        DEBUG_ASSERT_MSG(t%comm_env%rank == 0, *, "this kind of message is only expected at rank 0!")
        if (tree_comm_debug) then
          DEBUG_INFO('("PE", I6, " has been told that PE", I6, " has finished walking")', t%comm_env%rank, ipe_sender)
          DEBUG_INFO(*, 'comm_finished = ', comm_finished)
        end if

        if (.not. comm_finished(ipe_sender + 1)) then
          comm_finished(ipe_sender + 1) = .true.

          if (all(comm_finished)) then
            call broadcast_comm_finished(t)
            if (tree_comm_debug) then
              DEBUG_INFO(*, 'BCWF: comm_finished = ', comm_finished)
            end if
          end if
        else
          DEBUG_WARNING_ALL('("PE", I6, " has been told that PE", I6, " has finished walking, but already knew that. Obviously received duplicate TAG_FINISHED_PE, ignoring.")', t%comm_env%rank, ipe_sender)
        end if

      ! all PEs have finished their walk
      case (TREE_COMM_TAG_FINISHED_ALL)
        call MPI_RECV(tree_comm_dummy, 1, MPI_INTEGER, ipe_sender, TREE_COMM_TAG_FINISHED_ALL, &
                    t%comm_env%comm, MPI_STATUS_IGNORE, ierr) ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication

        if (tree_comm_debug) then
          DEBUG_INFO('("PE", I6, " has been told to terminate by PE", I6, " since all walks on all PEs are finished")', t%comm_env%rank, ipe_sender)
        end if

        call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_STOPPED)
          
      ! on some rank something went wrong - we shall dump our tree  
      case (TREE_COMM_TAG_DUMP_TREE_AND_ABORT)
            
        call MPI_RECV(tree_comm_dummy, 1, MPI_INTEGER, ipe_sender, TREE_COMM_TAG_DUMP_TREE_AND_ABORT, &
                    t%comm_env%comm, MPI_STATUS_IGNORE, ierr) ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication

        if (tree_comm_debug) then
          DEBUG_INFO('("PE", I6, " has been told to dump its tree and abort by PE", I6)', t%comm_env%rank, ipe_sender)
        end if

        call dump_tree_and_abort(t, ipe_sender)

    end select

    t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE) = t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE) +  (MPI_WTIME() - tcomm)
  end subroutine run_communication_loop_inner
end module module_tree_communicator
