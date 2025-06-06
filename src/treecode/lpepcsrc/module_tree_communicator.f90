! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2019 Juelich Supercomputing Centre,
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
   use module_tree, only: t_request_queue_entry, t_tree_communicator, TREE_COMM_REQUEST_QUEUE_LENGTH, TREE_COMM_ANSWER_BUFF_LENGTH, &
                          TREE_COMM_THREAD_STATUS_STOPPED, TREE_COMM_THREAD_STATUS_STARTING, TREE_COMM_THREAD_STATUS_STARTED, &
                          TREE_COMM_THREAD_STATUS_STOPPING, TREE_COMM_THREAD_STATUS_WAITING
   use module_pepc_kinds
   use module_pepc_types
   implicit none
   private

   !> debug flags - cannot be modified at runtime due to performance reasons
   logical, public, parameter :: tree_comm_debug = .false.

   integer, private :: tree_comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)
   integer, private, save :: tree_comm_thread_counter = 0

   ! tags to be used in communication
   integer, private, parameter :: TREE_COMM_TAG_REQUEST_KEY = 1257 !< message requesting child data for a certain key
   integer, private, parameter :: TREE_COMM_TAG_REQUEST_KEY_EAGER = 1258 !< message requesting child and grnadchild data for a certain key using the eager send algorithm
   integer, private, parameter :: TREE_COMM_TAG_REQUESTED_DATA = 1259 !< message that contains requested child data
   integer, private, parameter :: TREE_COMM_TAG_FINISHED_PE = 1260 !< message to rank 0, to announce requested stop of communicator
   integer, private, parameter :: TREE_COMM_TAG_FINISHED_ALL = 1261 !< message from rank 0, that communication has finished
   integer, private, parameter :: TREE_COMM_TAG_DUMP_TREE_AND_ABORT = 1262 !< message from some rank telling us, that we shall dump the tree and abort
   integer, private, parameter :: mintag = TREE_COMM_TAG_REQUEST_KEY
   integer, private, parameter :: maxtag = TREE_COMM_TAG_DUMP_TREE_AND_ABORT

   ! IDs for internal timing measurement
   integer, public, parameter :: TREE_COMM_TIMING_COMMLOOP = 1
   integer, public, parameter :: TREE_COMM_TIMING_RECEIVE = 2
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
      use module_pepc_types, only: MPI_TYPE_tree_node_package_vec
      use mpi
      implicit none

      integer(kind_default) :: msg_size_request, msg_size_data, buffsize, ierr

      if (.not. allocated(tree_comm_bsend_buffer)) then
         ! compute upper bounds for request and data message size
         call MPI_PACK_SIZE(1, MPI_KIND_KEY, mpi_comm_lpepc, msg_size_request, ierr)
         msg_size_request = msg_size_request + MPI_BSEND_OVERHEAD

         call MPI_PACK_SIZE(8, MPI_TYPE_tree_node_package_vec, mpi_comm_lpepc, msg_size_data, ierr)
         msg_size_data = msg_size_data + MPI_BSEND_OVERHEAD

         buffsize = (TREE_COMM_REQUEST_QUEUE_LENGTH * msg_size_request + TREE_COMM_ANSWER_BUFF_LENGTH * msg_size_data)
         ! reserve memory for buffered mpi communication
         allocate (tree_comm_bsend_buffer(buffsize))
         ! and tell mpi , where it can be found
         call MPI_BUFFER_ATTACH(tree_comm_bsend_buffer, buffsize, ierr)
      end if
   end subroutine tree_communicator_prepare

   !>
   !> Deallocate MPI communication buffers in pepc_finalize.
   !>
   subroutine tree_communicator_finalize()
      use mpi
      implicit none

      integer(kind_default) :: dummy, buffsize, ierr

      if (allocated(tree_comm_bsend_buffer)) then
         ! free our buffer that was reserved for buffered communication
         call MPI_BUFFER_DETACH(dummy, buffsize, ierr)
         deallocate (tree_comm_bsend_buffer)
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

      DEBUG_ASSERT(atomic_load_int(t%communicator%thread_status) .eq. TREE_COMM_THREAD_STATUS_STOPPED)
      if (tree_comm_debug) then
         DEBUG_INFO('("PE", I6, " run_communication_loop start.")', t%comm_env%rank)
      end if
      ! TODO: in future, need to handle multiple communicators.
      call timer_reset(t_comm_total)
      call timer_reset(t_comm_recv)
      call timer_reset(t_comm_sendreqs)

      tp = c_loc(t)
      call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_STARTING)
      tree_comm_thread_counter = tree_comm_thread_counter + 1
      ERROR_ON_FAIL(pthreads_createthread(t%communicator%comm_thread, c_funloc(run_communication_loop), tp, thread_type=THREAD_TYPE_COMMUNICATOR, counter=tree_comm_thread_counter))

      ! we have to wait here until the communicator has really started to find out its processor id
      do while (atomic_load_int(t%communicator%thread_status) .ne. TREE_COMM_THREAD_STATUS_STARTED)
         ERROR_ON_FAIL(pthreads_sched_yield())
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
      use mpi
      implicit none

      type(t_tree), intent(inout) :: t

      ! prevent multiple calls to this function
      if ((atomic_load_int(t%communicator%thread_status) .eq. TREE_COMM_THREAD_STATUS_STARTED) &
          .or. (atomic_load_int(t%communicator%thread_status) .eq. TREE_COMM_THREAD_STATUS_STARTING)) then

         ! notify rank the communicator that we are finished with our walk
         call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_STOPPING)

         ERROR_ON_FAIL(pthreads_jointhread(t%communicator%comm_thread))
         tree_comm_thread_counter = tree_comm_thread_counter - 1

         call timer_add(t_comm_total, t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP))
         call timer_add(t_comm_recv, t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE))
         call timer_add(t_comm_sendreqs, t%communicator%timings_comm(TREE_COMM_TIMING_SENDREQS))

         if (tree_comm_debug) then
            DEBUG_INFO('("PE", I6, " run_communication_loop end.")', t%comm_env%rank)
         end if
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
   !> Information from `particle` will be used for the eager send algorithm.
   !> If `particle` is omitted, no eager sending will be requested.
   !> `pos` can be used to explicitly override the particle position that
   !> is reported to the eager sending algorithm (eg for periodic boxes)
   !>
   subroutine tree_node_fetch_children(t, n, nidx, particle, pos)
      use module_tree, only: t_tree
      use module_pepc_types, only: t_tree_node, t_particle
      use module_atomic_ops, only: atomic_mod_increment_and_fetch_int, &
                                   atomic_write_barrier, atomic_read_barrier, &
                                   atomic_load_int
      use module_debug
      use mpi
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_tree_node), target, intent(inout) :: n
      integer(kind_node), intent(in) :: nidx
      type(t_particle), intent(in), optional :: particle
      real(kind_physics), intent(in), optional :: pos(3)

      integer :: local_queue_bottom
!!!#define MPI_MULTIPLE   <-- use CPP flag instead of hard-coded choice
#ifdef MPI_MULTIPLE
      logical :: flag
      type(t_request_queue_entry) :: req
#endif

      ! check whether the node has already been requested
      ! this if-construct should be secured against synchronous invocation (together with the modification while receiving data)
      ! otherwise it will be possible that two walk threads can synchronously post a particle to the request queue
      ! (this _may_ be fixed upon receiving data, ignoring already present data)
#ifdef MPI_MULTIPLE
!$OMP atomic seq_cst read
      flag = n%request_posted
!$OMP end atomic
!$OMP flush
      if (flag) then
         return
      end if
#else
      if (n%request_posted) then
         return
      end if
#endif

      ! we first flag the particle as having been already requested to prevent other threads from doing it while
      ! we are inside this function
#ifdef __PGI
      n%request_posted = .true. ! Set requested flag
#else
!$OMP atomic seq_cst write
      n%request_posted = .true. ! Set requested flag
!$OMP end atomic
!$OMP flush
#endif

#ifdef MPI_MULTIPLE
      ! we issue the request right here right now from the worker thread instead of using the message queue_top
      req%request%node = n%first_child
      req%request%parent = nidx
      req%node => n

      if (present(particle)) then
         ! eager request
         req%eager_request = .true.
         req%request%particle = particle

         if (present(pos)) then
            req%request%particle%x = pos
         end if
      else
         ! simple request
         req%eager_request = .false.
      end if

      ! we do not have to do this here, keep it anyway...
      req%entry_valid = .true.

      if (send_request(req, t%comm_env)) then
         req%entry_valid = .false.
      else
         DEBUG_ERROR(*, "Issue with sending request, send_request returned error")
      end if
#else
      ! thread-safe way of reserving storage for our request
      local_queue_bottom = atomic_mod_increment_and_fetch_int(t%communicator%req_queue_bottom, TREE_COMM_REQUEST_QUEUE_LENGTH)

      if (local_queue_bottom .eq. atomic_load_int(t%communicator%req_queue_top)) then
         call atomic_read_barrier() ! ensure all information is compared correctly (queue_top and entry_valid)
         DEBUG_ERROR(*, "Issue with request sending queue: TREE_COMM_REQUEST_QUEUE_LENGTH is too small: ", TREE_COMM_REQUEST_QUEUE_LENGTH)
      end if

      ! the communicator will check validity of the request and will only proceed as soon as the entry is valid
      !     -- this actually only work in debugging more, serialisation may be safer
      DEBUG_ASSERT(.not. t%communicator%req_queue(local_queue_bottom)%entry_valid)

      t%communicator%req_queue(local_queue_bottom)%request%node = n%first_child
      t%communicator%req_queue(local_queue_bottom)%request%parent = nidx
      t%communicator%req_queue(local_queue_bottom)%node => n

      if (present(particle)) then
         ! eager request
         t%communicator%req_queue(local_queue_bottom)%eager_request = .true.
         t%communicator%req_queue(local_queue_bottom)%request%particle = particle

         if (present(pos)) then
            t%communicator%req_queue(local_queue_bottom)%request%particle%x = pos
         end if
      else
         ! simple request
         t%communicator%req_queue(local_queue_bottom)%eager_request = .false.
      end if

      call atomic_write_barrier() ! make sure the above information is actually written before flagging the entry valid by writing the owner
      t%communicator%req_queue(local_queue_bottom)%entry_valid = .true.
#endif

#ifdef MPI_MULTIPLE
      ! we double the function for the moment being, simpler than changing the module structure
   contains

      function send_request(req, comm_env)
         use module_tree_node
         use module_pepc_types, only: MPI_TYPE_request_eager_sca
         use module_comm_env, only: t_comm_env
         use mpi
         implicit none

         logical :: send_request
         type(t_request_queue_entry), volatile, intent(inout) :: req
         type(t_comm_env), intent(inout) :: comm_env
         integer(kind_node) :: req_simple(2)

         integer(kind_default) :: ierr

         if (.not. btest(req%node%flags_local, TREE_NODE_FLAG_LOCAL_REQUEST_SENT)) then
            ! send a request to PE req_queue_owners(req_queue_top)
            ! telling, that we need child data for particle request_key(req_queue_top)

            if (req%eager_request) then
               call MPI_BSEND(req%request, 1, MPI_TYPE_request_eager_sca, req%node%owner, TREE_COMM_TAG_REQUEST_KEY_EAGER, &
                              comm_env%comm, ierr)
            else
               req_simple(1) = req%request%node
               req_simple(2) = req%request%parent
               call MPI_BSEND(req_simple, 2, MPI_KIND_NODE, req%node%owner, TREE_COMM_TAG_REQUEST_KEY, &
                              comm_env%comm, ierr)
            end if

            req%node%flags_local = ibset(req%node%flags_local, TREE_NODE_FLAG_LOCAL_REQUEST_SENT)
            send_request = .true.
         else
            send_request = .false.
         end if
      end function
#endif
   end subroutine tree_node_fetch_children

   !>
   !> Notify all ranks about a communication shutdown.
   !>
   subroutine broadcast_comm_finished(t)
      use module_tree, only: t_tree
      use module_debug
      use mpi
      implicit none

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
      use mpi
      implicit none

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
      use module_pepc_types, only: t_tree_node, t_tree_node_package, MPI_TYPE_tree_node_package_vec, t_request_eager
      use module_debug
      use mpi
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_tree_node_package), contiguous, intent(in) :: nodes(:)
      integer, intent(in) :: numnodes
      integer(kind_pe), intent(in) :: adressee
      integer(kind_default) :: ierr

      !write(*,'(a,i0,a,i0)') 'snd:', numnodes, '->', adressee
      DEBUG_ASSERT(nodes(1)%parent .ne. 0)
      ! Ship child data back to PE that requested it
      call MPI_BSEND(nodes, numnodes, MPI_TYPE_tree_node_package_vec, &
                     adressee, TREE_COMM_TAG_REQUESTED_DATA, t%comm_env%comm, ierr)

      ! statistics on number of sent children-packages
      t%communicator%sum_ships = t%communicator%sum_ships + 1

   end subroutine send_data

   !>
   !> Simply collect node and all its siblings
   !>
   subroutine answer_request_simple(t, req, ipe_sender)
      use module_tree, only: t_tree
      use module_pepc_types, only: t_tree_node, t_tree_node_package
      use module_tree_node
      use module_debug
      implicit none
      type(t_tree), intent(inout) :: t
      integer(kind_node), intent(in) :: req(2)
      integer(kind_pe), intent(in) :: ipe_sender

      type(t_tree_node_package) :: children_to_send(8) ! for an octtree, there will never be more than 8 direct children - no need for counting in advance
      integer :: nchild

      integer(kind_node) :: n

      if (tree_comm_debug) then
         DEBUG_INFO('("PE", I6, " answering request.                         node=", I22, ",        sender=", I6)', t%comm_env%rank, req(1), ipe_sender)
      end if

      if (req(1) .eq. NODE_INVALID) then
         DEBUG_WARNING_ALL('("Received request with node == NODE_INVALID from pe ", I0, ", Its tree data might be damaged. Dumping all trees and aborting.")', ipe_sender)
         call broadcast_dump_tree_and_abort(t)
         return
      end if

      nchild = 0
      n = req(1)
      DEBUG_ASSERT(req(2) .ne. 0)

      do
         nchild = nchild + 1
         call tree_node_pack(t%nodes(n), children_to_send(nchild))
         ! we set the parent for the entries correctly for insertion on the receiver side. request%parent is a valid node index there.
         children_to_send(nchild)%parent = req(2)
         n = tree_node_get_next_sibling(t%nodes(n))
         if (n .eq. NODE_INVALID) exit
      end do

      call send_data(t, children_to_send, nchild, ipe_sender)

   end subroutine answer_request_simple

   !>
   !> Collect node, child and grandchild nodes for `node`
   !> that are needed to fulfill the mac() from given position `pos`
   !> same for `node`s siblings
   !>
   subroutine answer_request_eager(t, request, ipe_sender)
      use module_tree, only: t_tree
      use module_pepc_types, only: t_tree_node, &
                                   t_tree_node_package, t_request_eager
      use module_tree_node
      use module_debug
      implicit none
      type(t_tree), intent(inout) :: t
      type(t_request_eager), intent(in) :: request
      integer(kind_pe), intent(in) :: ipe_sender

      type(t_tree_node), pointer :: parent
      type(t_tree_node_package), allocatable :: children_to_send(:)
      integer :: nchild, i

      integer(kind_node), allocatable :: child_nodes(:)

      if (tree_comm_debug) then
         DEBUG_INFO('("PE", I6, " answering eager request.                   node=", O22, ",        sender=", I6)', t%comm_env%rank, request%node, ipe_sender)
      end if

      nchild = 0

      ! first, we only collect pointers to all nodes that have to be sent and count them - therefore we need the maximum number of child nodes to be expected
      parent => t%nodes(t%nodes(request%node)%parent)
      allocate (child_nodes(parent%descendants)) ! enough space to keep all descendants in this array

      call eager_collect_traverse(request%node)

      ! collect real data from pointers in child_nodes
      if (nchild .gt. 0) then
         allocate (children_to_send(nchild))

         do i = 1, nchild
            call tree_node_pack(t%nodes(child_nodes(i)), children_to_send(i))
         end do

         ! we set the parent for the first entry correctly for insertion on the receiver side. request%parent is a valid node index there.
         children_to_send(1)%parent = request%parent
         DEBUG_ASSERT(request%parent .ne. 0)

         call send_data(t, children_to_send, nchild, ipe_sender)

         deallocate (children_to_send)
      end if

      deallocate (child_nodes)

   contains

      recursive subroutine eager_collect_traverse(node)
         use module_interaction_specific, only: mac
         use module_pepc_types, only: t_tree_node
         implicit none
         integer(kind_node), intent(in) :: node

         integer(kind_node) :: n, s
         type(t_tree_node), pointer :: n_ptr
         real(kind_physics) :: dist2
         real(kind_physics) :: delta(3)

         n = node

         do
            nchild = nchild + 1
            child_nodes(nchild) = n

            n_ptr => t%nodes(n)

            delta = request%particle%x - n_ptr%interaction_data%coc  ! Separation vector
            dist2 = dot_product(delta, delta)

            ! check MAC
            if (.not. (mac(IF_MAC_NEEDS_PARTICLE(request%particle) n_ptr%interaction_data, dist2, t%boxlength2(n_ptr%level)) .or. tree_node_is_leaf(n_ptr))) then
               ! resolve the node
               s = tree_node_get_first_child(n_ptr)
               if (s .ne. NODE_INVALID) call eager_collect_traverse(s)
            else
               ! according to the MAC, the node does not have to be resolved any further (or it is a leaf) - there is nothing to do here
            end if

            ! traverse to next sibling
            n = tree_node_get_next_sibling(n_ptr)
            if (n .eq. NODE_INVALID) exit
         end do
      end subroutine

   end subroutine answer_request_eager

   !>
   !> Insert incoming data into the tree.
   !>
   subroutine unpack_data(t, child_data, num_children, ipe_sender)
      use module_tree, only: t_tree
      use module_pepc_types, only: t_tree_node, t_tree_node_package
      use module_tree_node
      use module_atomic_ops, only: atomic_write_barrier
      use module_debug
      use mpi
      implicit none

      type(t_tree), intent(inout) :: t
      type(t_tree_node_package) :: child_data(num_children) !< child data that has been received
      integer :: num_children !< actual number of valid children in dataset
      integer(kind_pe), intent(in) :: ipe_sender

      integer(kind_node) :: parent_node
      integer :: ic

      DEBUG_ASSERT(num_children .gt. 0)

      parent_node = child_data(1)%parent
      DEBUG_ASSERT(parent_node .ne. 0)

      if (tree_node_children_available(t%nodes(parent_node))) then
         DEBUG_WARNING_ALL(*, 'Received some node data but parent flags indicate that respective children are already present. Ignoring these nodes.')
      else
         ic = 1
         call unpack_children(parent_node, .true.)

         ! count number of fetched nodes
         t%communicator%sum_fetches = t%communicator%sum_fetches + num_children

         DEBUG_ASSERT(num_children .eq. ic - 1) ! otherwise, the received list of children was not in appropriate traversal order
      end if

   contains

      recursive subroutine unpack_children(parent_idx, toplvl)
         use module_tree_node, only: NODE_INVALID
         use module_tree, only: tree_node_connect_children, tree_provision_node, tree_count_node
         implicit none
         integer(kind_node), intent(in) :: parent_idx
         logical, intent(in) :: toplvl

         integer :: lvl
         type(t_tree_node), pointer :: parent
         integer(kind_node) :: newnode
         integer(kind_node) :: child_nodes(1:8)
         integer :: nchild

         if (ic .gt. num_children) return

         lvl = child_data(ic)%level
         nchild = 0
         newnode = NODE_INVALID ! in case of an algorithmic error, this invalid pointer should be catched
         parent => t%nodes(parent_idx)

         do
            if (child_data(ic)%level .lt. lvl) then
               ! all nodes on this level (and below) have been inserted
               exit
            else if (child_data(ic)%level .gt. lvl) then
               call unpack_children(newnode, .false.)
            else
               newnode = tree_provision_node(t)
               call tree_node_unpack(child_data(ic), t%nodes(newnode))
               ! tree nodes coming from remote PEs are flagged for easier identification
               t%nodes(newnode)%flags_local = ibset(t%nodes(newnode)%flags_local, TREE_NODE_FLAG_LOCAL_HAS_REMOTE_CONTRIBUTIONS)
               call tree_count_node(t, t%nodes(newnode))

               ic = ic + 1
               nchild = nchild + 1
               child_nodes(nchild) = newnode

               if (tree_comm_debug) then
                  DEBUG_INFO('("PE", I6, " received answer. parent_key=", O22, ",  sender=", I6, ",  owner=", I6, ",  kchild=", O22)', t%comm_env%rank, parent%key, ipe_sender, t%nodes(newnode)%owner, t%nodes(newnode)%key)
               end if
            end if

            if (ic .gt. num_children) then; exit; end if
         end do

         call tree_node_connect_children(t, parent_idx, child_nodes(1:nchild))

         if (toplvl) then
            ! make sure children are actually inserted before indicating their presence
            ! in fact, this is only relevant for the topmost parent node as the others cannot be traversed
            ! before its direct children are present
            call atomic_write_barrier()
         end if
         ! set 'children-here'-flag for all parent addresses
         ! may only be done *after inserting all* children, hence not(!) during the loop above
         parent%flags_local = ibset(parent%flags_local, TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node

      end subroutine

   end subroutine unpack_data

   !>
   !> send all requests from our thread-safe list until we find an invalid one
   !>
   subroutine send_requests(q, b, t, rb, tsend, comm_env)
      use module_tree, only: t_request_queue_entry
      use module_comm_env, only: t_comm_env
      use module_atomic_ops, only: t_atomic_int, atomic_load_int, atomic_store_int, atomic_read_barrier, atomic_write_barrier
      use module_debug
      use mpi
      implicit none

      type(t_request_queue_entry), volatile, intent(inout) :: q(TREE_COMM_REQUEST_QUEUE_LENGTH) !< request queue
      type(t_atomic_int), intent(inout) :: b !< queue bottom
      type(t_atomic_int), intent(inout) :: t !< queue top
      integer, intent(inout) :: rb !< request balance
      real*8, intent(inout) :: tsend !< time measurement
      type(t_comm_env), intent(inout) :: comm_env !< communication environment for sending

      integer :: tmp_top

      tsend = tsend - MPI_WTIME()

      do while (atomic_load_int(t) .ne. atomic_load_int(b))
         tmp_top = mod(atomic_load_int(t), TREE_COMM_REQUEST_QUEUE_LENGTH) + 1

         ! first check whether the entry is actually valid
         if (q(tmp_top)%entry_valid) then
            call atomic_read_barrier() ! make sure that reads of parts of the queue entry occurr in the correct order

            if (tree_comm_debug) then
               DEBUG_INFO('("PE", I6, " sending request.      req_queue_top=", I5, ", request_key=", O22, ", request_owner=", I6)', comm_env%rank, tmp_top, q(tmp_top)%node%key, q(tmp_top)%node%owner)
            end if

            if (send_request(q(tmp_top), comm_env)) then
               rb = rb + 1
            end if

            ! we have to invalidate this request queue entry. this shows that we actually processed it and prevents it from accidentially being resent after the req_queue wrapped around
            q(tmp_top)%entry_valid = .false.
            call atomic_write_barrier() ! try to ensure entry_valid is stored along with 't' and known to all threads
            call atomic_store_int(t, tmp_top)
         else
            ! the next entry is not valid (obviously it has not been stored completely until now -> we abort here and try again later
            exit
         end if
      end do

      tsend = tsend + MPI_WTIME()

   contains

      function send_request(req, comm_env)
         use module_tree_node
         use module_pepc_types, only: MPI_TYPE_request_eager_sca
         use mpi
         implicit none

         logical :: send_request
         type(t_request_queue_entry), volatile, intent(inout) :: req
         type(t_comm_env), intent(inout) :: comm_env
         integer(kind_node) :: req_simple(2)

         integer(kind_default) :: ierr

         if (.not. btest(req%node%flags_local, TREE_NODE_FLAG_LOCAL_REQUEST_SENT)) then
            ! send a request to PE req_queue_owners(req_queue_top)
            ! telling, that we need child data for particle request_key(req_queue_top)

            if (req%eager_request) then
               call MPI_BSEND(req%request, 1, MPI_TYPE_request_eager_sca, req%node%owner, TREE_COMM_TAG_REQUEST_KEY_EAGER, &
                              comm_env%comm, ierr)
            else
               req_simple(1) = req%request%node
               req_simple(2) = req%request%parent
               call MPI_BSEND(req_simple, 2, MPI_KIND_NODE, req%node%owner, TREE_COMM_TAG_REQUEST_KEY, &
                              comm_env%comm, ierr)
            end if

            req%node%flags_local = ibset(req%node%flags_local, TREE_NODE_FLAG_LOCAL_REQUEST_SENT)
            send_request = .true.
         else
            send_request = .false.
         end if
      end function
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
      use pthreads_stuff, only: pthreads_sched_yield, get_my_core, pthreads_exitthread
      use module_atomic_ops, only: atomic_load_int, atomic_store_int, atomic_write_barrier
      use module_debug
      use mpi
      implicit none

      type(c_ptr) :: run_communication_loop
      type(c_ptr), value :: arg

      type(t_tree), pointer :: t
      !integer, intent(in) :: max_particles_per_thread
      integer, dimension(mintag:maxtag) :: nummessages
      integer :: messages_per_iteration !< tracks current number of received and transmitted messages per commloop iteration for adjusting particles_per_yield
      integer(kind_default) :: ierr
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

      nummessages = 0
      messages_per_iteration = 0

      allocate (comm_finished(t%comm_env%size))
      comm_finished = .false.

      t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP) = MPI_WTIME()

      do while (atomic_load_int(t%communicator%thread_status) .ne. TREE_COMM_THREAD_STATUS_STOPPED)
         t%communicator%comm_loop_iterations(1) = t%communicator%comm_loop_iterations(1) + 1

         ! send our requested keys
         call send_requests(t%communicator%req_queue, &
                            t%communicator%req_queue_bottom, &
                            t%communicator%req_queue_top, &
                            t%communicator%request_balance, &
                            t%communicator%timings_comm(TREE_COMM_TIMING_SENDREQS), &
                            t%comm_env)

         if (atomic_load_int(t%communicator%thread_status) .eq. TREE_COMM_THREAD_STATUS_STOPPING) then
            ! notify rank 0 that we are finished with our walk
            call MPI_BSEND(tree_comm_dummy, 1, MPI_INTEGER, 0, TREE_COMM_TAG_FINISHED_PE, &
                           t%comm_env%comm, ierr)
            call atomic_store_int(t%communicator%thread_status, TREE_COMM_THREAD_STATUS_WAITING)
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
         !else if ((particles_per_yield < max_particles_per_thread) .and. (messages_per_iteration < MIN_MESSAGES_PER_ITERATION)) then
         !  particles_per_yield = int(min(1.5 * particles_per_yield, 1. * max_particles_per_thread))
         !end if

         ! currently, there is no further communication request --> other threads may do something interesting
         ERROR_ON_FAIL(pthreads_sched_yield())

      end do ! while (.not. t%communicator%comm_thread_stopping)

      deallocate (comm_finished)

      t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP) = MPI_WTIME() - t%communicator%timings_comm(TREE_COMM_TIMING_COMMLOOP)

      run_communication_loop = c_null_ptr
      ERROR_ON_FAIL(pthreads_exitthread())
   end function run_communication_loop

   !>
   !> Checks for incoming MPI messages and acts on them.
   !>
   subroutine run_communication_loop_inner(t, comm_finished, nummessages)
      use module_tree, only: t_tree
      use module_pepc_types, only: t_tree_node_package, MPI_TYPE_tree_node_package_vec, t_request_eager, MPI_TYPE_request_eager_sca
      use module_atomic_ops, only: atomic_store_int
      use module_debug
      use mpi
      implicit none

      type(t_tree), intent(inout) :: t
      logical, intent(inout) :: comm_finished(:)
      integer, intent(inout), dimension(mintag:maxtag) :: nummessages

      integer(kind_default) :: ierr
      integer :: stat(MPI_STATUS_SIZE)
      type(t_request_eager) :: request
      integer(kind_node) :: req_simple(2)
      type(t_tree_node_package), allocatable :: child_data(:) ! child data to be received
      integer :: num_children
      integer(kind_pe) :: ipe_sender
      integer ::msg_tag
      real*8 :: tcomm
      logical :: msg_avail

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
            msg_tag = stat(MPI_TAG)

            ! the functions returns the number of received messages of any tag
            nummessages(msg_tag) = nummessages(msg_tag) + 1

            select case (msg_tag)
               ! another PE requested child data for a certain key
            case (TREE_COMM_TAG_REQUEST_KEY)
               ! actually receive this request...
               ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
               call MPI_RECV(req_simple, 2, MPI_KIND_NODE, ipe_sender, TREE_COMM_TAG_REQUEST_KEY, &
                             t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
               ! ... and answer it
               call answer_request_simple(t, req_simple, ipe_sender)

            case (TREE_COMM_TAG_REQUEST_KEY_EAGER)
               ! actually receive this request...
               ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
               call MPI_RECV(request, 1, MPI_TYPE_request_eager_sca, ipe_sender, TREE_COMM_TAG_REQUEST_KEY_EAGER, &
                             t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
               ! ... and answer it
               call answer_request_eager(t, request, ipe_sender)

               ! some PE answered our request and sends
            case (TREE_COMM_TAG_REQUESTED_DATA)
               ! actually receive the data...
               ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
               call MPI_GET_COUNT(stat, MPI_TYPE_tree_node_package_vec, num_children, ierr)
               allocate (child_data(num_children), stat=ierr)
               !write(*,'(a,i0,a,i0)') 'recv:', num_children, '     <-', ipe_sender
               DEBUG_ASSERT(ierr .eq. 0)
               DEBUG_ASSERT(num_children .gt. 0)
               child_data(:)%parent = -666 ! if a quick glance at the code is
               ! correct, the first value will be meaningful after the RECV, all
               ! others should be -1?
               call MPI_RECV(child_data, num_children, MPI_TYPE_tree_node_package_vec, ipe_sender, TREE_COMM_TAG_REQUESTED_DATA, &
                             t%comm_env%comm, MPI_STATUS_IGNORE, ierr)
               ! ... and put it into the tree and all other data structures
               !if (any(child_data(:)%parent.eq.0)) write(*,*) child_data(:)%parent
               !if(child_data(1)%parent.eq.0) write(*,*) child_data(1)
               DEBUG_ASSERT(all(child_data(:)%parent .gt. -2))
               call unpack_data(t, child_data, num_children, ipe_sender)
               deallocate (child_data)

               ! rank 0 does bookkeeping about which PE is already finished with its walk
               ! no one else will ever receive this message tag
            case (TREE_COMM_TAG_FINISHED_PE)
               ! actually receive the data (however, we are not interested in it here)
               ! TODO: use MPI_RECV_INIT(), MPI_START() and colleagues for faster communication
               call MPI_RECV(tree_comm_dummy, 1, MPI_INTEGER, ipe_sender, TREE_COMM_TAG_FINISHED_PE, &
                             t%comm_env%comm, MPI_STATUS_IGNORE, ierr)

               DEBUG_ASSERT_MSG(t%comm_env%rank .eq. 0, *, "this kind of message is only expected at rank 0!")
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

            call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, t%comm_env%comm, msg_avail, stat, ierr)
         end do ! while (msg_avail)

         t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE) = t%communicator%timings_comm(TREE_COMM_TIMING_RECEIVE) + (MPI_WTIME() - tcomm)
      end if ! msg_avail
   end subroutine run_communication_loop_inner
end module module_tree_communicator
