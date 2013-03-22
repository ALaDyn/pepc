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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates helper functions that simplify communication during tree traversal
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_walk_communicator
  use module_pepc_types, only: t_tree_node
  implicit none

  private

    integer, public :: walk_status
    integer*8, public :: comm_loop_iterations(3) !< number of comm loop iterations (total, sending, receiving)

    !> debug flags - cannot be modified at runtime due to performance reasons
    logical, public, parameter :: walk_comm_debug = .false.

    ! tags to be used in communication
    integer, public, parameter :: TAG_REQUEST_KEY    = 1257 !< message tag for walk communication: message for requesting child data for a certain key
    integer, public, parameter :: TAG_REQUESTED_DATA = 1258 !< message tag for walk communication: message that contains requested child data
    integer, public, parameter :: TAG_FINISHED_PE    = 1259 !< message tag for walk communication: message to rank 0, that the sender has finished with walking, i.e. its walk status == WALK_IAM_FINISHED
    integer, public, parameter :: TAG_FINISHED_ALL   = 1260 !< message tag for walk communication: message from rank 0, that all PEs are finished with walking, i.e. we can set walk_status := WALK_ALL_FINISHED
    integer, public, parameter :: mintag = TAG_REQUEST_KEY
    integer, public, parameter :: maxtag = TAG_FINISHED_ALL

    ! internal status
    integer, public, parameter :: WALK_STILL_RUNNING = 0 !< value of walk_status: there are still particles left to be processed
    integer, public, parameter :: WALK_IAM_FINISHED  = 1 !< value of walk_status: all particles on this PE have finished their walk
    integer, public, parameter :: WALK_ALL_MSG_DONE  = 2 !< value of walk_status: there are no more pending requests or answers in the message queue
    integer, public, parameter :: WALK_I_NOTIFIED_0  = 3 !< value of walk_status: rank 0 has been informed that this PE is finished
    integer, public, parameter :: WALK_ALL_FINISHED  = 4 !< value of walk_status: rank 0 told this pe, that every PE is finished with walking, i.e. no more communication is necessary

    ! internal communication variables - not to be touched from outside the module
    integer*1, private, allocatable, target :: bsend_buffer(:) !< buffer for bsend-alls
    integer, private :: comm_dummy = 123456 !< dummy variable for sending "empty" messages (those, where we are only interested in the tag)
    real*8, public :: timings_comm(3) !< array for storing internal timing information

    ! IDs for internal timing measurement
    integer, public, parameter :: TIMING_COMMLOOP = 1
    integer, public, parameter :: TIMING_RECEIVE  = 2
    integer, public, parameter :: TIMING_SENDREQS = 3

    !> data type for internal request queue
    type, public :: t_request_queue_entry
      integer*8 :: key
      integer   :: owner
      type(t_tree_node), pointer :: node
    end type


    public init_comm_data
    public uninit_comm_data
    public notify_walk_finished
    public send_request
    public send_data
    public unpack_data
    public send_walk_finished
    public broadcast_walk_finished

  contains


      subroutine init_comm_data(t, REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH)
        use module_tree, only: t_tree
        use module_pepc_types, only: mpi_type_tree_node
        implicit none
        include 'mpif.h'

        type(t_tree), intent(in) :: t
        integer, intent(in) :: REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH

        integer :: msg_size_request, msg_size_data
        integer :: buffsize !< size of bsend buffer in bytes
        integer :: ierr

        ! compute upper bounds for request and data message size
        call MPI_PACK_SIZE(1, MPI_INTEGER8, t%comm_env%comm, msg_size_request, ierr)
        msg_size_request = msg_size_request + MPI_BSEND_OVERHEAD

        call MPI_PACK_SIZE(8, MPI_TYPE_tree_node, t%comm_env%comm, msg_size_data, ierr)
        msg_size_data = msg_size_data + MPI_BSEND_OVERHEAD

        buffsize = (REQUEST_QUEUE_LENGTH * msg_size_request + ANSWER_BUFF_LENGTH * msg_size_data)
        ! reserve memory for buffered mpi communication
        allocate(bsend_buffer(buffsize))
        ! and tell mpi , where it can be found
        call MPI_BUFFER_ATTACH(bsend_buffer, buffsize, ierr)

        walk_status = WALK_STILL_RUNNING

        timings_comm = 0.

    end subroutine init_comm_data



      subroutine uninit_comm_data()
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer :: buffsize
        integer :: dummy

        ! free our buffer that was reserved for buffered communication
        call MPI_BUFFER_DETACH(dummy, buffsize, ierr) ! FIXME: what is the dummy thought for?
        deallocate(bsend_buffer)
      end subroutine uninit_comm_data


      subroutine notify_walk_finished()
        use treevars
        implicit none

        walk_status = max(walk_status, WALK_IAM_FINISHED)
      end subroutine notify_walk_finished


      subroutine send_walk_finished(t)
        use module_tree, only: t_tree
        implicit none
        include 'mpif.h'

        type(t_tree), intent(in) :: t

        integer :: ierr

        ! notify rank 0 that we are finished with our walk
        call MPI_BSEND(comm_dummy, 1, MPI_INTEGER, 0, TAG_FINISHED_PE, &
          t%comm_env%comm, ierr)

        walk_status = WALK_I_NOTIFIED_0
      end subroutine send_walk_finished


      subroutine broadcast_walk_finished(t)
        use module_tree, only: t_tree
        use module_debug
        implicit none
        include 'mpif.h'

        type(t_tree), intent(in) :: t
        integer :: i, ierr

         ! all PEs have to be informed
         ! TODO: need better idea here...
         if (walk_comm_debug) then
           DEBUG_INFO('("PE", I6, " has found out that all PEs have finished walking - telling them to exit now")', t%comm_env%rank)
         end if

         do i = 0, t%comm_env%size - 1
           call MPI_BSEND(comm_dummy, 1, MPI_INTEGER, i, TAG_FINISHED_ALL, &
             t%comm_env%comm, ierr)
         end do
      end subroutine


    subroutine send_data(t, requested_key, ipe_sender)
      use module_tree, only: t_tree, tree_lookup_node_critical
      use module_pepc_types, only: t_tree_node, MPI_TYPE_tree_node
      use module_debug
      use module_tree_node
      implicit none
      include 'mpif.h'

      type(t_tree), intent(inout) :: t
      integer*8, intent(in) :: requested_key
      integer, intent(in) :: ipe_sender

      type(t_tree_node), target :: children_to_send(8)
      type(t_tree_node), pointer :: n
      integer*8, dimension(8) :: key_child
      integer :: j, ic, ierr, nchild

      if (walk_comm_debug) then
        DEBUG_INFO('("PE", I6, " answering request.                         request_key=", O22, ",        sender=", I6)',
                       t%comm_env%rank, requested_key, ipe_sender )
      end if

      j = 0
      call tree_lookup_node_critical(t, requested_key, n, 'WALK:send_data:parentkey')
      call tree_node_get_childkeys(n, nchild, key_child)

      do ic = 1,nchild
        call tree_lookup_node_critical(t, key_child(ic), n, 'WALK:send_data:childkey')
        children_to_send(ic) = n
        children_to_send(ic)%flags = int(iand( children_to_send(ic)%flags, TREE_NODE_CHILDBYTE ))! Catch lowest 8 bits of childbyte - filter off requested and here flags
      end do
      
      ! Ship child data back to PE that requested it
      call MPI_BSEND( children_to_send(1:nchild), nchild, MPI_TYPE_tree_node, &
        ipe_sender, TAG_REQUESTED_DATA, t%comm_env%comm, ierr)

      ! statistics on number of sent children-packages
      t%sum_ships = t%sum_ships + 1

    end subroutine send_data


    function send_request(t, req)
      use module_tree, only: t_tree
      use module_tree_node
      implicit none
      include 'mpif.h'
      
      logical :: send_request
      type(t_tree), intent(inout) :: t
      type(t_request_queue_entry), intent(in) :: req

      integer :: ierr

      if (.not. btest( req%node%flags, TREE_NODE_FLAG_REQUEST_SENT ) ) then
        ! send a request to PE req_queue_owners(req_queue_top)
        ! telling, that we need child data for particle request_key(req_queue_top)
        call MPI_BSEND(req%key, 1, MPI_INTEGER8, req%owner, TAG_REQUEST_KEY, &
          t%comm_env%comm, ierr)

        req%node%flags = ibset(req%node%flags, TREE_NODE_FLAG_REQUEST_SENT )

        send_request = .true.
      else
        send_request = .false.
      end if

    end function


    subroutine unpack_data(t, child_data, num_children, ipe_sender)
      use module_tree, only: t_tree, tree_insert_or_update_node, tree_lookup_node_critical
      use module_tree_node
      use module_spacefilling
      use module_debug
      use module_atomic_ops
      implicit none
      include 'mpif.h'

      type(t_tree), intent(inout) :: t
      type(t_tree_node) :: child_data(num_children) !< child data that has been received
      integer :: num_children !< actual number of valid children in dataset
      integer, intent(in) :: ipe_sender

      type(t_tree_node), pointer :: parent
      integer*8 :: parent_key(0:num_children)
      integer :: num_parents
      integer :: ic

      num_parents = 0
      parent_key(0) = 0_8 ! TODO: use named constant

      do ic = 1, num_children
        ! save parent key - after (!) inserting all (!) children we can flag it: it`s children are then accessible
        parent_key(num_parents + 1) = parent_key_from_key( child_data(ic)%key )
        !parent_addr(num_parents + 1) = key2addr( kparent, 'WALK:unpack_data() - get parent address' )
        if (parent_key(num_parents) .ne. parent_key(num_parents + 1)) then
          num_parents = num_parents + 1
        end if

        if (walk_comm_debug) then
          DEBUG_INFO('("PE", I6, " received answer.                            parent_key=", O22, ",        sender=", I6, ",        owner=", I6, ", kchild=", O22)',
                         t%comm_env%rank, parent_key(num_parents), ipe_sender, child_data(ic)%owner, child_data(ic)%key)
        end if

        ! tree nodes coming from remote PEs are flagged for easier identification
        child_data(ic)%flags = ibset(child_data(ic)%flags, TREE_NODE_FLAG_HAS_REMOTE_CONTRIBUTIONS)

        ! TODO: print message if node exists, i.e., use tree_insert_node()
        call tree_insert_or_update_node(t, child_data(ic))
        ! count number of fetched nodes
        t%sum_fetches = t%sum_fetches+1
     end do

     call atomic_write_barrier()

     ! set 'children-here'-flag for all parent addresses
     ! may only be done *after inserting all* children, hence not(!) during the loop above
     do ic=1,num_parents
         call tree_lookup_node_critical(t, parent_key(ic), parent, 'WALK:unpack_data() - get parent node')
         parent%flags = ibset(parent%flags, TREE_NODE_FLAG_CHILDREN_AVAILABLE) ! Set children_HERE flag for parent node
     end do


    end subroutine unpack_data

end module module_walk_communicator
