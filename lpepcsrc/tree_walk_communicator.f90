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
      integer   :: addr
      integer   :: owner
    end type


    public init_comm_data
    public run_communication_loop_inner
    public uninit_comm_data
    public notify_walk_finished
    public send_request

  contains


      ! initializes bsend buffer and rwlock objects
      ! returns size of bsend buffer in buffsize in bytes
      subroutine init_comm_data(REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH
        integer :: msg_size_request, msg_size_data
        integer :: buffsize !< size of bsend buffer in bytes
        integer :: ierr

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

        walk_status = WALK_STILL_RUNNING

        timings_comm = 0.

    end subroutine init_comm_data



    subroutine run_communication_loop_inner(walk_finished, nummessages)
        implicit none
        include 'mpif.h'
        logical :: msg_avail
        integer :: ierr
        integer :: stat(MPI_STATUS_SIZE)
        integer :: reqhandle
        integer*8 :: requested_key
        type (multipole), allocatable :: child_data(:) ! child data to be received - maximum up to eight children per particle
        integer :: num_children
        integer :: ipe_sender, msg_tag
        integer, intent(inout), dimension(mintag:maxtag) :: nummessages
        logical, intent(inout) :: walk_finished(num_pe)
        integer :: dummy, i
        real*8 :: tcomm

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

          ipe_sender = stat(MPI_SOURCE)
          msg_tag    = stat(MPI_TAG)

          ! the functions returns the number of received messages of any tag
          nummessages(msg_tag) = nummessages(msg_tag) + 1

          select case (msg_tag)
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


    subroutine send_data(requested_key, ipe_sender)
      use module_htable
      implicit none
      include 'mpif.h'
      integer*8, intent(in) :: requested_key
      integer, intent(in) :: ipe_sender
      integer :: process_addr
      type(multipole), target :: children_to_send(8)
      type(multipole), pointer :: c
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
         c=>children_to_send(ic)
           c        = tree_nodes(node_child(ic))
           c%key    = key_child(ic)   ! TODO: this data is maybe not consistently stored in tree_nodes array
           c%byte   = byte_child(ic)  ! therefore, we have to take it directly form the htable --> repair this
           c%leaves = leaves_child(ic)
           c%owner  = owner_child(ic)
      end do

      ! Ship child data back to PE that requested it
      call MPI_IBSEND( children_to_send(1:nchild), nchild, MPI_type_multipole, ipe_sender, TAG_REQUESTED_DATA, &
              MPI_COMM_WORLD, reqhandle, ierr )
      call MPI_REQUEST_FREE(reqhandle, ierr)

      ! statistics on number of sent children-packages
      sum_ships = sum_ships + 1

    end subroutine send_data


    function send_request(req)
      use module_htable
      implicit none
      include 'mpif.h'
      logical :: send_request
      type(t_request_queue_entry), intent(in) :: req
      integer :: reqhandle, ierr

      if (.not. BTEST( htable(req%addr)%childcode, CHILDCODE_BIT_REQUEST_SENT ) ) then
        ! send a request to PE req_queue_owners(req_queue_top)
        ! telling, that we need child data for particle request_key(req_queue_top)
        call MPI_IBSEND(req%key, 1, MPI_INTEGER8, req%owner, TAG_REQUEST_KEY, MPI_COMM_WORLD, reqhandle, ierr )
        call MPI_REQUEST_FREE( reqhandle, ierr)

        htable(req%addr)%childcode = ibset(htable(req%addr)%childcode, CHILDCODE_BIT_REQUEST_SENT )

        send_request = .true.
      else
        send_request = .false.
      end if

    end function



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