module module_tree_walk_smpss_utils

  use treetypes
  use module_tree_walk_communicator

  integer :: nparticles
  type(t_particle), pointer, dimension(:) :: particles
  type(t_calc_force_params) :: cf_par
  real*8 :: vbox(3) !< real space shift vector of box to be processed

  !!!!! simiplification, max_rank = num_pe-1
  integer :: max_rank

  ! inside central box?
  logical :: fcentral

  !!!!! chunk variables
  integer :: chunk_number
  integer :: chunk_size_default
  integer, dimension(:), allocatable :: chunk_sizes
  integer*8, dimension(:,:), allocatable :: chunk_status
  integer*8, dimension(:,:), allocatable :: chunk_requests
  integer,   dimension(:,:), allocatable :: chunk_particles
  integer*8, dimension(:,:), allocatable :: chunk_rlvl1


  integer :: nfetches
  real    :: sleep_dummy

  !!!!! size of a node, needed for mac evaluation
  real*8,    dimension(:),   allocatable :: boxlength2

  ! constants
  integer, public, parameter :: CHILDCODE_BIT_CHILDREN_AVAILABLE =  9

  integer, parameter :: ANSWER_BUFF_LENGTH   = 10000 !< amount of possible entries in the BSend buffer for shipping child data
  integer, parameter :: REQUEST_QUEUE_LENGTH = 400000 !< maximum length of request queue

  
  logical,   dimension(:), allocatable         :: walk_finished 
  integer*1, dimension(:), allocatable, target :: bsend_buffer
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tree_walk_smpss_setconst(nparticles_, particles_, cf_par_, vbox_)

    use treevars, only: num_pe
    
    implicit none
    
    integer, intent(in) :: nparticles_
    type(t_particle), target, intent(in) :: particles_(:)
    type(t_calc_force_params), intent(in) :: cf_par_
    real*8, intent(in) :: vbox_(3) 
    
    nparticles       =  nparticles_
    particles        => particles_
    cf_par           =  cf_par_
    vbox             =  vbox_

    !cf_par%theta2 = cf_par%theta**2

    fcentral    = (vbox(1)**2 + vbox(2)**2 + vbox(3)**2 ) .eq. 0
 
    chunk_size_default = 500
   
    max_rank = num_pe-1

  end subroutine tree_walk_smpss_setconst
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tree_walk_smpss_setup_chunks
    
    use treevars, only: nlev, boxsize

    implicit none
    include 'mpif.h'
    
    integer :: full_chunks, remainder_chunk
    integer :: ccnt
    integer :: pcnt, lcnt

    integer :: tchunk, tposition

    integer :: ierr
    integer :: msg_size_request, msg_size_data, buffsize

    ! number of full chunks
    full_chunks     = nparticles / chunk_size_default
    ! last not fully populated chunk
    remainder_chunk = MOD(nparticles, chunk_size_default)
    ! total number of chunks (full + one semi populated)
    chunk_number    = full_chunks
    if(remainder_chunk .gt. 0) chunk_number = chunk_number + 1
    
    write(*,*) "number of full chunks: ", full_chunks, "; remaining particles: ", remainder_chunk
    
    ! allocate memory
    allocate(chunk_sizes    (chunk_number                    ))
    allocate(chunk_status   (chunk_size_default, chunk_number))
    allocate(chunk_particles(chunk_size_default, chunk_number))
    allocate(chunk_requests (chunk_size_default, chunk_number))
    allocate(chunk_rlvl1    (chunk_size_default, chunk_number))
    
    ! set the chunk sizes
    do ccnt=1,full_chunks
       chunk_sizes(ccnt) = chunk_size_default
    end do
    chunk_sizes(chunk_number) = remainder_chunk
    
    ! distribute particles on chunks
    ! clear status and requests list
    chunk_status = -1_8
    do pcnt=1, nparticles

       tchunk    = (pcnt-1) / chunk_size_default + 1
       tposition = pcnt - (tchunk-1)*chunk_size_default 

       chunk_particles(tposition, tchunk) = pcnt
       chunk_status   (tposition, tchunk) = 1_8
       chunk_requests (tposition, tchunk) = 0_8
       chunk_rlvl1    (tposition, tchunk) = 0_8
    end do

    ! setup box lengths for the mac evaluation
    allocate(boxlength2(0:nlev))

    boxlength2(0) = boxsize**2
    do lcnt=1, nlev
       boxlength2(lcnt) =  boxlength2(lcnt-1)/4._8
    end do

    allocate(walk_finished(max_rank+1))

    walk_finished = .false.

    call init_comm_data(REQUEST_QUEUE_LENGTH, ANSWER_BUFF_LENGTH)

    nfetches = 0
    sleep_dummy = 0

  end subroutine tree_walk_smpss_setup_chunks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sleep_loop(sd)

    implicit none

    real, intent(inout) :: sd

    real    :: a, b,c
    integer :: lcnt

    a = 1.045
    c = 0.998
    do lcnt=1, 1000000
       
       b = b + a*c
       
    end do

    sd = sd + b

  end subroutine sleep_loop

  subroutine tree_walk_smpss_comm_send_requests(short_requests, short_req_owner, short_size)
    
    implicit none
    include 'mpif.h'
    
    integer,   intent(in)                        :: short_size
    integer*8, intent(in), dimension(short_size) :: short_requests
    integer,   intent(in), dimension(short_size) :: short_req_owner

    integer :: reqhandle, ierr

    integer :: rcnt

    do rcnt=1, short_size

       call MPI_IBSEND(short_requests(rcnt), 1, MPI_INTEGER8, &
            short_req_owner(rcnt), TAG_REQUEST_KEY, &
            MPI_COMM_WORLD, reqhandle, ierr )
       call MPI_REQUEST_FREE( reqhandle, ierr)

    end do

  end subroutine tree_walk_smpss_comm_send_requests


  subroutine tree_walk_smpss_comm_loop_inner(short_requests, short_req_owner, short_size)

    use treevars

    implicit none
    include 'mpif.h'

    integer,   intent(in)                        :: short_size
    integer*8, intent(in), dimension(short_size) :: short_requests
    integer,   intent(in), dimension(short_size) :: short_req_owner

    integer :: received_answers

    logical :: msg_avail
    integer :: ierr
    integer :: stat(MPI_STATUS_SIZE)
    integer :: reqhandle
    integer*8 :: requested_key
    type (t_tree_node) :: child_data(8) ! child data to be received - maximum up to eight children per particle
    integer :: num_children
    integer :: ipe_sender

    integer :: comm_dummy

    integer :: i
    
    real*8 :: ta, tb

    nfetches = MAX(nfetches,short_size)

    ta = MPI_WTIME()

    ! notify rank 0 if we completed our traversal
    if (walk_status == WALK_IAM_FINISHED) call send_walk_finished()

    received_answers = 0

    ! send our requested keys
    call tree_walk_smpss_comm_send_requests(short_requests, short_req_owner, short_size)

    do ! exit statement at bottom of "infinite"-do-loop

       !call sleep_loop(sleep_dummy)

       ! probe for incoming messages
       ! non blocking for "just checking" the inbox, and potentially exit routine
       if(short_size .eq. received_answers) then
          call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_avail, stat, ierr)
       ! blocking, as there must/will be communication, due to outstanding answers
       else
          call MPI_PROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, stat, ierr)
          msg_avail = .true.
       end if

       !write(*,*) "comm wait on rank", me, "for ", received_answers, " out of ", short_size, " status:", walk_status

       if (msg_avail) then

          ipe_sender = stat(MPI_SOURCE)

          !call Extrae_eventandcounters(1001, 0) 

          select case (stat(MPI_TAG))
             ! another PE requested child data for a certain key
          case (TAG_REQUEST_KEY)

             !call Extrae_eventandcounters(1001, 1) 
             ! actually receive this request...
             call MPI_RECV( requested_key, 1, MPI_INTEGER8, ipe_sender, TAG_REQUEST_KEY, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
             ! ... and answer it
             call send_data(requested_key, ipe_sender)

             ! some PE answered our request and sends
          case (TAG_REQUESTED_DATA)
             !call Extrae_eventandcounters(1001, 2) 
             ! actually receive the data... 
             call MPI_GET_COUNT(stat, MPI_TYPE_tree_node, num_children, ierr)
             call MPI_RECV( child_data, num_children, MPI_TYPE_tree_node, ipe_sender, TAG_REQUESTED_DATA, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
             !call Extrae_eventandcounters(1001, 3) 
             ! ... and put it into the tree and all other data structures
             call unpack_data(child_data, num_children, ipe_sender)

             received_answers = received_answers + 1

             ! rank 0 does bookkeeping about which PE is already finished with its walk
             ! no one else will ever receive this message tag
          case (TAG_FINISHED_PE)
             !call Extrae_eventandcounters(1001, 4) 
             ! actually receive the data (however, we are not interested in it here)
             call MPI_RECV( comm_dummy, 1, MPI_INTEGER, ipe_sender, TAG_FINISHED_PE, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

             walk_finished(ipe_sender+1) = .true.

             write(*,*) "rank send finish: ", ipe_sender

             if ( all(walk_finished) ) then
                ! all PEs have to be informed

                do i=0,num_pe-1
                   call MPI_IBSEND(comm_dummy, 1, MPI_INTEGER, i, TAG_FINISHED_ALL, &
                        MPI_COMM_WORLD, reqhandle, ierr )
                   call MPI_REQUEST_FREE( reqhandle, ierr)
                end do
             end if


             ! all PEs have finished their walk
          case (TAG_FINISHED_ALL)
             !call Extrae_eventandcounters(1001, 5) 
             call MPI_RECV( comm_dummy, 1, MPI_INTEGER, ipe_sender, TAG_FINISHED_ALL, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) 

             walk_status = WALK_ALL_FINISHED

          end select

          !call Extrae_eventandcounters(1001, 0) 
       end if
       !else
       ! exit statement for do-loop
       if (received_answers .eq. short_size) exit
       !end if
    end do



    tb = MPI_WTIME()

    !if(me.eq.0) write(*,*) "time in comm: ", tb-ta

  end subroutine tree_walk_smpss_comm_loop_inner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tree_walk_smpss_gather_request_numbers(reqs, reqs_size, rnumber, snumber)

    use treevars, only: me
    use module_htable
    implicit none
    include 'mpif.h'

    integer,   intent(in)                          :: reqs_size
    integer*8, intent(in),  dimension(reqs_size)   :: reqs
    integer,   intent(out), dimension(0:max_rank)  :: rnumber, snumber

    integer :: rcnt, pcnt
    
    integer :: ierr

    integer :: owner

    !!! clear my request number list
    do pcnt=0, max_rank
       rnumber(pcnt) = 0
    end do

    !!! gather my request number list
    do rcnt=1, reqs_size
       owner = htable(key2addr(reqs(rcnt), "gather reqs number")) % owner  
       rnumber(owner) = rnumber(owner) + 1
    end do

    !!! exchange req numbers -> snumber
    call MPI_ALLTOALL(rnumber, 1, MPI_INTEGER, snumber, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    do pcnt=0, max_rank
       write(*,*) "rank ", me, " sends ", snumber(pcnt), " to ", pcnt
    end do
    

  end subroutine tree_walk_smpss_gather_request_numbers


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tree_walk_smpss_print_keylist(list, lsize, name)

    use module_htable

    implicit none

    integer*8,        intent(in), dimension(lsize)  :: list
    integer,          intent(in)                    :: lsize
    character(LEN=*), intent(in)                    :: name

    integer :: cnt
    integer*8 :: tkey

    do cnt=1, lsize
       tkey = list(cnt)
       write(*,'(3a,i5,o21,i5)') "key list ", name, "; index, key, owner", &
            cnt, tkey, htable(key2addr(tkey, "print keylist")) % owner
    end do

  end subroutine tree_walk_smpss_print_keylist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tree_walk_smpss_unique_keylist(llist, lsize, slist, olist, ssize)

    use module_htable

    implicit none

    integer*8, intent(in), dimension(lsize)  :: llist
    integer,   intent(in)                    :: lsize

    integer*8, intent(out), dimension(lsize) :: slist
    integer,   intent(out), dimension(lsize) :: olist
    integer,   intent(out)                   :: ssize

    !!! counter
    integer :: scnt, lcnt

    !!! tmp
    integer*8 :: tkey

    !!! flags
    logical :: ffound
    logical :: falreadyfetched

    ssize = 0
    do lcnt=1, lsize
       
       tkey = llist(lcnt)

       if (tkey.gt.0) then
          
          ffound = .false.
          do scnt=1, ssize
             if(tkey .eq. slist(scnt)) ffound = .true.
          end do

          falreadyfetched = btest( htable(key2addr(tkey, "unique key list child test"))%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE)
          if(.not. ffound .and. .not.falreadyfetched) then
             ssize = ssize + 1
             slist(ssize) = tkey
             olist(ssize) = htable(key2addr(tkey, "unique key list get owner")) % owner
          end if
       end if
    end do

  end subroutine tree_walk_smpss_unique_keylist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tree_walk_smpss_global_finish(local_finished, global_finished)

    implicit none

    include 'mpif.h'

    logical, intent(in) :: local_finished
    logical, intent(out) :: global_finished

    integer :: ierr

    call MPI_ALLREDUCE(local_finished, global_finished, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  end subroutine tree_walk_smpss_global_finish

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine tree_walk_smpss_free_chunks

    implicit none
    include 'mpif.h'

    integer :: ierr, buffsize, dummy

    deallocate(chunk_sizes)
    deallocate(chunk_status)
    deallocate(chunk_requests)
    deallocate(chunk_particles)
    deallocate(chunk_rlvl1)

    deallocate(boxlength2)

    deallocate(walk_finished)

    call uninit_comm_data()

  end subroutine tree_walk_smpss_free_chunks

  subroutine tree_walk_smpss_comm_serve_inner
    use treevars
    
    implicit none
    include 'mpif.h'
    
    logical :: msg_avail
    integer :: ierr
    integer :: stat(MPI_STATUS_SIZE)
    integer :: reqhandle
    integer*8 :: requested_key
    type (t_tree_node) :: child_data(8) ! child data to be received - maximum up to eight children per particle
    integer :: num_children
    integer :: ipe_sender
    
    integer :: comm_dummy
    
    integer :: i

    msg_avail = .true.
    do while(msg_avail .eqv. .true.)
       call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, msg_avail, stat, ierr)
     
       if (msg_avail) then
          
          ipe_sender = stat(MPI_SOURCE)
          
          select case (stat(MPI_TAG))
             ! another PE requested child data for a certain key
          case (TAG_REQUEST_KEY)
             
             !call Extrae_eventandcounters(1001, 1) 
             ! actually receive this request...
             call MPI_RECV( requested_key, 1, MPI_INTEGER8, ipe_sender, TAG_REQUEST_KEY, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
             ! ... and answer it
             call send_data(requested_key, ipe_sender)
             
          case (TAG_FINISHED_PE)
             !call Extrae_eventandcounters(1001, 4) 
             ! actually receive the data (however, we are not interested in it here)
             call MPI_RECV( comm_dummy, 1, MPI_INTEGER, ipe_sender, TAG_FINISHED_PE, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
             
             walk_finished(ipe_sender+1) = .true.
             
             write(*,*) "rank send finish: ", ipe_sender
             
             if ( all(walk_finished) ) then
                ! all PEs have to be informed
                
                do i=0,num_pe-1
                   call MPI_IBSEND(comm_dummy, 1, MPI_INTEGER, i, TAG_FINISHED_ALL, &
                        MPI_COMM_WORLD, reqhandle, ierr )
                   call MPI_REQUEST_FREE( reqhandle, ierr)
                end do
             end if
          
             
             ! all PEs have finished their walk
          case (TAG_FINISHED_ALL)
             !call Extrae_eventandcounters(1001, 5) 
             call MPI_RECV( comm_dummy, 1, MPI_INTEGER, ipe_sender, TAG_FINISHED_ALL, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) 
             
             walk_status = WALK_ALL_FINISHED
             
          end select
          
          !call Extrae_eventandcounters(1001, 0) 
       end if
    end do

  end subroutine tree_walk_smpss_comm_serve_inner
  
end module module_tree_walk_smpss_utils


