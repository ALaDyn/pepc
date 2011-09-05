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
module tree_walk_pthreads
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

      allocate(threaddata(num_walk_threads))

      ! start the worker threads...
      do ith = 1,num_walk_threads
        threaddata(ith)%id = ith
        threaddata(ith)%finished = .false.
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

      deallocate(threaddata)

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
      use tree_walk_communicator
      implicit none
      deallocate(boxlength2)
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
            comm_on_shared_processor = .true.
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
      call notify_walk_finished()

      if (same_core_as_communicator) comm_on_shared_processor = .false.

      walk_worker_thread = c_null_ptr

      my_threaddata%runtime_seconds = MPI_WTIME() - my_threaddata%runtime_seconds

      call retval(pthreads_exitthread(), "walk_worker_thread:pthread_exit")

    end function walk_worker_thread




   function walk_single_particle(myidx, nodeidx, defer_list, defer_list_entries, listlengths, partner_leaves, my_threaddata)
      use tree_walk_communicator
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

