!>
!> Accelerator module, in this case for a GPU.
!>
module module_accelerator
   use module_interaction_specific_types
   use module_atomic_ops

   implicit none

   integer, public, parameter :: ACC_THREAD_STATUS_WAITING  = 1
   integer, public, parameter :: ACC_THREAD_STATUS_STARTING = 2
   integer, public, parameter :: ACC_THREAD_STATUS_STARTED  = 3
   integer, public, parameter :: ACC_THREAD_STATUS_STOPPING = 4
   integer, public, parameter :: ACC_THREAD_STATUS_STOPPED  = 5

   type(t_critical_section), pointer :: queue_lock

   contains
   !>
   !> main routine to drive GPU
   !>
   !> this is a separate thread since multi-threaded GPU support is ... sth else rather
   !>
   function acc_loop() bind(c)
      use, intrinsic :: iso_c_binding
      use pthreads_stuff, only: pthreads_sched_yield, get_my_core, pthreads_exitthread
      use module_debug
#ifdef __OPENACC
      use openacc
#endif
      implicit none
      include 'mpif.h'
   
      type(c_ptr) :: acc_loop
      type(t_particle_thread) :: tmp_particle
      integer :: tmp_top, q_tmp
#ifdef __OPENACC
      integer :: strt, stp, tck
#endif
      ! the following lines in case we want to submit an argument...
      ! right now the data-strucutre is 'global' in module_interaction_specific_types
      !type(c_ptr), value :: arg
      !type(acc_driver), pointer :: acc
      !
      !gpu => null()
      !call c_f_pointer(arg, acc)
      !DEBUG_ASSERT(associated(acc))
   
      ! store ID of comm-thread processor
      acc%processor_id = get_my_core()
      call atomic_write_barrier()

      write(*,*) 'GPU thread on core ', acc%processor_id

      ! signal we're starting...
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTING)

      ! reset all queues
      call atomic_allocate_int(acc%q_top)
      call atomic_store_int(acc%q_top, 0)
      call atomic_allocate_int(acc%q_bottom)
      call atomic_store_int(acc%q_bottom, 0)
      call atomic_allocate_int(acc%q_len)
      call atomic_store_int(acc%q_len, ACC_QUEUE_LENGTH)

      ! init critical lock
      call critical_section_allocate(queue_lock)

      ! start GPU
#ifdef __OPENACC
      ! start GPU stuff here to reduce time spent in walk
      write(*,*) '[',42,'] Starting GPU Context...' ! me
      call system_clock(strt)
      call acc_init(acc_device_nvidia)
      call system_clock(stp,tck)
      write(*,*) '[',42,'] ...done', real(stp-strt)/tck,'secs' ! me
#endif

      ! signal successfull start
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)

      ! create GPU memory
#ifdef __OPENACC
      !$acc data create(gpu_l, gpu)
#endif

      do while (atomic_load_int(acc%thread_status) .ne. ACC_THREAD_STATUS_STOPPED)
   
         ! check list for available data
         do while (atomic_load_int(acc%q_top) .ne. atomic_load_int(acc%q_bottom))
            tmp_top = mod(atomic_load_int(acc%q_top), ACC_QUEUE_LENGTH) + 1
   
            ! first check whether the entry is actually valid	  
            if (acc%acc_queue(tmp_top)%entry_valid) then
               call atomic_read_barrier() ! make sure that reads of parts of the queue entry occurr in the correct order

               !! move list, copy data
               !acc%acc_queue(tmp_top)%partner ...
               !acc%acc_queue(tmp_top)%queued ...

               !! run GPU kernel

               !! kill list

               ! copy particle info to temporal local copy
               tmp_particle = acc%acc_queue(tmp_top)%particle
               tmp_particle%queued = acc%acc_queue(tmp_top)%queued
               tmp_particle%partner => acc%acc_queue(tmp_top)%partner

               ! call GPU
               call kernel_node(tmp_particle, acc%acc_queue(tmp_top)%eps, acc%acc_queue(tmp_top)%pen)

               call critical_section_enter(queue_lock)

                  ! free lists
                  deallocate(acc%acc_queue(tmp_top)%partner)
                  nullify(acc%acc_queue(tmp_top)%partner)
   
                  call atomic_store_int(acc%q_top, tmp_top)
                  q_tmp = atomic_fetch_and_increment_int(acc%q_len)
   
                  ! we have to invalidate this request queue entry.
                  ! this shows that we actually processed it and prevents it from accidentially being resent after the queue wrapped around
                  acc%acc_queue(tmp_top)%entry_valid = .false.

               call critical_section_leave(queue_lock)

            else
               ! the next entry is not valid (obviously it has not been stored completely until now -> we abort here and try again later
               exit
            end if
         end do
   
      end do

      ! free GPU memory
#ifdef __OPENACC
      !$acc end data
#endif

      ! free queues and lock
      call critical_section_deallocate(queue_lock)
      call atomic_deallocate_int(acc%q_top)
      call atomic_deallocate_int(acc%q_bottom)
      call atomic_deallocate_int(acc%q_len)

      write(*,*) 'GPU thread terminating'
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STOPPED)
   
      acc_loop = c_null_ptr
      ERROR_ON_FAIL(pthreads_exitthread())
   
   end function acc_loop

   subroutine dispatch_list(particle, eps2, penalty)
      implicit none

      type(t_particle_thread), target, intent(in) :: particle
      real*8, intent(in) :: eps2, penalty
      integer :: local_queue_bottom, gpu_status, q_tmp

      ! safety check if the GPU is up and running
      gpu_status = atomic_load_int(acc%thread_status)
      do while (gpu_status .ne. ACC_THREAD_STATUS_STARTED)
         if (gpu_status .lt. ACC_THREAD_STATUS_STARTED) then
            write(*,*) '..waiting for GPU thread..'
            call system('sleep 1s')
            gpu_status = atomic_load_int(acc%thread_status)
         else
            DEBUG_ERROR(*, "NO GPU!")
         endif
      end do

      call critical_section_enter(queue_lock)

         do while( atomic_load_int(acc%q_len) .le. 0 )
            ! busy loop while the queue is processed
            ERROR_ON_FAIL(pthreads_sched_yield())
         end do
   
         ! thread-safe way of reserving storage for our request
         local_queue_bottom = atomic_mod_increment_and_fetch_int(acc%q_bottom, ACC_QUEUE_LENGTH)
   
         if (local_queue_bottom == atomic_load_int(acc%q_top)) then
            DEBUG_ERROR(*, "ACC_QUEUE_LENGTH is too small: ", ACC_QUEUE_LENGTH)
            write(*,*) 'ACC queue exhausted...'
         end if
   
         ! the communicator will check validity of the request and will only proceed as soon as the entry is valid -- this actually serializes the requests
         DEBUG_ASSERT(.not. acc%acc_queue(local_queue_bottom)%entry_valid)
   
         ! store link to particle, so we know what forces to update
         acc%acc_queue(local_queue_bottom)%particle = particle
         ! move list information to queue, list info in particle will be deleter by calling subroutine...
         acc%acc_queue(local_queue_bottom)%queued = particle%queued
         acc%acc_queue(local_queue_bottom)%partner => particle%partner
         ! store epsilon and work_penalty to preserve original calc_force functionality
         acc%acc_queue(local_queue_bottom)%eps = eps2
         acc%acc_queue(local_queue_bottom)%pen = penalty
   
         q_tmp = atomic_fetch_and_decrement_int(acc%q_len)
   
         call atomic_write_barrier() ! make sure the above information is actually written before flagging the entry valid by writing the owner
         acc%acc_queue(local_queue_bottom)%entry_valid    = .true.

      call critical_section_leave(queue_lock)

      return
   end subroutine dispatch_list

end module module_accelerator
