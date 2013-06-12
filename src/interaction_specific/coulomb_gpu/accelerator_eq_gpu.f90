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
   integer, public, parameter :: ACC_THREAD_STATUS_FLUSH    = 4
   integer, public, parameter :: ACC_THREAD_STATUS_STOPPING = 5
   integer, public, parameter :: ACC_THREAD_STATUS_STOPPED  = 6

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
      use cudadevice
      use cudafor
#endif
      implicit none
      include 'mpif.h'
   
      type(c_ptr) :: acc_loop
      type(t_particle_thread) :: tmp_particle
      integer :: tmp_top, q_tmp
#ifdef __OPENACC
      integer :: strt, stp, tck
#endif
      ! kernel data
      real*8, dimension(MAX_IACT_PARTNERS, size(gpu)) :: e_1, e_2, e_3, pot
      type point
         type(t_particle_results), pointer :: results
         real*8, pointer :: work
      end type point
      type(point), dimension(size(gpu)) :: ptr
      real*8 :: e_1_, e_2_, e_3_, pot_

      ! GPU kernel stuff... move to a function again?
      integer :: idx, idx_, queued(size(gpu))
      real*8 :: dist2, eps2, WORKLOAD_PENALTY_INTERACTION

      real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6

      integer :: flushy
      logical :: eps_update

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
      gpu_id = 0
      flushy = 0
      eps_update = .true.

      ! signal successfull start
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)

      ! create GPU memory
!!!      !$acc data create(gpu_l, gpu, eps2)
#ifdef __OPENACC
!!!#define ONE
#ifdef ONE
      !$acc data create(gpu_l, gpu, e_1, e_2, e_3, pot)
#else
      !$acc data create(gpu_l, gpu)
#endif
#endif
      
      do while (atomic_load_int(acc%thread_status) .lt. ACC_THREAD_STATUS_STOPPING)
   
         ! check list for available data
         do while (atomic_load_int(acc%q_top) .ne. atomic_load_int(acc%q_bottom))
            tmp_top = mod(atomic_load_int(acc%q_top), ACC_QUEUE_LENGTH) + 1
   
            ! first check whether the entry is actually valid	  
            if (acc%acc_queue(tmp_top)%entry_valid) then
               call atomic_read_barrier() ! make sure that reads of parts of the queue entry occurr in the correct order

#define INLINEKERNEL
#ifdef INLINEKERNEL
               ! find a stream
               gpu_id = mod(gpu_id,size(gpu)) + 1
!!!               !$acc wait(gpu_id)

#ifdef __OPENACC
               flushy = mod(flushy+1,5000)
               if ( flushy .eq. 0 ) write(*,*) 'sync ',cudaDeviceSynchronize()
#endif

               if(gpu_id .lt. 0 .or. gpu_id .gt. size(gpu)) write(*,*) 'BUGGER'
               ! move list, copy data
               queued(gpu_id) = acc%acc_queue(tmp_top)%particle%queued
               ptr(gpu_id)%results => acc%acc_queue(tmp_top)%particle%results
               ptr(gpu_id)%work => acc%acc_queue(tmp_top)%particle%work
               do idx = 1, queued(gpu_id)
                  gpu(gpu_id)%delta1(idx) = acc%acc_queue(tmp_top)%partner(idx)%delta(1)
                  gpu(gpu_id)%delta2(idx) = acc%acc_queue(tmp_top)%partner(idx)%delta(2)
                  gpu(gpu_id)%delta3(idx) = acc%acc_queue(tmp_top)%partner(idx)%delta(3)
                  gpu(gpu_id)%charge(idx) = acc%acc_queue(tmp_top)%partner(idx)%node%charge
                  gpu(gpu_id)%dip1(idx)   = acc%acc_queue(tmp_top)%partner(idx)%node%dip(1)
                  gpu(gpu_id)%dip2(idx)   = acc%acc_queue(tmp_top)%partner(idx)%node%dip(2)
                  gpu(gpu_id)%dip3(idx)   = acc%acc_queue(tmp_top)%partner(idx)%node%dip(3)
                  gpu(gpu_id)%quad1(idx)  = acc%acc_queue(tmp_top)%partner(idx)%node%quad(1)
                  gpu(gpu_id)%quad2(idx)  = acc%acc_queue(tmp_top)%partner(idx)%node%quad(2)
                  gpu(gpu_id)%quad3(idx)  = acc%acc_queue(tmp_top)%partner(idx)%node%quad(3)
                  gpu(gpu_id)%xyquad(idx) = acc%acc_queue(tmp_top)%partner(idx)%node%xyquad
                  gpu(gpu_id)%yzquad(idx) = acc%acc_queue(tmp_top)%partner(idx)%node%yzquad
                  gpu(gpu_id)%zxquad(idx) = acc%acc_queue(tmp_top)%partner(idx)%node%zxquad
               enddo

               ! update GPU data
#ifdef __OPENACC
#ifdef ONE
               !$acc update device(gpu(gpu_id:gpu_id)) async(gpu_id)
#else
               !$acc update device(gpu(gpu_id:gpu_id))
#endif
#endif

#ifdef ONE
               ! really, this is not necessary. we write to all of those entries anyway...
!!               e_1(:, gpu_id) = 0.d0
!!               e_2(:, gpu_id) = 0.d0
!!               e_3(:, gpu_id) = 0.d0
!!               pot(:, gpu_id) = 0.d0
#ifdef __OPENACC
!!               !$acc update device(e_1(:,gpu_id:gpu_id), e_2(:,gpu_id:gpu_id), e_3(:,gpu_id:gpu_id), pot(:,gpu_id:gpu_id)) async(gpu_id)
#endif
#else
               ! flush reduction variables (should not be necessary)
               e_1_ = 0.d0
               e_2_ = 0.d0
               e_3_ = 0.d0
               pot_ = 0.d0
#endif

               ! read penalty and epsilon
               eps2 = acc%acc_queue(tmp_top)%eps
!!!               !$acc update device(eps2) if(eps_update)   !<<< this should WORK!
               eps_update = .false.

               WORKLOAD_PENALTY_INTERACTION = acc%acc_queue(tmp_top)%pen

               ! run GPU kernel
#ifdef __OPENACC
#ifdef ONE
               !$acc parallel loop                                                                                 &       
               !$acc present(gpu(gpu_id:gpu_id))                                                                   &
               !$acc present(e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), pot(:,gpu_id))                           &
               !$acc private(dist2,rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6)  &
               !$acc async(gpu_id)
#else
               !$acc parallel loop                                                                                 &
               !$acc present(gpu(gpu_id:gpu_id))                                                                   &
               !$acc reduction(+: e_1_, e_2_, e_3_, pot_)                                                          &
               !$acc private(dist2,rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6)
#endif
#endif
               do idx = 1, queued(gpu_id)
             
                  dist2     =         gpu(gpu_id)%delta1(idx) * gpu(gpu_id)%delta1(idx)
                  dist2     = dist2 + gpu(gpu_id)%delta2(idx) * gpu(gpu_id)%delta2(idx)
                  dist2     = dist2 + gpu(gpu_id)%delta3(idx) * gpu(gpu_id)%delta3(idx)
            
                  dx = gpu(gpu_id)%delta1(idx)
                  dy = gpu(gpu_id)%delta2(idx)
                  dz = gpu(gpu_id)%delta3(idx)
            
                  r  = sqrt(dist2+1e-4)!eps2) ! eps2 is added in calling routine to have plummer instead of coulomb here
                  rd = 1.d0/r
                  rd2 = rd *rd
                  rd3 = rd *rd2
                  rd5 = rd3*rd2
                  rd7 = rd5*rd2
            
                  dx2 = dx*dx
                  dy2 = dy*dy
                  dz2 = dz*dz
                  dx3 = dx*dx2
                  dy3 = dy*dy2
                  dz3 = dz*dz2
            
                  fd1 = 3.d0*dx2*rd5 - rd3
                  fd2 = 3.d0*dy2*rd5 - rd3
                  fd3 = 3.d0*dz2*rd5 - rd3
                  fd4 = 3.d0*dx*dy*rd5
                  fd5 = 3.d0*dy*dz*rd5
                  fd6 = 3.d0*dx*dz*rd5
            
#ifdef ONE
                  pot(idx, gpu_id) = gpu(gpu_id)%charge(idx)*rd                                                          &  !  monopole term
#else
                  pot_ = pot_ + gpu(gpu_id)%charge(idx)*rd                                                               &  !  monopole term
#endif
                        + (dx*gpu(gpu_id)%dip1(idx) + dy*gpu(gpu_id)%dip2(idx) + dz*gpu(gpu_id)%dip3(idx))*rd3           &  !  dipole
                        + 0.5d0*(fd1*gpu(gpu_id)%quad1(idx)  + fd2*gpu(gpu_id)%quad2(idx)  + fd3*gpu(gpu_id)%quad3(idx)) &  !  quadrupole
                        +        fd4*gpu(gpu_id)%xyquad(idx) + fd5*gpu(gpu_id)%yzquad(idx) + fd6*gpu(gpu_id)%zxquad(idx)
            
#ifdef ONE
                  e_1(idx, gpu_id) = gpu(gpu_id)%charge(idx)*dx*rd3                                                      &  ! monopole term
#else
                  e_1_ = e_1_ + gpu(gpu_id)%charge(idx)*dx*rd3                                                                  &  ! monopole term
#endif
                            + fd1*gpu(gpu_id)%dip1(idx) + fd4*gpu(gpu_id)%dip2(idx) + fd6*gpu(gpu_id)%dip3(idx)          &  ! dipole term
                            + 3.d0   * (                                                                                 &  ! quadrupole term
                               0.5d0 * (                                                                                 &
                                   ( 5.d0*dx3   *rd7 - 3.d0*dx*rd5 )*gpu(gpu_id)%quad1(idx)                              &
                                 + ( 5.d0*dx*dy2*rd7 -      dx*rd5 )*gpu(gpu_id)%quad2(idx)                              &
                                 + ( 5.d0*dx*dz2*rd7 -      dx*rd5 )*gpu(gpu_id)%quad3(idx)                              &
                               )                                                                                         &
                               + ( 5.d0*dy*dx2  *rd7 - dy*rd5 )*gpu(gpu_id)%xyquad(idx)                                  &
                               + ( 5.d0*dz*dx2  *rd7 - dz*rd5 )*gpu(gpu_id)%zxquad(idx)                                  &
                               + ( 5.d0*dx*dy*dz*rd7          )*gpu(gpu_id)%yzquad(idx)                                  &
                              )
            
#ifdef ONE
                  e_2(idx, gpu_id) = gpu(gpu_id)%charge(idx)*dy*rd3                                                      &
#else
                  e_2_= e_2_ + gpu(gpu_id)%charge(idx)*dy*rd3                                                                   &
#endif
                            + fd2*gpu(gpu_id)%dip2(idx) + fd4*gpu(gpu_id)%dip1(idx) + fd5*gpu(gpu_id)%dip3(idx)          &
                            + 3 * (                                                                                      &
                               0.5d0 * (                                                                                 &
                                   ( 5*dy3*rd7    - 3*dy*rd5 )*gpu(gpu_id)%quad2(idx)                                    &
                                 + ( 5*dy*dx2*rd7 -   dy*rd5 )*gpu(gpu_id)%quad1(idx)                                    &
                                 + ( 5*dy*dz2*rd7 -   dy*rd5 )*gpu(gpu_id)%quad3(idx)                                    &
                               )                                                                                         &
                               + ( 5*dx*dy2  *rd7 - dx*rd5 )*gpu(gpu_id)%xyquad(idx)                                     &
                               + ( 5*dz*dy2  *rd7 - dz*rd5 )*gpu(gpu_id)%yzquad(idx)                                     &
                               + ( 5*dx*dy*dz*rd7          )*gpu(gpu_id)%zxquad(idx)                                     &
                              )
            
#ifdef ONE
                  e_3(idx, gpu_id) = gpu(gpu_id)%charge(idx)*dz*rd3                                                      &
#else
                  e_3_ = e_3_ + gpu(gpu_id)%charge(idx)*dz*rd3                                                                  &
#endif
                            + fd3*gpu(gpu_id)%dip3(idx) + fd5*gpu(gpu_id)%dip2(idx) + fd6*gpu(gpu_id)%dip1(idx)          &
                            + 3 * (                                                                                      &
                               0.5d0 * (                                                                                 &
                                 + ( 5*dz3   *rd7 - 3*dz*rd5 )*gpu(gpu_id)%quad3(idx)                                    &
                                 + ( 5*dz*dy2*rd7 -   dz*rd5 )*gpu(gpu_id)%quad2(idx)                                    &
                                 + ( 5*dz*dx2*rd7 -   dz*rd5 )*gpu(gpu_id)%quad1(idx)                                    &
                                              )                                                                          &
                               + ( 5*dx*dz2  *rd7 - dx*rd5 )*gpu(gpu_id)%zxquad(idx)                                     &
                               + ( 5*dy*dz2  *rd7 - dy*rd5 )*gpu(gpu_id)%yzquad(idx)                                     &
                               + ( 5*dx*dy*dz*rd7          )*gpu(gpu_id)%xyquad(idx)                                     &
                              )
               end do
#ifdef __OPENACC
               !$acc end parallel loop
#endif

!!               ! perform the reduction...
!!               !$acc parallel loop present(e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), pot(:,gpu_id)) async(gpu_id)
!!               do idx =  2, queued
!!                  e_1(1, gpu_id) = e_1(1, gpu_id) + e_1(idx, gpu_id)
!!                  e_2(1, gpu_id) = e_2(1, gpu_id) + e_2(idx, gpu_id)
!!                  e_3(1, gpu_id) = e_3(1, gpu_id) + e_3(idx, gpu_id)
!!                  pot(1, gpu_id) = pot(1, gpu_id) + pot(idx, gpu_id)
!!               end do
!!               !$acc end parallel loop
! this was too naive. the threadblocks are not synchronised. the sum will fail...

               ! get data from GPU
#ifdef ONE
#ifdef __OPENACC
!!!               !$acc update host(e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), pot(:,gpu_id)) async(gpu_id)
#endif

!!!#define FAILSRED
#ifdef FAILSRED
               !$acc update host(e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), pot(:,gpu_id))
               !$acc wait
               acc%acc_queue(tmp_top)%particle%results%e(1) = acc%acc_queue(tmp_top)%particle%results%e(1) + sum(e_1(1:queued(gpu_id),gpu_id))
               acc%acc_queue(tmp_top)%particle%results%e(2) = acc%acc_queue(tmp_top)%particle%results%e(2) + sum(e_2(1:queued(gpu_id),gpu_id))
               acc%acc_queue(tmp_top)%particle%results%e(3) = acc%acc_queue(tmp_top)%particle%results%e(3) + sum(e_3(1:queued(gpu_id),gpu_id))
               acc%acc_queue(tmp_top)%particle%results%pot  = acc%acc_queue(tmp_top)%particle%results%pot  + sum(pot(1:queued(gpu_id),gpu_id))
               acc%acc_queue(tmp_top)%particle%work         = acc%acc_queue(tmp_top)%particle%work + queued(gpu_id) * WORKLOAD_PENALTY_INTERACTION

#else
               if (gpu_id .eq. size(gpu)) then
  !             if (.true.) then
#define CPURED
#ifdef CPURED
  !                !$acc update host(e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), pot(:,gpu_id))
!!                  !$acc wait
                  !$acc update host(e_1, e_2, e_3, pot)
                  do idx = 1,size(gpu)
  !                do idx = gpu_id,gpu_id
                     ptr(idx)%results%e(1) = ptr(idx)%results%e(1) + sum(e_1(1:queued(idx),idx))
                     ptr(idx)%results%e(2) = ptr(idx)%results%e(2) + sum(e_2(1:queued(idx),idx))
                     ptr(idx)%results%e(3) = ptr(idx)%results%e(3) + sum(e_3(1:queued(idx),idx))
                     ptr(idx)%results%pot  = ptr(idx)%results%pot  + sum(pot(1:queued(idx),idx))
                     ptr(idx)%work         = ptr(idx)%work + queued(idx) * WORKLOAD_PENALTY_INTERACTION
                  enddo
!                  do idx = 1,size(gpu)
!                     write(*,'(i4,1x,z10,1x,f15.10)') gpu_id, loc(ptr(idx)%results), e_1(1,idx)
!                  enddo
#else
                  !$acc wait
                  do idx = 1,size(gpu)
  !                do idx = gpu_id, gpu_id
                     !$acc parallel loop                                        &
                     !$acc reduction(+: e_1_, e_2_, e_3_, pot_)                 &
                     !$acc present(e_1(:,idx), e_2(:,idx), e_3(:,idx), pot(:,idx))
                     do idx_ = 1, queued(idx)
                        e_1_ = e_1_ + e_1(idx_, idx)
                        e_2_ = e_2_ + e_2(idx_, idx)
                        e_3_ = e_3_ + e_3(idx_, idx)
                        pot_ = pot_ + pot(idx_, idx)
                     enddo
                     !$acc end parallel loop
                     ptr(idx)%results%e(1) = ptr(idx)%results%e(1) + e_1_ 
                     ptr(idx)%results%e(2) = ptr(idx)%results%e(2) + e_2_
                     ptr(idx)%results%e(3) = ptr(idx)%results%e(3) + e_3_
                     ptr(idx)%results%pot  = ptr(idx)%results%pot  + pot_
                     ptr(idx)%work         = ptr(idx)%work + queued(idx) * WORKLOAD_PENALTY_INTERACTION
                  enddo
#endif
               endif
#endif
#else
               acc%acc_queue(tmp_top)%particle%results%e   = acc%acc_queue(tmp_top)%particle%results%e + [e_1_, e_2_, e_3_]
               acc%acc_queue(tmp_top)%particle%results%pot = acc%acc_queue(tmp_top)%particle%results%pot + pot_
               acc%acc_queue(tmp_top)%particle%work        = acc%acc_queue(tmp_top)%particle%work + &
                                                             queued(gpu_id) * WORKLOAD_PENALTY_INTERACTION
#endif
#else
               ! copy particle info to temporal local copy
               call kernel_node(acc%acc_queue(tmp_top)%particle, acc%acc_queue(tmp_top)%eps, acc%acc_queue(tmp_top)%pen)
#endif

               ! kill list
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

         ! flush GPU buffers at end of timestep
         if (atomic_load_int(acc%thread_status) .eq. ACC_THREAD_STATUS_FLUSH) then
            ! got signal to flush, so check if there is data to flush...
#ifdef ONE
            if ( .not. (atomic_load_int(acc%q_top) .ne. atomic_load_int(acc%q_bottom)) .and. &
                 .not. (gpu_id .eq. size(gpu))                                               &
               ) then
               write(*,*) 'flushing GPU - ', gpu_id,' entries, ',sum(queued(1:gpu_id)),' interactions'
#ifdef CPURED
!               !$acc wait
               !$acc update host(e_1, e_2, e_3, pot)
               do idx = 1,gpu_id
                  ptr(idx)%results%e(1) = ptr(idx)%results%e(1) + sum(e_1(1:queued(idx),idx))
                  ptr(idx)%results%e(2) = ptr(idx)%results%e(2) + sum(e_2(1:queued(idx),idx))
                  ptr(idx)%results%e(3) = ptr(idx)%results%e(3) + sum(e_3(1:queued(idx),idx))
                  ptr(idx)%results%pot  = ptr(idx)%results%pot  + sum(pot(1:queued(idx),idx))
                  ptr(idx)%work         = ptr(idx)%work + queued(idx) * WORKLOAD_PENALTY_INTERACTION
               enddo
#else
               !$acc wait
               do idx = 1,gpu_id
                  !$acc parallel loop                                        &
                  !$acc reduction(+: e_1_, e_2_, e_3_, pot_)                 &
                  !$acc present(e_1(:,idx), e_2(:,idx), e_3(:,idx), pot(:,idx))
                  do idx_ = 1, queued(idx)
                     e_1_ = e_1_ + e_1(idx_, idx)
                     e_2_ = e_2_ + e_2(idx_, idx)
                     e_3_ = e_3_ + e_3(idx_, idx)
                     pot_ = pot_ + pot(idx_, idx)
                  enddo
                  !$acc end parallel loop
                  ptr(idx)%results%e(1) = ptr(idx)%results%e(1) + e_1_ 
                  ptr(idx)%results%e(2) = ptr(idx)%results%e(2) + e_2_
                  ptr(idx)%results%e(3) = ptr(idx)%results%e(3) + e_3_
                  ptr(idx)%results%pot  = ptr(idx)%results%pot  + pot_
                  ptr(idx)%work         = ptr(idx)%work + queued(idx) * WORKLOAD_PENALTY_INTERACTION
               enddo
#endif
            endif
#endif
            ! reset queue
            gpu_id = 0
            ! tell others we're ready again...
            call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)
         endif
   
      end do

      ! free GPU memory
#ifdef __OPENACC
      !$acc wait
      !$acc end data
      write(*,*) 'sync ',cudaDeviceSynchronize()
      write(*,*) 'reset ',cudaDeviceReset()
      !call acc_shutdown(acc_device_nvidia)
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

      ! when this is called, the accelerator thread is up and running
      ! (we check its status before we continue from calc_force_prepare)

      do while( atomic_load_int(acc%q_len) .le. 5 )
         ! busy loop while the queue is processed
         ERROR_ON_FAIL(pthreads_sched_yield())
      end do
   
      ! have a critical section to make sure only one thread enters data into queue
      call critical_section_enter(queue_lock)

         ! thread-safe way of reserving storage for our request
         local_queue_bottom = atomic_mod_increment_and_fetch_int(acc%q_bottom, ACC_QUEUE_LENGTH)
   
         if (local_queue_bottom == atomic_load_int(acc%q_top)) then
            ! we should not get here because of the busy-wait right above...
            write(*,*) 'ACC queue exhausted...'
            DEBUG_ERROR(*, "ACC_QUEUE_LENGTH is too small: ", ACC_QUEUE_LENGTH)
         end if
   
         ! the accelerator thread will check validity of the request and will only proceed as soon as the entry is valid -- this actually serializes the requests
         DEBUG_ASSERT(.not. acc%acc_queue(local_queue_bottom)%entry_valid)
   
         ! store link to particle, so we know what forces to update
         acc%acc_queue(local_queue_bottom)%particle = particle
         ! move list information to queue, list info in particle will be deleted by calling subroutine...
         acc%acc_queue(local_queue_bottom)%queued = particle%queued
         acc%acc_queue(local_queue_bottom)%partner => particle%partner
         ! store epsilon and work_penalty to preserve original calc_force functionality
         acc%acc_queue(local_queue_bottom)%eps = eps2
         acc%acc_queue(local_queue_bottom)%pen = penalty
         ! change queue indicator 
         q_tmp = atomic_fetch_and_decrement_int(acc%q_len)
   
         call atomic_write_barrier() ! make sure the above information is actually written before flagging the entry valid by writing the owner
         acc%acc_queue(local_queue_bottom)%entry_valid    = .true.

      call critical_section_leave(queue_lock)

      return
   end subroutine dispatch_list

end module module_accelerator
