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

#define DELTA1 (MAX_IACT_PARTNERS * (1-1))
#define DELTA2 (MAX_IACT_PARTNERS * (2-1))
#define DELTA3 (MAX_IACT_PARTNERS * (3-1))
#define CHARGE (MAX_IACT_PARTNERS * (4-1))
#define DIP1   (MAX_IACT_PARTNERS * (5-1))
#define DIP2   (MAX_IACT_PARTNERS * (6-1))
#define DIP3   (MAX_IACT_PARTNERS * (7-1))
#define QUAD1  (MAX_IACT_PARTNERS * (8-1))
#define QUAD2  (MAX_IACT_PARTNERS * (9-1))
#define QUAD3  (MAX_IACT_PARTNERS * (10-1))
#define XYQUAD (MAX_IACT_PARTNERS * (11-1))
#define YZQUAD (MAX_IACT_PARTNERS * (12-1))
#define ZXQUAD (MAX_IACT_PARTNERS * (13-1))

!>
!> Accelerator module, in this case for a GPU using OmpSs with tasks (for the GPU).
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

   !> Data structures to be fed to the GPU
   type :: mpdelta
      real*8 :: delta1(1:MAX_IACT_PARTNERS)
      real*8 :: delta2(1:MAX_IACT_PARTNERS)
      real*8 :: delta3(1:MAX_IACT_PARTNERS)
      real*8 :: charge(1:MAX_IACT_PARTNERS)
      real*8 :: dip1(1:MAX_IACT_PARTNERS)
      real*8 :: dip2(1:MAX_IACT_PARTNERS)
      real*8 :: dip3(1:MAX_IACT_PARTNERS)
      real*8 :: quad1(1:MAX_IACT_PARTNERS)
      real*8 :: quad2(1:MAX_IACT_PARTNERS)
      real*8 :: quad3(1:MAX_IACT_PARTNERS)
      real*8 :: xyquad(1:MAX_IACT_PARTNERS)
      real*8 :: yzquad(1:MAX_IACT_PARTNERS)
      real*8 :: zxquad(1:MAX_IACT_PARTNERS)
   end type mpdelta

   ! kernel data
   type point
      type(t_particle_results), pointer :: results
      real*8, pointer :: work
   end type point

   integer :: zer(1)

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
      implicit none
      include 'mpif.h'
   
      type(c_ptr) :: acc_loop

      integer :: gpu_id          ! to keep track of streams

      type(t_particle_thread), target :: tmp_particle
      integer :: tmp_top, q_tmp

      ! kernel data
      real*8, dimension(:,:), allocatable :: e_1, e_2, e_3, pot
real*8, dimension(:,:), allocatable :: pot_, e_1_, e_2_, e_3_
      type(point), dimension(GPU_STREAMS) :: ptr

      ! GPU kernel stuff... move to a function again?
      integer :: idx, idx_, queued(GPU_STREAMS), ccc
integer :: lloc
      real*8 :: dist2, eps2, WORKLOAD_PENALTY_INTERACTION

      real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6

#define DBG_ARR 32
real*8 :: debug_data(DBG_ARR)
      ! kernel interface
      interface
         !$omp target device(opencl) ndrange(1, queued, 128) file(ocl_kernel.cl) copy_deps
         !$omp task in(queued, eps2, partner, dummy) inout(pot, e_1, e_2, e_3, debug_data)
         subroutine ocl_gpu_kernel(queued, eps2, partner, pot, e_1, e_2, e_3, dummy, debug_data)
            use module_interaction_specific_types
            implicit none
            integer, value :: queued, dummy
            real*8, value :: eps2
            real*8 :: partner(13*MAX_IACT_PARTNERS)
            real*8 :: pot(MAX_IACT_PARTNERS), e_1(MAX_IACT_PARTNERS), e_2(MAX_IACT_PARTNERS), e_3(MAX_IACT_PARTNERS)
real*8 :: debug_data(DBG_ARR)
         end subroutine ocl_gpu_kernel
      end interface

#ifdef MONITOR
      external :: Extrae_event
! define the Event types
#define COPYBACK_FORCE_NO 66666
#define COPYBACK_NO 66667
#define WORK 66668
#define FILL 66669
#define COPYBACK 66670
#define DISPATCH 66671
#define CRIT_LOCK 66672
#define LIST_LEN 66673
#define WORK_NO 66674
#endif

#ifndef OMPSS_TASKS
#ifndef NO_NANOS
      external :: nanos_admit_current_thread, nanos_expel_current_thread
      call nanos_admit_current_thread()
#endif
#endif

      ! store ID of comm-thread processor
      acc%processor_id = get_my_core()
      call atomic_write_barrier()

      write(*,*) 'GPU thread on core ', acc%processor_id

      ! allocate GPU structures on heap...
      allocate(e_1(MAX_IACT_PARTNERS, GPU_STREAMS)) ! this corresponds to iact+(strm-1)*MAX_IACT_PARTNERS
      allocate(e_2(MAX_IACT_PARTNERS, GPU_STREAMS))
      allocate(e_3(MAX_IACT_PARTNERS, GPU_STREAMS))
      allocate(pot(MAX_IACT_PARTNERS, GPU_STREAMS))
      allocate(pot_(MAX_IACT_PARTNERS, GPU_STREAMS))
      allocate(e_1_(MAX_IACT_PARTNERS, GPU_STREAMS))
      allocate(e_2_(MAX_IACT_PARTNERS, GPU_STREAMS))
      allocate(e_3_(MAX_IACT_PARTNERS, GPU_STREAMS))

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

      ! initialise GPU variables
      gpu_id = 0

      ! signal successfull start
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)

      do while (atomic_load_int(acc%thread_status) .lt. ACC_THREAD_STATUS_STOPPING)
   
         ! check list for available data
         do while (atomic_load_int(acc%q_top) .ne. atomic_load_int(acc%q_bottom))
            tmp_top = mod(atomic_load_int(acc%q_top), ACC_QUEUE_LENGTH) + 1
   
            ! first check whether the entry is actually valid	  
            if (acc%acc_queue(tmp_top)%entry_valid) then
               call atomic_read_barrier() ! make sure that reads of parts of the queue entry occurr in the correct order

               ! find a stream
               gpu_id = mod(gpu_id,GPU_STREAMS) + 1
               if(gpu_id .lt. 0 .or. gpu_id .gt. GPU_STREAMS) write(*,*) 'BUGGER'
               ! wait for the stream in a task
         !$OMP taskwait
               !$OMP target device(smp) copy_deps
               !$OMP task firstprivate(gpu_id, tmp_top, eps2) &
               !$OMP private(tmp_particle, idx, ccc, q_tmp, zer) &
               !$OMP inout(ptr(gpu_id), acc%acc_queue(tmp_top), e_1(1,gpu_id), e_2(1,gpu_id), e_3(1,gpu_id), pot(1,gpu_id), queued(gpu_id)) &
               !$OMP shared(acc)

#ifdef MONITOR
               call Extrae_event(FILL, gpu_id)
#endif
               ! copy list
               tmp_particle = acc%acc_queue(tmp_top)%particle
               tmp_particle%queued = acc%acc_queue(tmp_top)%queued
               tmp_particle%partner => acc%acc_queue(tmp_top)%partner
               eps2 = acc%acc_queue(tmp_top)%eps
               WORKLOAD_PENALTY_INTERACTION = acc%acc_queue(tmp_top)%pen

               ! free list entry in acc_queue
               ! (do this now so that the loop can issue more tasks)
#ifdef MONITOR
               call Extrae_event(CRIT_LOCK, gpu_id)
#endif
               call critical_section_enter(queue_lock)

               nullify(acc%acc_queue(tmp_top)%partner)
   
               call atomic_store_int(acc%q_top, tmp_top)
               q_tmp = atomic_fetch_and_increment_int(acc%q_len)

               call critical_section_leave(queue_lock)
#ifdef MONITOR
               q_tmp = atomic_load_int(acc%q_len)
               zer = q_tmp
               call Extrae_event(LIST_LEN, zer(1))
               zer = 0
               call Extrae_event(CRIT_LOCK, zer(1))
#endif
               
               ! copy data in GPU structures
               ! only get positions of 'individual arrays'
               queued(gpu_id) = tmp_particle%queued
               ptr(gpu_id)%results => tmp_particle%results
               ptr(gpu_id)%work => tmp_particle%work

#ifdef MONITOR
               zer = 0
               call Extrae_event(FILL, zer(1))
#endif
pot(:,gpu_id) = 0.d0
pot_(:,gpu_id) = 0.d0
e_1(:,gpu_id) = 0.d0
e_1_(:,gpu_id) = 0.d0
e_2(:,gpu_id) = 0.d0
e_2_(:,gpu_id) = 0.d0
e_3(:,gpu_id) = 0.d0
e_3_(:,gpu_id) = 0.d0
               ! run (GPU) kernel
#define OCL_KERNEL 1
#ifdef OCL_KERNEL
               call ocl_gpu_kernel(queued(gpu_id), eps2, tmp_particle%partner(:), pot(:,gpu_id), e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), gpu_id, debug_data)
#else
               call smp_kernel(queued(gpu_id), eps2, tmp_particle%partner(:), pot(:,gpu_id), e_1(:,gpu_id), e_2(:,gpu_id), e_3(:,gpu_id), gpu_id)
#endif

               ! wait for the task to finish. a global one will do here since we only 'posted' 1
               write(*,*) pot(1,gpu_id), pot(queued(gpu_id),gpu_id), gpu_id, loc(pot(1,gpu_id)), ' bef'
               !$OMP taskwait
!               write(*,*) debug_data(1:4)
!               write(*,*) tmp_particle%partner(DELTA1+1), tmp_particle%partner(DELTA2+1), tmp_particle%partner(DELTA3+1), tmp_particle%partner(CHARGE+1)
               write(*,*) pot(1,gpu_id), pot(queued(gpu_id),gpu_id), gpu_id, loc(pot(1,gpu_id)), '    GPU'
               call smp_kernel(queued(gpu_id), eps2, tmp_particle%partner(:), pot_(:,gpu_id), e_1_(:,gpu_id), e_2_(:,gpu_id), e_3_(:,gpu_id), gpu_id)
               !$OMP taskwait
               write(*,*) pot_(1,gpu_id), pot_(queued(gpu_id),gpu_id), gpu_id, loc(pot_(1,gpu_id)), '       CPU'
               !diff the data
               lloc = sum(maxloc(pot(:,gpu_id)-pot_(:,gpu_id)))
               write(*,*) queued(gpu_id), lloc, pot(lloc,gpu_id), pot_(lloc,gpu_id), pot(lloc,gpu_id)-pot_(lloc,gpu_id) 
               lloc = sum(maxloc(e_1(:,gpu_id)-e_1_(:,gpu_id)))
               write(*,*) queued(gpu_id), lloc, e_1(lloc,gpu_id), e_1_(lloc,gpu_id), e_1(lloc,gpu_id)-e_1_(lloc,gpu_id) 
               lloc = sum(maxloc(e_2(:,gpu_id)-e_2_(:,gpu_id)))
               write(*,*) queued(gpu_id), lloc, e_2(lloc,gpu_id), e_2_(lloc,gpu_id), e_2(lloc,gpu_id)-e_2_(lloc,gpu_id) 
               lloc = sum(maxloc(e_3(:,gpu_id)-e_3_(:,gpu_id)))
               write(*,*) queued(gpu_id), lloc, e_3(lloc,gpu_id), e_3_(lloc,gpu_id), e_3(lloc,gpu_id)-e_3_(lloc,gpu_id) 
               ! get data from GPU
               ! now free memory
               deallocate(tmp_particle%partner)
               nullify(tmp_particle%partner)
#ifdef MONITOR
               call Extrae_event(COPYBACK_NO, queued(gpu_id))
               call Extrae_event(COPYBACK, gpu_id)
#endif
               ptr(gpu_id)%results%e(1) = ptr(gpu_id)%results%e(1) + sum(e_1(1:queued(gpu_id),gpu_id))
               ptr(gpu_id)%results%e(2) = ptr(gpu_id)%results%e(2) + sum(e_2(1:queued(gpu_id),gpu_id))
               ptr(gpu_id)%results%e(3) = ptr(gpu_id)%results%e(3) + sum(e_3(1:queued(gpu_id),gpu_id))
               ptr(gpu_id)%results%pot  = ptr(gpu_id)%results%pot  + sum(pot(1:queued(gpu_id),gpu_id))
               ptr(gpu_id)%work         = ptr(gpu_id)%work + queued(gpu_id) * WORKLOAD_PENALTY_INTERACTION
#ifdef MONITOR
               zer = 0
               call Extrae_event(COPYBACK, zer(1))
               call Extrae_event(COPYBACK_NO, zer(1))
#endif
               !$OMP end task
   
               ! we have to invalidate this request queue entry.
               ! this shows that we actually processed it and prevents it from accidentially being resent after the queue wrapped around
               acc%acc_queue(tmp_top)%entry_valid = .false.
               ! need to do this here, since looping the queue may be quicker than starting a task on tmp_top

            else
               ! the next entry is not valid (obviously it has not been stored completely until now -> we abort here and try again later
               exit
            end if
         end do

         ! flush GPU buffers at end of timestep
         if (atomic_load_int(acc%thread_status) .eq. ACC_THREAD_STATUS_FLUSH) then
#ifndef NO_NANOS
            ! need to wait until all tasks are finished (should be the case anyway)
            !$OMP taskwait
#endif
            ! reset queue
            gpu_id = 0
            ! tell others we're ready again...
            call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)
         endif

      end do

      ! should not hurt to wait for child-tasks
      !$OMP taskwait

      ! free queues and lock
      call critical_section_deallocate(queue_lock)
      call atomic_deallocate_int(acc%q_top)
      call atomic_deallocate_int(acc%q_bottom)
      call atomic_deallocate_int(acc%q_len)

      ! free GPU structures
      deallocate(e_1)
      deallocate(e_2)
      deallocate(e_3)
      deallocate(pot)
      deallocate(pot_)
      deallocate(e_1_)
      deallocate(e_2_)
      deallocate(e_3_)

      write(*,*) 'GPU thread terminating'
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STOPPED)
   
      acc_loop = c_null_ptr
#ifndef OMPSS_TASKS
#ifndef NO_NANOS
      call nanos_expel_current_thread()
#endif
      ERROR_ON_FAIL(pthreads_exitthread())
#endif
   
   end function acc_loop

   subroutine dispatch_list(particle, eps2, penalty)
      implicit none

      type(t_particle_thread), target, intent(in) :: particle
      real*8, intent(in) :: eps2, penalty
      integer :: local_queue_bottom, gpu_status, q_tmp

      ! when this is called, the accelerator thread is up and running
      ! (we check its status before we continue from calc_force_prepare)

#ifdef MONITOR
      call Extrae_event(DISPATCH, 1)
#endif

      do
         if (atomic_fetch_and_decrement_int(acc%q_len) .gt. 1) then
            ! decrease queue indicator and check if we can add in one got

            ! have a critical section to make sure only one thread enters data into queue
            call critical_section_enter(queue_lock)

            ! thread-safe way of reserving storage for our request
            local_queue_bottom = atomic_mod_increment_and_fetch_int(acc%q_bottom, ACC_QUEUE_LENGTH)
   
            if (local_queue_bottom .eq. atomic_load_int(acc%q_top)) then
               ! we should not get here because of the busy-wait for free entries in the queue...
               write(*,*) 'ACC queue exhausted...'
               DEBUG_ERROR(*, "ACC_QUEUE_LENGTH is too small: ", ACC_QUEUE_LENGTH)
            end if

            ! the accelerator thread will check validity of the request and will only proceed as soon as the entry is valid -- this actually serializes the requests
            DEBUG_ASSERT(.not. acc%acc_queue(local_queue_bottom)%entry_valid)

            ! store link to particle, so we know what forces to update
            acc%acc_queue(local_queue_bottom)%particle = particle
            ! move list information to queue, list info in particle will be deleted by calling subroutine...
            !                                 \_ compute_iact_list will nullify particle%partner
            acc%acc_queue(local_queue_bottom)%queued = particle%queued
            acc%acc_queue(local_queue_bottom)%partner => particle%partner
            ! store epsilon and work_penalty to preserve original calc_force functionality
            acc%acc_queue(local_queue_bottom)%eps = eps2
            acc%acc_queue(local_queue_bottom)%pen = penalty

            call atomic_write_barrier() ! make sure the above information is actually written before flagging the entry valid by writing the owner
            acc%acc_queue(local_queue_bottom)%entry_valid    = .true.

            call critical_section_leave(queue_lock)

#ifdef MONITOR
            q_tmp = atomic_load_int(acc%q_len)
            call Extrae_event(LIST_LEN, q_tmp)
            call Extrae_event(DISPATCH, 0)
#endif

            return

         else
            ! add to queue indicator again
            q_tmp = atomic_fetch_and_increment_int(acc%q_len)

            ! busy loop while the queue is processed
            ERROR_ON_FAIL(pthreads_sched_yield())
         end if
      end do

#ifdef MONITOR
            call Extrae_event(DISPATCH, 0)
#endif

      return

   end subroutine dispatch_list

   !$omp target device(smp) copy_deps
   !$omp task in(queued, eps2, partner, id) inout(pot, e_1, e_2, e_3)
   subroutine smp_kernel(queued, eps2, partner, pot, e_1, e_2, e_3, id)
      use module_interaction_specific_types
      implicit none
      integer, value :: queued, id
      real*8, value :: eps2
      real*8 :: partner(13*MAX_IACT_PARTNERS)
      real*8 :: pot(MAX_IACT_PARTNERS), e_1(MAX_IACT_PARTNERS), e_2(MAX_IACT_PARTNERS), e_3(MAX_IACT_PARTNERS)

      type(mpdelta), dimension(:), allocatable :: gpu
      integer :: idx, gpu_id, zer(1)
      real*8 :: dist2
      real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
#ifdef MONITOR
      external :: Extrae_event
#endif

      gpu_id = 1
      zer = id
      allocate(gpu(gpu_id))

#ifdef MONITOR
               call Extrae_event(FILL, zer(1))
#endif
      do idx = 1, queued
         gpu(gpu_id)%delta1(idx) = partner(DELTA1+idx)
         gpu(gpu_id)%delta2(idx) = partner(DELTA2+idx)
         gpu(gpu_id)%delta3(idx) = partner(DELTA3+idx)
         gpu(gpu_id)%charge(idx) = partner(CHARGE+idx)
         gpu(gpu_id)%dip1(idx)   = partner(DIP1  +idx)
         gpu(gpu_id)%dip2(idx)   = partner(DIP2  +idx)
         gpu(gpu_id)%dip3(idx)   = partner(DIP3  +idx)
         gpu(gpu_id)%quad1(idx)  = partner(QUAD1 +idx)
         gpu(gpu_id)%quad2(idx)  = partner(QUAD2 +idx)
         gpu(gpu_id)%quad3(idx)  = partner(QUAD3 +idx)
         gpu(gpu_id)%xyquad(idx) = partner(XYQUAD+idx)
         gpu(gpu_id)%yzquad(idx) = partner(YZQUAD+idx)
         gpu(gpu_id)%zxquad(idx) = partner(ZXQUAD+idx)

      enddo
#ifdef MONITOR
      zer = 0
      call Extrae_event(FILL, zer(1))

      zer = id
      call Extrae_event(WORK, zer(1))
      zer = queued
      call Extrae_event(WORK_NO, zer(1))
#endif

      do idx = 1, queued

         dx = gpu(gpu_id)%delta1(idx)
         dy = gpu(gpu_id)%delta2(idx)
         dz = gpu(gpu_id)%delta3(idx)
         dx2 = dx*dx
         dy2 = dy*dy
         dz2 = dz*dz
         dx3 = dx*dx2
         dy3 = dy*dy2
         dz3 = dz*dz2

         dist2     =         dx2
         dist2     = dist2 + dy2
         dist2     = dist2 + dz2

         r  = sqrt(dist2+eps2) ! eps2 is added in calling routine to have plummer instead of coulomb here
         rd = 1.d0/r
         rd2 = rd *rd
         rd3 = rd *rd2
         rd5 = rd3*rd2
         rd7 = rd5*rd2

         fd1 = 3.d0*dx2*rd5 - rd3
         fd2 = 3.d0*dy2*rd5 - rd3
         fd3 = 3.d0*dz2*rd5 - rd3
         fd4 = 3.d0*dx*dy*rd5
         fd5 = 3.d0*dy*dz*rd5
         fd6 = 3.d0*dx*dz*rd5

         pot(idx) = gpu(gpu_id)%charge(idx)*rd                                                          &  !  monopole term
            + (dx*gpu(gpu_id)%dip1(idx) + dy*gpu(gpu_id)%dip2(idx) + dz*gpu(gpu_id)%dip3(idx))*rd3           &  !  dipole
            + 0.5d0*(fd1*gpu(gpu_id)%quad1(idx)  + fd2*gpu(gpu_id)%quad2(idx)  + fd3*gpu(gpu_id)%quad3(idx)) &  !  quadrupole
            +        fd4*gpu(gpu_id)%xyquad(idx) + fd5*gpu(gpu_id)%yzquad(idx) + fd6*gpu(gpu_id)%zxquad(idx)

         e_1(idx) = gpu(gpu_id)%charge(idx)*dx*rd3                                                      &  ! monopole term
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

         e_2(idx) = gpu(gpu_id)%charge(idx)*dy*rd3                                                      &
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

         e_3(idx) = gpu(gpu_id)%charge(idx)*dz*rd3                                                      &
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
#ifdef MONITOR
      zer = 0
      call Extrae_event(WORK, zer(1))
      call Extrae_event(WORK_NO, zer(1))
#endif

      deallocate(gpu)

      return

   end subroutine smp_kernel

end module module_accelerator
