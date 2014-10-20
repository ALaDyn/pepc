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

#define POT    (MAX_IACT_PARTNERS * (1-1))
#define E_1    (MAX_IACT_PARTNERS * (2-1))
#define E_2    (MAX_IACT_PARTNERS * (3-1))
#define E_3    (MAX_IACT_PARTNERS * (4-1))

#define BLN 128
! parameter for GPU kernel size

! check we have picked a kernel type
#if !defined ( OCL_KERNEL ) && !defined ( SMP_ALTERNATIVE ) && !defined ( SMP_KERNEL )
#error DEFINE AT LEAST ONE KERNEL
#endif

!>
!> Accelerator module, in this case for a GPU using OmpSs with tasks (for the GPU).
!>
module module_accelerator
   use module_interaction_specific_types
   use module_atomic_ops
   use treevars, only: num_threads

   implicit none

   integer, public, parameter :: ACC_THREAD_STATUS_WAITING  = 1
   integer, public, parameter :: ACC_THREAD_STATUS_STARTING = 2
   integer, public, parameter :: ACC_THREAD_STATUS_STARTED  = 3
   integer, public, parameter :: ACC_THREAD_STATUS_FLUSH    = 4
   integer, public, parameter :: ACC_THREAD_STATUS_STOPPING = 5
   integer, public, parameter :: ACC_THREAD_STATUS_STOPPED  = 6

   type(t_critical_section), pointer :: queue_lock

   !> Data structures to be fed to the GPU
   type mpdelta
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

   ! kernel data, ACC thread local data, etc
   type t_acc_l
      integer :: acc_id                                                   ! to keep track of used stream
      type(t_atomic_int), pointer :: q_len                                ! how many streams are available (max = ACC_QUEUE_LENGTH)
      type(t_acc_queue_entry) :: qentry(ACC_QUEUE_LENGTH)                 ! stores interaction data per stream, including results
      real*8, dimension(4*MAX_IACT_PARTNERS, ACC_QUEUE_LENGTH) :: results ! store intermediate results within task per stream
                                                                          ! we have 4 results (pot, E_) per iact partner
   end type t_acc_l
   type(t_acc_l), allocatable :: tl_acc(:)

   integer :: zer(1)

   contains

   subroutine acc_start()
      ! this will initialise and allocate ACC (thread local) fields
      ! called after we know how many worker tasks/threads we will have during the tree walk
      use, intrinsic :: iso_c_binding
      use pthreads_stuff, only: pthreads_sched_yield, get_my_core, pthreads_exitthread
      use module_debug
      implicit none
      include 'mpif.h'

      integer :: i

      ! signal we're starting...
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTING)
      write(*,*) '== initialising ACC structures'

      ! allocate ACC/GPU structures (on heap...)
      allocate(tl_acc(num_threads))

      do i = 1, num_threads
         ! reset all queues
         call atomic_allocate_int(tl_acc(i)%q_len)
         call atomic_store_int(tl_acc(i)%q_len, ACC_QUEUE_LENGTH)
         ! initialise ACC/GPU variables
         tl_acc(i)%acc_id = 0
      enddo

      ! init critical lock
      call critical_section_allocate(queue_lock)

      ! signal we're done
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)

   end subroutine acc_start

   subroutine acc_flush()
      ! this is a tricky one - we have to wait for all ACC tasks to be finished...
      implicit none

      integer :: i, j

      ! flush GPU buffers at end of timestep

      ! need to wait until all tasks are finished (should be the case anyway)
      do j = 1, num_threads
         do i = 1, ACC_QUEUE_LENGTH
            ! wait for kernels
            !$OMP taskwait on (tl_acc(j)%results(1,i))
         enddo
      enddo

      ! tell others we're ready again...
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STARTED)

   end subroutine acc_flush

   subroutine acc_stop()
      ! get rid of all ACC fields cleanly and signal we're done
      use, intrinsic :: iso_c_binding
      use pthreads_stuff, only: pthreads_sched_yield, get_my_core, pthreads_exitthread
      use module_debug
      implicit none
      include 'mpif.h'
   
      ! free queues and lock
      call critical_section_deallocate(queue_lock)

      ! free ACC/GPU structures
      deallocate(tl_acc)

      write(*,*) '== removing ACC structures'
      call atomic_store_int(acc%thread_status, ACC_THREAD_STATUS_STOPPED)
   
   end subroutine acc_stop

   subroutine dispatch_list(particle, eps2, penalty)
      implicit none

      type(t_particle_thread), target, intent(in) :: particle
      real*8, intent(in) :: eps2, penalty
      integer :: local_queue_bottom, gpu_status, q_tmp, acc_id, t_id
      real*8 :: dist2, WORKLOAD_PENALTY_INTERACTION

      ! kernel interface
#ifdef OCL_KERNEL
      interface
         !$omp target device(opencl) ndrange(1, queued, BLN) file(ocl_kernel.cl) copy_deps
         !$omp task in(queued, eps2, id) inout(partner, results) label(OCL-kernel)
         subroutine ocl_gpu_kernel(queued, eps2, partner, results, id)
            use module_interaction_specific_types
            implicit none
            integer, value :: queued, id
            real*8, value :: eps2
            real*8 :: partner(13*MAX_IACT_PARTNERS)
            real*8 :: results(4*((queued-1)/BLN + 1 ))
         end subroutine ocl_gpu_kernel
      end interface
#endif
#ifdef SMP_ALTERNATIVE
      interface
#ifdef OCL_KERNEL
         !$omp target device(smp) implements(ocl_gpu_kernel) copy_deps
#else
         !$omp target device(smp) copy_deps
#endif
         !$omp task in(queued, eps2, id) inout(partner, results) label(SMP-alternative)
         subroutine ocl_smp_kernel(queued, eps2, partner, results, id)
            use module_interaction_specific_types
            implicit none
            integer, value :: queued, id
            real*8, value :: eps2
            real*8 :: partner(13*MAX_IACT_PARTNERS)
            real*8 :: results(4*((queued-1)/BLN + 1 ))
         end subroutine ocl_smp_kernel
      end interface
#endif

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

      call Extrae_event(DISPATCH, 1)
#endif

      ! get the id of the thread we're on, to keep all data thread- (i.e. task-) local
      ! (thanks to being able to track global dependencies...)
      t_id = particle%thread_id

      do
         q_tmp = atomic_fetch_and_decrement_int(tl_acc(t_id)%q_len) 
         if (q_tmp .gt. 0) then
            ! decrease queue indicator and check if we can add in one task

            ! get a stream id
            tl_acc(t_id)%acc_id = mod(tl_acc(t_id)%acc_id,ACC_QUEUE_LENGTH) + 1
            ! take a shorthand copy
            acc_id = tl_acc(t_id)%acc_id

#ifdef MONITOR
            call Extrae_event(FILL, acc_id)
#endif

            ! copy subroutine data to our thread-/task-local ACC queue to keep it alive
            tl_acc(t_id)%qentry(acc_id)%particle = particle
            ! move list information to queue, list info in particle will be deleted by calling subroutine...
            !                                 \_ compute_iact_list will nullify particle%partner
            tl_acc(t_id)%qentry(acc_id)%queued = particle%queued
            tl_acc(t_id)%qentry(acc_id)%partner => particle%partner
            ! store epsilon and work_penalty to preserve original calc_force functionality
            tl_acc(t_id)%qentry(acc_id)%eps = eps2
            tl_acc(t_id)%qentry(acc_id)%pen = penalty

            WORKLOAD_PENALTY_INTERACTION = tl_acc(t_id)%qentry(acc_id)%pen

#ifdef MONITOR
            q_tmp = atomic_load_int(tl_acc%q_len(t_id))
            zer = q_tmp
            call Extrae_event(LIST_LEN, zer(1))
            zer = 0
            call Extrae_event(FILL, zer(1))
#endif

            ! clean results buffer from last iteration - possibly drop this useless line?
            t_acc(t_id)%results(:,acc_id) = 0.d0

            ! run (GPU) kernel
write(*,'(a,i3,2(" 0x",z16.16))') 'submitting ', acc_id, loc(tl_acc(t_id)%qentry(acc_id)%particle%partner), loc(tl_acc(t_id)%results(1, acc_id))
#if defined OCL_KERNEL || defined SMP_ALTERNATIVE
#ifdef OCL_KERNEL
            call ocl_gpu_kernel(tl_acc(t_id)%qentry(acc_id)%queued, eps2, tl_acc(t_id)%qentry(acc_id)%particle%partner(:), tl_acc(t_id)%results(1:4*( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ),acc_id), acc_id)
            !                                                                                                                |------- no of blocks -------|
#else
            call ocl_smp_kernel(tl_acc(t_id)%qentry(acc_id)%queued, eps2, tl_acc(t_id)%qentry(acc_id)%particle%partner(:), tl_acc(t_id)%results(1:4*( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ), acc_id), acc_id)
#endif
#endif
#ifdef SMP_KERNEL
            call smp_kernel(tl_acc(t_id)%qentry(acc_id)%queued, eps2, tl_acc(t_id)%qentry(acc_id)%particle%partner(:), tl_acc(t_id)%results(:, acc_id), acc_id)
#endif

!$OMP taskwait
            ! wait for the task to finish. a global one will do here since we only 'posted' 1
!$$!            !$OMP target device(SMP) copy_deps
!$$!            !$OMP task firstprivate(t_id, acc_id) private(zer) &
!$$!            !$OMP inout(tl_acc(t_id)%results(1:4*( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ), acc_id), tl_acc(t_id)%qentry(acc_id)%queued) &
!$$!            !$OMP inout(tl_acc(t_id)%qentry(acc_id), tl_acc(t_id)%qentry(acc_id)%particle%partner) &
!$$!            !$OMP label(reduction)
            ! have a task to post-process data
            ! get data from GPU
#ifdef MONITOR
            call Extrae_event(COPYBACK_NO, tl_acc(t_id)%qentry(acc_id)%queued)
            call Extrae_event(COPYBACK, acc_id)
#endif

#ifdef SMP_KERNEL
            tl_acc(t_id)%qentry(acc_id)%particle%results%e(1) = tl_acc(t_id)%qentry(acc_id)%particle%results%e(1) + sum(tl_acc(t_id)%results((E_1+1):(E_1+tl_acc(t_id)%qentry(acc_id)%queued), acc_id))
            tl_acc(t_id)%qentry(acc_id)%particle%results%e(2) = tl_acc(t_id)%qentry(acc_id)%particle%results%e(2) + sum(tl_acc(t_id)%results((E_2+1):(E_2+tl_acc(t_id)%qentry(acc_id)%queued), acc_id))
            tl_acc(t_id)%qentry(acc_id)%particle%results%e(3) = tl_acc(t_id)%qentry(acc_id)%particle%results%e(3) + sum(tl_acc(t_id)%results((E_3+1):(E_3+tl_acc(t_id)%qentry(acc_id)%queued), acc_id))
            tl_acc(t_id)%qentry(acc_id)%particle%results%pot  = tl_acc(t_id)%qentry(acc_id)%particle%results%pot  + sum(tl_acc(t_id)%results((POT+1):(POT+tl_acc(t_id)%qentry(acc_id)%queued), acc_id))
#else
            tl_acc(t_id)%qentry(acc_id)%particle%results%e(1) = tl_acc(t_id)%qentry(acc_id)%particle%results%e(1) + sum(tl_acc(t_id)%results(((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (2-1))+1):((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (2-1))+( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 )), acc_id))
            tl_acc(t_id)%qentry(acc_id)%particle%results%e(2) = tl_acc(t_id)%qentry(acc_id)%particle%results%e(2) + sum(tl_acc(t_id)%results(((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (3-1))+1):((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (3-1))+( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 )), acc_id))
            tl_acc(t_id)%qentry(acc_id)%particle%results%e(3) = tl_acc(t_id)%qentry(acc_id)%particle%results%e(3) + sum(tl_acc(t_id)%results(((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (4-1))+1):((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (4-1))+( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 )), acc_id))
            tl_acc(t_id)%qentry(acc_id)%particle%results%pot  = tl_acc(t_id)%qentry(acc_id)%particle%results%pot  + sum(tl_acc(t_id)%results(((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (1-1))+1):((((tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 ) * (1-1))+( (tl_acc(t_id)%qentry(acc_id)%queued-1)/BLN + 1 )), acc_id))
#endif
            tl_acc(t_id)%qentry(acc_id)%particle%work         = tl_acc(t_id)%qentry(acc_id)%particle%work + tl_acc(t_id)%qentry(acc_id)%queued * WORKLOAD_PENALTY_INTERACTION
#ifdef MONITOR
            zer = 0
            call Extrae_event(COPYBACK, zer(1))
            call Extrae_event(COPYBACK_NO, zer(1))
#endif
            ! now free memory
write(*,'(a,i3,2(" 0x",z16.16))') 'dealloc          ', acc_id, loc(tl_acc(t_id)%qentry(acc_id)%particle%partner), loc(tl_acc(t_id)%results(1, acc_id))
            deallocate(tl_acc(t_id)%qentry(acc_id)%particle%partner)
            nullify(tl_acc(t_id)%qentry(acc_id)%particle%partner)

            ! make room for other tasks
            q_tmp = atomic_fetch_and_increment_int(tl_acc(t_id)%q_len) 
!$$!            !$OMP end task

#ifdef MONITOR
            call Extrae_event(LIST_LEN, q_tmp)
            call Extrae_event(DISPATCH, 0)
#endif

            return

         else
            ! add to queue indicator again
            q_tmp = atomic_fetch_and_increment_int(tl_acc(t_id)%q_len) 

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
   !$omp task in(queued, eps2, partner, id) out(results) label(SMP-kernel)
   subroutine smp_kernel(queued, eps2, partner, results, id)
      use module_interaction_specific_types
      implicit none
      integer, value :: queued, id
      real*8, value :: eps2
      real*8 :: partner(13*MAX_IACT_PARTNERS)
      real*8 :: results(4*MAX_IACT_PARTNERS)

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

         results(POT+idx) = gpu(gpu_id)%charge(idx)*rd                                                       &  !  monopole term
            + (dx*gpu(gpu_id)%dip1(idx) + dy*gpu(gpu_id)%dip2(idx) + dz*gpu(gpu_id)%dip3(idx))*rd3           &  !  dipole
            + 0.5d0*(fd1*gpu(gpu_id)%quad1(idx)  + fd2*gpu(gpu_id)%quad2(idx)  + fd3*gpu(gpu_id)%quad3(idx)) &  !  quadrupole
            +        fd4*gpu(gpu_id)%xyquad(idx) + fd5*gpu(gpu_id)%yzquad(idx) + fd6*gpu(gpu_id)%zxquad(idx)

         results(E_1+idx) = gpu(gpu_id)%charge(idx)*dx*rd3                                                   &  ! monopole term
            + fd1*gpu(gpu_id)%dip1(idx) + fd4*gpu(gpu_id)%dip2(idx) + fd6*gpu(gpu_id)%dip3(idx)              &  ! dipole term
            + 3.d0   * (                                                                                     &  ! quadrupole term
            0.5d0 * (                                                                                        &
            ( 5.d0*dx3   *rd7 - 3.d0*dx*rd5 )*gpu(gpu_id)%quad1(idx)                                         &
            + ( 5.d0*dx*dy2*rd7 -      dx*rd5 )*gpu(gpu_id)%quad2(idx)                                       &
            + ( 5.d0*dx*dz2*rd7 -      dx*rd5 )*gpu(gpu_id)%quad3(idx)                                       &
            )                                                                                                &
            + ( 5.d0*dy*dx2  *rd7 - dy*rd5 )*gpu(gpu_id)%xyquad(idx)                                         &
            + ( 5.d0*dz*dx2  *rd7 - dz*rd5 )*gpu(gpu_id)%zxquad(idx)                                         &
            + ( 5.d0*dx*dy*dz*rd7          )*gpu(gpu_id)%yzquad(idx)                                         &
            )

         results(E_2+idx) = gpu(gpu_id)%charge(idx)*dy*rd3                                                   &
            + fd2*gpu(gpu_id)%dip2(idx) + fd4*gpu(gpu_id)%dip1(idx) + fd5*gpu(gpu_id)%dip3(idx)              &
            + 3 * (                                                                                          &
            0.5d0 * (                                                                                        &
            ( 5*dy3*rd7    - 3*dy*rd5 )*gpu(gpu_id)%quad2(idx)                                               &
            + ( 5*dy*dx2*rd7 -   dy*rd5 )*gpu(gpu_id)%quad1(idx)                                             &
            + ( 5*dy*dz2*rd7 -   dy*rd5 )*gpu(gpu_id)%quad3(idx)                                             &
            )                                                                                                &
            + ( 5*dx*dy2  *rd7 - dx*rd5 )*gpu(gpu_id)%xyquad(idx)                                            &
            + ( 5*dz*dy2  *rd7 - dz*rd5 )*gpu(gpu_id)%yzquad(idx)                                            &
            + ( 5*dx*dy*dz*rd7          )*gpu(gpu_id)%zxquad(idx)                                            &
            )

         results(E_3+idx) = gpu(gpu_id)%charge(idx)*dz*rd3                                                   &
            + fd3*gpu(gpu_id)%dip3(idx) + fd5*gpu(gpu_id)%dip2(idx) + fd6*gpu(gpu_id)%dip1(idx)              &
            + 3 * (                                                                                          &
            0.5d0 * (                                                                                        &
            + ( 5*dz3   *rd7 - 3*dz*rd5 )*gpu(gpu_id)%quad3(idx)                                             &
            + ( 5*dz*dy2*rd7 -   dz*rd5 )*gpu(gpu_id)%quad2(idx)                                             &
            + ( 5*dz*dx2*rd7 -   dz*rd5 )*gpu(gpu_id)%quad1(idx)                                             &
            )                                                                                                &
            + ( 5*dx*dz2  *rd7 - dx*rd5 )*gpu(gpu_id)%zxquad(idx)                                            &
            + ( 5*dy*dz2  *rd7 - dy*rd5 )*gpu(gpu_id)%yzquad(idx)                                            &
            + ( 5*dx*dy*dz*rd7          )*gpu(gpu_id)%xyquad(idx)                                            &
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

subroutine ocl_smp_kernel(queued, eps2, partner, results, id)
   use module_accelerator, only: mpdelta
   use module_interaction_specific_types
   implicit none
   integer, value :: queued, id
   real*8, value :: eps2
   real*8 :: partner(13*MAX_IACT_PARTNERS)
   real*8 :: results(4*((queued-1)/BLN + 1 ))

   real*8, allocatable :: l_results(:)
   type(mpdelta), dimension(:), allocatable :: gpu
   integer :: idx, gpu_id, zer(1)
   real*8 :: dist2
   real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
#ifdef MONITOR
   external :: Extrae_event
#endif

   ! allocate local temp workspace...
   allocate(l_results(4*MAX_IACT_PARTNERS))
   l_results = 0.d0

   gpu_id = 1
   zer = id
   allocate(gpu(gpu_id))

write(*,'(a," 0x",z16.16)') 'working on     ', loc(partner)
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

      l_results(POT+idx) = gpu(gpu_id)%charge(idx)*rd                                                     &  !  monopole term
         + (dx*gpu(gpu_id)%dip1(idx) + dy*gpu(gpu_id)%dip2(idx) + dz*gpu(gpu_id)%dip3(idx))*rd3           &  !  dipole
         + 0.5d0*(fd1*gpu(gpu_id)%quad1(idx)  + fd2*gpu(gpu_id)%quad2(idx)  + fd3*gpu(gpu_id)%quad3(idx)) &  !  quadrupole
         +        fd4*gpu(gpu_id)%xyquad(idx) + fd5*gpu(gpu_id)%yzquad(idx) + fd6*gpu(gpu_id)%zxquad(idx)

      l_results(E_1+idx) = gpu(gpu_id)%charge(idx)*dx*rd3                                                 &  ! monopole term
         + fd1*gpu(gpu_id)%dip1(idx) + fd4*gpu(gpu_id)%dip2(idx) + fd6*gpu(gpu_id)%dip3(idx)              &  ! dipole term
         + 3.d0   * (                                                                                     &  ! quadrupole term
         0.5d0 * (                                                                                        &
         ( 5.d0*dx3   *rd7 - 3.d0*dx*rd5 )*gpu(gpu_id)%quad1(idx)                                         &
         + ( 5.d0*dx*dy2*rd7 -      dx*rd5 )*gpu(gpu_id)%quad2(idx)                                       &
         + ( 5.d0*dx*dz2*rd7 -      dx*rd5 )*gpu(gpu_id)%quad3(idx)                                       &
         )                                                                                                &
         + ( 5.d0*dy*dx2  *rd7 - dy*rd5 )*gpu(gpu_id)%xyquad(idx)                                         &
         + ( 5.d0*dz*dx2  *rd7 - dz*rd5 )*gpu(gpu_id)%zxquad(idx)                                         &
         + ( 5.d0*dx*dy*dz*rd7          )*gpu(gpu_id)%yzquad(idx)                                         &
         )

      l_results(E_2+idx) = gpu(gpu_id)%charge(idx)*dy*rd3                                                 &
         + fd2*gpu(gpu_id)%dip2(idx) + fd4*gpu(gpu_id)%dip1(idx) + fd5*gpu(gpu_id)%dip3(idx)              &
         + 3 * (                                                                                          &
         0.5d0 * (                                                                                        &
         ( 5*dy3*rd7    - 3*dy*rd5 )*gpu(gpu_id)%quad2(idx)                                               &
         + ( 5*dy*dx2*rd7 -   dy*rd5 )*gpu(gpu_id)%quad1(idx)                                             &
         + ( 5*dy*dz2*rd7 -   dy*rd5 )*gpu(gpu_id)%quad3(idx)                                             &
         )                                                                                                &
         + ( 5*dx*dy2  *rd7 - dx*rd5 )*gpu(gpu_id)%xyquad(idx)                                            &
         + ( 5*dz*dy2  *rd7 - dz*rd5 )*gpu(gpu_id)%yzquad(idx)                                            &
         + ( 5*dx*dy*dz*rd7          )*gpu(gpu_id)%zxquad(idx)                                            &
         )

      l_results(E_3+idx) = gpu(gpu_id)%charge(idx)*dz*rd3                                                 &
         + fd3*gpu(gpu_id)%dip3(idx) + fd5*gpu(gpu_id)%dip2(idx) + fd6*gpu(gpu_id)%dip1(idx)              &
         + 3 * (                                                                                          &
         0.5d0 * (                                                                                        &
         + ( 5*dz3   *rd7 - 3*dz*rd5 )*gpu(gpu_id)%quad3(idx)                                             &
         + ( 5*dz*dy2*rd7 -   dz*rd5 )*gpu(gpu_id)%quad2(idx)                                             &
         + ( 5*dz*dx2*rd7 -   dz*rd5 )*gpu(gpu_id)%quad1(idx)                                             &
         )                                                                                                &
         + ( 5*dx*dz2  *rd7 - dx*rd5 )*gpu(gpu_id)%zxquad(idx)                                            &
         + ( 5*dy*dz2  *rd7 - dy*rd5 )*gpu(gpu_id)%yzquad(idx)                                            &
         + ( 5*dx*dy*dz*rd7          )*gpu(gpu_id)%xyquad(idx)                                            &
         )
   end do
#ifdef MONITOR
   zer = 0
   call Extrae_event(WORK, zer(1))
   call Extrae_event(WORK_NO, zer(1))
#endif

   ! now perform the (partial-)reduction the GPU is doing as well... MIND THE SIZE BLN
   do idx = 1, queued, BLN
      results((((queued-1)/BLN + 1 ) * (1-1))+(idx-1)/BLN+1) = sum(l_results(POT+idx:min(POT+idx+BLN-1,POT+queued)))
      results((((queued-1)/BLN + 1 ) * (2-1))+(idx-1)/BLN+1) = sum(l_results(E_1+idx:min(E_1+idx+BLN-1,E_1+queued)))
      results((((queued-1)/BLN + 1 ) * (3-1))+(idx-1)/BLN+1) = sum(l_results(E_2+idx:min(E_2+idx+BLN-1,E_2+queued)))
      results((((queued-1)/BLN + 1 ) * (4-1))+(idx-1)/BLN+1) = sum(l_results(E_3+idx:min(E_3+idx+BLN-1,E_3+queued)))
   enddo

   deallocate(gpu)

   return

end subroutine ocl_smp_kernel
