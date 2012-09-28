! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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

module pthreads_stuff
  use, intrinsic :: iso_c_binding
  implicit none

  integer(kind = c_int), bind(C, name='RWLOCKS_BUSY') :: RWLOCKS_BUSY
  integer(kind = c_int), bind(C, name='BARRIER_THEONE') :: BARRIER_THEONE

  type t_barrier
    type(c_ptr) :: p
  end type t_barrier

  interface
    integer(c_int) function get_my_core() bind(C, name='get_my_core')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface


  interface

    integer(c_int) function pthreads_init(numthreads) bind(C, name='pthreads_init')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: numthreads
    end function

     integer(c_int) function pthreads_uninit() bind(C, name='pthreads_uninit')
      use, intrinsic :: iso_c_binding
      implicit none
    end function
  end interface

  interface
    integer(c_int) function pthreads_createthread(id, start_routine, arg) bind(C, name='pthreads_createthread')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int),   intent(in), value :: id
      type( c_funptr ), intent(in), value :: start_routine
      type( c_ptr ),    intent(in), value :: arg
    end function

    integer(c_int) function pthreads_jointhread(id) bind(C, name='pthreads_jointhread')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function pthreads_exitthread() bind(C, name='pthreads_exitthread')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function pthreads_sched_yield() bind(C, name='pthreads_sched_yield')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

  end interface




  interface

    integer(c_int) function rwlocks_init(numlocks) bind(C, name='rwlocks_init')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: numlocks
    end function

    integer(c_int) function rwlocks_uninit() bind(C, name='rwlocks_uninit')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function rwlocks_wrlock(id) bind(C, name='rwlocks_wrlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function rwlocks_trywrlock(id) bind(C, name='rwlocks_trywrlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function rwlocks_rdlock(id) bind(C, name='rwlocks_rdlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function rwlocks_tryrdlock(id) bind(C, name='rwlocks_tryrdlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

    integer(c_int) function rwlocks_unlock(id) bind(C, name='rwlocks_unlock')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: id
    end function

  end interface


  interface

    type(c_ptr) function c_barrier_alloc() bind(C, name='_barrier_alloc')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function c_barrier_init(barrier, numthreads) bind(C, name='_barrier_init')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: barrier
      integer(c_int), intent(in), value :: numthreads
    end function

    integer(c_int) function c_barrier_destroy_and_free(barrier) bind(C, name='_barrier_destroy_and_free')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: barrier
    end function

    subroutine c_barrier_free(barrier) bind(C, name='_barrier_free')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: barrier
    end subroutine

    integer(c_int) function c_barrier_wait(barrier) bind(C, name='_barrier_wait')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: barrier
    end function

  end interface


  interface

    integer(c_int) function get_my_tid() bind(C, name='get_my_tid')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

    integer(c_int) function get_my_pid() bind(C, name='get_my_pid')
      use, intrinsic :: iso_c_binding
      implicit none
    end function

  end interface


  contains

    function getfullid()
      implicit none
      character(20) :: getfullid

      write(getfullid,'("{", I8, ".", I8, "}")') get_my_pid(), get_my_tid()
   end function

  subroutine barrier_allocate(barrier)
    implicit none

    type(t_barrier), pointer, intent(out) :: barrier

    type(t_barrier), pointer :: tmp_f
    type(c_ptr) :: tmp_c

    allocate(tmp_f)
    tmp_c = c_barrier_alloc()

    if (.not. (associated(tmp_f) .and. c_associated(tmp_c))) then
      call c_barrier_free(tmp_c)
      deallocate(tmp_f)
      barrier => null()
    else
      tmp_f%p = tmp_c
      barrier => tmp_f
    end if
  end subroutine

  integer function barrier_init(barrier, numthreads)
    implicit none

    type(t_barrier), intent(in) :: barrier
    integer, intent(in) :: numthreads

    barrier_init = c_barrier_init(barrier%p, numthreads)
  end function

  integer function barrier_destroy_and_deallocate(barrier)
    implicit none

    type(t_barrier), pointer, intent(inout) :: barrier

    integer :: ret

    ret = c_barrier_destroy_and_free(barrier%p)
    
    if (ret == 0) then
      deallocate(barrier)
      barrier => null()
    end if

    barrier_destroy_and_deallocate = ret
  end function

  integer function barrier_wait(barrier)
    implicit none

    type(t_barrier), intent(in) :: barrier

    barrier_wait = c_barrier_wait(barrier%p)
  end function

end module pthreads_stuff
