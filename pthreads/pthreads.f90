module pthreads_stuff

  interface
    function get_my_core()
      implicit none
      integer :: get_my_core
    end function
  end interface


  interface

    function pthreads_init(numthreads)
      implicit none
      integer :: pthreads_init
      integer, intent(in), value :: numthreads
    end function

    function pthreads_uninit()
      implicit none
      integer :: pthreads_uninit
    end function
  end interface

  interface
    function pthreads_createthread(id, start_routine, arg)
      use iso_c_binding
      implicit none
      integer :: pthreads_createthread
      integer, intent(in), value :: id
      type( c_funptr ), value :: start_routine
      type( c_ptr ), value :: arg
    end function

    function pthreads_jointhread(id)
      use iso_c_binding
      implicit none
      integer :: pthreads_jointhread
      integer, intent(in), value :: id
    end function

    function pthreads_exitthread()
      use iso_c_binding
      implicit none
      integer :: pthreads_exitthread
    end function

    function pthreads_sched_yield()
      use iso_c_binding
      implicit none
      integer :: pthreads_sched_yield
    end function

  end interface




  interface

    function rwlocks_init(numlocks)
      implicit none
      integer :: rwlocks_init
      integer, intent(in), value :: numlocks
    end function

    function rwlocks_uninit()
      implicit none
      integer :: rwlocks_uninit
    end function

    function rwlocks_wrlock(id)
      implicit none
      integer :: rwlocks_wrlock
      integer, intent(in), value :: id
    end function

    function rwlocks_rdlock(id)
      implicit none
      integer :: rwlocks_rdlock
      integer, intent(in), value :: id
    end function

    function rwlocks_unlock(id)
      implicit none
      integer :: rwlocks_unlock
      integer, intent(in), value :: id
    end function

  end interface



  interface

    function pthreads_conds_init(numconds)
      implicit none
      integer :: pthreads_conds_init
      integer, intent(in), value :: numconds
    end function

    function pthreads_conds_uninit()
      implicit none
      integer :: pthreads_conds_uninit
    end function

    function pthreads_conds_signal(id)
      implicit none
      integer :: pthreads_conds_signal
      integer, intent(in), value :: id
    end function

    function pthreads_conds_broadcast(id)
      implicit none
      integer :: pthreads_conds_broadcast
      integer, intent(in), value :: id
    end function

    function pthreads_conds_wait(id)
      implicit none
      integer :: pthreads_conds_wait
      integer, intent(in), value :: id
    end function

    function pthreads_conds_timedwait(id, microseconds)
      implicit none
      integer :: pthreads_conds_timedwait
      integer, intent(in), value :: id
      integer, intent(in), value :: microseconds
    end function

    function pthreads_conds_mutex_lock(id)
      implicit none
      integer :: pthreads_conds_mutex_lock
      integer, intent(in), value :: id
    end function

    function pthreads_conds_mutex_unlock(id)
      implicit none
      integer :: pthreads_conds_mutex_unlock
      integer, intent(in), value :: id
    end function

     function pthreads_conds_mutex_timedlock(id, microseconds)
      implicit none
      integer :: pthreads_conds_mutex_timedlock
      integer, intent(in), value :: id
      integer, intent(in), value :: microseconds
    end function

    function pthreads_nanosleep(microseconds)
      implicit none
      integer :: pthreads_nanosleep
      integer, intent(in), value :: microseconds
    end function

   end interface



  interface

    function get_my_tid()
      implicit none
      integer :: get_my_tid
    end function

    function get_my_pid()
      implicit none
      integer :: get_my_pid
    end function

  end interface


  contains

    function getfullid()
      implicit none
      character(20) :: getfullid

      write(getfullid,'("{", I8, ".", I8, "}")') get_my_pid(), get_my_tid()
   end function

end module pthreads_stuff