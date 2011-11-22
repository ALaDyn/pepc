!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> All stuff concerning timing: timings are contained in a single
!> array. certain entries therein are addressed via integer
!> parameters, eg.   tim(t_allocate) = 0.1234
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module timings
  implicit none

    !> number of timing entries - dont forget to adjust if you add some timing variables
    integer, private, parameter :: numtimings = 53

    ! global timings
    integer, parameter :: t_domains            =  1
    integer, parameter :: t_allocate           =  2
    integer, parameter :: t_build              =  3
    integer, parameter :: t_exchange_branches_pack        =  4
    integer, parameter :: t_exchange_branches_allgatherv  =  5
    integer, parameter :: t_exchange_branches_integrate   =  6
    integer, parameter :: t_restore            =  7
    integer, parameter :: t_walk               =  8
    integer, parameter :: t_walk_local         =  9
    integer, parameter :: t_unused10           = 10
    integer, parameter :: t_deallocate         = 11
    integer, parameter :: t_all                = 12
    integer, parameter :: t_local              = 13
    integer, parameter :: t_exchange_branches  = 14
    integer, parameter :: t_global             = 15
    integer, parameter :: t_lattice            = 16
    ! fields internal
    integer, parameter :: t_fields_begin       = 17
    integer, parameter :: t_fields_tree        = 18
    integer, parameter :: t_unused19           = 19
    integer, parameter :: t_fields_passes      = 20
    integer, parameter :: t_fields_stats       = 21
    integer, parameter :: t_unused22           = 22
    ! tree_domains
    integer, parameter :: t_domains_keys       = 23
    integer, parameter :: t_domains_sort       = 24
    integer, parameter :: t_domains_sort_pure  = 25
    integer, parameter :: t_domains_ship       = 26
    integer, parameter :: t_domains_bound      = 27
    ! tree_allocate
    integer, parameter :: t_unused28           = 28
    ! tree_build
    integer, parameter :: t_build_neigh        = 29
    integer, parameter :: t_build_part         = 30
    integer, parameter :: t_build_byte         = 31
    ! tree_branches
    integer, parameter :: t_branches_find      = 32
    integer, parameter :: t_exchange_branches_admininstrative = 33
    integer, parameter :: t_build_pure         = 34
    ! tree_fill
    integer, parameter :: t_unused_35          = 35
    integer, parameter :: t_unused_36          = 36
    ! tree_props
    integer, parameter :: t_props_leafs        = 37
    integer, parameter :: t_props_twigs        = 38
    integer, parameter :: t_unused39           = 39
    integer, parameter :: t_unused40           = 40
    ! timings for outside fields()
    integer, parameter :: t_tot                = 41
    ! timings for tree_walk_communicator
    integer, parameter :: t_comm_total         = 42
    integer, parameter :: t_comm_recv          = 43
    integer, parameter :: t_comm_sendreqs      = 44
    ! tree_domains, additional timings
    integer, parameter :: t_domains_add_sort   = 45
    integer, parameter :: t_domains_add_pack   = 46
    integer, parameter :: t_domains_add_unpack = 47
    integer, parameter :: t_domains_add_alltoallv = 48

    integer, parameter :: t_unused_49          = 49
    integer, parameter :: t_unused_50          = 50
    integer, parameter :: t_unused_51          = 51
    integer, parameter :: t_unused_52          = 52
    integer, parameter :: t_unused_53          = 53

    !> array for local timings
    real*8, private, dimension(1:numtimings) :: tim = 0.

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Resets a certain timer to zero
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function timer_read(id)
      implicit none
      integer, intent(in) :: id !< the affected timer address
      real*8 :: timer_read

      timer_read = tim(id)

    end function timer_read

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Resets a certain timer to zero
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_reset(id)
      implicit none
      integer, intent(in) :: id !< the affected timer address

      tim(id) = 0.

    end subroutine timer_reset

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Resets all timers to zero
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_reset_all()
      implicit none

      tim(1:numtimings) = 0.

    end subroutine timer_reset_all

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Logs current time, i.e. sets
    !>      tim(id) = MPI_WTIME()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_stamp(id)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: id !< the affected timer address

      tim(id) = MPI_WTIME()

    end subroutine timer_stamp


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Starts a timer, i.e. sets
    !>      tim(id) = MPI_WTIME()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_start(id)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: id !< the affected timer address

      tim(id) = MPI_WTIME()

    end subroutine timer_start


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Stops a timer, i.e. sets
    !>      tim(id) = MPI_WTIME() - tim(id)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_stop(id)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: id !< the affected timer address

      tim(id) = MPI_WTIME() - tim(id)

    end subroutine timer_stop


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Adds the given value to a timer for cumulative measurements
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timer_add(id, val)
      implicit none
      integer, intent(in) :: id !< the affected timer address
      real*8, intent(in) :: val !< value to be added

      tim(id) = tim(id) + val

    end subroutine timer_add


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Outputs the given timing array to a file with the given filename
    !> @todo: add file header
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timings_ToFile(itime, tdata, filename)
      implicit none
      integer, intent(in) :: itime !< current timestep
      character(*), intent(in) :: filename !< output filename
      real*8, dimension(1:numtimings), intent(in) :: tdata !< array with timing data

      open  (60, file=filename,STATUS='UNKNOWN', POSITION = 'APPEND')
      write (60,*) itime, tdata
      close (60)

    end subroutine timings_ToFile




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Outputs all local timing data to timing_XXXX.dat
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timings_LocalOutput(itime)
      use treevars

      implicit none
      integer, intent(in) :: itime !< current timestep
      character(30) :: cfile

      if ( timing_file_debug ) then
         call system("mkdir -p " // "timing")
         write(cfile,'("timing/timing_",i6.6,".dat")') me
         call timings_ToFile(itime, tim, cfile)
      end if

    end subroutine timings_LocalOutput


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Gathers global timing data and outputs to
    !> timing_avg.dat, timing_min.dat, timing_max.dat
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine timings_GatherAndOutput(itime)
      use treevars
      implicit none
      include 'mpif.h'
      integer, intent(in) :: itime !< current timestep
      integer :: ierr

      real*8, dimension(1:numtimings) :: tim_max
      real*8, dimension(1:numtimings) :: tim_avg
      real*8, dimension(1:numtimings) :: tim_min
      real*8, dimension(1:numtimings) :: tim_dev

      call MPI_REDUCE(tim, tim_max, numtimings, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD,ierr);
      call MPI_REDUCE(tim, tim_min, numtimings, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD,ierr);
      call MPI_REDUCE(tim, tim_avg, numtimings, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,ierr);

     if (me==0) then
        tim_avg = tim_avg / num_pe
        tim_dev = tim_max - tim_min
        call timings_ToFile(itime, tim_max, 'timing_max.dat')
        call timings_ToFile(itime, tim_avg, 'timing_avg.dat')
        call timings_ToFile(itime, tim_min, 'timing_min.dat')
        call timings_ToFile(itime, tim_dev, 'timing_dev_abs.dat')
        tim_dev = tim_dev / tim_min
        call timings_ToFile(itime, tim_dev, 'timing_dev_rel.dat')


        write(*,'(a20,f16.10," s")') "t_all = ",       tim(t_all)
        write(*,'(a20,f16.10," s")') "t_tot = ",       tim(t_tot)
     endif

    end subroutine timings_GatherAndOutput

end module timings
