! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

!>
!>  Encapsulates some debugging and i/o specific routines
!>
module module_debug
   use module_pepc_kinds
   implicit none
   save
   private

   integer, parameter, public :: debug_ipefile = 21
   integer, public :: debug_my_rank = -1
   integer, parameter, public :: debug_stdout = 6

   integer, public :: debug_level = 0 !< or-combination of the bitmasks below
   ! set to the following values to get the old behaviour:
   !             db_level = 0      --> debug_level = 0
   !             db_level = 1      --> debug_level = 1 + 2 + 32 = 35
   !             db_level = 2      --> debug_level = 1 + 2 + 32 + 64 = 99
   !             db_level = 3      --> debug_level = 1 + 2 + 32 + 64 + 8 + 2048 = 2155
   !             db_level = 4      --> debug_level = 1 + 2 + 32 + 64 + 8 + 2048 + 4 + 16 = 2175
   !             db_level = 5      --> debug_level = 128
   !             db_level = 6      --> debug_level = 1 + 2 + 32 + 64 + 8 + 2048 + 4 + 16 + 128 + 256 + 512 = 3071
   !
   integer, parameter, public :: DBG_STATUS      = int(B'0000000000000001')    !&  1
   integer, parameter, public :: DBG_TREE        = int(B'0000000000000010')    !&  2
   integer, parameter, public :: DBG_BUILD       = int(B'0000000000000100')    !&  4
   integer, parameter, public :: DBG_DOMAIN      = int(B'0000000000001000')    !&  8
   integer, parameter, public :: DBG_BRANCH      = int(B'0000000000010000')    !&  16
   integer, parameter, public :: DBG_STATS       = int(B'0000000000100000')    !&  32
   integer, parameter, public :: DBG_WALKSUMMARY = int(B'0000000001000000')    !&  64
   integer, parameter, public :: DBG_DUMPTREE    = int(B'0000000010000000')    !&  128
   integer, parameter, public :: DBG_TIMINGFILE  = int(B'0000000100000000')    !&  256
   ! deprecated: integer, parameter, public :: DBG_LOADFILE    = int(B'0000001000000000')    ! 512
   integer, parameter, public :: DBG_WALK        = int(B'0000010000000000')    !&  1024
   integer, parameter, public :: DBG_PERIODIC    = int(B'0000100000000000')    !&  2048

   character(30), private :: debug_ipefile_name

   public dbg
   public pepc_status
   public debug_mpi_abort
   public debug_barrier
   public debug_print_timestamp
   public debug_initialize
   public debug_finalize

contains

   !>
   !>  lpepc status output
   !>
   subroutine pepc_status(stat)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use treevars, only: me
      implicit none
      character(*), intent(in) :: stat

      if (dbg(DBG_STATUS)) then
         call debug_print_timestamp(debug_ipefile)
         write (debug_ipefile, '(" LPEPC | ", a)') stat
         flush (debug_ipefile)

         if (me .eq. 0) then
            call debug_print_timestamp(output_unit)
            write (*, '(" LPEPC | ", a)') stat
         end if
      end if
   end subroutine

   !>
   !>  module initialization (is automatically called when PEPC is initialized)
   !>
   subroutine debug_initialize()
      use treevars
      use module_utils, only: create_directory
      use mpi
      implicit none

      character(MPI_MAX_PROCESSOR_NAME) :: procname
      integer(kind_default) :: resultlen, ierr

      debug_my_rank = me
      call MPI_GET_PROCESSOR_NAME(procname, resultlen, ierr)

      write (debug_ipefile_name, '("diag/diag_",i6.6,".dat")') me
      call create_directory("diag")

      open (debug_ipefile, file=trim(debug_ipefile_name), STATUS='UNKNOWN', POSITION='REWIND')
      call debug_print_timestamp(debug_ipefile)
      write (debug_ipefile, '(a)') " PEPC on ["//procname(1:resultlen)//"]"
      flush (debug_ipefile)
   end subroutine debug_initialize

   !>
   !> module finalization (is automatically called when PEPC is finalized)
   !>
   subroutine debug_finalize()
      implicit none

      close (debug_ipefile)
   end subroutine

   !>
   !> Writes the current date and time into the argument in the format YYYY-MM-DD hh:mm:ss
   !>
   subroutine debug_timestamp(str)
      implicit none
      character(*), intent(out) :: str

      ! [
      !  year,
      !  month of year,
      !  day of month,
      !  difference to UTC in minutes,
      !  hour of day,
      !  minutes of hour,
      !  seconds of minute,
      !  milliseconds of second
      ! ]
      integer :: v(8)

      if (len(str) .ge. 19) then
         call date_and_time(values=v)
         write (str, '(i4, "-", i2.2, "-", i2.2, " ", i2.2, ":", i2.2, ":", i2.2)') v(1), v(2), v(3), v(5), v(6), v(7)
      end if
   end subroutine debug_timestamp

   !>
   !> Writes the output of `debug_timestamp` surrounded by angled brackets to the I/O unit specified in the argument
   !> without advancing to the next line.
   !>
   subroutine debug_print_timestamp(iounit)
      implicit none

      integer(kind_default), intent(in) :: iounit

      character(19) :: str

      call debug_timestamp(str)

      write (iounit, '(".lt.", a, ".gt.")', advance='no') str
   end subroutine debug_print_timestamp

   !>
   !>  calls MPI_ABORT(MPI_COMM_lpepc, 1, ierr)
   !>
   subroutine debug_mpi_abort()
      use treevars, only: MPI_COMM_lpepc
#if defined(__ICC) || defined(__INTEL_COMPILER)
      use ifcore
#endif
      use mpi
      implicit none
      integer(kind_default) :: ierr

#if defined(__ICC) || defined(__INTEL_COMPILER)
      ! http://software.intel.com/sites/products/documentation/hpc/composerxe/en-us/2011Update/fortran/lin/lref_for/source_files/rftrace.htm
      call tracebackqq("stacktrace", -1)
#elif defined(__IBMC__) || defined(__IBMCPP__)
      ! http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp?topic=%2Fcom.ibm.xlf101a.doc%2Fxlflr%2Fsup-xltrbk.htm
      call xl__trbk()
#elif defined(__GNUC__) && !defined(__PGI)
      ! starting from GCC version 4.8, a backtrace() subroutine is provided by gfortran
#define GCC_VERSION (__GNUC__ * 10000 \
      +__GNUC_MINOR__ * 100\
      +__GNUC_PATCHLEVEL__)
#if GCC_VERSION >= 40800
      ! http://gcc.gnu.org/onlinedocs/gfortran/BACKTRACE.html
      call backtrace()
#endif
#endif

      call MPI_ABORT(MPI_COMM_lpepc, 1, ierr)
   end subroutine

   !>
   !>  calls MPI_BARRIER(MPI_COMM_lpepc, ierr)
   !>
   subroutine debug_barrier()
      use treevars, only: MPI_COMM_lpepc
      use mpi
      implicit none
      integer(kind_default) :: ierr

      call MPI_BARRIER(MPI_COMM_lpepc, ierr)
   end subroutine

   !>
   !>  Debug flag query function
   !>
   !>  Usage:
   !>      if (dbg(DBG_DOMAIN)) call do_some_debug()
   !>
   function dbg(flag)
      implicit none
      logical :: dbg
      integer, intent(in) :: flag

      dbg = (iand(debug_level, flag) .ne. 0)
   end function
end module module_debug
