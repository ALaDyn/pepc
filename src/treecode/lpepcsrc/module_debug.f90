!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates some debugging and i/o specific routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_debug
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character(len=255), public :: debug !< string that will be output by the print_debug()-routine

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
      integer, public :: ipefile        = 20 !< local output stream, TODO: open this file on demand only (for all frontends)

      integer, parameter, public :: DBG_STATUS      = B'0000000000000001'    ! 1
      integer, parameter, public :: DBG_TREE        = B'0000000000000010'    ! 2
      integer, parameter, public :: DBG_BUILD       = B'0000000000000100'    ! 4
      integer, parameter, public :: DBG_DOMAIN      = B'0000000000001000'    ! 8
      integer, parameter, public :: DBG_BRANCH      = B'0000000000010000'    ! 16
      integer, parameter, public :: DBG_STATS       = B'0000000000100000'    ! 32
      integer, parameter, public :: DBG_WALKSUMMARY = B'0000000001000000'    ! 64
      integer, parameter, public :: DBG_DUMPTREE    = B'0000000010000000'    ! 128
      integer, parameter, public :: DBG_TIMINGFILE  = B'0000000100000000'    ! 256
      integer, parameter, public :: DBG_LOADFILE    = B'0000001000000000'    ! 512
      integer, parameter, public :: DBG_WALK        = B'0000010000000000'    ! 1024
      integer, parameter, public :: DBG_PERIODIC    = B'0000100000000000'    ! 2048

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      public print_debug
      public dbg
      public pepc_status


   contains

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  lpepc status output
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine pepc_status(stat)
       use treevars, only : me
       implicit none
       character(*), intent(in) :: stat

       if (dbg(DBG_STATUS)) then
          write(ipefile,'("LPEPC | ", a)') stat
          if (me==0) write(*,'("LPEPC | ", a)') stat
       endif

     end subroutine


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  Debug flag query function
     !>
     !>  Usage:
     !>      if (dbg(DBG_DOMAIN)) call do_some_debug()
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     function dbg(flag)
       implicit none
       logical :: dbg
       integer, intent(in) :: flag

       dbg = (iand(debug_level, flag) .ne. 0)
     end function



     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !>
     !>  Debug output helper routine
     !>
     !>  Usage:
     !>      write (debug, *) 'Debug output text and', number
     !>      call print_debug(general_condition_logical, 2, [ipefile, 6], [.true., myrank.eq.0])
     !>
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine print_debug(condition, nstreams, streams, conditions)
       logical, intent(in) :: condition !< debug output is only prnted if condition==.true.
       integer*4, intent(in) :: nstreams !< number of output streams given in
       integer*4, intent(in) :: streams(nstreams) !< output-stream numbers to be printed to
       logical, intent(in) :: conditions(nstreams)!< debug output will only be printed if condition(i) is true for stream(i)

       integer*4 :: idx

       if (condition) then
         do idx = 1,nstreams
           if (conditions(idx) .eqv. .true.) then
             write(streams(idx),*) trim(debug)
           end if
         end do
       end if

     end subroutine print_debug

end module module_debug
