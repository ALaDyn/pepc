!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Tries to find out the maximum of available memory by successive allocation
!> starting from megabytes. The search is only performed downwards to be
!> able to limit the total amount of used memory
!>
!> If something seriously goes wrong, the routine returns 0 and throws
!> an error message
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NegotiateFreeMem(megabytes)
  implicit none
  integer, intent(inout) :: megabytes
  integer :: res

    ! interface to C-routine
    interface
       integer function internal_negotiatefreemem(mb)
         integer, intent(in) :: mb
       end function internal_negotiatefreemem
    end interface

  res = internal_negotiatefreemem(megabytes)

  if (res > 0) then
    megabytes = res
  else
    write (*,'("Evil error: allocating", I6, "MB of memory failed and even no smaller block was available")') megabytes
    megabytes = 0
  end if

end subroutine NegotiateFreeMem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Returns actually occupied memory in bytes
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetMemUsage(bytes)
  implicit none
  integer*8, intent(out) :: bytes !< used memory in bytes

    ! interface to C-routine
    interface
       integer*8 function internal_getmemusage()
       end function internal_getmemusage
    end interface


  bytes = internal_getmemusage()

end subroutine GetMemUsage



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Prints some more or less informative memory information
!>  Parameter text can be used for arbitrary information that
!>  shall be included in output
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintMemUsage(text)
  implicit none
  character(LEN=*), intent(in) :: text !< arbitrary text that will be included in the output

  integer*8 :: bytes

  call GetMemUsage(bytes)

  write(*,'("using ",I16," bytes (",I9,"KB | ",I5,"MB) @ ",a)') bytes, NINT(real(bytes)/1024.),  NINT(real(bytes)/1024./1024.), text

end subroutine PrintMemUsage



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Outputs current memory usage to the file with handle file_mem_debug
!> if debug_mem == .true.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OutputMemUsage(ID, text, debug_memory, file_mem_debug)
  implicit none
  integer :: ID !< ID of position in program (for easier identification in plots)
  character(LEN=*), intent(in) :: text !< arbitrary text that will be included in the output
  logical :: debug_memory
  integer :: file_mem_debug
  integer*8 :: bytes

  if (debug_memory) then
    call GetMemUsage(bytes)
    write(file_mem_debug,'(I5,I20,I20,I20,"        """,a,"""")') ID, bytes, NINT(real(bytes)/1024.),  NINT(real(bytes)/1024./1024.), text
  end if

end subroutine OutputMemUsage


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Installs a signal Handler that tries to create a stack trace when
!> receiving sigterm, sigkill, sigsegv, sigint
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitSignalHandler()
  implicit none
    ! interface to C-routine
    interface
       subroutine init_signal_handler()
       end subroutine init_signal_handler
    end interface

    call init_signal_handler()

end subroutine InitSignalHandler
