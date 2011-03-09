!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Test and demonstration program for memwatch.f90 and memwatch.c
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params
  integer*8, parameter :: k = 500000 !< number of integer variables that will be used as test case
  integer :: globaltest(1:k)
end module params

subroutine allocatemem(n)
  implicit none
  integer*8, intent(in) :: n
  integer, allocatable :: s(:)

  call PrintMemUsage("[before allocation]")
  allocate(s(n))
  s(1)=5
  call PrintMemUsage("[after allocation]")
  deallocate(s)
  call PrintMemUsage("[after deallocation]")

end subroutine allocatemem


subroutine allocatedyn(n)
  implicit none
  integer*8, intent(in) :: n
  integer :: s(1:n)

  call PrintMemUsage("[inside routine]")
  s(1)=5

end subroutine allocatedyn


subroutine allocatestat()
  use params
  implicit none
  integer :: s(1:k)

  call PrintMemUsage("[inside routine]")
  s(1)=5

end subroutine allocatestat

subroutine evilsigsegv()
  implicit none

  ! in this call the first parameter is missing --> SIGSEGV
  call OutputMemUsage("[evilsigsegv]", .true., 60)

end subroutine evilsigsegv




program memwatchtest
    use params
    implicit none

    integer :: i = 500
    integer*8 :: bytes

    call InitSignalHandler()

    call OutputMemUsage(0, "[at start of program]", .true., 60)
    call PrintMemUsage("[at program startup]")
    globaltest(1)=5

    call GetMemUsage(bytes)
    write(*,*) "Some user defined output after fetching our information with GetMemUsage(..):", bytes

    write(*,*) "************* Allocating", 4*k, "bytes of memory with allocate(..) statement *************"
    call allocatemem(k)

    write(*,*) "********* Allocating", 4*k, "bytes of memory with dynamically-sized local variable *******"
    call PrintMemUsage("[before routine]")
    call allocatedyn(k)
    call PrintMemUsage("[after routine]")

    write(*,*) "********* Allocating", 4*k, "bytes of memory with statically-sized local variable *******"
    call PrintMemUsage("[before routine]")
    call allocatestat()
    call PrintMemUsage("[after routine]")

    write(*,*) "Input desired amount of memory in MB (0 to abort)"
    read (*,*) i
    do while (i .ne. 0)
      write(*,*) "----------------------------- Trying to allocate ", i, "MB of memory"
      call NegotiateFreeMem(i)
      write (*,'("You can allocate", I6, "MB of memory --------------------------------")') i
      read (*,*) i
    end do

    call OutputMemUsage(1, "[at end of program]", .true., 60)

    call evilsigsegv()

end program memwatchtest
