subroutine abort
  implicit none
  include 'mpif.h'
  integer :: ierr
  call closefiles
  call MPI_FINALIZE(ierr)

end subroutine abort
