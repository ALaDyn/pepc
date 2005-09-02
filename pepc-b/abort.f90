subroutine cleanup
  use physvars
  implicit none
  include 'mpif.h'
  integer :: ierr
  if (my_rank==0) write(0,*) 'User-abort of program: cleaning up ...'
  call closefiles
  call MPI_FINALIZE(ierr)
  stop
end subroutine cleanup
