!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Frees user-defined MPI types
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine libpepc_finalize()
  use treevars
  use tree_utils
  implicit none
  include 'mpif.h'

  integer :: ierr

  call MPI_TYPE_FREE( mpi_type_particle, ierr)
  call MPI_TYPE_FREE( mpi_type_results, ierr)
  call MPI_TYPE_FREE( mpi_type_multipole, ierr)

end subroutine libpepc_finalize







