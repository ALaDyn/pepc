!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Frees user-defined MPI types
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine libpepc_finalize()
  use treevars
  use module_fmm_framework
  implicit none
  include 'mpif.h'

  integer :: ierr

  call MPI_TYPE_FREE( mpi_type_particle, ierr)
  call MPI_TYPE_FREE( mpi_type_results, ierr)
  call MPI_TYPE_FREE( mpi_type_multipole, ierr)

  ! finalize framework for lattice contributions
  call fmm_framework_finalize()

end subroutine libpepc_finalize






