!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Frees user-defined MPI types
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine libpepc_finalize()
  use treevars
  use module_fmm_framework
  use module_branching
  implicit none

  call status('FINALIZE')

  ! deregister mpi types
  call free_lpepc_mpi_types()

  ! finalize framework for lattice contributions
  call fmm_framework_finalize()

  ! finalize data structures in module_branches
  call branches_finalize()


end subroutine libpepc_finalize







