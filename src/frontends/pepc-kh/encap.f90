module encap
   use iso_c_binding
   implicit none

    ! variables for MPI within pepc
   type, bind(c) :: pepc_comm_t
      integer(c_int) :: mpi_size, mpi_rank, mpi_comm
   end type pepc_comm_t

   type, bind(c) :: pepc_pars_t
      integer(c_int) :: np, npp, pdump, fdump
      real(c_double) :: theta
      type(pepc_comm_t) :: pepc_comm
   end type pepc_pars_t

   type :: time_pars_t
      real(kind=8) :: te, dt
      integer :: nsteps
   end type time_pars_t

end module encap
