module encap
  use module_pepc_kinds
  use module_pepc_types
  use iso_c_binding
  implicit none

  ! variables for MPI within pepc
  type, bind(c) :: pepc_comm_t
    integer(c_int) :: mpi_size, mpi_rank, mpi_comm
  end type pepc_comm_t

  type, bind(c) :: pepc_pars_t
    integer(c_int) :: pdump, fdump, cdump
    integer(c_int64_t) :: np
    type(pepc_comm_t) :: pepc_comm
  end type pepc_pars_t

  type :: physics_pars_t
    real(kind=8) :: B0, vte, vti, qe, qi, me, mi, shear_halfwidth, shear_velocity
    real(kind=8), dimension(3) :: l_plasma
    integer(kind = kind_particle) :: ni
  end type physics_pars_t

  type :: time_pars_t
    real(kind=8) :: tresume, dt
    integer :: nsteps, nresume
  end type time_pars_t

  type field_grid_t
    integer(kind = kind_particle), dimension(2) :: n
    integer(kind = kind_particle) :: ntot, nl
    real(kind = 8), dimension(2) :: offset, extent, dx
    type(t_particle), dimension(:), allocatable :: p
    real(kind = 8), dimension(:,:), allocatable :: ne, ni, vex, vey, vix, viy, ne_from_left, ni_from_left
  end type field_grid_t
end module encap
