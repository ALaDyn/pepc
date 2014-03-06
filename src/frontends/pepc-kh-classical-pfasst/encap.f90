module mpi
  implicit none
  include 'mpif.h'
end module

module encap
   use module_pepc_types
   use iso_c_binding
   implicit none

   integer, parameter :: WM_BORIS_SDC   = 1
   integer, parameter :: WM_BORIS_MLSDC = 2
   integer, parameter :: WM_BORIS       = 3
   integer, parameter :: WM_BENEDIKT    = 4

   integer, parameter :: IFILE_SUMMAND_NFEVAL = 50
   integer, parameter :: IFILE_SUMMAND_NITER  = 70

    ! variables for MPI within pepc
   type :: pepc_comm_t
      ! MPI variables
      integer(kind_pe) :: rank_space, nrank_space
      integer(kind_pe) :: rank_world, nrank_world !< rank/num_ranks in MPI_COMM_WORLD (used for globally unique output)
      integer(kind_pe) :: rank_time,  nrank_time
      integer(kind_default) :: comm_space, comm_time
      !
      logical :: root_stdio
      logical :: root_file
   end type pepc_comm_t

   type :: pepc_pars_t
      integer(kind_default) :: pdump, fdump, cdump
      integer(kind_particle) :: np
      type(pepc_comm_t) :: pepc_comm
      integer :: workingmode
   end type pepc_pars_t

  type :: physics_pars_t  ! FIXME * m/Q in B0 ??
    real(kind=8) :: B0, vte, vti, qe, qi, me, mi, shear_halfwidth, shear_velocity
    real(kind=8), dimension(3) :: l_plasma
    integer(kind = kind_particle) :: ni
  end type physics_pars_t

   type :: time_pars_t
      real(kind=8) :: te, dt
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
