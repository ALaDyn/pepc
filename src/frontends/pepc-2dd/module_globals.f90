! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
!                         Forschungszentrum Juelich GmbH,
!                         Germany
!
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!>
!> global declarations for pepc-andreev
!>
module module_globals
  use module_pepc_kinds
  use module_pepc_types
  implicit none
  save

  !integer(kind_dim), parameter :: dim = 2

  ! MPI variables
  integer :: my_rank, n_ranks,communicator
  logical :: root

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt               ! number of timesteps
  integer(kind_particle) :: tnp ! total number of particles
  integer(kind_particle) :: np  ! local number of particles
  real(kind_particle)    :: eps2,theta2
  integer(kind_particle) :: initial_setup ! initial setup
  logical :: particle_output  ! turn vtk output on/off
  logical :: domain_output    ! turn vtk output on/off
  logical :: particle_test    ! turn direct summation on/off
  integer :: diag_interval    ! number of timesteps between all diagnostics and IO
  integer :: restart_step    ! number of timesteps between all diagnostics and IO
  real(kind_particle) :: Lx                       ! box size - x direction
  real(kind_particle) :: Ly                       ! box size - y direction
  real(kind_particle) :: Lz                       ! box size - z direction

!  real(kind_particle), parameter :: pi = acos(-1.0)   !< ludolfs number :-)


  integer(kind_physics), dimension(2) :: n = [ 0, 0 ]
  real(kind_physics)   , dimension(2) :: offset = [ 0., 0. ]
  real(kind_physics)   , dimension(2) :: extent = [ 0., 0. ]

  character(*), parameter :: FRONTEND_NAME = "pepc-2dd"
  character(255)          :: ischeme,restart_file


end module
