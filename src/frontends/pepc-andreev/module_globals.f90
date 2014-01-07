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
module pepca_globals
  use module_pepc_types
  implicit none
  save
  
  integer(kind_dim), parameter :: dim = 2
  
  ! grid for density output
  integer :: Ngrid(1:3) = [500, 1000, 1]

  ! MPI variables
  integer(kind_pe) :: my_rank, n_ranks
  logical :: root

  ! time variables
  real*8 :: dt       !< timestep (initially in fs, later in simunits)
  integer :: nt      !< number of timesteps

  ! interaction cutoff parameter
  real*8 :: eps

  ! control variables
  integer :: particle_output_interval      !< turn vtk output on/off
  integer :: domain_output_interval        !< turn vtk output on/off

end module
