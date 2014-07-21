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
!> Contains all lpepc-specific types and routines for registering them to MPI
!>
module module_pepc_kinds
  implicit none
  
  private

  include 'mpif.h'

  ! ATTENTION: if kind_particle, kind_default, or kind_key are modified, respective adaptions have to be performed in the first 50 lines of sl_pepckeys.h (just grep for the modified kind_XX values)
  integer, public, parameter :: kind_particle     = kind(1_8) ! ATTENTION, see above
  integer, public, parameter :: MPI_KIND_PARTICLE = MPI_INTEGER8
  integer, public, parameter :: kind_node         = kind_particle
  integer, public, parameter :: MPI_KIND_NODE     = MPI_KIND_PARTICLE
  integer, public, parameter :: kind_key          = kind(1_8) ! ATTENTION, see above
  integer, public, parameter :: MPI_KIND_KEY      = MPI_INTEGER8
  integer, public, parameter :: kind_byte         = kind(1_1)
  integer, public, parameter :: MPI_KIND_BYTE     = MPI_BYTE
  integer, public, parameter :: kind_level        = kind(1_1)
  integer, public, parameter :: MPI_KIND_LEVEL    = MPI_BYTE
  
  integer, public, parameter :: kind_default      = kind(1_4) !< default integer kind as the MPI standard calls it. This should be kind(1). We use kind(1_4) instead to allow for switching the default integer kind with certain compiler command line parameters  ! ATTENTION, see above
  integer, public, parameter :: MPI_KIND_DEFAULT  = MPI_INTEGER
  integer, public, parameter :: kind_pe           = kind_default !< this has to be of default integer kind - otherwise MPI gets angry if we use an owner as target rank of an MPI operation
  integer, public, parameter :: MPI_KIND_PE       = MPI_INTEGER
  integer, public, parameter :: kind_dim          = kind_level
  
  integer, public, parameter :: kind_physics      = kind(1._8)
  integer, public, parameter :: MPI_KIND_PHYSICS  = MPI_REAL8

end module module_pepc_kinds
