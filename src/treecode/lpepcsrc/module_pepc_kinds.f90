! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2019 Juelich Supercomputing Centre,
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
!> Contains all lpepc-specific integer constants to be used as kind type parameters and associated MPI type tags
!>
module module_pepc_kinds
   ! This module should use the portable kind constants defined as part of the iso_fortran_env intrinsic module in the 2008
   ! standard. However, as of 2014, support for this feature is still considered bleeding edge by some vendors. Therefore,
   ! the constants are redefined below for integer types to the best of our knowledge.
   !use, intrinsic :: iso_fortran_env
   implicit none

   private

   include 'mpif.h'

   !> an 8-bit integer
   integer, parameter :: int8 = selected_int_kind(0)
   !> a 16-bit integer
   integer, parameter :: int16 = selected_int_kind(3)
   !> a 32-bit integer
   integer, parameter :: int32 = selected_int_kind(6)
   !> a 64-bit integer
   integer, parameter :: int64 = selected_int_kind(10)

   ! These constants should correspond to the IEEE precision classes.
   ! Note: Not all of these are supported on all platforms. E.g.: XLF on BG/Q does not have quad_precision.
   ! Note: Precision classes cannot be mapped to storage sizes in a straightforward way. This makes specifying the MPI datatype
   ! for a given precision class somewhat tricky. Usually precision classes shourd correspond to MPI datatypes like this:
   !   single precision -> MPI_REAL4
   !   double precision -> MPI_REAL8
   !   quad precision -> MPI_REAL16
   ! Note: Nonetheless, for the purpose of numerical analysis, these should be preferred to real32, real64 and real128 as defined
   ! in iso_fortran_env. Storage size does not define precision. E.g., on x86_64 using GNU Fortran, real128 is an 80-bit extended
   ! precision floating-point number, on BG/Q with XLF, real128 implements double double arithmetic.
   !> a floating point number offering precision and exponent range of IEEE single precision
   integer, parameter :: single_precision = selected_real_kind(6, 37)
   !> a floating point number offering precision and exponent range of IEEE double precision
   integer, parameter :: double_precision = selected_real_kind(15, 307)
   !> a floating point number offering precision and exponent range of IEEE quadruple precision
   integer, parameter :: quad_precision = selected_real_kind(33, 4931)

   ! ATTENTION: if kind_particle, kind_default, or kind_key are modified, respective adaptions have to be performed in the first 50
   ! lines of sl_pepckeys.h (just grep for the modified kind_XX values)
   !> an integer identifying a particle, e.g. an index into an array of particles
   integer, public, parameter :: kind_particle     = int64 !&  ATTENTION, see above
   integer, public, parameter :: MPI_KIND_PARTICLE = MPI_INTEGER8       !&
   !> an integer identifying a tree node, e.g. return value of `tree_node_get_first_child()`
   integer, public, parameter :: kind_node         = kind_particle      !&
   integer, public, parameter :: MPI_KIND_NODE     = MPI_KIND_PARTICLE  !&
   !> a hashed tree scheme binary key
   integer, public, parameter :: kind_key          = int64 !&  ATTENTION, see above
   integer, public, parameter :: MPI_KIND_KEY      = MPI_INTEGER8       !&
   !> a byte, i.e. an 8-bit integer
   integer, public, parameter :: kind_byte         = int8               !&
   integer, public, parameter :: MPI_KIND_BYTE     = MPI_BYTE           !&
   !> an integer identifying a level in the tree
   integer, public, parameter :: kind_level        = int8               !&
   integer, public, parameter :: MPI_KIND_LEVEL    = MPI_BYTE           !&

   !> default integer kind as the MPI standard calls it. This should be kind(1). We use int32 instead to allow for switching
   !> the default integer kind with certain compiler command line parameters
   integer, public, parameter :: kind_default      = int32 !&  ATTENTION, see above
   integer, public, parameter :: MPI_KIND_DEFAULT  = MPI_INTEGER        !&
   !> this has to be of default integer kind - otherwise MPI gets angry if we use an owner as target rank of an MPI operation
   integer, public, parameter :: kind_pe           = kind_default       !&
   integer, public, parameter :: MPI_KIND_PE       = MPI_INTEGER        !&
   !> an integer identifying a number of dimensions or one specific dimension
   integer, public, parameter :: kind_dim          = kind_level         !&

   !> the floating point kind for physical quantities, e.g. positions, interaction results, ...
   integer, public, parameter :: kind_physics      = double_precision   !&
   integer, public, parameter :: MPI_KIND_PHYSICS  = MPI_REAL8          !&
end module module_pepc_kinds
