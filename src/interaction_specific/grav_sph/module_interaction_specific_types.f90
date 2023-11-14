! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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
!> Contains all types that are specific to a certain interaction
!> all subroutines and types within this module are obligatory
!>
module module_interaction_specific_types
   use module_pepc_kinds
   implicit none

   integer, parameter :: max_neighbour_particles = 50

   integer :: num_neighbour_particles = 12

   !> Data structure for storing interaction-specific particle data
   type t_particle_data
      real*8 :: q               !< charge (equals mass in case of gravity)
      real*8 :: v(3)            !< velocity (same time as x)
      real*8 :: v_minus_half(3) !< velocity (1/2 time step after x (t-1/2), for leap frog integrator)
      real*8 :: temperature
      integer :: type           !< a bitfield for storing particle properties
   end type t_particle_data
   integer, private, parameter :: nprops_particle_data = 5

   !> Data structure for results
   type t_particle_results
      real*8 :: maxdist2        !< maxval(dist2)
      integer :: maxidx         !< maxloc(dist2)
      integer*8:: neighbour_nodes(max_neighbour_particles) ! FIXME: this should be integer(kind_node) and MPI_KIND_NODE
      real*8 :: dist2(max_neighbour_particles)
      real*8 :: dist_vector(3, max_neighbour_particles)                           ! distance_vectors from particle to neighbour with respact to periodic shift vector
      real*8 :: rho             !< density for sph
      real*8 :: h               !< smoothing-length for sph
      real*8 :: sph_force(1:3)
      real*8 :: temperature_change
      real*8, dimension(3) :: e !< gravitational field (variable name taken from coulomb specific)
      real*8 :: pot             !< potential energy (in gravitational field)
   end type t_particle_results
   integer, private, parameter :: nprops_particle_results = 11

   !> Data structure for storing multiple moments of tree nodes
   type t_tree_node_interaction_data
      real*8 :: coc(3)     !< center of charge
      real*8 :: charge     !< net charge sum
      real*8 :: dip(3)     !< dipole moment
      real*8 :: quad(3)    !< diagonal quadrupole moments
      real*8 :: xyquad     !< other quadrupole moments
      real*8 :: yzquad
      real*8 :: zxquad
      real*8 :: bmax       !< TODO: what is this variable?
      real*8 :: v(1:3)     !< velocity
      real*8 :: temperature
      real*8 :: rho        !< sph density
      real*8 :: h          !< sph smoothing-length
   end type t_tree_node_interaction_data
   integer, private, parameter :: nprops_tree_node_interaction_data = 12

   ! bit switches for particles types. use only powers of 2, combine with IOR, eg.: ior(PARTICLE_TYPE_FIXED, PARTICLE_TYPE_NONGAS)
   integer, parameter :: PARTICLE_TYPE_DEFAULT = 0    !< for setting all bits to 0, default values: moving, sph, ...
   integer, parameter :: PARTICLE_TYPE_FIXED = 1      !< fixed particles are not moved
   integer, parameter :: PARTICLE_TYPE_NONGAS = 2     !< treat particle as non-gas and compute no sph force
   ! TODO: implement nongas particles in neighbour search and force summation

contains

   !>
   !> Creates and registers interaction-specific MPI-types
   !> is automatically called from register_libpepc_mpi_types()
   !>
   subroutine register_interaction_specific_mpi_types(mpi_type_particle_data, MPI_TYPE_tree_node_interaction_data, mpi_type_particle_results)
      use mpi
      implicit none
      integer, intent(out) :: mpi_type_particle_data, MPI_TYPE_tree_node_interaction_data, mpi_type_particle_results

      integer, parameter :: max_props = nprops_particle_data + nprops_particle_results + nprops_tree_node_interaction_data ! maxval([..]) would be enough, but ifort does notlike that

      integer :: ierr
      ! address calculation
      integer, dimension(1:max_props) :: blocklengths, displacements, types
      integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
      ! dummies for address calculation
      type(t_particle_data)    :: dummy_particle_data
      type(t_particle_results) :: dummy_particle_results
      type(t_tree_node_interaction_data)   :: dummy_tree_node_interaction_data

      ! register particle data type
      blocklengths(1:nprops_particle_data) = [1, 3, 3, 1, 1]
      types(1:nprops_particle_data)        = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_INTEGER] !&
      call MPI_GET_ADDRESS(dummy_particle_data,              address(0), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_data%q,            address(1), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_data%v,            address(2), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_data%v_minus_half, address(3), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_data%temperature,  address(4), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_data%type,         address(5), ierr) !&
      displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
      call MPI_TYPE_STRUCT(nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr)
      call MPI_TYPE_COMMIT(mpi_type_particle_data, ierr)

      ! register results data type
      blocklengths(1:nprops_particle_results) = [1, 1, max_neighbour_particles, max_neighbour_particles, 3 * max_neighbour_particles, 1, 1, 3, 1, 3, 1]
      types(1:nprops_particle_results)        = [MPI_REAL8, MPI_INTEGER, MPI_INTEGER8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8] !&
      call MPI_GET_ADDRESS(dummy_particle_results,                    address(0), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%maxdist2,           address(1), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%maxidx,             address(2), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%neighbour_nodes,    address(3), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%dist2,              address(4), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%dist_vector,        address(5), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%rho,                address(6), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%h,                  address(7), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%sph_force,          address(8), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%temperature_change, address(9), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%e,                  address(10), ierr) !&
      call MPI_GET_ADDRESS(dummy_particle_results%pot,                address(11), ierr) !&
      displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
      call MPI_TYPE_STRUCT(nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr)
      call MPI_TYPE_COMMIT(mpi_type_particle_results, ierr)

      ! register multipole data type
      blocklengths(1:nprops_tree_node_interaction_data) = [3, 1, 3, 3, 1, 1, 1, 1, 3, 1, 1, 1]
      types(1:nprops_tree_node_interaction_data)        = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8] !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data,             address(0), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%coc,         address(1), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%charge,      address(2), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%dip,         address(3), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quad,        address(4), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%xyquad,      address(5), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%yzquad,      address(6), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%zxquad,      address(7), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%bmax,        address(8), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%v,           address(9), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%temperature, address(10), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%rho,         address(11), ierr) !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%h,           address(12), ierr) !&
      displacements(1:nprops_tree_node_interaction_data) = int(address(1:nprops_tree_node_interaction_data) - address(0))
      call MPI_TYPE_STRUCT(nprops_tree_node_interaction_data, blocklengths, displacements, types, MPI_TYPE_tree_node_interaction_data, ierr)
      call MPI_TYPE_COMMIT(MPI_TYPE_tree_node_interaction_data, ierr)
   end subroutine register_interaction_specific_mpi_types
end module module_interaction_specific_types
