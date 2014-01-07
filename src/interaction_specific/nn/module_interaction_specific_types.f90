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
!> Contains all types that are specific to a certain interaction
!> all subroutines and types within this module are obligatory
!>
module module_interaction_specific_types
      implicit none

      integer, parameter :: kind_key_i = kind(1_8)

      integer, parameter :: max_neighbour_particles = 50

      integer :: num_neighbour_particles = 12


      !> Data structure for storing interaction-specific particle data
      type t_particle_data
        integer*8 :: particle_id !< frontend-defined ID for recognizing the particle inside the tree
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 1

      !> Data structure for storing multiple moments of tree nodes
      type t_tree_node_interaction_data
        real*8 :: coc(3)     !< center of charge
        integer*8 :: particle_id !< corresponding t_particle_data%particle_id (only valid for leaves)
      end type t_tree_node_interaction_data
      integer, private, parameter :: nprops_tree_node_interaction_data = 2

      type p_tree_node_interaction_data
        type(t_tree_node_interaction_data), pointer :: p
      end type

      !> Data structure for results
      type t_particle_results
         real*8 :: maxdist2       !< maxval(dist2)
         integer :: maxidx        !< maxloc(dist2)
         type(p_tree_node_interaction_data) :: neighbour_nodes(max_neighbour_particles) ! FIXME: this should be integer(kind_node) and MPI_KIND_NODE
         real*8 :: dist2(max_neighbour_particles)
         real*8 :: dist_vector(3,max_neighbour_particles) ! distance_vectors from particle to neighbour with respect to periodic shift vector
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 5

      type t_particle_pack
         real*8, allocatable :: maxdist2(:)       !< maxval(dist2)
         integer, allocatable :: maxidx(:)        !< maxloc(dist2)
         type(p_tree_node_interaction_data), allocatable :: neighbour_nodes(:,:) ! FIXME: this should be integer(kind_node) and MPI_KIND_NODE
         real*8, allocatable :: dist2(:,:)
         real*8, allocatable :: dist_vector(:,:,:) ! distance_vectors from particle to neighbour with respect to periodic shift vector
      end type t_particle_pack

      contains

      !>
      !> Creates and registers interaction-specific MPI-types
      !> is automatically called from register_libpepc_mpi_types()
      !>
      subroutine register_interaction_specific_mpi_types(mpi_type_particle_data, MPI_TYPE_tree_node_interaction_data, mpi_type_particle_results)
        implicit none
        include 'mpif.h'
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
        blocklengths(1:nprops_particle_data)  = [1]
        types(1:nprops_particle_data)         = [MPI_INTEGER8]
        call MPI_GET_ADDRESS( dummy_particle_data,              address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%particle_id,  address(1), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [1, 1, max_neighbour_particles, max_neighbour_particles, 3*max_neighbour_particles]
        types(1:nprops_particle_results)         = [MPI_REAL8, MPI_INTEGER, MPI_INTEGER8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_results,                    address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%maxdist2,           address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%maxidx,             address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%neighbour_nodes,    address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%dist2,              address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%dist_vector,        address(5), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_tree_node_interaction_data)  = [3, 1]
        types(1:nprops_tree_node_interaction_data)         = [MPI_REAL8, MPI_INTEGER8]
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data,             address(0), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%coc,         address(1), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%particle_id, address(2), ierr )
        displacements(1:nprops_tree_node_interaction_data) = int(address(1:nprops_tree_node_interaction_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_tree_node_interaction_data, blocklengths, displacements, types, MPI_TYPE_tree_node_interaction_data, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_interaction_data, ierr)
      end subroutine register_interaction_specific_mpi_types
end module module_interaction_specific_types
