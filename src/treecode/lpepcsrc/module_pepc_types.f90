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
!> Contains all lpepc-specific types and routines for registering them to MPI
!>
module module_pepc_types
   use module_pepc_kinds
   use module_interaction_specific_types
   use mpi
   implicit none

   private


   public :: t_particle_data
   public :: t_particle_results
   public :: t_tree_node_interaction_data

   integer :: MPI_TYPE_particle_data_sca,              & ! we have these private for now to better control where they are used
              MPI_TYPE_tree_node_interaction_data_sca, &
              MPI_TYPE_particle_results_sca
   integer, public :: MPI_TYPE_particle_sca,          & ! we introduce those in a scalar and vector variant
                      MPI_TYPE_particle_vec,          & ! make sure to use the correct ones
                      MPI_TYPE_tree_node_sca,         &
                      MPI_TYPE_tree_node_vec,         &
                      MPI_TYPE_tree_node_package_sca, &
                      MPI_TYPE_tree_node_package_vec, &
                      MPI_TYPE_request_eager_sca,     &
                      MPI_TYPE_request_eager_vec

   !> Data structure for shipping single particles
   integer, private, parameter :: nprops_particle = 7 ! # particle properties to ship
   type, public :: t_particle
      real(kind_physics) :: x(1:3)      !< coordinates
      real*8 :: work        !< work load from force sum, ATTENTION: the sorting library relies on this being a real*8
      integer(kind_key) :: key      !< particle key, i.e. key on highest tree level
      integer(kind_node) :: node_leaf !< node index of corresponding leaf (tree node)
      integer(kind_particle) :: label     !< particle label (only for diagnostic purposes, can be used freely by the frontend
      type(t_particle_data) :: data       !< real physics (charge, etc.)
      type(t_particle_results) :: results !< results of calc_force_etc and companions
   end type t_particle

   !> Data structure for tree nodes
   type, public :: t_tree_node
      integer(kind_key) :: key
      integer(kind_byte) :: flags_global !< flags which are globally valid (and have to be shipped to other ranks)
      integer(kind_byte) :: flags_local  !< flags which are only locally valid (may not be shipped)
      logical(kind_byte) :: request_posted !< is set to .true. after a request for child data has been put onto the request list
      integer(kind_level) :: level
      integer(kind_pe) :: owner
      integer(kind_node) :: leaves       !< total number of leaf nodes below this node
      integer(kind_node) :: descendants  !< total number of descendants (tree nodes and leaves) below this node
      integer(kind_node) :: parent
      integer(kind_node) :: first_child
      integer(kind_node) :: next_sibling
      type(t_tree_node_interaction_data) :: interaction_data
   end type t_tree_node

   !> Data structure for shipping tree nodes
   integer, private, parameter :: nprops_tree_node_package = 10
   type, public :: t_tree_node_package
      integer(kind_key) :: key
      integer(kind_byte) :: flags_global
      integer(kind_level) :: level ! an integer*1 is sufficient. we place it here, to avoid excessive padding
      integer(kind_byte) :: dummy ! manual padding - so we know what exactly is happening here
      integer(kind_pe) :: owner
      integer(kind_node) :: leaves !< total number of leaf nodes below this node
      integer(kind_node) :: descendants  !< total number of descendants (tree nodes and leaves) below this node
      integer(kind_node) :: parent
      integer(kind_node) :: first_child
      type(t_tree_node_interaction_data) :: interaction_data
   end type t_tree_node_package

   !> Data structure for requesting tree nodes with eager send algorithm
   integer, private, parameter :: nprops_request_eager = 3
   type, public :: t_request_eager
      integer(kind_node) :: node
      integer(kind_node) :: parent
      type(t_particle) :: particle
   end type t_request_eager

   public register_lpepc_mpi_types
   public free_lpepc_mpi_types

contains

   !>
   !> Creates and registers lpepc-MPI types
   !>
   subroutine register_lpepc_mpi_types()
      use module_interaction_specific_types
      use treevars, only : me
      implicit none

      integer, parameter :: max_props = nprops_particle + nprops_tree_node_package + nprops_request_eager

      integer(kind_default) :: ierr
      ! address calculation
      integer, dimension(1:max_props) :: blocklengths, types
      integer(KIND=MPI_ADDRESS_KIND), dimension(1:max_props) :: displacements
      integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
      integer(KIND=MPI_ADDRESS_KIND) :: extent !< to store the extent to the next type in arrays
      ! dummies for address calculation
      type(t_particle)  :: dummy_particle(2)
      type(t_tree_node_package) :: dummy_tree_node_package(2)
      type(t_request_eager) :: dummy_request(2)

      ! first register the interaction-specific MPI types since they are embedded into the lpepc-datatypes
      ! the types returned by this call are for scalar use only
      call register_interaction_specific_mpi_types(MPI_TYPE_particle_data_sca, MPI_TYPE_tree_node_interaction_data_sca, MPI_TYPE_particle_results_sca)

      ! register particle type
      blocklengths(1:nprops_particle)  = [3, 1, 1, 1, 1, 1, 1]
      types(1:nprops_particle)         = [MPI_KIND_PHYSICS, MPI_REAL8, MPI_KIND_KEY, MPI_KIND_NODE, MPI_KIND_PARTICLE, &
         MPI_TYPE_particle_data_sca, MPI_TYPE_particle_results_sca]
      call MPI_GET_ADDRESS( dummy_particle(2),           extent, ierr )
      call MPI_GET_ADDRESS( dummy_particle(1),           address(0), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%x,         address(1), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%work,      address(2), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%key,       address(3), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%node_leaf, address(4), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%label,     address(5), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%data,      address(6), ierr )
      call MPI_GET_ADDRESS( dummy_particle(1)%results,   address(7), ierr )
      displacements(1:nprops_particle) = address(1:nprops_particle) - address(0)
      extent = extent - address(0)
      call MPI_TYPE_CREATE_STRUCT( nprops_particle, blocklengths, displacements, types, MPI_TYPE_particle_sca, ierr )
      call MPI_TYPE_COMMIT( MPI_TYPE_particle_sca, ierr)
      call MPI_TYPE_CREATE_RESIZED( MPI_TYPE_particle_sca, 0_MPI_ADDRESS_KIND, extent, MPI_TYPE_particle_vec, ierr )
      call MPI_TYPE_COMMIT( MPI_TYPE_particle_vec, ierr)

      ! register tree_node type
      blocklengths(1:nprops_tree_node_package)  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
      types(1:nprops_tree_node_package)         = [MPI_KIND_KEY, MPI_KIND_BYTE, MPI_KIND_LEVEL, MPI_KIND_BYTE, &
         MPI_KIND_PE, MPI_KIND_NODE, MPI_KIND_NODE, MPI_KIND_NODE, MPI_KIND_NODE, MPI_TYPE_tree_node_interaction_data_sca]
      call MPI_GET_ADDRESS( dummy_tree_node_package(2),                  extent, ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1),                  address(0), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%key,              address(1), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%flags_global,     address(2), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%level,            address(3), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%dummy,            address(4), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%owner,            address(5), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%leaves,           address(6), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%descendants,      address(7), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%parent,           address(8), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%first_child,      address(9), ierr )
      call MPI_GET_ADDRESS( dummy_tree_node_package(1)%interaction_data, address(10), ierr )
      displacements(1:nprops_tree_node_package) = address(1:nprops_tree_node_package) - address(0)
      extent = extent - address(0)
      call MPI_TYPE_CREATE_STRUCT( nprops_tree_node_package, blocklengths, displacements, types, MPI_TYPE_tree_node_package_sca, ierr )
      call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_package_sca, ierr)
      call MPI_TYPE_CREATE_RESIZED( MPI_TYPE_tree_node_package_sca, 0_MPI_ADDRESS_KIND, extent, MPI_TYPE_tree_node_package_vec, ierr )
      call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_package_vec, ierr )

      ! register request type
      blocklengths(1:nprops_request_eager)  = [1, 1, 1]
      types(1:nprops_request_eager)         = [MPI_KIND_NODE, MPI_KIND_NODE, MPI_TYPE_particle_sca]
      call MPI_GET_ADDRESS( dummy_request(2),                  extent, ierr )
      call MPI_GET_ADDRESS( dummy_request(1),                  address(0), ierr )
      call MPI_GET_ADDRESS( dummy_request(1)%node,             address(1), ierr )
      call MPI_GET_ADDRESS( dummy_request(1)%parent,           address(2), ierr )
      call MPI_GET_ADDRESS( dummy_request(1)%particle,         address(3), ierr )
      displacements(1:nprops_request_eager) = address(1:nprops_request_eager) - address(0)
      extent = extent - address(0)
      call MPI_TYPE_CREATE_STRUCT( nprops_request_eager, blocklengths, displacements, types, MPI_TYPE_request_eager_sca, ierr )
      call MPI_TYPE_COMMIT( MPI_TYPE_request_eager_sca, ierr)
      call MPI_TYPE_CREATE_RESIZED( MPI_TYPE_request_eager_sca, 0_MPI_ADDRESS_KIND, extent, MPI_TYPE_request_eager_vec, ierr )
      call MPI_TYPE_COMMIT( MPI_TYPE_request_eager_vec, ierr )

      block
         integer :: sze
         integer(kind=MPI_ADDRESS_KIND) :: lb, ex

         if (me == 0) then
            write(*,*) 'MPI_TYPE_request_eager'
            call MPI_TYPE_SIZE(MPI_TYPE_request_eager_sca, sze, ierr)
            write(*,*) 'MPI_TYPE_SIZE            sca', sze
            call MPI_TYPE_SIZE(MPI_TYPE_request_eager_vec, sze, ierr)
            write(*,*) '                         vec', sze
            call MPI_TYPE_GET_EXTENT(MPI_TYPE_request_eager_sca, lb, ex, ierr)
            write(*,*) 'MPI_TYPE_GET_EXTENT      sca', lb, ex
            call MPI_TYPE_GET_EXTENT(MPI_TYPE_request_eager_vec, lb, ex, ierr)
            write(*,*) '                         vec', lb, ex
            call MPI_TYPE_GET_TRUE_EXTENT(MPI_TYPE_request_eager_sca, lb, ex, ierr)
            write(*,*) 'MPI_TYPE_GET_TRUE_EXTENT sca', lb, ex
            call MPI_TYPE_GET_TRUE_EXTENT(MPI_TYPE_request_eager_vec, lb, ex, ierr)
            write(*,*) '                         vec', lb, ex
            write(*,*) 'MPI_TYPE_particle'
            call MPI_TYPE_SIZE(MPI_TYPE_particle_sca, sze, ierr)
            write(*,*) 'MPI_TYPE_SIZE            sca', sze
            call MPI_TYPE_SIZE(MPI_TYPE_particle_vec, sze, ierr)
            write(*,*) '                         vec', sze
            call MPI_TYPE_GET_EXTENT(MPI_TYPE_particle_sca, lb, ex, ierr)
            write(*,*) 'MPI_TYPE_GET_EXTENT      sca', lb, ex
            call MPI_TYPE_GET_EXTENT(MPI_TYPE_particle_vec, lb, ex, ierr)
            write(*,*) '                         vec', lb, ex
            call MPI_TYPE_GET_TRUE_EXTENT(MPI_TYPE_particle_sca, lb, ex, ierr)
            write(*,*) 'MPI_TYPE_GET_TRUE_EXTENT sca', lb, ex
            call MPI_TYPE_GET_TRUE_EXTENT(MPI_TYPE_particle_vec, lb, ex, ierr)
            write(*,*) '                         vec', lb, ex
            write(*,*) 'MPI_TYPE_tree_node_package'
            call MPI_TYPE_SIZE(MPI_TYPE_tree_node_package_sca, sze, ierr)
            write(*,*) 'MPI_TYPE_SIZE            sca', sze
            call MPI_TYPE_SIZE(MPI_TYPE_tree_node_package_vec, sze, ierr)
            write(*,*) '                         vec', sze
            call MPI_TYPE_GET_EXTENT(MPI_TYPE_tree_node_package_sca, lb, ex, ierr)
            write(*,*) 'MPI_TYPE_GET_EXTENT      sca', lb, ex
            call MPI_TYPE_GET_EXTENT(MPI_TYPE_tree_node_package_vec, lb, ex, ierr)
            write(*,*) '                         vec', lb, ex
            call MPI_TYPE_GET_TRUE_EXTENT(MPI_TYPE_tree_node_package_sca, lb, ex, ierr)
            write(*,*) 'MPI_TYPE_GET_TRUE_EXTENT sca', lb, ex
            call MPI_TYPE_GET_TRUE_EXTENT(MPI_TYPE_tree_node_package_vec, lb, ex, ierr)
            write(*,*) '                         vec', lb, ex
         end if
      end block


   end subroutine register_lpepc_mpi_types


   !>
   !> Deregisters lpepc- and interaction-specific MPI types
   !>
   subroutine free_lpepc_mpi_types()
      implicit none
      integer(kind_default) :: ierr

      call MPI_TYPE_FREE( MPI_TYPE_tree_node_package_sca,            ierr)
      call MPI_TYPE_FREE( MPI_TYPE_tree_node_package_vec,            ierr)
      call MPI_TYPE_FREE( MPI_TYPE_particle_sca,                     ierr)
      call MPI_TYPE_FREE( MPI_TYPE_particle_vec,                     ierr)
      call MPI_TYPE_FREE( MPI_TYPE_particle_results_sca,             ierr)
      call MPI_TYPE_FREE( MPI_TYPE_request_eager_sca,                ierr)
      call MPI_TYPE_FREE( MPI_TYPE_request_eager_vec,                ierr)
      call MPI_TYPE_FREE( MPI_TYPE_tree_node_interaction_data_sca,   ierr)
      call MPI_TYPE_FREE( MPI_TYPE_particle_data_sca,                ierr)

   end subroutine free_lpepc_mpi_types

end module module_pepc_types
