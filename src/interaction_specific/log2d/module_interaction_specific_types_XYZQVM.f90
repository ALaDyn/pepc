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

      integer, parameter :: pMultipole = 10

      integer, parameter :: CALC_FORCE_SOURCE_KIND_PARTICULAR = 1
      integer, parameter :: CALC_FORCE_SOURCE_KIND_DIRICHLET = 2
      integer, parameter :: CALC_FORCE_SOURCE_KIND_NEUMANN = 3

      !> Data structure for storing interaction-specific particle data
      type t_particle_data
         real(kind_physics) :: q
         real(kind_physics) :: v(3)
         real(kind_physics) :: m
#ifdef BEM2D
         real(kind_physics) :: phi
         real(kind_physics) :: ra(2)
         real(kind_physics) :: rb(2)
         integer :: source_kind
#endif
      end type t_particle_data
#ifdef BEM2D
      integer, private, parameter :: nprops_particle_data = 7
#else
      integer, private, parameter :: nprops_particle_data = 3
#endif

      !> Data structure for shipping results
      type t_particle_results
         real(kind_physics), dimension(2) :: e
         real(kind_physics) :: pot
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 2

      type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results([0., 0.], 0.)

      !> Data structure for storing multiple moments of tree nodes
      type t_tree_node_interaction_data
        real(kind_physics) :: coc(3)     ! centre of charge
        real(kind_physics) :: charge     ! net charge sum
        real(kind_physics) :: abs_charge !  absolute charge sum
        complex(kind_physics) :: omega(pMultipole)  ! multipole moments
        real(kind_physics) :: bmax
      end type t_tree_node_interaction_data
      integer, private, parameter :: nprops_tree_node_interaction_data = 5

      contains

      !>
      !> Writes particle interaction data and results to a VTK file.
      !>
      subroutine vtk_write_particle_data_results(d, r, vtkf)
        use module_vtk
        implicit none

        type(t_particle_data), intent(in) :: d(:)
        type(t_particle_results), intent(in) :: r(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf

        call vtkf%write_data_array("q", d(:)%q)
        call vtkf%write_data_array("m", d(:)%m)
        call vtkf%write_data_array("v", d(:)%v(1), d(:)%v(2), d(:)%v(3))

        call vtkf%write_data_array("pot", r(:)%pot)
        call vtkf%write_data_array("field", r(:)%e(1), r(:)%e(2))
      end subroutine vtk_write_particle_data_results


      !>
      !> Writes (a sensible subset of) tree node interaction data to a VTK file.
      !>
      subroutine vtk_write_node_interaction_data(d, vtkf)
        use module_vtk
        implicit none

        type(t_tree_node_interaction_data), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf

        call vtkf%write_data_array("charge", d(:)%charge)
        call vtkf%write_data_array("abs_charge", d(:)%abs_charge)
      end subroutine vtk_write_node_interaction_data


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
        integer :: MPI_KIND_PHYSICS_COMPLEX

        ! register particle data type
#ifdef BEM2D
        blocklengths(1:nprops_particle_data) = [1, 3, 1, 1, 2, 2, 1]
        types(1:nprops_particle_data)        = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
          MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_INTEGER]
#else
        blocklengths(1:nprops_particle_data) = [1, 3, 1]
        types(1:nprops_particle_data)        = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS]
#endif
        call MPI_GET_ADDRESS( dummy_particle_data, address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%q, address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%v, address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%m, address(3), ierr )
#ifdef BEM2D
        call MPI_GET_ADDRESS( dummy_particle_data%phi, address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%ra, address(5), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%rb, address(6), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%source_kind, address(7), ierr )
#endif
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [2, 1]
        types(1:nprops_particle_results)         = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS]
        call MPI_GET_ADDRESS( dummy_particle_results,      address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%e,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%pot,  address(2), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        select case(MPI_KIND_PHYSICS)
          ! all of these are optionally defined and if missing cause compilation error
          ! MPI_REAL2 and MPI_COMPLEX4 seem to be rare, disable for now
          !case(MPI_REAL2)
          !  MPI_KIND_PHYSICS_COMPLEX = MPI_COMPLEX4
          case(MPI_REAL4)
            MPI_KIND_PHYSICS_COMPLEX = MPI_COMPLEX8
          case(MPI_REAL8)
            MPI_KIND_PHYSICS_COMPLEX = MPI_COMPLEX16
          case(MPI_REAL16)
            MPI_KIND_PHYSICS_COMPLEX = MPI_COMPLEX32
          case default
            stop
        end select

        ! register multipole data type
        blocklengths(1:nprops_tree_node_interaction_data)  = [3, 1, 1, pMultipole, 1]
        types(1:nprops_tree_node_interaction_data)         = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
          MPI_KIND_PHYSICS_COMPLEX, MPI_KIND_PHYSICS]
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data,            address(0), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%coc,        address(1), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%charge,     address(2), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%abs_charge, address(3), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%omega,      address(4), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%bmax,       address(5), ierr )
        displacements(1:nprops_tree_node_interaction_data) = int(address(1:nprops_tree_node_interaction_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_tree_node_interaction_data, blocklengths, displacements, types, MPI_TYPE_tree_node_interaction_data, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_interaction_data, ierr)

      end subroutine register_interaction_specific_mpi_types
end module module_interaction_specific_types
