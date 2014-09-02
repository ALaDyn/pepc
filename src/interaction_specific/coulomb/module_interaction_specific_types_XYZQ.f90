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
      use module_pepc_kinds
      implicit none

      !> Data structure for storing interaction-specific particle data
      type t_particle_data
         real(kind_physics) :: q
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 1

      !> Data structure for shipping results
      type t_particle_results
         real(kind_physics), dimension(3) :: e
         real(kind_physics) :: pot
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 2

      type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results([0., 0., 0.], 0.)

      !> Data structure for storing multiple moments of tree nodes
      type t_multipole_moments
        real(kind_physics) :: charge     ! net charge sum
        real(kind_physics) :: abs_charge !  absolute charge sum
        real(kind_physics) :: dip(3)     ! dipole moment
        real(kind_physics) :: quad(3)    ! diagonal quadrupole moments
        real(kind_physics) :: xyquad     ! other quadrupole moments
        real(kind_physics) :: yzquad
        real(kind_physics) :: zxquad
        real(kind_physics) :: bmax
      end type t_multipole_moments
      integer, private, parameter :: nprops_multipole_moments = 8

      type t_local_coefficients
        real(kind_physics) :: f     
        real(kind_physics) :: fx     
        real(kind_physics) :: fy     
        real(kind_physics) :: fz     
        real(kind_physics) :: fxx
        real(kind_physics) :: fyy
        real(kind_physics) :: fzz
        real(kind_physics) :: fxy
        real(kind_physics) :: fyz
        real(kind_physics) :: fzx
      end type t_local_coefficients  

      type(t_local_coefficients), parameter :: EMPTY_LOCAL_COEFFICIENTS  = t_local_coefficients(0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)

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

        call vtkf%write_data_array("pot", r(:)%pot)
        call vtkf%write_data_array("field", r(:)%e(1), r(:)%e(2), r(:)%e(3))
      end subroutine vtk_write_particle_data_results


      !>
      !> Writes (a sensible subset of) tree node interaction data to a VTK file.
      !>
      subroutine vtk_write_multipole_moments(d, vtkf)
        use module_vtk
        implicit none

        type(t_multipole_moments), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf

        call vtkf%write_data_array("charge", d(:)%charge)
        call vtkf%write_data_array("abs_charge", d(:)%abs_charge)
      end subroutine vtk_write_multipole_moments


      !>
      !> Creates and registers interaction-specific MPI-types
      !> is automatically called from register_libpepc_mpi_types()
      !>
      subroutine register_interaction_specific_mpi_types(mpi_type_particle_data, MPI_TYPE_multipole_moments, mpi_type_particle_results)
        implicit none
        include 'mpif.h'
        integer, intent(out) :: mpi_type_particle_data, MPI_TYPE_multipole_moments, mpi_type_particle_results

        integer, parameter :: max_props = nprops_particle_data + nprops_particle_results + nprops_multipole_moments ! maxval([..]) would be enough, but ifort does notlike that

        integer :: ierr
        ! address calculation
        integer, dimension(1:max_props) :: blocklengths, displacements, types
        integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
        ! dummies for address calculation
        type(t_particle_data)    :: dummy_particle_data
        type(t_particle_results) :: dummy_particle_results
        type(t_multipole_moments)   :: dummy_multipole_moments

        ! register particle data type
        blocklengths(1:nprops_particle_data)  = [1]
        types(1:nprops_particle_data)         = [MPI_KIND_PHYSICS]
        call MPI_GET_ADDRESS( dummy_particle_data,   address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%q, address(1), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [3, 1]
        types(1:nprops_particle_results)         = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS]
        call MPI_GET_ADDRESS( dummy_particle_results,      address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%e,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%pot,  address(2), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_multipole_moments)  = [1, 1, 3, 3, 1, 1, 1, 1]
        types(1:nprops_multipole_moments)         = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS]
        call MPI_GET_ADDRESS( dummy_multipole_moments,            address(0), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%charge,     address(1), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%abs_charge, address(2), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%dip,        address(3), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%quad,       address(4), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xyquad,     address(5), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yzquad,     address(6), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zxquad,     address(7), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%bmax,       address(8), ierr )
        displacements(1:nprops_multipole_moments) = int(address(1:nprops_multipole_moments) - address(0))
        call MPI_TYPE_STRUCT( nprops_multipole_moments, blocklengths, displacements, types, MPI_TYPE_multipole_moments, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_multipole_moments, ierr)
      end subroutine register_interaction_specific_mpi_types
end module module_interaction_specific_types
