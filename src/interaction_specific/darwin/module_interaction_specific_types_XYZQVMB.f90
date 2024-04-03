! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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
!>
!>
!> Backend for two dimensional darwin approximation
!>
!>
module module_interaction_specific_types
   use module_pepc_kinds
   implicit none

   !> Data structure for storing interaction-specific particle data
   type t_particle_data
      real(kind_physics) :: q             !< electric charge
      real(kind_physics) :: v(3)          !< velocity
      real(kind_physics) :: m             !< mass
      real(kind_physics) :: g             !< mass
      !real(kind_physics) :: b(3)          !< magnetic filed at particle position (due to external fields applied in frontend)
   end type t_particle_data
   integer, private, parameter :: nprops_particle_data = 4

   !> Data structure for shipping results
   type t_particle_results
      real(kind_physics)               :: pot    !< Scalar potential
      real(kind_physics), dimension(3) :: E      !< tranversal part electric field - irrotational part of Electric Field
      real(kind_physics), dimension(3) :: A      !< Vector potential
      real(kind_physics), dimension(3) :: B      !< Magnetic Field - it has only Bz, orthogonal to xy plane
      real(kind_physics), dimension(3) :: J      !< Total Current density
      real(kind_physics), dimension(3) :: Jirr   !< Irrotational Current density
      real(kind_physics), dimension(3) :: dxA    !< Gradiend x Vector potential
      real(kind_physics), dimension(3) :: dyA    !< Gradiend y Vector potential
      real(kind_physics), dimension(3) :: dzA    !< Gradiend z Vector potential
!         real(kind_physics), dimension(3) :: dxE    !< Gradiend x tranversal part electric field
!         real(kind_physics), dimension(3) :: dyE    !< Gradiend y tranversal part electric field
!         real(kind_physics), dimension(3) :: dxxA   !< Hessian xx Vector potential
!         real(kind_physics), dimension(3) :: dxyA   !< Hessian xy Vector potential
!         real(kind_physics), dimension(3) :: dyyA   !< Hessian yy Vector potential
!         real(kind_physics)               :: rho    !< Charge density

   end type t_particle_results
   integer, private, parameter :: nprops_particle_results = 9

   type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results(0., [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], &
                                                                                      [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.])
   !type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results(0.,[0., 0., 0.])

   !> Data structure for storing multiple moments of tree nodes
   type t_tree_node_interaction_data
      real(kind_physics) :: coc(3)     ! centre of charge
      real(kind_physics) :: charge     ! net charge sum
      real(kind_physics) :: abs_charge !  absolute charge sum
      real(kind_physics) :: dip(3)     ! dipole moment
      real(kind_physics) :: quad(3)    ! diagonal quadrupole moments
      real(kind_physics) :: xyquad     ! other quadrupole moments
      real(kind_physics) :: yzquad
      real(kind_physics) :: zxquad
      real(kind_physics) :: monoj(3)   ! current density monopole - qj*vj
      real(kind_physics) :: dipjx(3)   ! current density dipole - qj*vxj*x - qj*vyj*x - qj*vzj*x
      real(kind_physics) :: dipjy(3)   ! current density dipole - qj*vxj*y - qj*vyj*y - qj*vzj*y
      real(kind_physics) :: dipjz(3)   ! current density dipole - qj*vxj*z - qj*vyj*z - qj*vzj*z
      real(kind_physics) :: quadjx(3)  ! current density quadrupole - qj*vxj*x**2 - qj*vyj*x**2 - qj*vzj*x**2
      real(kind_physics) :: quadjy(3)  ! current density quadrupole - qj*vxj*y**2 - qj*vyj*y**2 - qj*vzj*y**2
      real(kind_physics) :: quadjz(3)  ! current density quadrupole - qj*vxj*z**2 - qj*vyj*z**2 - qj*vzj*z**2
      real(kind_physics) :: quadjxy(3) ! current density quadrupole - qj*vxj*x*y - qj*vyj*x*y - qj*vzj*x*y
      real(kind_physics) :: quadjyz(3) ! current density quadrupole - qj*vxj*z*y - qj*vyj*z*y - qj*vzj*z*y
      real(kind_physics) :: quadjzx(3) ! current density quadrupole - qj*vxj*x*z - qj*vyj*x*z - qj*vzj*x*z
      real(kind_physics) :: current(3) ! net current density sum
!        real(kind_physics) :: g
      real(kind_physics) :: bmax
   end type t_tree_node_interaction_data
   integer, private, parameter :: nprops_tree_node_interaction_data = 20

contains

   !>
   !> Writes particle interaction data and results to a VTK file.
   !>
   subroutine vtk_write_particle_data_results(d, r, vtkf)
      use module_vtk
      implicit none

      type(t_particle_data), intent(in)    :: d(:)
      type(t_particle_results), intent(in)    :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf

      call vtkf%write_data_array("q", d(:)%q)
      call vtkf%write_data_array("v", d(:)%v(1), d(:)%v(2), d(:)%v(3))
      call vtkf%write_data_array("m", d(:)%m)
      call vtkf%write_data_array("gamma", d(:)%g)

      call vtkf%write_data_array("pot", r(:)%pot)
!        call vtkf%write_data_array("rho"    , r(:)%rho                                  )
      call vtkf%write_data_array("E", r(:)%e(1), r(:)%e(2), r(:)%e(3))
      call vtkf%write_data_array("A", r(:)%A(1), r(:)%A(2), r(:)%A(3))
!        call vtkf%write_data_array("dxA"    , r(:)%dxA(1) , r(:)%dxA(2) , r(:)%dxA(3)   )
!        call vtkf%write_data_array("dyA"    , r(:)%dyA(1) , r(:)%dyA(2) , r(:)%dyA(3)   )
!        call vtkf%write_data_array("dzA"    , r(:)%dzA(1) , r(:)%dzA(2) , r(:)%dzA(3)   )
      call vtkf%write_data_array("B", r(:)%B(1), r(:)%B(2), r(:)%B(3))
      call vtkf%write_data_array("J", r(:)%J(1), r(:)%J(2), r(:)%J(3))
!        call vtkf%write_data_array("Jirr"   , r(:)%Jirr(1), r(:)%Jirr(2), r(:)%Jirr(3)  )

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
!        call vtkf%write_data_array("current density", d(:)%current   )
   end subroutine vtk_write_node_interaction_data

   !>
   !> Creates and registers interaction-specific MPI-types
   !> is automatically called from register_libpepc_mpi_types()
   !>
   subroutine register_interaction_specific_mpi_types(MPI_TYPE_particle_data, MPI_TYPE_tree_node_interaction_data, MPI_TYPE_particle_results)
      use mpi
      implicit none
      integer, intent(out) :: MPI_TYPE_particle_data, MPI_TYPE_tree_node_interaction_data, MPI_TYPE_particle_results

      integer, parameter :: max_props = nprops_particle_data + nprops_particle_results + nprops_tree_node_interaction_data ! maxval([..]) would be enough, but ifort does notlike that

      integer :: ierr
      ! address calculation
      integer, dimension(1:max_props) :: blocklengths, types
      integer(KIND=MPI_ADDRESS_KIND), dimension(1:max_props) :: displacements
      integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
      ! dummies for address calculation
      type(t_particle_data)    :: dummy_particle_data
      type(t_particle_results) :: dummy_particle_results
      type(t_tree_node_interaction_data)   :: dummy_tree_node_interaction_data

      ! register particle data type
      blocklengths(1:nprops_particle_data) = [1, 3, 1, 1]
      types(1:nprops_particle_data) = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS]
      call MPI_GET_ADDRESS(dummy_particle_data,   address(0), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_data%q, address(1), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_data%v, address(2), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_data%m, address(3), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_data%g, address(4), ierr)  !&
      displacements(1:nprops_particle_data) = address(1:nprops_particle_data) - address(0)
      call MPI_TYPE_CREATE_STRUCT(nprops_particle_data, blocklengths, displacements, types, MPI_TYPE_particle_data, ierr)
      call MPI_TYPE_COMMIT(MPI_TYPE_particle_data, ierr)

      ! register results data type
      blocklengths(1:nprops_particle_results) = [1, 3, 3, 3, 3, 3, 3, 3, 3]
      types(1:nprops_particle_results) = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
                                          MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
                                          MPI_KIND_PHYSICS]

      call MPI_GET_ADDRESS(dummy_particle_results     , address(0), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%pot , address(1), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%E   , address(2), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%A   , address(3), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%B   , address(4), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%J   , address(5), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%Jirr, address(6), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%dxA , address(7), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%dyA , address(8), ierr)  !&
      call MPI_GET_ADDRESS(dummy_particle_results%dzA , address(9), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_particle_results%dxE , address(7), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_particle_results%dyE , address(8), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_particle_results%dxxA, address(12), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_particle_results%dxyA, address(13), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_particle_results%dyyA, address(14), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_particle_results%rho , address(9), ierr)  !&

      displacements(1:nprops_particle_results) = address(1:nprops_particle_results) - address(0)
      call MPI_TYPE_CREATE_STRUCT(nprops_particle_results, blocklengths, displacements, types, MPI_TYPE_particle_results, ierr)
      call MPI_TYPE_COMMIT(MPI_TYPE_particle_results, ierr)

      ! register multipole data type
      blocklengths(1:nprops_tree_node_interaction_data) = [3, 1, 1, 3, 3, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1]
      types(1:nprops_tree_node_interaction_data) = [MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
                                                    MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
                                                    MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, &
                                                    MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS, MPI_KIND_PHYSICS]
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data,            address(0),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%coc,        address(1),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%charge,     address(2),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%abs_charge, address(3),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%dip,        address(4),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quad,       address(5),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%xyquad,     address(6),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%yzquad,     address(7),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%zxquad,     address(8),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%monoj,      address(9),  ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%dipjx,      address(10), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%dipjy,      address(11), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%dipjz,      address(12), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quadjx,     address(13), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quadjy,     address(14), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quadjz,     address(15), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quadjxy,    address(16), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quadjyz,    address(17), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%quadjzx,    address(18), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%g,          address(19), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%current,    address(19), ierr)  !&
      call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%bmax,       address(20), ierr)  !&
!        call MPI_GET_ADDRESS(dummy_tree_node_interaction_data%jj,         address(9), ierr)  !&

      displacements(1:nprops_tree_node_interaction_data) = address(1:nprops_tree_node_interaction_data) - address(0)
      call MPI_TYPE_CREATE_STRUCT(nprops_tree_node_interaction_data, blocklengths, displacements, types, MPI_TYPE_tree_node_interaction_data, ierr)
      call MPI_TYPE_COMMIT(MPI_TYPE_tree_node_interaction_data, ierr)
   end subroutine register_interaction_specific_mpi_types
end module module_interaction_specific_types
