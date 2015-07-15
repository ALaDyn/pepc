! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2015 Juelich Supercomputing Centre, 
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
!> Contains all routines that are specific to a certain multipole expansion
!> i.e. shifting along the tree etc.
!>
module module_interaction_specific_types
      implicit none

      !> Data structure for storing interaction-specific particle data
      type t_particle_data
         real*8 :: alpha(3)    ! vorticity or better: alpha = vorticity * volume
         real*8 :: x_rk(3)     ! position temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
         real*8 :: alpha_rk(3) ! vorticity temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
         real*8 :: u_rk(3)     ! velocity temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
         real*8 :: af_rk(3)    ! vorticity RHS temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 5

      !> Data structure for shipping results
      type t_particle_results
         real*8, dimension(3) :: u   ! velocities
         real*8, dimension(3) :: af  ! RHS for vorticity/alpha ODE
         real*8 :: div               ! divergence
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 3

      type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results([0., 0., 0.], [0., 0., 0.], 0.)

      !> Data structure for storing multiple moments of tree nodes
      type t_multipole_moments
        real*8 :: abs_charge ! absolute charge sum
        real*8 :: chargex    ! 3D monopole = 3 entries
        real*8 :: chargey
        real*8 :: chargez
        real*8 :: xdip1      ! 3D dipole = 3*3 entries
        real*8 :: ydip1
        real*8 :: zdip1
        real*8 :: xdip2
        real*8 :: ydip2
        real*8 :: zdip2
        real*8 :: xdip3
        real*8 :: ydip3
        real*8 :: zdip3
        real*8 :: xxquad1    ! 3D quadrupole = 3*3*3 entries (minus symmetries)
        real*8 :: xyquad1
        real*8 :: xzquad1
        real*8 :: yzquad1
        real*8 :: yyquad1
        real*8 :: zzquad1
        real*8 :: xxquad2
        real*8 :: xyquad2
        real*8 :: xzquad2
        real*8 :: yzquad2
        real*8 :: yyquad2
        real*8 :: zzquad2
        real*8 :: xxquad3
        real*8 :: xyquad3
        real*8 :: xzquad3
        real*8 :: yzquad3
        real*8 :: yyquad3
        real*8 :: zzquad3
        real*8 :: bmax
      end type t_multipole_moments
      integer, private, parameter :: nprops_multipole_moments = 32

      type t_local_coefficients
        ! empty for now
      end type t_local_coefficients

      type(t_local_coefficients), parameter :: EMPTY_LOCAL_COEFFICIENTS = t_local_coefficients()

      contains

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
        blocklengths(1:nprops_particle_data)  = [3,3,3,3,3]
        types(1:nprops_particle_data)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_data,       address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%alpha, address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%x_rk, address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%alpha_rk, address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%u_rk, address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%af_rk, address(5), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [3, 3, 1]
        types(1:nprops_particle_results)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_results,      address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%u,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%af,   address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%div,  address(3), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_multipole_moments)  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        types(1:nprops_multipole_moments)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_multipole_moments,            address(0), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%abs_charge, address(1), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%chargex,    address(2), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%chargey,    address(3), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%chargez,    address(4), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xdip1,      address(5), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%ydip1,      address(6), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zdip1,      address(7), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xdip2,      address(8), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%ydip2,      address(9), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zdip2,      address(10), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xdip3,      address(11), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%ydip3,      address(12), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zdip3,      address(13), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xxquad1,    address(14), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xyquad1,    address(15), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xzquad1,    address(16), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yzquad1,    address(17), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yyquad1,    address(18), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zzquad1,    address(19), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xxquad2,    address(20), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xyquad2,    address(21), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xzquad2,    address(22), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yzquad2,    address(23), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yyquad2,    address(24), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zzquad2,    address(25), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xxquad3,    address(26), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xyquad3,    address(27), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%xzquad3,    address(28), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yzquad3,    address(29), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%yyquad3,    address(30), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%zzquad3,    address(31), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_moments%bmax,       address(32), ierr )
        displacements(1:nprops_multipole_moments) = int(address(1:nprops_multipole_moments) - address(0))
        call MPI_TYPE_STRUCT( nprops_multipole_moments, blocklengths, displacements, types, MPI_TYPE_multipole_moments, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_multipole_moments, ierr)
      end subroutine register_interaction_specific_mpi_types
end module module_interaction_specific_types
