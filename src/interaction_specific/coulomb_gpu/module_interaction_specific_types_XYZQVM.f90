! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all types that are specific to a certain interaction
!> all subroutines and types within this module are obligatory
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_interaction_specific_types
      use pthreads_stuff
      use module_atomic_ops, only: t_atomic_int
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !> Data structure for storing multiple moments of tree nodes
      type t_tree_node_interaction_data
        real*8 :: coc(3)     ! centre of charge
        real*8 :: charge     ! net charge sum
        real*8 :: abs_charge ! absolute charge sum
        real*8 :: dip(3)     ! dipole moment
        real*8 :: quad(3)    ! diagonal quadrupole moments
        real*8 :: xyquad     ! other quadrupole moments
        real*8 :: yzquad
        real*8 :: zxquad
        real*8 :: bmax
      end type t_tree_node_interaction_data
      integer, private, parameter :: nprops_tree_node_interaction_data = 9

      !> Data structure for storing interaction-specific particle data
      type t_particle_data
         real*8 :: q
         real*8 :: v(3)
         real*8 :: m
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 3

      !> Data structure for shipping results
      type t_particle_results
         real*8, dimension(3) :: e
         real*8 :: pot
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 2

      type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results([0., 0., 0.], 0.)

      !> Data structure for thread local storage of single particles
      !> This includes lists of the interaction partners
      integer, public, parameter :: MAX_IACT_PARTNERS = 2 * 256 * 2 * 2 ! length of vectors for GPU, multiples of 256 because of gang size?
      integer, public, parameter :: ACC_QUEUE_LENGTH  = 16              ! how may lists we accept from the workers at one time
      integer, public, parameter :: GPU_STREAMS = 8                     ! number of GPU lists to prepare, multiple streams may require than one

      !> Thread local data structure to store extra interaction information
      ! thread local, since we do not want to ship this via MPI
      type t_particle_thread
         real*8 :: x(1:3)                             !< coordinates
         real*8, pointer :: work                      !< work load from force sum
         integer*8 :: key                             !< particle key, i.e. key on highgest tree level
         integer*8 :: node_leaf                       !< key of corresponding leaf (tree node)
         integer*8 :: label                           !< particle label (only for diagnostic purposes, can be used freely by the frontend
         integer :: pid                               !< particle owner
         type(t_particle_data) :: data                !< real physics (charge, etc.)
         type(t_particle_results), pointer :: results !< results of calc_force_etc and companions
         integer :: queued = -1                       !< counter for actually queued particles
         real*8, pointer :: partner(:)                !< 1-D array to store all info about interaction partner (delta, tree_node_interaction_data)
         logical, pointer :: leaf(:)                  !< 1-D array to store particle_is_leaf flag
         integer*8 :: my_idx = -1
      end type t_particle_thread

      type, public :: t_acc_queue_entry
         type(t_particle_thread) :: particle
         type(t_particle_results), pointer :: results
         real*8 :: eps, pen
         logical :: entry_valid
         integer :: queued
         real*8, pointer :: partner(:)
      end type

      type :: t_acc
         type(t_atomic_int), pointer :: thread_status
         integer :: processor_id
         type(t_pthread_with_type) :: acc_thread
         type(t_acc_queue_entry) :: acc_queue(ACC_QUEUE_LENGTH)
         type(t_atomic_int), pointer :: q_top, q_bottom, q_len
      end type
      type(t_acc) :: acc

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Creates and registers interaction-specific MPI-types
      !> is automatically called from register_libpepc_mpi_types()
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        blocklengths(1:nprops_particle_data)  = [1, 3, 1]
        types(1:nprops_particle_data)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_data,         address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%q,       address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%v,       address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%m,       address(3), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [3, 1]
        types(1:nprops_particle_results)         = [MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_results,      address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%e,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%pot,  address(2), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_tree_node_interaction_data)  = [3, 1, 1, 3, 3, 1, 1, 1, 1]
        types(1:nprops_tree_node_interaction_data)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data,            address(0), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%coc,        address(1), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%charge,     address(2), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%abs_charge, address(3), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%dip,        address(4), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%quad,       address(5), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%xyquad,     address(6), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%yzquad,     address(7), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%zxquad,     address(8), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%bmax,       address(9), ierr )
        displacements(1:nprops_tree_node_interaction_data) = int(address(1:nprops_tree_node_interaction_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_tree_node_interaction_data, blocklengths, displacements, types, MPI_TYPE_tree_node_interaction_data, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_interaction_data, ierr)

      end subroutine register_interaction_specific_mpi_types


end module module_interaction_specific_types