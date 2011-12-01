!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all lpepc-specific types and routines for registering them to MPI
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module treetypes
  use module_interaction_specific
  implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: MPI_TYPE_particle_data,   &
                 MPI_TYPE_multipole_data,   &
                 MPI_TYPE_particle_results, &
                 MPI_TYPE_particle,         &
                 MPI_TYPE_tree_node

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, parameter, private :: real8_dummy = 1._8

      !> Data structure for user-defined variables that are directly involved into the force calculation
      type t_calc_force_params
        !TODO: make this eps2
        real    :: eps          = 0.0    !< short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
        real    :: force_const  = 1.0    !< force constant
        integer :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
        logical :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
        integer :: mac          = 0      !< selector for multipole acceptance criterion, mac==0: Barnes-Hut, currently unused
        real    :: theta        = 0.6    !< multipole opening angle
        real    :: theta2        = 0.6**2. !< (multipole opening angle)^2 - is set automatically by tree_aswalk_pthreads
        real*8  :: spatial_interaction_cutoff(3) = huge(real8_dummy) * [1., 1., 1.] !< all nodes, where any(abs(coc(1:3)-particle_position(1:3)) > spatial_interaction_cutoff(1:3) are ignored when calculating interactions
      end type t_calc_force_params

      !> Data structure for shipping single particles
      integer, private, parameter :: nprops_particle = 7 ! # particle properties to ship
      type t_particle
         real*8 :: x(1:3)    ! coords
         real*8 :: work  ! work load from force sum
         integer*8 :: key           ! Key
         integer :: label    ! label
         integer :: pid      ! owner
         type(t_particle_data) :: data ! real physics (charge, etc.)
         type(t_particle_results) :: results ! results of calc_force_etc companions
      end type t_particle

      ! Data structure for shipping multiple moments of child nodes
      integer, private, parameter :: nprops_tree_node = 5 ! Number of multipole properties to ship
      type t_tree_node
         integer*8 :: key     ! key
         integer   :: byte    ! byte code
         integer   :: leaves  ! # leaves contained
         integer   :: owner   ! owner where multipole resides
         type(t_multipole_data) :: m ! real physics
      end type t_tree_node


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
      !> Creates and registers lpepc-MPI types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine register_lpepc_mpi_types()
        use module_interaction_specific
        implicit none
        include 'mpif.h'
        integer, parameter :: max_props = maxval([nprops_particle, nprops_tree_node])

        integer :: ierr
        ! address calculation
        integer, dimension(1:max_props) :: blocklengths, displacements, types
        integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
        ! dummies for address calculation
        type(t_particle)  :: dummy_particle
        type(t_tree_node) :: dummy_tree_node

        ! first register the interaction-specific MPI types since they are embedded into the lpepc-datatypes
        call register_interaction_specific_mpi_types(MPI_TYPE_particle_data, MPI_TYPE_multipole_data, MPI_TYPE_particle_results)

        ! register particle type
        blocklengths(1:nprops_particle)  = [3, 1, 1, 1, 1, 1, 1]
        types(1:nprops_particle)         = [MPI_REAL8, MPI_REAL8, MPI_INTEGER8, MPI_INTEGER, MPI_INTEGER, MPI_TYPE_particle_data, MPI_TYPE_particle_results]
        call MPI_GET_ADDRESS( dummy_particle,          address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle%x,        address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle%work,     address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle%key,      address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle%label,    address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle%pid,      address(5), ierr )
        call MPI_GET_ADDRESS( dummy_particle%data,     address(6), ierr )
        call MPI_GET_ADDRESS( dummy_particle%results,  address(7), ierr )
        displacements(1:nprops_particle) = int(address(1:nprops_particle) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, MPI_TYPE_particle, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_particle, ierr)

        ! register tree_node type
        blocklengths(1:nprops_tree_node)  = [1, 1, 1, 1, 1]
        types(1:nprops_tree_node)         = [MPI_INTEGER8, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_TYPE_multipole_data]
        call MPI_GET_ADDRESS( dummy_tree_node,        address(0), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%key,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%byte,   address(2), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%leaves, address(3), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%owner,  address(4), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node%m    ,  address(5), ierr )
        displacements(1:nprops_tree_node) = int(address(1:nprops_tree_node) - address(0))
        call MPI_TYPE_STRUCT( nprops_tree_node, blocklengths, displacements, types, MPI_TYPE_tree_node, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_tree_node, ierr)

    end subroutine register_lpepc_mpi_types


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Deregisters lpepc- and interaction-specific MPI types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine free_lpepc_mpi_types()
        implicit none
        include 'mpif.h'
        integer :: ierr

        call MPI_TYPE_FREE( MPI_TYPE_tree_node,        ierr)
        call MPI_TYPE_FREE( MPI_TYPE_particle,         ierr)
        call MPI_TYPE_FREE( MPI_TYPE_particle_results, ierr)
        call MPI_TYPE_FREE( MPI_TYPE_multipole_data,   ierr)
        call MPI_TYPE_FREE( MPI_TYPE_particle_data,    ierr)

      end subroutine free_lpepc_mpi_types

end module treetypes
