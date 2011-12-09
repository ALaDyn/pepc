!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all routines that are specific to a certain multipole expansion
!> i.e. shifting along the tree etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_interaction_specific
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

      integer, parameter :: num_neighbour_particles = 50

      !> Data structure for storing interaction-specific particle data
      type t_particle_data
         real*8 :: q                 !< charge
         real*8 :: v(3)              !< velocity (same time as x)
         real*8 :: v_and_half(3)     !< velocity (1/2 time step after x (t+1/2), for leap frog integrator)
         real*8 :: temperature
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 4

      !> Data structure for results
      type t_particle_results
         real*8 :: maxdist2       !< maxval(dist2)
         integer :: maxidx        !< maxloc(dist2)
         integer*8:: neighbour_nodes(num_neighbour_particles)
         real*8 :: dist2(num_neighbour_particles)
         real*8 :: dist_vector(3,num_neighbour_particles)                           ! distance_vectors from particle to neighbour with respact to periodic shift vector
         real*8 :: rho            !< density for sph
         real*8 :: h              !< smoothing-length for sph
         real*8 :: sph_force(1:3)
         real*8 :: temperature_change
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 9

      !> Data structure for storing multiple moments of tree nodes
      type t_tree_node_interaction_data
        real*8 :: coc(3)     !< center of charge
        real*8 :: q          !< charge (for particles)
        real*8 :: v(1:3)     !< velocity
        real*8 :: temperature
        real*8 :: rho        !< sph density
        real*8 :: h          !< sph smoothing-length
      end type t_tree_node_interaction_data
      integer, private, parameter :: nprops_tree_node_interaction_data = 6

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

        integer, parameter :: max_props = maxval([nprops_particle_data, nprops_particle_results, nprops_tree_node_interaction_data])

        integer :: ierr
        ! address calculation
        integer, dimension(1:max_props) :: blocklengths, displacements, types
        integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
        ! dummies for address calculation
        type(t_particle_data)    :: dummy_particle_data
        type(t_particle_results) :: dummy_particle_results
        type(t_tree_node_interaction_data)   :: dummy_tree_node_interaction_data

        ! register particle data type
        blocklengths(1:nprops_particle_data)  = [1, 3, 3, 1]
        types(1:nprops_particle_data)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_data,             address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%q,           address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%v,           address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%v_and_half,  address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%temperature, address(4), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [1, 1, num_neighbour_particles, num_neighbour_particles, 3*num_neighbour_particles, 1, 1, 3, 1]
        types(1:nprops_particle_results)         = [MPI_REAL8, MPI_INTEGER, MPI_INTEGER8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_results,                    address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%maxdist2,           address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%maxidx,             address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%neighbour_nodes,    address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%dist2,              address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%dist_vector,        address(5), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%rho,                address(6), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%h,                  address(7), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%sph_force,          address(8), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%temperature_change, address(9), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_tree_node_interaction_data)  = [3, 1, 3, 1, 1, 1]
        types(1:nprops_tree_node_interaction_data)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data,             address(0), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%coc,         address(1), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%q,           address(2), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%v,           address(3), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%temperature, address(4), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%rho,         address(5), ierr )
        call MPI_GET_ADDRESS( dummy_tree_node_interaction_data%h,           address(6), ierr )
        displacements(1:nprops_tree_node_interaction_data) = int(address(1:nprops_tree_node_interaction_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_tree_node_interaction_data, blocklengths, displacements, types, MPI_TYPE_tree_node_interaction_data, ierr )
        call MPI_TYPE_COMMIT( MPI_TYPE_tree_node_interaction_data, ierr)

      end subroutine register_interaction_specific_mpi_types

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Computes multipole properties of a single particle
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        implicit none
        real*8, intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole

        ! use velocity (v) at same time step as coordinate, not v_and_half
        multipole = t_tree_node_interaction_data(particle_pos,particle%q, particle%v, particle%temperature, -13._8, -13._8 )
        ! set rho to -13 as dummy.
        ! TODO: find a better place to store rho

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: nchild, i

        nchild = size(children)

        parent%coc = [0._8, 0._8, 0._8]
        parent%q = 0._8

        do i=1,nchild
          parent%coc = parent%coc + children(i)%coc
          parent%q   = parent%q   + children(i)%q
        end do

        parent%coc = parent%coc / nchild

        ! set velocity, temperature and rho for tree node to a nonsense value which may be recognised as nonsense when one tries to use them
        parent%v = [-13._8,-13._8,-13._8]
        parent%temperature = -13._8
        parent%rho = -13._8
        parent%h = -13._8

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> adds res2 to res1
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in) :: res2

        ! TODONN: ist das wirklich alles?
        res1 = res2

      end subroutine


end module module_interaction_specific
