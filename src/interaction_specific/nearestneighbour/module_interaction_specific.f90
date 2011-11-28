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
         real*8 :: q
         integer :: label !< TODONN: remove this nonsense
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 2

      !> Data structure for results
      type t_particle_results
         real*8 :: maxdist2      !< maxval(dist2)
         integer :: maxidx       !< maxloc(dist2)
         integer*8:: neighbour_keys(num_neighbour_particles)
         real*8 :: dist2(num_neighbour_particles)
         real*8 :: work
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 5

      !> Data structure for storing multiple moments of tree nodes
      type t_multipole_data
        real*8 :: coc(3)     !< center of charge
        integer :: label     !< particle label for leaf nodes (nonsense)
      end type t_multipole_data
      integer, private, parameter :: nprops_multipole_data = 2

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
      subroutine register_interaction_specific_mpi_types(mpi_type_particle_data, mpi_type_multipole_data, mpi_type_particle_results)
        implicit none
        include 'mpif.h'
        integer, intent(out) :: mpi_type_particle_data, mpi_type_multipole_data, mpi_type_particle_results

        integer, parameter :: max_props = maxval([nprops_particle_data, nprops_particle_results, nprops_multipole_data])

        integer :: ierr
        ! address calculation
        integer, dimension(1:max_props) :: blocklengths, displacements, types
        integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
        ! dummies for address calculation
        type(t_particle_data)    :: dummy_particle_data
        type(t_particle_results) :: dummy_particle_results
        type(t_multipole_data)   :: dummy_multipole_data

        ! register particle data type
        blocklengths(1:nprops_particle_data)  = [1, 1]
        types(1:nprops_particle_data)         = [MPI_REAL8, MPI_INTEGER]
        call MPI_GET_ADDRESS( dummy_particle_data,       address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%q,     address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%label, address(2), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [1, 1, num_neighbour_particles, num_neighbour_particles, 1]
        types(1:nprops_particle_results)         = [MPI_REAL8, MPI_INTEGER, MPI_INTEGER8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_results,      address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%maxdist2,       address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%maxidx,         address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%neighbour_keys, address(3), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%dist2,          address(4), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%work,           address(5), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_multipole_data)  = [3, 1]
        types(1:nprops_multipole_data)         = [MPI_REAL8, MPI_INTEGER]
        call MPI_GET_ADDRESS( dummy_multipole_data,            address(0), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%coc,        address(1), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%label,      address(2), ierr )
        displacements(1:nprops_multipole_data) = int(address(1:nprops_multipole_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_multipole_data, blocklengths, displacements, types, mpi_type_multipole_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_multipole_data, ierr)

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
        type(t_multipole_data), intent(out) :: multipole

        multipole = t_multipole_data(particle_pos, particle%label)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_multipole_data), intent(out) :: parent
        type(t_multipole_data), intent(in) :: children(:)

        integer :: nchild

        nchild = size(children)

        parent%label = nchild

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

        ! TODONN: muss hier was getan werden?

      end subroutine



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> clears results datatype
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elemental subroutine results_clear(res)
        implicit none
        type(t_particle_results), intent(out) :: res
        real*8 :: realdummy

        res%maxdist2       = huge(realdummy)
        res%maxidx         = 1
        res%neighbour_keys = 0
        res%dist2          = huge(realdummy)
        res%work           = 1.

      end subroutine


end module module_interaction_specific
