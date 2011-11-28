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


      !> Data structure for storing interaction-specific particle data
      type t_particle_data
         real*8 :: alpha(3)   ! vorticity or better: alpha = vorticity * volume
      end type t_particle_data
      integer, private, parameter :: nprops_particle_data = 1

      !> Data structure for shipping results
      type t_particle_results
         real*8, dimension(3) :: u   ! velocities
         real*8, dimension(3) :: af  ! RHS for vorticity/alpha ODE
         real*8 :: work
      end type t_particle_results
      integer, private, parameter :: nprops_particle_results = 3

      type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results([0., 0., 0.], [0., 0., 0.], 1.)

      !> Data structure for storing multiple moments of tree nodes
      type t_multipole_data
        real*8 :: coc(3)     ! centre of charge
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
      end type t_multipole_data
      integer, private, parameter :: nprops_multipole_data = 32

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
        blocklengths(1:nprops_particle_data)  = [3]
        types(1:nprops_particle_data)         = [MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_data,       address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_data%alpha, address(1), ierr )
        displacements(1:nprops_particle_data) = int(address(1:nprops_particle_data) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_data, blocklengths, displacements, types, mpi_type_particle_data, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_data, ierr)

        ! register results data type
        blocklengths(1:nprops_particle_results)  = [3, 3, 1]
        types(1:nprops_particle_results)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_particle_results,      address(0), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%u,    address(1), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%af,   address(2), ierr )
        call MPI_GET_ADDRESS( dummy_particle_results%work, address(3), ierr )
        displacements(1:nprops_particle_results) = int(address(1:nprops_particle_results) - address(0))
        call MPI_TYPE_STRUCT( nprops_particle_results, blocklengths, displacements, types, mpi_type_particle_results, ierr )
        call MPI_TYPE_COMMIT( mpi_type_particle_results, ierr)

        ! register multipole data type
        blocklengths(1:nprops_multipole_data)  = [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        types(1:nprops_multipole_data)         = [MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, &
                                                  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
        call MPI_GET_ADDRESS( dummy_multipole_data,            address(0), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%coc,        address(1), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%abs_charge, address(2), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%chargex,    address(3), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%chargey,    address(4), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%chargez,    address(5), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xdip1,      address(6), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%ydip1,      address(7), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%zdip1,      address(8), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xdip2,      address(9), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%ydip2,      address(10), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%zdip2,      address(11), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xdip3,      address(12), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%ydip3,      address(13), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%zdip3,      address(14), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xxquad1,    address(15), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xyquad1,    address(16), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xzquad1,    address(17), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%yzquad1,    address(18), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%yyquad1,    address(19), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%zzquad1,    address(20), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xxquad2,    address(21), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xyquad2,    address(22), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xzquad2,    address(23), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%yzquad2,    address(24), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%yyquad2,    address(25), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%zzquad2,    address(26), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xxquad3,    address(27), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xyquad3,    address(28), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%xzquad3,    address(29), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%yzquad3,    address(30), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%yyquad3,    address(31), ierr )
        call MPI_GET_ADDRESS( dummy_multipole_data%zzquad3,    address(32), ierr )
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

        multipole = t_multipole_data(particle_pos, &
                                     particle%q,   &
                                 abs(particle%q),  &
                                     [0., 0., 0.], &
                                     [0., 0., 0.], &
                                      0., 0., 0.    )
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

        integer :: nchild, j

        real*8 :: shift(1:3)

        nchild = size(children)

        parent%charge     = SUM( children(1:nchild)%charge )
        parent%abs_charge = SUM( children(1:nchild)%abs_charge )

        ! centre of charge
        parent%coc        = [0., 0., 0.]

        do j=1,nchild
          parent%coc(1:3) = parent%coc(1:3) + ( children(j)%coc(1:3) * children(j)%abs_charge )
        end do

        parent%coc(1:3) = parent%coc(1:3) / parent%abs_charge

        ! multipole properties
        parent%dip    = [0., 0., 0.]
        parent%quad   = [0., 0., 0.]
        parent%xyquad = 0.
        parent%yzquad = 0.
        parent%zxquad = 0.

        do j=1,nchild
          shift(1:3) = parent%coc(1:3) - children(j)%coc

          ! dipole moment
          parent%dip = parent%dip + children(j)%dip - children(j)%charge*shift(1:3)

          ! quadrupole moment
          parent%quad(1:3) = parent%quad(1:3) + children(j)%quad(1:3) - 2*children(j)%dip(1:3)*shift(1:3) + children(j)%charge*shift(1:3)**2

          parent%xyquad = parent%xyquad + children(j)%xyquad - children(j)%dip(1)*shift(2) - children(j)%dip(2)*shift(1) + children(j)%charge*shift(1)*shift(2)
          parent%yzquad = parent%yzquad + children(j)%yzquad - children(j)%dip(2)*shift(3) - children(j)%dip(3)*shift(2) + children(j)%charge*shift(2)*shift(3)
          parent%zxquad = parent%zxquad + children(j)%zxquad - children(j)%dip(3)*shift(1) - children(j)%dip(1)*shift(3) + children(j)%charge*shift(3)*shift(1)
        end do
      end subroutine


end module module_interaction_specific
