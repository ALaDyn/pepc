!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates functions for accessing, manipulating, and verifying hash table data
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_mirror_boxes
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! lattice basis vectors
      real*8, public :: t_lattice_1(3) = [1, 0, 0] !< 1st vector of lattice basis
      real*8, public :: t_lattice_2(3) = [0, 1, 0] !< 2nd vector of lattice basis
      real*8, public :: t_lattice_3(3) = [0, 0, 1] !< 3rd vector of lattice basis
      real*8, public :: Lattice(3,3) !< holds the lattice transformation matrix
      real*8, public :: LatticeInv(3,3) !< holds the inverse lattice transformation matrix
      real*8, public :: LatticeCenter(3) !< holds the central point of the lattice box
      logical, public :: periodicity(3) = [.false., .false., .false.]  !< boolean switches for determining periodicity directions
      !> variables that should not be written to
      integer, public :: num_neighbour_boxes = 1
      integer, dimension(:,:), allocatable, public :: neighbour_boxes !dimensions in 3D case (3,27)
      !> is set to OR(periodicity), to be able to exit from all procedures as fast as possible
      !> if nothing is to be done
      logical, public :: do_periodic
      ! number of boxes to include into each direction
      integer, public :: periodicity_switches(3)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public calc_neighbour_boxes
      public init_movement_constraint
      public constrain_periodic
      public lattice_vect

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! true for lattices with axes parallel to cartesian system
      logical :: simplelattice

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates vector with respect to lattice base vectors
        !> @param[in] ijk lattice indices
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function lattice_vect(ijk)
          implicit none
          integer, intent(in) :: ijk(3)
          real*8 :: lattice_vect(3)

          lattice_vect = ijk(1)*t_lattice_1 + ijk(2)*t_lattice_2 + ijk(3)*t_lattice_3

        end function lattice_vect


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Prepares the list of neighbour boxes within ws
        !> stores their number in num_neighbour_boxes and their logical
        !> indices/coordinates in neighbour_boxes
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_neighbour_boxes(ws)
        implicit none
        integer :: i,j,k,idx
        integer, intent(in) :: ws

          if (allocated(neighbour_boxes)) deallocate(neighbour_boxes)

          do i = 1,3
            if (periodicity(i)) then
               periodicity_switches(i) = ws
             else
               periodicity_switches(i) = 0
             end if
          end do

          idx = 0

          num_neighbour_boxes = product(2*periodicity_switches+1)
          allocate(neighbour_boxes(3,num_neighbour_boxes))

          ! the boxes of the shells that surround the central box
          ! are ordered in a way, that a box and its central-symmetric counterpart are grouped together
          do i = -periodicity_switches(1),-1
            do j = -periodicity_switches(2),periodicity_switches(2)
              do k = -periodicity_switches(3),periodicity_switches(3)
                  idx = idx + 1
                  neighbour_boxes(:,idx) = [ i,  j ,  k]
                  idx = idx + 1
                  neighbour_boxes(:,idx) = [-i, -j , -k]
              end do
            end do
          end do

          i = 0
          do j = -periodicity_switches(2),-1
            do k = -periodicity_switches(3),periodicity_switches(3)
                idx = idx + 1
                neighbour_boxes(:,idx) = [ i,  j ,  k]
                idx = idx + 1
                neighbour_boxes(:,idx) = [-i, -j , -k]
            end do
          end do

          i = 0
          j = 0
          do k = -periodicity_switches(3),-1
            idx = idx + 1
            neighbour_boxes(:,idx) = [ i,  j ,  k]
            idx = idx + 1
            neighbour_boxes(:,idx) = [-i, -j , -k]
          end do

          neighbour_boxes(:,idx+1) = [0, 0, 0] ! center box is put to the back of the boxlist for easier iteration

        end subroutine calc_neighbour_boxes


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Initializes transformation matrices between cartesian system and lattice basis
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_movement_constraint
          use module_math_tools
          implicit none

          integer :: i,j

          Lattice(1,:) = t_lattice_1
          Lattice(2,:) = t_lattice_2
          Lattice(3,:) = t_lattice_3
          LatticeInv = Inverse3(Lattice)

          ! simplify the movement constraint if the lattice is really simple
          simplelattice = .true.
          do i = 1,3
            do j = 1,3
              if (i.ne.j) then
                simplelattice = simplelattice .and. (Lattice(i,j) == 0)
              else
                simplelattice = simplelattice .and. (Lattice(i,j) >  0)
              end if
            end do
          end do


        end subroutine init_movement_constraint


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Shifts all particles that left the original box
        !> back into it
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine constrain_periodic(x, y, z, np_local)
            use module_math_tools

            implicit none

            real*8, intent(inout), dimension(1:np_local) :: x,y,z
            integer, intent(in) :: np_local

            integer :: p !< loop variable
            real*8 :: lattice_coord(3), real_coord(3)

            if (simplelattice) then

                do p = 1,np_local

                  if (periodicity(1)) x(p) = modulo(x(p), t_lattice_1(1))
                  if (periodicity(2)) y(p) = modulo(y(p), t_lattice_2(2))
                  if (periodicity(3)) z(p) = modulo(z(p), t_lattice_3(3))

                end do

            else
              ! the lattice axes are not parallel to the outer cartesian axes
                do p = 1,np_local

                    lattice_coord = matmul([x(p), y(p), z(p)], LatticeInv)

                    if (any(lattice_coord > 1.D+0) .or. any(lattice_coord < 0.D+0)) then

                      lattice_coord = modulo(lattice_coord, 1.D+0)

                      real_coord = matmul(lattice_coord, Lattice)

                      if (periodicity(1)) x(p) = real_coord(1)
                      if (periodicity(2)) y(p) = real_coord(2)
                      if (periodicity(3)) z(p) = real_coord(3)

                    end if
                end do
            end if

        end subroutine

end module module_mirror_boxes
