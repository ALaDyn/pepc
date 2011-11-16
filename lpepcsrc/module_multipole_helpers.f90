!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all routines that are specific to a certain multipole expansion
!> i.e. shifting along the tree etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_multipole_helpers
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Data structure for storing multiple moments of tree nodes
      type t_multipole_data
        real*8 :: charge     ! net charge sum
        real*8 :: abs_charge !  absolute charge sum
        real*8 :: coc(3)     ! centre of charge
        real*8 :: dip(3)     ! dipole moment
        real*8 :: quad(3)    ! diagonal quadrupole moments
        real*8 :: xyquad     ! other quadrupole moments
        real*8 :: yzquad
        real*8 :: zxquad
      end type t_multipole_data

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

      !> TODO: remove this function
      subroutine multipole_to_data(m,d)
        use treetypes
        implicit none
        type(t_tree_node), intent(in) :: m
        type(t_multipole_data), intent(out) :: d

        d%charge     =  m%charge
        d%abs_charge =  m%abs_charge
        d%coc        =  m%coc
        d%dip        = [m%xdip,   m%ydip,   m%zdip]
        d%quad       = [m%xxquad, m%yyquad, m%zzquad]

        d%xyquad     =  m%xyquad
        d%yzquad     =  m%yzquad
        d%zxquad     =  m%zxquad

      end subroutine


      !> TODO: remove this function
      subroutine data_to_multipole(d,m)
        use treetypes
        implicit none
        type(t_tree_node), intent(inout) :: m
        type(t_multipole_data), intent(in) :: d

        m%charge     = d%charge
        m%abs_charge = d%abs_charge
        m%coc        = d%coc
        m%xdip       = d%dip(1)
        m%ydip       = d%dip(2)
        m%zdip       = d%dip(3)
        m%xxquad     = d%quad(1)
        m%yyquad     = d%quad(2)
        m%zzquad     = d%quad(3)

        m%xyquad     = d%xyquad
        m%yzquad     = d%yzquad
        m%zxquad     = d%zxquad

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Allocates all tree- and multipole-specific arrays,
      !> Initializes several estimations on maximum used memory
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_nodes_up(parent, children)
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


end module module_multipole_helpers
