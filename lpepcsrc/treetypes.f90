
module treetypes
  use module_multipole_helpers
  implicit none

  real*8, parameter, private :: real8_dummy = 1.

  !> Data structure for user-defined variables that are directly involved into the force calculation
  type t_calc_force_params
    real    :: eps          = 0.0    !< short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)
    real    :: force_const  = 1.0    !< force constant
    integer :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
    logical :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
    integer :: mac          = 0      !< selector for multipole acceptance criterion, mac==0: Barnes-Hut, currently unused
    real    :: theta        = 0.6    !< multipole opening angle
    real*8  :: spatial_interaction_cutoff(3) = huge(real8_dummy) * [1., 1., 1.] !< all nodes, where any(abs(coc(1:3)-particle_position(1:3)) > spatial_interaction_cutoff(1:3) are ignored when calculating interactions
  end type t_calc_force_params

  ! Data structure for shipping multiple moments of child nodes
  integer, parameter :: nprops_tree_node = 18 ! Number of multipole properties to ship
  type t_tree_node
     integer*8 :: key     ! key
     integer   :: byte    ! byte code
     integer   :: leaves  ! # leaves contained
     integer :: owner     ! owner where multipole resides
     type(t_multipole_data) :: m ! real physics
  end type t_tree_node


end module treetypes
