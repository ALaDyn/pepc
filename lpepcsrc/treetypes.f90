
module treetypes

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

  ! Data structure for shipping results
  integer, parameter :: nprops_particle_results = 5       ! # results to ship
  type t_particle_results
     real*8, dimension(3) :: e
     real*8 :: pot
     real*8 :: work
  end type t_particle_results

  type(t_particle_results), parameter :: EMPTY_PARTICLE_RESULTS = t_particle_results([0., 0., 0.], 0., 0.)


 ! Data structure for shipping single particles
  integer, parameter :: nprops_particle = 11 ! # particle properties to ship
  type t_particle
     real*8 :: x(1:3)    ! coords
     real*8 :: u(1:3)    ! momenta
     real*8 :: q     ! charge
     real*8 :: work  ! work load from force sum
     integer*8 :: key           ! Key
     integer :: label    ! label
     integer :: pid      ! owner
  end type t_particle


  ! Data structure for shipping multiple moments of child nodes
  integer, parameter :: nprops_multipole = 22 ! Number of multipole properties to ship
  type t_multipole
     integer*8 :: key     ! key
     integer   :: byte    ! byte code
     integer   :: leaves  ! # leaves contained
     integer :: owner     ! owner where multipole resides
     real*8 :: charge     ! net charge sum
     real*8 :: abs_charge !  absolute charge sum
     real*8 :: coc(3)       ! centre of charge
     integer :: level     ! level of tree node (for fast lookup instead of calculation) - stored in close proximity to coc-coordinates for better caching posibilities in walk
     real*8 :: xdip     ! dipole moment
     real*8 :: ydip
     real*8 :: zdip
     real*8 :: xxquad   ! quadrupole moment
     real*8 :: yyquad
     real*8 :: zzquad
     real*8 :: xyquad
     real*8 :: yzquad
     real*8 :: zxquad
     real*8 :: xshift
     real*8 :: yshift
     real*8 :: zshift
  end type t_multipole


end module treetypes
