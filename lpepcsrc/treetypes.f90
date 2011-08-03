
module treetypes

  ! Data structure for user-defined variables that are directly involved into the force calculation
  type calc_force_params
    real    :: eps
    real    :: force_const
    integer :: force_law
  end type calc_force_params

  ! Data structure for shipping results
  type results
     real*8 :: Ex
     real*8 :: Ey
     real*8 :: Ez
     real*8 :: pot
     real*8 :: work
     integer :: label
  end type results


  ! Hash table datatype - 36 bytes per entry
  type hash
     integer   :: node          ! Address of particle/pseudoparticle data
     integer*8 :: key           ! Key
     integer   :: link          ! Pointer to next empty address in table in case of collision
     integer   :: leaves        ! # leaves contained within twig (=1 for leaf, npart for root)
     integer   :: childcode     ! Byte code indicating position of children (twig node); particle label (leaf node)
     integer*8 :: next          ! Pointer to next key to examine in tree-walk
     integer   :: owner         ! Node owner (for branches)
  end type hash


 ! Data structure for shipping single particles
  type particle
     real*8 :: x    ! coords
     real*8 :: y
     real*8 :: z
     real*8 :: ux    ! momenta
     real*8 :: uy
     real*8 :: uz 
     real*8 :: q     ! charge
     real*8 :: work  ! work load from force sum
     integer*8 :: key           ! Key
     integer :: label    ! label
     integer :: pid      ! owner
  end type particle


  ! Data structure for shipping multiple moments of child nodes
  type multipole
     integer*8 :: key     ! key
     integer   :: byte    ! byte code
     integer   :: leaves  ! # leaves contained
     integer :: owner    ! owner where multipole resides
     integer*8 :: next    ! next key on walk
     real*8 :: q        ! net charge sum
     real*8 :: absq     !  absolute charge sum
     real*8 :: xcoc     ! centre of charge
     real*8 :: ycoc
     real*8 :: zcoc
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
  end type multipole


end module treetypes
