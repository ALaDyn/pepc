! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars

  implicit none
  include 'mpif.h'

  config: select case(system_config)

  case(1)              ! Set up particles according to geometry

  case(2)
     call special_start(ispecial)

  end select config

end subroutine configure

