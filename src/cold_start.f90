! ==============================================
!
!                COLD_START
!
!  Initialise set of particles with zero velocity
!
! ==============================================

subroutine cold_start(i1,n)

  use treevars
  use utils
  implicit none


  integer :: i1,n

  ux(i1:i1+n) = 0.
  uy(i1:i1+n) = 0.
  uz(i1:i1+n) = 0.

end subroutine cold_start
