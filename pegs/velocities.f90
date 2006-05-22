!  ===================================================================
!
!                              VELOCITIES
!
!   $Revision: 1.7 $
!
!   Calculate velocities from accelerations
! 
!
!
!
!  ===================================================================


subroutine velocities(p_start,p_finish,delta_t)

  use treevars
  use utils

  implicit none
  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer p, i, ne_loc

     ! unconstrained motion by default

     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * ex(p)
	uy(p) = uy(p) + delta_t * ey(p)
	uz(p) = uz(p) + delta_t * ez(p)
     end do


end subroutine velocities
