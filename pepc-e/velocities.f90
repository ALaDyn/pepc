!  ===================================================================
!
!                              VELOCITIES
!
!   $Revision$
!
!   Calculate velocities from accelerations
! 
!   apply thermodynamic constraint according to ensemble
!
!
!
!  ===================================================================


subroutine velocities(p_start,p_finish,delta_t)


  use physvars
  use utils
  implicit none
  include 'mpif.h'

  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer p, i, ne_loc, ierr
  real, dimension(pe_capacity) :: accx, accy, accz
  real :: sum_vxe, sum_vye, sum_vze, sum_v2e, sum_2ve, Te0, Te_uncor, Ti0, Ti_uncor, chie, chii
  real :: sum_vxi, sum_vyi, sum_vzi, sum_v2i, sum_2vi, mass_eqm
  real :: global_v2e, global_v2i, gammah, delta_Te, delta_Ti, Te_loc

  !  Available ensemble modes
  !      1 = NVE - total energy conserved

! Accelerations
  do i=1,npp
     accx(i) = q(i)*ex(i)/m(i)
     accy(i) = q(i)*ey(i)/m(i)
     accz(i) = q(i)*ez(i)/m(i)
  end do

  ! unconstrained motion by default (scheme=1)

  do p = p_start, p_finish
     ux(p) = ux(p) + delta_t * accx(p)
     uy(p) = uy(p) + delta_t * accy(p)
     uz(p) = uz(p) + delta_t * accz(p)
  end do


end subroutine velocities




