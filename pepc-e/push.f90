 
!  ===============================================================
!
!                           PMOVE
!
!   Update particle positions - used with leap-frog scheme
!
!  ===============================================================

subroutine push_x(ips,ipf,delt)

  use physvars
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: i,p
  real :: gamma

  !  relativistic particle push in space

  do p=ips,ipf
     x(p)=x(p)+ux(p)*delt
     y(p)=y(p)+uy(p)*delt
     z(p)=z(p)+uz(p)*delt
  end do

end subroutine push_x
