 
!  ===============================================================
!
!                           PMOVE
!
!   Update particle positions - used with leap-frog scheme
!
!  ===============================================================

subroutine push(ips,ipf,delt)

  use treevars
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: i,p
  real :: gamma

  !  relativistic particle push in space

  do p=ips,ipf
     gamma = sqrt(1.0 + ux(p)**2 + uy(p)**2 + uz(p)**2)
     x(p)=x(p)+ux(p)/gamma*delt
     y(p)=y(p)+uy(p)/gamma*delt
     z(p)=z(p)+uz(p)/gamma*delt

  end do

end subroutine push
