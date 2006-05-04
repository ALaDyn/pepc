 
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

  !  particle push in space

  do p=ips,ipf
     x(p)=x(p)+ux(p)*delt
     y(p)=y(p)+uy(p)*delt
     z(p)=z(p)+uz(p)*delt
  end do

end subroutine push
