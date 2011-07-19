 
!  ===============================================================
!
!                           PUSH_NONREL
!
!   Nonrelativistic particle position update - used with leap-frog scheme
!
!  ===============================================================

subroutine push_nonrel(ips,ipf,delt)

  use physvars
  use treevars
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: p


  do p=ips,ipf

     x(p)=x(p)+ux(p)*delt
     y(p)=y(p)+uy(p)*delt
     z(p)=z(p)+uz(p)*delt

  end do

end subroutine push_nonrel
