 
!  ===============================================================
!
!                           PMOVE
!
!   Update particle positions - used with leap-frog scheme
!
!  ===============================================================

subroutine push_x(ips,ipf,delt)

  use physvars
  use treevars
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: p
  real :: gamma, r_glue2, dx, dy, dz, dr2

  r_glue2 = (max(xl,yl,zl)*glue_radius)**2

  !  relativistic particle push in space

  do p=ips,ipf
! find radius from plasma_centre
     dx = x(p)-plasma_centre(1)
     dy = y(p)-plasma_centre(2)
     dz = z(p)-plasma_centre(3)
     dr2 = dx**2+dy**2+dz**2
     if (dr2 < r_glue2 ) then 
       gamma = sqrt(1.0 + ux(p)**2 + uy(p)**2 + uz(p)**2)
	gamma=1.0
       x(p)=x(p)+ux(p)/gamma*delt
       if (idim > 1) y(p)=y(p)+uy(p)/gamma*delt
       if (idim == 3) z(p)=z(p)+uz(p)/gamma*delt
     else
!	if (mod(me,200).eq.0) write(*,*) "particle ",p," glued at ",sqrt(dr2)
	! leave particle where it is (should flag it to remove from force lists)
     endif
  end do

end subroutine push_x