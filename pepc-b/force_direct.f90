!  ====================================================================
!
!                              FORCE_DIRECT
!
!   Direct PP force sum for error estimates
!
!  ====================================================================

subroutine force_direct(n,ntest,x,y,z,q, list, eps, const, ex, ey, ez, pot)
  implicit none
  real, intent(in) :: eps, const
  integer, intent(in) :: n, ntest
  integer, intent(in), dimension(n) :: list 
  real, intent(in), dimension(n) :: x, y, z, q  ! coords and charge 
  real,  dimension(ntest) :: ex, ey, ez, pot  ! fields and potential to return

  real :: eps2, d, d3, dx, dy, dz
  integer :: i,j,k
!  write (*,*) 'Direct sum params: ',eps,const
!  write (*,'(a10,a20/(i6,4f15.3))') 'DIRECT | ','Initial buffers: ',(i, x(i), y(i), z(i), q(i),i=1,n) 
 eps2=eps**2
  ex(1:ntest) = 0.
  ey(1:ntest) = 0.
  ez(1:ntest) = 0.
  pot(1:ntest) = 0.

  !  PP contribution from simulation volume

  do  k=1,ntest
     i=list(k)
!     write (*,*) i,'x_i=',x(i)
     do  j=1,n
        if (j.ne.i) then
           dx=x(i)-x(j)
           dy=y(i)-y(j)
           dz=z(i)-z(j)
           d=sqrt(dx**2+dy**2+dz**2+eps2)
           d3=d**3
           ex(k) = ex(k) + const*q(j)*dx/d3
           ey(k) = ey(k) + const*q(j)*dy/d3
           ez(k) = ez(k) + const*q(j)*dz/d3
           pot(k) = pot(k) + const*q(j)/d
!          write (*,'(i5,a5,f12.3,a5,f12.3,a5,f12.3)') j,'q_j=',q(j),' x_j=',x(j), ' d=',d
        endif
     end do
  end do


end subroutine force_direct
