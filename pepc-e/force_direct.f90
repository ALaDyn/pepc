!  ====================================================================
!
!                              FORCE_DIRECT
!
!   Direct PP force sum for error estimates
!
!  ====================================================================

subroutine force_direct(n,x,y,z,q, eps, const, ex, ey, ez, pot)
  implicit none
  real, intent(in) :: eps, const
  integer, intent(in) :: n
  real*8, intent(in), dimension(n) :: x, y, z, q  ! coords and charge 
  real*8,  dimension(n) :: ex, ey, ez, pot  ! fields and potential to return

  real*8 :: eps2, d, d3, dx, dy, dz
  integer :: i,j
!  write (*,*) 'Direct sum params: ',eps,const
!  write (*,'(a10,a20/(i6,4f15.3))') 'DIRECT | ','Initial buffers: ',(i, x(i), y(i), z(i), q(i),i=1,n) 
 eps2=eps**2
  ex(1:n) = 0.
  ey(1:n) = 0.
  ez(1:n) = 0.
  pot(1:n) = 0.

  !  PP contribution from simulation volume

  do  i=1,n
!     write (*,*) i,'x_i=',x(i)
     do  j=1,n
        if (j.ne.i) then
           dx=x(i)-x(j)
           dy=y(i)-y(j)
           dz=z(i)-z(j)
           d=sqrt(dx**2+dy**2+dz**2+eps2)
           d3=d**3
           ex(i) = ex(i) + const*q(j)*dx/d3
           ey(i) = ey(i) + const*q(j)*dy/d3
           ez(i) = ez(i) + const*q(j)*dz/d3
           pot(i) = pot(i) + const*q(j)/d
!          write (*,'(i5,a5,f12.3,a5,f12.3,a5,f12.3)') j,'q_j=',q(j),' x_j=',x(j), ' d=',d
        endif
     end do
  end do


end subroutine force_direct



