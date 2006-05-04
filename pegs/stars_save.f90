subroutine stars(delta_t)

  use treevars
  use utils
  use physvars

  implicit none
  real, intent(in) :: delta_t
  real :: dxstar,dystar,dzstar,dr
  real, dimension(2) :: ax_star, ay_star, az_star
  integer :: i,j,p


  ! Distance between stars
               dxstar = x_star(1)-x_star(2)
               dystar = y_star(1)-y_star(2)
               dzstar = z_star(1)-z_star(2)
               dr = dxstar*dxstar + dystar*dystar + dzstar*dzstar
               dr = sqrt(dr)

  !        call intlistind(npart+i,xs(1,i),xs(2,i),xs(3,i))
  !        call gravforcind(npart+i,xs(1,i),xs(2,i),xs(3,i),mstar(i))
 

  !   Accelerations:
             do i=1,ni
              do j=1,ni              
               if(i.ne.j) then
               ax_star(i) = (-1)**i*gamma*(m_star(j)*dxstar/(dr**3))
               ay_star(i) = (-1)**i*gamma*(m_star(j)*dystar/(dr**3)) 
               az_star(i) = (-1)**i*gamma*(m_star(j)*dzstar/(dr**3))          
 ! potential  
               pot_star(i) = -gamma*m_star(j)/dr
              end if
              end do
             end do

  !  Velocities
 
         do i = 1,ni
	   ux_star(i) = ux_star(i) + delta_t * ax_star(i)
           uy_star(i) = uy_star(i) + delta_t * ay_star(i)
	   uz_star(i) = uz_star(I) + delta_t * az_star(i)
         end do

   !  star push in space

  do p=1,ni
     x_star(p) = x_star(p) + ux_star(p)*delta_t
     y_star(p) = y_star(p) + uy_star(p)*delta_t
     z_star(p) = z_star(p) + uz_star(p)*delta_t
  end do
  write(16,*)x_star(1),y_star(1),z_star(1),x_star(2),y_star(2),z_star(2)


end subroutine stars
