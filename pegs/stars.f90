subroutine stars(delta_t)

  use treevars
  use utils
  use physvars

  implicit none
  include 'mpif.h'


  real, intent(in) :: delta_t
  real*8 :: dxstar,dystar,dzstar,dr,dxpart,dypart,dzpart,drpart
  real*8 :: distx, disty, distz, dists
  real*8, dimension(2) :: ax_star, ay_star, az_star
  real*8 :: ax_partial, ay_partial, az_partial, pot_partial  ! partial sums
  real*8 :: ax_sum, ay_sum, az_sum, pot_sum  ! global sums
  integer :: i,j,p, ierr


  ! Forces from star on dust particles: add to field sums from disc self-forces
    
     do p = 1,npp     
      do j=1,nstar
       distx=x(p)-x_star(j)
       disty=y(p)-y_star(j)
       distz=z(p)-z_star(j)
       dists=sqrt(distx*distx+disty*disty+distz*distz+eps**2)
       ex(p)=ex(p)-gamma*m_star(j)*distx/(dists**3)
       ey(p)=ey(p)-gamma*m_star(j)*disty/(dists**3)
       ez(p)=ez(p)-gamma*m_star(j)*distz/(dists**3)
       pot(p) = pot(p)-gamma*m_star(j)/dists  ! potential contribution
      end do
     end do


  !
  ! STAR  STAR  STAR  STAR  STAR  STAR  STAR  STAR  STAR  STAR 
  ! 
  ! Distance between stars
  dxstar = x_star(1)-x_star(2)
  dystar = y_star(1)-y_star(2)
  dzstar = z_star(1)-z_star(2)
  dr = dxstar*dxstar + dystar*dystar + dzstar*dzstar+eps**2
  dr = sqrt(dr)

!  write (*,*) 'CPU ',me,' star sep',dr
  !   Accelerations:
  do i=1,nstar
     do j=1,nstar              
        if(i.ne.j) then
           ax_star(i) = (-1)**i*gamma*(m_star(j)*dxstar/(dr**3))
           ay_star(i) = (-1)**i*gamma*(m_star(j)*dystar/(dr**3)) 
           az_star(i) = (-1)**i*gamma*(m_star(j)*dzstar/(dr**3))          
           ! potential  
           pot_star(i) = -gamma*m_star(j)/dr
        end if
     end do
  end do

 !         write(6,*)"star1 without",ax_star(1),ay_star(1),az_star(1),pot_star(1)
 !         write(6,*)"star2 without",ax_star(2),ay_star(2),az_star(2),pot_star(2)


  !
  ! DISC  DISC  DISC  DISC  DISC  DISC  DISC  DISC  DISC  DISC  DISC 
  !
  ! Distance between star and particle
  ! initialise partial sums over own disc particles: 1..npp

  do i=1,nstar
     ax_partial = 0.
     ay_partial = 0.
     az_partial = 0.
     pot_partial = 0.

     do j= 1,npp
        dxpart = x_star(i)-x(j)
        dypart = y_star(i)-y(j)
        dzpart = z_star(i)-z(j)
        drpart = dxpart*dxpart + dypart*dypart + dzpart*dzpart+eps**2
        drpart = sqrt(drpart)
        !               write(6,*)"star pos",i,x_star(i),y_star(i),z_star(i)
        !               write(6,*)"part dist",i,dxpart,dypart,dzpart,drpart
        !   Accelerations:
        ax_partial = ax_partial - gamma*m(j)*dxpart/drpart**3
        ay_partial = ay_partial - gamma*m(j)*dypart/drpart**3 
        az_partial = az_partial - gamma*m(j)*dzpart/drpart**3          
        ! potential  
        pot_partial = pot_partial - gamma*m(j)/drpart
     end do

     ! add partial sums from each processor together 
     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE(ax_partial, ax_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(ay_partial, ay_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(az_partial, az_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(pot_partial, pot_sum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD,ierr)

     ! add results to star arrays
     ax_star(i) = ax_star(i) + ax_sum
     ay_star(i) = ay_star(i) + ay_sum
     az_star(i) = az_star(i) + az_sum
     pot_star(i) = pot_star(i) + pot_sum

  end do

  !        write(6,*)"star1",ax_star(1),ay_star(1),az_star(1),pot_star(1)

  !        write(6,*)"star2",ax_star(2),ay_star(2),az_star(2),pot_star(2)
  !        stop

  !  Star velocities

  do i = 1,nstar
     ux_star(i) = ux_star(i) + delta_t * ax_star(i)
     uy_star(i) = uy_star(i) + delta_t * ay_star(i)
     uz_star(i) = uz_star(I) + delta_t * az_star(i)
  end do

  !  star push in space

  do p=1,nstar
     x_star(p) = x_star(p) + ux_star(p)*delta_t
     y_star(p) = y_star(p) + uy_star(p)*delta_t
     z_star(p) = z_star(p) + uz_star(p)*delta_t
  end do
  !  write(16,*)x_star(1),y_star(1),z_star(1),x_star(2),y_star(2),z_star(2)


end subroutine stars
