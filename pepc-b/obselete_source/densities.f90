!  =================================
!
!    3D Density gather for rhoi
!
!  =================================

subroutine densities

  use physvars
  use treevars

  implicit none
  include 'mpif.h'

  real :: rdx, rdy, rdz, dx, dy, dz, cweight
  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za
  integer :: i, ng, i1, i2, j1, j2, k1, k2
  integer :: ierr

  dx = xl/ngx
  dy = yl/ngy
  dz = zl/ngz
  rdx = 1./dx
  rdy = 1./dy
  rdz = 1./dz


  !  field box limits: (0-xl, 0-yl, 0-zl)
  !  Any particle outside gets put in ghost cells 0, ngx+1

  !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

  rhoi_loc(0:ngx+1,0:ngy+1,0:ngz+1) = 0.

  do i=1,npp

     xa=x(i)*rdx
     ya=y(i)*rdy
     za=z(i)*rdz

     !  indices
     i1=xa+1
     i2=i1+1
     j1=ya+1
     j2=j1+1
     k1=za+1
     k2=k1+1

     i1 = min(max(0,i1),ngx+1)
     i2 = min(max(0,i2),ngx+1)
     j1 = min(max(0,j1),ngy+1)
     j2 = min(max(0,j2),ngy+1)
     k1 = min(max(0,k1),ngz+1)
     k2 = min(max(0,k2),ngz+1)

     !  linear weighting
     fx2=min(max(i1-xa,0.),1.)  ! Prevent overflow/negative weighting for particles outside box
     fx1=1.-fx2
     fy2=min(max(j1-ya,0.),1.)
     fy1=1.-fy2
     fz2=min(max(k1-za,0.),1.)
     fz1=1.-fz2

     !  gather charge at nearest grid points
     if (q(i)>0) then
        cweight = q(i)*rdx*rdy*rdz       ! charge weighting factor

        rhoi_loc(i1,j1,k1)=rhoi_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
        rhoi_loc(i2,j1,k1)=rhoi_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
        rhoi_loc(i1,j2,k1)=rhoi_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
        rhoi_loc(i2,j2,k1)=rhoi_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
        rhoi_loc(i1,j1,k2)=rhoi_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
        rhoi_loc(i2,j1,k2)=rhoi_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
        rhoi_loc(i1,j2,k2)=rhoi_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
        rhoi_loc(i2,j2,k2)=rhoi_loc(i2,j2,k2) + cweight*fx2*fy2*fz2
     else

     endif
  end do

  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints
! gather on root
  call MPI_REDUCE(rhoi_loc, rhoi, ng, MPI_REAL, MPI_SUM, 0,  MPI_COMM_WORLD, ierr)
!  call MPI_REDUCE(rhoe_loc, rhoe, ng, MPI_REAL, MPI_SUM, 0,  MPI_COMM_WORLD, ierr)

end subroutine densities