!  =================================
!
!    3D Density gather for rhoi used in Helmholtz solver
!  - 1D grid with 3x3 y-z cells for dummy counts
!
!  =================================

subroutine density_helmholtz

  use physvars
  use treevars

  implicit none
  include 'mpif.h'

  real :: rdx, rdy, rdz, dx, dy, dz, cweight
  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za
  integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, jfoc, kfoc
  integer :: ierr
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis
  real*4 :: rho_loc(0:nxh+1,0:2,0:2), rho_glob(0:nxh+1,0:2,0:2)

! Helmholtz grid limits
  dxh = (xh_end-xh_start)/nxh

! tranverse resolution defined by laser spot size
  dy = sigma/5.
  dz = sigma/5.
  rdx = 1./dxh
  rdy = 1./dy
  rdz = 1./dz


  !  Any particle outside gets put in ghost cells 0, 2

  !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

  rho_loc(0:nxh+1,0:2,0:2) = 0.

  do i=1,npp

     xa=(x(i) - xh_start)*rdx
     ya=(y(i) - focus(2))*rdy
     za=(z(i) - focus(3))*rdz

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

        rho_loc(i1,j1,k1)=rho_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
        rho_loc(i2,j1,k1)=rho_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
        rho_loc(i1,j2,k1)=rho_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
        rho_loc(i2,j2,k1)=rho_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
        rho_loc(i1,j1,k2)=rho_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
        rho_loc(i2,j1,k2)=rho_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
        rho_loc(i1,j2,k2)=rho_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
        rho_loc(i2,j2,k2)=rho_loc(i2,j2,k2) + cweight*fx2*fy2*fz2
     else

     endif
  end do

  ng = (nxh+2)*9                         ! total # gridpoints
! gather on root
  call MPI_REDUCE(rho_loc, rho_glob, ng, MPI_REAL, MPI_SUM, 0,  MPI_COMM_WORLD, ierr)

! Store density lineout for Helmholtz solver
  rho_helm(0:nxh+1) = rho_glob(0:nxh+1,1,1) 

end subroutine density_helmholtz
