!  =================================
!
!    3D Density gather for rhoi, rhoe
!
!  =================================

subroutine densities

  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none

  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: rhoi_loc, rhoe_loc
  real, dimension(0:ngx+1) :: rho1d
  real :: rdx, rdy, rdz, dx, dy, dz, cweight
  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za
  integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, jfoc, kfoc
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis

!VAMPINST subroutine_start
       CALL VTENTER(IF_densities,VTNOSCL,VTIERR)
!      write(*,*) 'VT: densities S>',VTIERR,
!     *    IF_densities,ICLASSH
!
  dx = xl/ngx
  dy = yl/ngy
  dz = zl/ngz
  rdx = 1./dx
  rdy = 1./dy
  rdz = 1./dz


!  field box limits: (0-xl, 0-yl, 0-zl)
!  Any particle outside gets put in ghost cells 0, ngx+1

  !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

  rhoe_loc(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
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
     fx1=i1-xa
     fx2=1.-fx1
     fy1=j1-ya
     fy2=1.-fy1
     fz1=k1-za
     fz2=1.-fz1

     !  gather charge at nearest grid points

     cweight = q(i)*rdx*rdy*rdz       ! charge weighting factor
     if (q(i)<0) then
	rhoe_loc(i1,j1,k1)=rhoe_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
	rhoe_loc(i2,j1,k1)=rhoe_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
	rhoe_loc(i1,j2,k1)=rhoe_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
	rhoe_loc(i2,j2,k1)=rhoe_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
	rhoe_loc(i1,j1,k2)=rhoe_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
	rhoe_loc(i2,j1,k2)=rhoe_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
	rhoe_loc(i1,j2,k2)=rhoe_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
	rhoe_loc(i2,j2,k2)=rhoe_loc(i2,j2,k2) + cweight*fx2*fy2*fz2
     else
	rhoi_loc(i1,j1,k1)=rhoi_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
	rhoi_loc(i2,j1,k1)=rhoi_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
	rhoi_loc(i1,j2,k1)=rhoi_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
	rhoi_loc(i2,j2,k1)=rhoi_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
	rhoi_loc(i1,j1,k2)=rhoi_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
	rhoi_loc(i2,j1,k2)=rhoi_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
	rhoi_loc(i1,j2,k2)=rhoi_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
	rhoi_loc(i2,j2,k2)=rhoi_loc(i2,j2,k2) + cweight*fx2*fy2*fz2
     endif
  end do

  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints
  call MPI_ALLREDUCE(rhoi_loc, rhoi, ng, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_ALLREDUCE(rhoe_loc, rhoe, ng, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: densities S<',VTIERR,ICLASSH
!
end subroutine densities
