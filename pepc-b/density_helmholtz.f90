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
  real :: yh_start, zh_start ! Start of HH grid in transverse directions
  integer, parameter :: nyh=4 ! # additional points in transverse direction
  integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, jfoc, kfoc
  integer :: ierr
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis
  real*4 :: rho_loc(0:nxh+1,0:nyh,0:nyh), rho_glob(0:nxh+1,0:nyh,0:nyh)
  real*4 :: charge_sum, charge_tot

! Helmholtz grid limits
  dxh = (xh_end-xh_start)/nxh
! tranverse resolution defined by particle spacing 
  dy = max(dxh,5*a_ii)
  dz = max(dxh,5*a_ii)
!  dy = dxh
!  dz = dxh
  yh_start = focus(2)-nyh/2*dy
  zh_start = focus(3)-nyh/2*dz

  rdx = 1./dxh
  rdy = 1./dy
  rdz = 1./dz


  !  Any particle outside gets put in ghost cells 0, 2

  cweight = qi*rdx*rdy*rdz       ! charge weighting factor
!  if (me==0)   write(6,'(//a,3f12.3)') 'cw,dx,dy',cweight,dxh,dy

  rho_loc(0:nxh+1,0:nyh,0:nyh) = 0.

  do i=1,npp

     xa=(x(i) - xh_start)*rdx
     ya=(y(i) - yh_start)*rdy
     za=(z(i) - zh_start)*rdz

     !  indices
     i1=xa+1
     i2=i1+1
     j1=ya
     j2=j1+1
     k1=za
     k2=k1+1

     i1 = min(max(0,i1),nxh+1)
     i2 = min(max(0,i2),nxh+1)
     j1 = min(max(0,j1),nyh)
     j2 = min(max(0,j2),nyh)
     k1 = min(max(0,k1),nyh)
     k2 = min(max(0,k2),nyh)

     !  linear weighting
     fx2=i1-xa  ! Prevent overflow/negative weighting for particles outside box
     fx1=1.-fx2
     fy2=ya-j1
     fy1=1.-fy2
     fz2=za-k1
     fz1=1.-fz2

     !  gather ion charge at nearest grid points
     if (q(i)>0) then

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

  ng = (nxh+2)*(nyh+1)**2                         ! total # gridpoints
! gather on all
  call MPI_ALLREDUCE(rho_loc, rho_glob, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

! Store density lineout for Helmholtz solver
  do i=0,nxh+1
    rho_helm(i) = rho_glob(i,nyh/2,nyh/2)
  end do 
 charge_sum = SUM(rho_glob(1:nxh,nyh/2,nyh/2))/cweight
 charge_tot = SUM(rho_glob)/cweight
if (me==0 .and. current_step==1) then
  write (6,*) "Charge sum on HH grid lineout/total:",charge_sum,charge_tot
!  write(6,*) "rho20",(rho_glob(i,nyh/2,0),i=1,nxh)
!  write(6,*) "rho21",(rho_glob(i,nyh/2,1),i=1,nxh)
!  write(6,*) "rho22",(rho_glob(i,nyh/2,2),i=1,nxh)
!  write(6,*) "rho23",(rho_glob(i,nyh/2,3),i=1,nxh)
!  write(6,*) "rho24",(rho_glob(i,nyh/2,4),i=1,nxh)
endif
end subroutine density_helmholtz
