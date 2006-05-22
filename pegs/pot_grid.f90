!  =================================
!
!    3D gather for potential and E-field
!
!  =================================

subroutine pot_grid

  use treevars
  use physvars

  implicit none
  include 'mpif.h'

  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: phi_loc, Egx_loc, Egy_loc, Egz_loc, &
       g_w, g_w2

  real :: rdx, rdy, rdz, dx, dy, dz, cweight, c_pot
  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za, fr1
  real :: xpmin, ypmin, zpmin, xd,yd,zd
  integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, jfoc, kfoc, ierr, root=0
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis

  dx = xl/ngx
  dy = yl/ngy
  dz = zl/ngz
  rdx = 1./dx
  rdy = 1./dy
  rdz = 1./dz


  !  field box limits: (0-xl, 0-yl, 0-zl)
  !  Any particle outside gets put in ghost cells 0, ngx+1

  !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

  phi_loc(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Egx_loc(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Egy_loc(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Egz_loc(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  g_w(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  g_w2(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  xpmin = 0.
  ypmin = 0.
  zpmin = 0.

  do i=1,npp
     xd=(x(i)-xpmin)
     yd=(y(i)-ypmin)
     zd=(z(i)-zpmin)
     xa=xd*rdx
     ya=yd*rdy
     za=zd*rdz


     !  indices
     i1=xa+1
     j1=ya+1
     k1=za+1
     i1 = min(max(0,i1),ngx+1)
     j1 = min(max(0,j1),ngy+1)
     k1 = min(max(0,k1),ngz+1)

     fr1 = sqrt(xd**2+yd**2+zd**2+eps**2)
     phi_loc(i1,j1,k1)=phi_loc(i1,j1,k1) + pot(i)/fr1 ! NGP 1/r weighting
     Egx_loc(i1,j1,k1)=Egx_loc(i1,j1,k1) + ax(i)*m(i)/q(i)/fr1**2
     Egy_loc(i1,j1,k1)=Egy_loc(i1,j1,k1) + ay(i)*m(i)/q(i)/fr1**2
     Egz_loc(i1,j1,k1)=Egz_loc(i1,j1,k1) + az(i)*m(i)/q(i)/fr1**2
     g_w(i1,j1,k1) =  g_w(i1,j1,k1) + 1./fr1  ! Sum S(r) weights
     g_w2(i1,j1,k1) =  g_w(i1,j1,k1) + 1./fr1**2  ! Sum S(r) weights

     !  gather charge at nearest grid points


  end do
  g_w(1:ngx,1:ngy,1:ngz) = max(1.,g_w(1:ngx,1:ngy,1:ngz))
  g_w2(1:ngx,1:ngy,1:ngz) = max(1.,g_w2(1:ngx,1:ngy,1:ngz))
  C_pot = 1.
  Egx_loc(1:ngx,1:ngy,1:ngz) = Egx_loc(1:ngx,1:ngy,1:ngz)/ &
       (g_w2(1:ngx,1:ngy,1:ngz))  ! Average potential in cell and convert to MV
  Egy_loc(1:ngx,1:ngy,1:ngz) = Egy_loc(1:ngx,1:ngy,1:ngz)/ &
       (g_w2(1:ngx,1:ngy,1:ngz))
  Egz_loc(1:ngx,1:ngy,1:ngz) = Egz_loc(1:ngx,1:ngy,1:ngz)/ &
       (g_w2(1:ngx,1:ngy,1:ngz))
  phi_loc(1:ngx,1:ngy,1:ngz) = C_pot*phi_loc(1:ngx,1:ngy,1:ngz)/ &
       (g_w(1:ngx,1:ngy,1:ngz))

  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints
  ! Gather sum on root
  call MPI_REDUCE(phi_loc, phi_g, ng, MPI_REAL, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(Egx_loc, Ex_g, ng, MPI_REAL, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(Egy_loc, Ey_g, ng, MPI_REAL, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(Egz_loc, Ez_g, ng, MPI_REAL, MPI_SUM, root, MPI_COMM_WORLD, ierr)


end subroutine pot_grid
