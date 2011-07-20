!  =================================
!
!    3D gather for DC fields
!  =================================

subroutine sum_fields

  use physvars
  use treevars

  implicit none

  real :: rdx, rdy, rdz, dx, dy, dz, cweight, jxweight, jyweight, jzweight
  real :: tweight

  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za, gamma, fr1, fr2
  integer :: i, j, k, i1, i2, j1, j2, k1, k2, nelecs, nions
  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: ex_w, ey_w, ez_w
  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: bx_w, by_w, bz_w
  real :: gmin=1.e-3

  dx = xl/ngx
  dy = yl/ngy
  dz = zl/ngz
  rdx = 1./dx
  rdy = 1./dy
  rdz = 1./dz

  do k=1,ngz
     do j=1,ngy
        do i=1,ngx
           ex_w(i,j,k) = 0.
           ey_w(i,j,k) = 0.
           ez_w(i,j,k) = 0.
           bx_w(i,j,k) = 0.
           by_w(i,j,k) = 0.
           bz_w(i,j,k) = 0.
           g_ion(i,j,k) = 0.
           g_ele(i,j,k) = 0.
           rhoi_loc(i,j,k) = 0.
           rhoe_loc(i,j,k) = 0.
           jxe_loc(i,j,k)=0.
           jye_loc(i,j,k)=0.
           jze_loc(i,j,k)=0.
           te_loc(i,j,k) = 0.
           ti_loc(i,j,k) = 0.
        end do
     end do
  end do

  !  field box limits: (0-xl, 0-yl, 0-zl)
  !  Any particle outside gets put in ghost cells 0, ngx+1

  !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy


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
     gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
     cweight = abs(q(i))*rdx*rdy*rdz       ! charge weighting factor

     jxweight = cweight*ux(i)/gamma
     jyweight = cweight*uy(i)/gamma
     jzweight = cweight*uz(i)/gamma
     tweight = (gamma-1.)  ! K.E. of particle in keV
     fr1 = sqrt(fx1**2+fy1**2+fz1**2+eps**2)
     fr2 = sqrt(fx2**2+fy2**2+fz2**2+eps**2)


     if (q(i)<0) then
       	g_ele(i1,j1,k1)=g_ele(i1,j1,k1) + fx1*fy1*fz1  ! weighted # electrons
	g_ele(i2,j1,k1)=g_ele(i2,j1,k1) + fx2*fy1*fz1
	g_ele(i1,j2,k1)=g_ele(i1,j2,k1) + fx1*fy2*fz1
	g_ele(i2,j2,k1)=g_ele(i2,j2,k1) + fx2*fy2*fz1
	g_ele(i1,j1,k2)=g_ele(i1,j1,k2) + fx1*fy1*fz2
	g_ele(i2,j1,k2)=g_ele(i2,j1,k2) + fx2*fy1*fz2
	g_ele(i1,j2,k2)=g_ele(i1,j2,k2) + fx1*fy2*fz2
	g_ele(i2,j2,k2)=g_ele(i2,j2,k2) + fx2*fy2*fz2

	rhoe_loc(i1,j1,k1)=rhoe_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
	rhoe_loc(i2,j1,k1)=rhoe_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
	rhoe_loc(i1,j2,k1)=rhoe_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
	rhoe_loc(i2,j2,k1)=rhoe_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
	rhoe_loc(i1,j1,k2)=rhoe_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
	rhoe_loc(i2,j1,k2)=rhoe_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
	rhoe_loc(i1,j2,k2)=rhoe_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
	rhoe_loc(i2,j2,k2)=rhoe_loc(i2,j2,k2) + cweight*fx2*fy2*fz2

        jxe_loc(i1,j1,k1)=jxe_loc(i1,j1,k1) + jxweight*fx1*fy1*fz1
        jxe_loc(i2,j1,k1)=jxe_loc(i2,j1,k1) + jxweight*fx2*fy1*fz1
        jxe_loc(i1,j2,k1)=jxe_loc(i1,j2,k1) + jxweight*fx1*fy2*fz1
        jxe_loc(i2,j2,k1)=jxe_loc(i2,j2,k1) + jxweight*fx2*fy2*fz1
        jxe_loc(i1,j1,k2)=jxe_loc(i1,j1,k2) + jxweight*fx1*fy1*fz2
        jxe_loc(i2,j1,k2)=jxe_loc(i2,j1,k2) + jxweight*fx2*fy1*fz2
        jxe_loc(i1,j2,k2)=jxe_loc(i1,j2,k2) + jxweight*fx1*fy2*fz2
        jxe_loc(i2,j2,k2)=jxe_loc(i2,j2,k2) + jxweight*fx2*fy2*fz2

        jye_loc(i1,j1,k1)=jye_loc(i1,j1,k1) + jyweight*fx1*fy1*fz1
        jye_loc(i2,j1,k1)=jye_loc(i2,j1,k1) + jyweight*fx2*fy1*fz1
        jye_loc(i1,j2,k1)=jye_loc(i1,j2,k1) + jyweight*fx1*fy2*fz1
        jye_loc(i2,j2,k1)=jye_loc(i2,j2,k1) + jyweight*fx2*fy2*fz1
        jye_loc(i1,j1,k2)=jye_loc(i1,j1,k2) + jyweight*fx1*fy1*fz2
        jye_loc(i2,j1,k2)=jye_loc(i2,j1,k2) + jyweight*fx2*fy1*fz2
        jye_loc(i1,j2,k2)=jye_loc(i1,j2,k2) + jyweight*fx1*fy2*fz2
        jye_loc(i2,j2,k2)=jye_loc(i2,j2,k2) + jyweight*fx2*fy2*fz2

        jze_loc(i1,j1,k1)=jze_loc(i1,j1,k1) + jzweight*fx1*fy1*fz1
        jze_loc(i2,j1,k1)=jze_loc(i2,j1,k1) + jzweight*fx2*fy1*fz1
        jze_loc(i1,j2,k1)=jze_loc(i1,j2,k1) + jzweight*fx1*fy2*fz1
        jze_loc(i2,j2,k1)=jze_loc(i2,j2,k1) + jzweight*fx2*fy2*fz1
        jze_loc(i1,j1,k2)=jze_loc(i1,j1,k2) + jzweight*fx1*fy1*fz2
        jze_loc(i2,j1,k2)=jze_loc(i2,j1,k2) + jzweight*fx2*fy1*fz2
        jze_loc(i1,j2,k2)=jze_loc(i1,j2,k2) + jzweight*fx1*fy2*fz2
        jze_loc(i2,j2,k2)=jze_loc(i2,j2,k2) + jzweight*fx2*fy2*fz2

        Te_loc(i1,j1,k1)=Te_loc(i1,j1,k1) + tweight*fx1*fy1*fz1
        Te_loc(i2,j1,k1)=Te_loc(i2,j1,k1) + tweight*fx2*fy1*fz1
        Te_loc(i1,j2,k1)=Te_loc(i1,j2,k1) + tweight*fx1*fy2*fz1
        Te_loc(i2,j2,k1)=Te_loc(i2,j2,k1) + tweight*fx2*fy2*fz1
        Te_loc(i1,j1,k2)=Te_loc(i1,j1,k2) + tweight*fx1*fy1*fz2
        Te_loc(i2,j1,k2)=Te_loc(i2,j1,k2) + tweight*fx2*fy1*fz2
        Te_loc(i1,j2,k2)=Te_loc(i1,j2,k2) + tweight*fx1*fy2*fz2
        Te_loc(i2,j2,k2)=Te_loc(i2,j2,k2) + tweight*fx2*fy2*fz2

     else
       	g_ion(i1,j1,k1)=g_ion(i1,j1,k1) + fx1*fy1*fz1  ! weighted # ions
	g_ion(i2,j1,k1)=g_ion(i2,j1,k1) + fx2*fy1*fz1
	g_ion(i1,j2,k1)=g_ion(i1,j2,k1) + fx1*fy2*fz1
	g_ion(i2,j2,k1)=g_ion(i2,j2,k1) + fx2*fy2*fz1
	g_ion(i1,j1,k2)=g_ion(i1,j1,k2) + fx1*fy1*fz2
	g_ion(i2,j1,k2)=g_ion(i2,j1,k2) + fx2*fy1*fz2
	g_ion(i1,j2,k2)=g_ion(i1,j2,k2) + fx1*fy2*fz2
	g_ion(i2,j2,k2)=g_ion(i2,j2,k2) + fx2*fy2*fz2

        rhoi_loc(i1,j1,k1)=rhoi_loc(i1,j1,k1) + cweight*fx1*fy1*fz1
        rhoi_loc(i2,j1,k1)=rhoi_loc(i2,j1,k1) + cweight*fx2*fy1*fz1
        rhoi_loc(i1,j2,k1)=rhoi_loc(i1,j2,k1) + cweight*fx1*fy2*fz1
        rhoi_loc(i2,j2,k1)=rhoi_loc(i2,j2,k1) + cweight*fx2*fy2*fz1
        rhoi_loc(i1,j1,k2)=rhoi_loc(i1,j1,k2) + cweight*fx1*fy1*fz2
        rhoi_loc(i2,j1,k2)=rhoi_loc(i2,j1,k2) + cweight*fx2*fy1*fz2
        rhoi_loc(i1,j2,k2)=rhoi_loc(i1,j2,k2) + cweight*fx1*fy2*fz2
        rhoi_loc(i2,j2,k2)=rhoi_loc(i2,j2,k2) + cweight*fx2*fy2*fz2
       ! ion temp
        Ti_loc(i1,j1,k1)=Ti_loc(i1,j1,k1) + tweight*fx1*fy1*fz1
        Ti_loc(i2,j1,k1)=Ti_loc(i2,j1,k1) + tweight*fx2*fy1*fz1
        Ti_loc(i1,j2,k1)=Ti_loc(i1,j2,k1) + tweight*fx1*fy2*fz1
        Ti_loc(i2,j2,k1)=Ti_loc(i2,j2,k1) + tweight*fx2*fy2*fz1
        Ti_loc(i1,j1,k2)=Ti_loc(i1,j1,k2) + tweight*fx1*fy1*fz2
        Ti_loc(i2,j1,k2)=Ti_loc(i2,j1,k2) + tweight*fx2*fy1*fz2
        Ti_loc(i1,j2,k2)=Ti_loc(i1,j2,k2) + tweight*fx1*fy2*fz2
        Ti_loc(i2,j2,k2)=Ti_loc(i2,j2,k2) + tweight*fx2*fy2*fz2

     endif
     ! Electric field - use ngp, softened 1/r^2 weight
     ex_w(i1,j1,k1)=ex_w(i1,j1,k1) + Ex(i)/fr1**2
     ey_w(i1,j1,k1)=ey_w(i1,j1,k1) + Ey(i)/fr1**2
     ez_w(i1,j1,k1)=ez_w(i1,j1,k1) + Ez(i)/fr1**2
     bz_w(i1,j1,k1)=bz_w(i1,j1,k1) + Bz(i)/fr1**2 ! B-field

     !        phig(i1,j1,k1)=phig(i1,j1,k1) + phi(i)/fr1

  end do

  ! normalise averaged quantities
  nelecs = SUM(g_ele(1:ngx,1:ngy,1:ngz))
  nions = SUM(g_ion(1:ngx,1:ngy,1:ngz))
  if (debug_level>1) write(ipefile,*) 'density integrals: ',nelecs, nions
  g_ele = max(gmin,g_ele)
  g_ion = max(gmin,g_ion)
  Te_loc = .511*Te_loc   ! Temperature in MeV (KE per particle)
  Ti_loc = .511*mass_ratio*Ti_loc   ! K.E. in MeV

  Ex_loc = Ex_loc + ex_w/(g_ion+g_ele)/navcycle    ! Accumulate normalised fields
  Ey_loc = Ey_loc + ey_w/(g_ion+g_ele)/navcycle
  Ez_loc = Ez_loc + ez_w/(g_ion+g_ele)/navcycle

  Bz_loc = bz_w/(g_ion+g_ele)  ! Instantaneous B

! laser fields

end subroutine sum_fields





