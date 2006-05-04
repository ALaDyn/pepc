!  =================================
!
!    Postprocessor for particle data
!
!  =================================

program ppfields


  implicit none
  real, parameter :: pi=3.141592654
  integer :: n ! # particles
  integer :: ngx, ngy, ngz  ! grid dimensions
  integer :: ngux, nguy, nguz  ! momentum grid dimensions
  integer :: npx, npy, npz  ! reduced grid dimensions
  integer :: nalpha
  real, allocatable :: rho_ion(:,:,:), rho_ele(:,:,:), jx_ele(:,:,:), jx_ion(:,:,:),&
       jy_ele(:,:,:), jy_ion(:,:,:), jz_ele(:,:,:), jz_ion(:,:,:)
  real, allocatable :: Egx(:,:,:), Egy(:,:,:), Egz(:,:,:), Te(:,:,:), Ti(:,:,:), &
       phig(:,:,:), phi_whole(:,:,:)
  real, allocatable :: g_ion(:,:,:), g_ele(:,:,:), g_w(:,:,:)  ! particle weights for averages

  real, allocatable :: fpxpy_ion(:,:), fpxpy_ele(:,:)
  real, allocatable :: fpxpz_ion(:,:), fpxpz_ele(:,:)
  real, allocatable :: fpxx_ion(:,:), fpxx_ele(:,:), fpzx_ion(:,:)
  real, allocatable :: fang_ion(:,:), fang_ion1(:,:), fang_ion2(:,:), fang_back(:,:)
  real, allocatable :: uang_ion(:,:)

  real, allocatable :: x(:),y(:),z(:),ux(:),uy(:),uz(:),q(:),m(:),ex(:),ey(:),ez(:),phi(:)
  integer, allocatable :: own(:), label(:)
  real, allocatable :: aveni(:),avene(:),avejxi(:),avejxe(:),aveex(:),aveti(:),avete(:),avephi(:)
  real, allocatable :: fhote(:), fhoti_for(:), fhoti_back(:), avenhot(:)

  real, dimension(10) :: xslice, yslice, zslice
  real :: shift(3)

  real :: xmin, xmax, ymin, ymax, zmin, zmax
  real :: rdx, rdy, rdz, dx, dy, dz, cweight, jxweight, jyweight, jzweight
  real :: omega, tweight, exweight, eyweight
  real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za, xl, yl, zl, fr1, fr2, xxa
  real :: uxmin, uxmax, uymin, uymax, uzmin, uzmax, uxl, uyl, uzl, dux, duy, duz, du
  real :: uximin, uximax, uyimin, uyimax, uzimin, uzimax, dalpha
  real :: temin, temax, timin, timax, jemax, jimax, emax, fpmax, rhomax, tcold
  real :: xbox, ybox, zbox, yshift, zshift, pshift, jevec, jivec, evec, epsilon
  real :: xtick, ytick, ztick, uxtick, uytick, uztick, uxitick, uyitick, uzitick
  real :: umevmax, uimevmax, umev, uimev1, uimev2, xrear
  real :: xparmin, xparmax, yparmin, yparmax, zparmin,zparmax, vxmin, vxmax, massmax

  integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, islice, jslice, kslice, nave, norm
  integer :: iu, ii1, kmax, nk
  character(30) :: cfile,cfile1, cinfo, csnap
  character(30), dimension(15) :: cfout

  character(5) :: cme
  character(6) :: cdump, cvis
  logical :: found
  real :: xc1, rho_track, gamma, mass_ratio, phimin,phimax, logrho, mdisc

  integer :: icm, icp, nout, nelecs, nions, nsel, nmerger

  integer ::   itime_start, npartr,  ndustr, nstarr, ischeme, ntr
  real :: xlr, ylr, zlr, boxr, epsr, thetar, dtr, trun, mdiscr, qs
  real :: aimax, aimin, alphaz, alphay
  real :: xpmin, xpmax, ypmin, ypmax, zpmin, zpmax, xptick, yptick, xd, yd, zd
  real :: sumjx, sumjy, sumjz, sumex, sumey, sumez
  integer :: iaz, iay, iskip3d

  ! Conversion factors
  real :: conv_length  ! length conversion
  real :: conv_time  ! time conversion
  real :: C_dens ! mass density scaling factor for plot
  real :: C_pot ! potential/mc2 -> MV

  ! Diagnostics
  integer :: ntnsa, nfront, nblow
  real :: utnsa, ufront, ublow, rear_edge, front_edge, ub_max, uf_max, ur_max

  ! Determine grid and particle dimensions for postprocessing 
  cfile = 'grid_defs.in'
  write(*,*) 'Reading grid definitions from ',cfile
  open(20,file=cfile)
  read(20,'(8(9x,i10/))') n, ngx, ngy, ngz, ngux, nguy, nguz, nalpha
  read(20,'(9(9x,f12.5/))') xmin, xmax, xtick, ymin, ymax, ytick, zmin, zmax, ztick
  read(20,'(10(9x,f12.5/))') uximin, uximax, uxitick, uyimin, uyimax, uyitick, uzimin, uzimax, uzitick, umevmax
  read(20,'(14(9x,f12.5/))') uxmin, uxmax, uxtick, uymin, uymax, uytick, uzmin, uzmax, uztick, uimevmax,aimin,aimax, uimev1, uimev2
  read(20,'(18(9x,f12.5/))') yslice(1), zslice(1), rear_edge, mass_ratio, rhomax, temin, temax, timin, timax, &
       jemax, jimax, emax, fpmax, jevec, jivec, evec, epsr, tcold
  read(20,'(6(9x,f12.5/))') xbox, ybox, zbox, yshift, zshift, pshift    ! plot dimensions and scale positions in inches
  read(20,'(9x,i6)') iskip3d  ! skip stride for writing out 3D xyz plots
  close(20)
  write(*,'(a/,7(a20,i8/))') ' Grid check:' &
	,'# particles: ',n &
	,'x points: ',ngx &
	,'y points: ',ngy &
	,'z points: ',ngz &
	,'vx points: ',ngux &
	,'vy points: ',nguy &
	,'vz points: ',nguz

  ! open GMT defs file 
  cfile = 'grid_defs.gmt'
  open(50,file=cfile)

  allocate ( x(n), y(n), z(n), ux(n), uy(n), uz(n), q(n), m(n), ex(n), ey(n), ez(n), phi(n), own(n), label(n) )
  allocate ( rho_ion(0:ngx+1,0:ngy+1,0:ngz+1), rho_ele(0:ngx+1,0:ngy+1,0:ngz+1), &
       jx_ele(0:ngx+1,0:ngy+1,0:ngz+1), jx_ion(0:ngx+1,0:ngy+1,0:ngz+1), &
       jz_ele(0:ngx+1,0:ngy+1,0:ngz+1), jz_ion(0:ngx+1,0:ngy+1,0:ngz+1), &
       jy_ele(0:ngx+1,0:ngy+1,0:ngz+1), jy_ion(0:ngx+1,0:ngy+1,0:ngz+1) )
  allocate ( Egx(0:ngx+1,0:ngy+1,0:ngz+1), Egy(0:ngx+1,0:ngy+1,0:ngz+1), &
       Egz(0:ngx+1,0:ngy+1,0:ngz+1), &
       Te(0:ngx+1,0:ngy+1,0:ngz+1), Ti(0:ngx+1,0:ngy+1,0:ngz+1),&
       phig(0:ngx+1,0:ngy+1,0:ngz+1), phi_whole(0:ngx+1,0:ngy+1,0:ngz+1)) 
  allocate ( g_ion(0:ngx+1,0:ngy+1,0:ngz+1), g_ele(0:ngx+1,0:ngy+1,0:ngz+1), &
       g_w(0:ngx+1,0:ngy+1,0:ngz+1) )
  allocate ( fpxpy_ion(0:ngux+1,0:nguy+1), fpxpy_ele(0:ngux+1,0:nguy+1), &
       fpxpz_ion(0:ngux+1,0:nguz+1), fpxpz_ele(0:ngux+1,0:nguz+1), &
       fpxx_ion(0:ngx+1,0:ngux+1), fpxx_ele(0:ngx+1,0:ngux+1), &
       fpzx_ion(0:ngx+1,0:nguz+1), &
       fang_ion(0:nalpha+1,0:nalpha+1),   fang_ion2(0:nalpha+1,0:nalpha+1),&
       fang_back(0:nalpha+1,0:nalpha+1), fang_ion1(0:nalpha+1,0:nalpha+1), &
       uang_ion(0:nalpha+1,0:nalpha+1), &
       fhote(0:ngux+1), fhoti_for(0:ngux+1), fhoti_back(0:ngux+1) )
  allocate ( aveni(ngx),avene(ngx),avejxi(ngx),avejxe(ngx),aveex(ngx),aveti(ngx), & 
       avete(ngx),avephi(ngx), avenhot(ngx) )

  ! input data file

  read(5,'(a)') cdump

  !  read info block in run directory to get copy of run parameters

  cinfo="dumps/parts_info."//cdump
  !  cinfo="dumps/parts_info.000000"
  open(80,file=cinfo)

  cinfo=cdump//"/pp.log"
  open(10,file=cinfo)

!   read(80,'(i12,x20)')  &    ! info block - skip variable names
!       itime_start  ! time stamp

  read(80,'(5(i12/),8(f12.3/),2(i12/),f12.3)')  &    ! info block - skip variable names
       itime_start, &  ! time stamp
       npartr, &  ! total # particles
       ndustr, &  ! # dust particles
       nstarr, &  ! # stars
       ischeme,  & ! # integrator scheme
       mdiscr, &  ! disc mass
       xlr, &  ! box dimensions
       ylr, &
       zlr, &
       boxr, &
       shift(1), &
       epsilon, &  ! Coulomb smoothing radius
       thetar, &   ! Force clumping parameter
       nmerger, &            ! merge switch
       ntr, &     ! # timesteps
       dtr    ! timestep

  close(80)
  trun = ntr*dtr

  ! Compute time, length, number conversion factors

!  conv_length = lambda/2./pi*omega
!  conv_time = 5./3./pi*lambda*omega
  conv_length = 1.
  conv_time = 1.

  C_dens = mass_ratio 
  C_pot = 0.511

  write(10,*) 'Conversion factors: '

  write(10,*) 'Length           ',conv_length
  write(10,*) 'Time             ',conv_time
  write(10,*) 'Mass density     ',C_dens
  write(10,*) 'Potential  ',C_pot
  write(10,*)

  cfile = 'dumps/parts_dump.'//cdump
  write(*,*) 'Reading data from ',cfile
  open(20,file=cfile)
  read(20,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),m(i),q(i),ex(i),ey(i),ez(i),phi(i),own(i),label(i),i=1,n)
  close(20)

! statistics
  xparmin=minval(x)
  xparmax=maxval(x)
  yparmin=minval(y)
  yparmax=maxval(y)
  zparmin=minval(z)
  zparmax = maxval(z)
  vxmin = minval(ux)
  vxmax = maxval(ux)
  massmax = maxval(m)
  write (*,'(/a/,9(a20,1pe12.3,1pe14.2/))') ' Data check:             particles      box ' &
	,'xmin :',xparmin, xmin &
	,'xmax :',xparmax, xmax &
	,'ymin :',yparmin, ymin &
	,'ymax :',yparmax, ymax &
	,'zmin :',zparmin, zmin &
	,'zmax :',zparmax, zmax &
	,'vxmin :',vxmin, uxmin &
	,'vxmax :',vxmax, uxmax &
	,'masses:',minval(m),massmax
  write(*,*) 'z-position of xy slice: ',zslice(1)
  write(*,*) 'y-position of xz slice: ',yslice(1)
  write(*,*) 'Density scaling factor: ', mass_ratio
  write(*,*) 'Disc mass :', mdiscr

  ! Dump for AVS-EXPRESS
  cfile = 'ions_dump.'//cdump
  cfile1 = 'dumps/elecs_dump.'//cdump

  ! write(*,*) 'Writing ions to ',cfile
!  write(*,*) 'Writing hot electrons to ',cfile1
  ! open(20,file=cfile)
  open(21,file=cfile1)
  nsel=0
  do i=1,n
     gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
     umev  = 511*(gamma-1.0)
     if (q(i)<0 .and. umev>tcold) then
        nsel=nsel+1
        write(21,'(3f12.3,3(1pe12.3),1pe12.3)') x(i),y(i),z(i),ux(i),uy(i),uz(i),umev  ! hot electrons
     else
        !      write(20,*) x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),ex(i),ey(i),ez(i),phi(i),own(i),label(i)
     endif
  end do
  ! close(20)
  close(21)
  ! write AVS field record
  open(60,file='elecs.fld',position='append')
  write(60,'(3a)') 'time file='//cfile1//' filetype=ascii'
  write(60,'(a,i8)') 'dim1=',nsel
  write(60,'(3a)') 'coord 1 file='//cfile1//' filetype=ascii offset=0 stride=7'
  write(60,'(3a)') 'coord 2 file='//cfile1//' filetype=ascii offset=1 stride=7'
  write(60,'(3a)') 'coord 3 file='//cfile1//' filetype=ascii offset=2 stride=7'
  write(60,'(3a)') 'variable 1 file='//cfile1//' filetype=ascii offset=3 stride=7'
  write(60,'(3a)') 'variable 2 file='//cfile1//' filetype=ascii offset=4 stride=7'
  write(60,'(3a)') 'variable 3 file='//cfile1//' filetype=ascii offset=5 stride=7'
  write(60,'(3a)') 'variable 4 file='//cfile1//' filetype=ascii offset=6 stride=7'
  write(60,'(a3)') 'EOT'
  close(60)
  ! stop
  ! COORDINATE SPACE
  ! box limits 
  xl = xmax-xmin
  yl = ymax-ymin
  zl = zmax-zmin

  dx = xl/ngx
  dy = yl/ngy
  dz = zl/ngz
  rdx = 1./dx
  rdy = 1./dy
  rdz = 1./dz

  ! grid defs for xy slices

  write(50,'(a8,a3,4(f12.3,a1))') 'PA','"',xmin,'/',xmax,'/',ymin,'/',ymax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.3,a1))') 'DR','"',xmin+dx,'/',xmax,'/',ymin+dy,'/',ymax,'"'  ! data region 
  write(50,'(a8,a3,2(f12.3,a1))') 'MESH','"',dx,'/',dy,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.3,a2))') 'BOX','"X',xbox,'i/',ybox,'i"'   ! plot size in inches 
  write(50,'(a8,a3,2(f12.3,a1))') 'AXES','"',xtick,'/',ytick,'"'   ! tick intervals

  ! grid defs for xz slices

  write(50,'(a8,a3,4(f12.3,a1))') 'ZPA','"',xmin,'/',xmax,'/',zmin,'/',zmax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.3,a1))') 'ZDR','"',xmin+dx,'/',xmax,'/',zmin+dz,'/',zmax,'"'  ! data region 
  write(50,'(a8,a3,2(f12.3,a1))') 'ZMESH','"',dx,'/',dz,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.3,a2))') 'ZBOX','"X',xbox,'i/',zbox,'i"'   ! plot size in inches 
  write(50,'(a8,a3,2(f12.3,a1))') 'ZAXES','"',xtick,'/',ztick,'"'   ! tick intervals

  write(50,'(a8,a3,3(f12.3,a1))') 'RHOMAP','"',rhomax-3.,'/',rhomax,'/',0.1,'"'  ! density map
  write(50,'(a8,a3,3(f12.3,a1))') 'JEMAP','"',0.0,'/',jemax,'/',jemax/20.,'"'  ! current map
  write(50,'(a8,a3,3(f12.3,a1))') 'JIMAP','"',0.0,'/',jimax,'/',jimax/20.,'"'  ! Ji map
  write(50,'(a8,a3,3(f12.3,a1))') 'EMAP','"',0.0,'/',emax,'/',emax/20.,'"'  ! field map
  write(50,'(a8,a3,3(f12.5,a1))') 'TEMAP','"',temin,'/',temax,'/',temax/50.,'"'  ! Te map
  write(50,'(a8,a3,2(f12.3,a1))') 'TERANGE','"',temin,'/',temax,'"'  ! Te range
  write(50,'(a8,a3,3(f12.5,a1))') 'TIMAP','"',timin,'/',timax,'/',timax/50.,'"'  ! Ti map
  write(50,'(a8,a3,2(f12.3,a1))') 'TIRANGE','"',timin,'/',timax,'"'  ! Ti range
  write(50,'(a8,a3,f12.3,a1)') 'SSRHO','"',rhomax/2,'"'  ! density scale ticks 
  write(50,'(a8,a3,f12.3,a1)') 'TCOLD','"',tcold,'"'  ! threshold for cold/hot plots
  write(50,'(a8,a3,f12.3,a1)') 'SSTE','"',temax/2,'"'  ! temp scale ticks 
  write(50,'(a8,a3,f12.3,a1)') 'SSTI','"',timax/2,'"'  ! Ti scale ticks 
  write(50,'(a8,a3,f12.3,a2)') 'YSHIFT','"',yshift,' "'  ! yshift for color scale 
  write(50,'(a8,a3,f12.3,a2)') 'ZSHIFT','"',zshift,' "'  ! yshift for color scale 
  write(50,'(a8,a3,f12.3,a2)') 'PSHIFT','"',pshift,' "'  ! yshift for color scale 
  write(50,'(a8,a3,f12.3,a2)') 'JEVEC','"',jevec,' "'  ! Je vector scale (units/inch) 
  write(50,'(a8,a3,f12.3,a2)') 'JIVEC','"',jivec,' "'  ! Ji scale 
  write(50,'(a8,a3,f12.3,a2)') 'EVEC','"',evec,' "'  ! E-field vector scale 


  !  Any particle outside gets put in ghost cells 0, ngx+1

  !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

  Te(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Ti(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  rho_ele(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  rho_ion(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  jx_ele(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  jx_ion(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  jy_ele(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  jy_ion(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  jz_ele(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  jz_ion(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  g_ele(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  g_ion(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Egx(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Egy(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  Egz(0:ngx+1,0:ngy+1,0:ngz+1) = 0.
  phig(0:ngx+1,0:ngy+1,0:ngz+1) = 0.

  do i=1,n

     xa=(x(i)-xmin)*rdx
     ya=(y(i)-ymin)*rdy
     za=(z(i)-zmin)*rdz

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

     cweight = C_dens*m(i)*rdx*rdy*rdz       ! mass weighting factor
!     write(*,*) cweight
     tweight = .5*m(i)*(ux(i)**2+uy(i)**2+uz(i)**2)  ! K.E. of particle in keV
     fr1 = sqrt(fx1**2+fy1**2+fz1**2+epsilon**2)
     fr2 = sqrt(fx2**2+fy2**2+fz2**2+epsilon**2)
     jxweight = cweight*ux(i)
     jyweight = cweight*uy(i)
     jzweight = cweight*uz(i)


     if (q(i)<0) then
        ! # electrons per cell
        g_ele(i1,j1,k1)=g_ele(i1,j1,k1) + fx1*fy1*fz1  
        g_ele(i2,j1,k1)=g_ele(i2,j1,k1) + fx2*fy1*fz1
        g_ele(i1,j2,k1)=g_ele(i1,j2,k1) + fx1*fy2*fz1
        g_ele(i2,j2,k1)=g_ele(i2,j2,k1) + fx2*fy2*fz1
        g_ele(i1,j1,k2)=g_ele(i1,j1,k2) + fx1*fy1*fz2
        g_ele(i2,j1,k2)=g_ele(i2,j1,k2) + fx2*fy1*fz2
        g_ele(i1,j2,k2)=g_ele(i1,j2,k2) + fx1*fy2*fz2
        g_ele(i2,j2,k2)=g_ele(i2,j2,k2) + fx2*fy2*fz2

 !  
           ! 'cold' electrons
           rho_ele(i1,j1,k1)=rho_ele(i1,j1,k1) + cweight*fx1*fy1*fz1
           rho_ele(i2,j1,k1)=rho_ele(i2,j1,k1) + cweight*fx2*fy1*fz1
           rho_ele(i1,j2,k1)=rho_ele(i1,j2,k1) + cweight*fx1*fy2*fz1
           rho_ele(i2,j2,k1)=rho_ele(i2,j2,k1) + cweight*fx2*fy2*fz1
           rho_ele(i1,j1,k2)=rho_ele(i1,j1,k2) + cweight*fx1*fy1*fz2
           rho_ele(i2,j1,k2)=rho_ele(i2,j1,k2) + cweight*fx2*fy1*fz2
           rho_ele(i1,j2,k2)=rho_ele(i1,j2,k2) + cweight*fx1*fy2*fz2
           rho_ele(i2,j2,k2)=rho_ele(i2,j2,k2) + cweight*fx2*fy2*fz2

           jx_ele(i1,j1,k1)=jx_ele(i1,j1,k1) + jxweight*fx1*fy1*fz1
           jx_ele(i2,j1,k1)=jx_ele(i2,j1,k1) + jxweight*fx2*fy1*fz1
           jx_ele(i1,j2,k1)=jx_ele(i1,j2,k1) + jxweight*fx1*fy2*fz1
           jx_ele(i2,j2,k1)=jx_ele(i2,j2,k1) + jxweight*fx2*fy2*fz1
           jx_ele(i1,j1,k2)=jx_ele(i1,j1,k2) + jxweight*fx1*fy1*fz2
           jx_ele(i2,j1,k2)=jx_ele(i2,j1,k2) + jxweight*fx2*fy1*fz2
           jx_ele(i1,j2,k2)=jx_ele(i1,j2,k2) + jxweight*fx1*fy2*fz2
           jx_ele(i2,j2,k2)=jx_ele(i2,j2,k2) + jxweight*fx2*fy2*fz2

           jy_ele(i1,j1,k1)=jy_ele(i1,j1,k1) + jyweight*fx1*fy1*fz1
           jy_ele(i2,j1,k1)=jy_ele(i2,j1,k1) + jyweight*fx2*fy1*fz1
           jy_ele(i1,j2,k1)=jy_ele(i1,j2,k1) + jyweight*fx1*fy2*fz1
           jy_ele(i2,j2,k1)=jy_ele(i2,j2,k1) + jyweight*fx2*fy2*fz1
           jy_ele(i1,j1,k2)=jy_ele(i1,j1,k2) + jyweight*fx1*fy1*fz2
           jy_ele(i2,j1,k2)=jy_ele(i2,j1,k2) + jyweight*fx2*fy1*fz2
           jy_ele(i1,j2,k2)=jy_ele(i1,j2,k2) + jyweight*fx1*fy2*fz2
           jy_ele(i2,j2,k2)=jy_ele(i2,j2,k2) + jyweight*fx2*fy2*fz2

           jz_ele(i1,j1,k1)=jz_ele(i1,j1,k1) + jzweight*fx1*fy1*fz1
           jz_ele(i2,j1,k1)=jz_ele(i2,j1,k1) + jzweight*fx2*fy1*fz1
           jz_ele(i1,j2,k1)=jz_ele(i1,j2,k1) + jzweight*fx1*fy2*fz1
           jz_ele(i2,j2,k1)=jz_ele(i2,j2,k1) + jzweight*fx2*fy2*fz1
           jz_ele(i1,j1,k2)=jz_ele(i1,j1,k2) + jzweight*fx1*fy1*fz2
           jz_ele(i2,j1,k2)=jz_ele(i2,j1,k2) + jzweight*fx2*fy1*fz2
           jz_ele(i1,j2,k2)=jz_ele(i1,j2,k2) + jzweight*fx1*fy2*fz2
           jz_ele(i2,j2,k2)=jz_ele(i2,j2,k2) + jzweight*fx2*fy2*fz2


     else
        ! # ions per cell
        g_ion(i1,j1,k1)=g_ion(i1,j1,k1) + fx1*fy1*fz1  
        g_ion(i2,j1,k1)=g_ion(i2,j1,k1) + fx2*fy1*fz1
        g_ion(i1,j2,k1)=g_ion(i1,j2,k1) + fx1*fy2*fz1
        g_ion(i2,j2,k1)=g_ion(i2,j2,k1) + fx2*fy2*fz1
        g_ion(i1,j1,k2)=g_ion(i1,j1,k2) + fx1*fy1*fz2
        g_ion(i2,j1,k2)=g_ion(i2,j1,k2) + fx2*fy1*fz2
        g_ion(i1,j2,k2)=g_ion(i1,j2,k2) + fx1*fy2*fz2
        g_ion(i2,j2,k2)=g_ion(i2,j2,k2) + fx2*fy2*fz2

        rho_ion(i1,j1,k1)=rho_ion(i1,j1,k1) + cweight*fx1*fy1*fz1
        rho_ion(i2,j1,k1)=rho_ion(i2,j1,k1) + cweight*fx2*fy1*fz1
        rho_ion(i1,j2,k1)=rho_ion(i1,j2,k1) + cweight*fx1*fy2*fz1
        rho_ion(i2,j2,k1)=rho_ion(i2,j2,k1) + cweight*fx2*fy2*fz1
        rho_ion(i1,j1,k2)=rho_ion(i1,j1,k2) + cweight*fx1*fy1*fz2
        rho_ion(i2,j1,k2)=rho_ion(i2,j1,k2) + cweight*fx2*fy1*fz2
        rho_ion(i1,j2,k2)=rho_ion(i1,j2,k2) + cweight*fx1*fy2*fz2
        rho_ion(i2,j2,k2)=rho_ion(i2,j2,k2) + cweight*fx2*fy2*fz2


        jx_ion(i1,j1,k1)=jx_ion(i1,j1,k1) + ux(i)/fr1**2

        jy_ion(i1,j1,k1)=jy_ion(i1,j1,k1) + uy(i)/fr1**2

        jz_ion(i1,j1,k1)=jz_ion(i1,j1,k1) + uz(i)/fr1**2

        ! ion temp
        Ti(i1,j1,k1)=Ti(i1,j1,k1) + tweight*fx1*fy1*fz1
        Ti(i2,j1,k1)=Ti(i2,j1,k1) + tweight*fx2*fy1*fz1
        Ti(i1,j2,k1)=Ti(i1,j2,k1) + tweight*fx1*fy2*fz1
        Ti(i2,j2,k1)=Ti(i2,j2,k1) + tweight*fx2*fy2*fz1
        Ti(i1,j1,k2)=Ti(i1,j1,k2) + tweight*fx1*fy1*fz2
        Ti(i2,j1,k2)=Ti(i2,j1,k2) + tweight*fx2*fy1*fz2
        Ti(i1,j2,k2)=Ti(i1,j2,k2) + tweight*fx1*fy2*fz2
        Ti(i2,j2,k2)=Ti(i2,j2,k2) + tweight*fx2*fy2*fz2

        ! electric fields
        Egx(i1,j1,k1)=Egx(i1,j1,k1) + ex(i)/fr1**2

        Egy(i1,j1,k1)=Egy(i1,j1,k1) + ey(i)/fr1**2
        Egz(i1,j1,k1)=Egz(i1,j1,k1) + ez(i)/fr1**2

        phig(i1,j1,k1)=phig(i1,j1,k1) + phi(i)/fr1

     endif
  end do

  ! normalise averaged quantities
  nelecs = SUM(g_ele(1:ngx,1:ngy,1:ngz))
  write(*,*) 'number density integral: ', SUM(g_ion(1:ngx,1:ngy,1:ngz))
  write(*,*) 'mass integral: ', SUM(rho_ion(1:ngx,1:ngy,1:ngz))
  g_ele = max(1.,g_ele)
  g_ion = max(1.,g_ion)
  Te = .511*Te/g_ele   ! Temperature in MeV (KE per particle)
  Ti = .511*mass_ratio*Ti/g_ion   ! K.E. in MeV

  jx_ion = jx_ion/g_ion  ! Velocity fields
  jy_ion = jy_ion/g_ion
  jz_ion = jz_ion/g_ion
  Egx = Egx/(g_ion)    ! Normalise fields/potential
  Egy = Egy/(g_ion)
  Egz = Egz/(g_ion)
  phig = phig/(g_ion)



  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints

  ! output data file
  csnap = cdump//'/snap_defs.gmt'
  open(80,file=csnap)

  write(*,'(/a,a)') 'Writing out snapshot defs to ',csnap
  write(80,'(a8,a3,f12.3,a2)') 'TRUN','"',trun,' "'  ! total simulation time 

  nout = 7
  cfout(1) = cdump//'/xy_slice_den'
  cfout(2) = cdump//'/xy_slice_vx'
  cfout(3) = cdump//'/xy_slice_vy'
  cfout(4) = cdump//'/xy_slice_temp'
  cfout(5) = cdump//'/xy_slice_ex'
  cfout(6) = cdump//'/xy_slice_ey'
  cfout(7) = cdump//'/xy_slice_phi'
  cfout(8) = cdump//'/xy_slice_eden'
  cfout(9) = cdump//'/xy_slice_vxe'
  cfout(10) = cdump//'/xy_slice_vye'
  cfout(11) = cdump//'/xy_slice_te'

  do i=0,nout-1
     write(*,'(2a)') 'Writing slice ',cfout(i+1)
     open(20+i,file=cfout(i+1))
  end do

  ! density slice in xy-plane along laser axis: converted to n/nc


  kslice = (zslice(1)-zmin)*rdz
  kmax = .5*zl*rdz    ! average over all z\
  nk = 2*kmax+1
  do j=1,ngy
     do i=1,ngx
  !      write(20,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(abs(rho_ele(i,j,kslice-kmax:kslice+kmax)))/nk
        logrho = log10(max(1.e-5,SUM(rho_ion(i,j,1:ngz))/ngz))
	sumjx = sum(jx_ion(i,j,kslice-1:kslice+1))/3  ! projection average in xy
	sumEx = sum(egx(i,j,kslice-1:kslice+1))/3  ! projection average in xy
	sumjy = sum(jy_ion(i,j,kslice-1:kslice+1))/3  ! projection average in xy
	sumEy = sum(egy(i,j,kslice-1:kslice+1))/3  ! projection average in xy
        write(20,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, logrho
        write(21,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sumjx
        write(22,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sumjy
        write(23,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(ti(i,j,kslice-1:kslice+1))/3
        write(24,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sumEx
        write(25,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sumEy
        write(26,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, phig(i,j,kslice)
  !      write(26,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(te(i,j,1:ngz))/ngz
  !      write(21,'(2f13.4,1pe12.3)') i*dx+xmin, j*dy+ymin, SUM(rho_ion(i,j,1:ngz))/ngz
  !      write(22,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, jx_ele(i,j,kslice)
  !      write(23,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, jy_ele(i,j,kslice)
  !      write(31,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(abs(rho_hot(i,j,1:ngz)))/ngz
  !      write(32,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, jx_hot(i,j,kslice)
  !      write(33,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, jy_hot(i,j,kslice)
  !      write(34,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, jz_hot(i,j,kslice)

     end do
  end do

  do i=0,nout-1
     close(20+i)
  end do


  ! output data file
  nout = 7
  cfout(1) = cdump//'/xz_slice_den'
  cfout(2) = cdump//'/xz_slice_vx'
  cfout(3) = cdump//'/xz_slice_vz'
  cfout(4) = cdump//'/xz_slice_temp'
  cfout(5) = cdump//'/xz_slice_ex'
  cfout(6) = cdump//'/xz_slice_ez'
  cfout(7) = cdump//'/xz_slice_phi'
  cfout(8) = cdump//'/xz_slice_eden'
  cfout(9) = cdump//'/xz_slice_jxe'
  cfout(10) = cdump//'/xz_slice_jze'
  cfout(11) = cdump//'/xz_slice_te'

  do i=0,nout-1
     write(*,'(2a)') 'Writing slice ',cfout(i+1)
     open(20+i,file=cfout(i+1))
 end do

  ! density slice in xz-plane along laser axis: converted to n/nc


  jslice = (yslice(1)-ymin)*rdy
  do k=1,ngz
     do i=1,ngx
        logrho = log10(max(1.e-5,SUM(rho_ion(i,1:ngy,k))/ngy))
	sumjx = sum(jx_ion(i,jslice-1:jslice+1,k))/3  ! projection average in xz
	sumEx = sum(egx(i,jslice-1:jslice+1,k))/3  ! projection average in xz
	sumjz = sum(jz_ion(i,jslice-1:jslice+1,k))/3  ! projection average in xz
	sumEz = sum(egz(i,jslice-1:jslice+1,k))/3  ! projection average in xz
        write(20,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, logrho
        write(21,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, sumjx
        write(22,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, sumjz
        write(23,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, sum(ti(i,jslice-1:jslice+1,k))/3
        write(24,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, sumEx
        write(25,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, sumEz
        write(36,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, phig(i,jslice,k)
 !       write(20,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, abs(rho_ele(i,jslice,k))
 !       write(22,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jx_ele(i,jslice,k)
 !       write(23,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jz_ele(i,jslice,k)
 !       write(26,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, te(i,jslice,k)
     end do
  end do

  do i=0,nout-1
     close(20+i)
  end do

  ! 1D lineouts

  cfout(1) = 'fields_pp/lineout.'//cdump


  write(*,'(2a)') 'Writing 1D lineout ',cfout(1)
  open(20,file=cfout(1))
  jslice = (yslice(1)-ymin)*rdy
  kslice = (zslice(1)-zmin)*rdz
  avene = 0.
  aveni = 0.
  avenhot = 0.
  avejxi = 0.
  avejxe = 0.
  aveti = 0.
  avete = 0.
  aveex = 0.
  nave = 5
  norm = (2*nave+1)**2
  do k=kslice-nave,kslice+nave
     do j=jslice-nave,jslice+nave
        avene(1:ngx) = avene(1:ngx) + rho_ele(1:ngx,j,k)/norm
        aveni(1:ngx) = aveni(1:ngx) + rho_ion(1:ngx,j,k)/norm
        avejxi(1:ngx) = avejxi(1:ngx) + jx_ion(1:ngx,j,k)/norm
        avejxe(1:ngx) = avejxe(1:ngx) + jx_ele(1:ngx,j,k)/norm
        aveex(1:ngx) = aveex(1:ngx) + egx(1:ngx,j,k)/norm
        avete(1:ngx) = avete(1:ngx) + te(1:ngx,j,k)/norm
        aveti(1:ngx) = aveti(1:ngx) + ti(1:ngx,j,k)/norm
        avephi(1:ngx) = avephi(1:ngx) + phig(1:ngx,j,k)/norm
     end do
  end do
  write(20,'(6(1pe16.4))') (i*dx+xmin,  aveni(i), avejxi(i), aveex(i), &
       aveti(i), avephi(i), i=1,ngx )
  close(20)

  cfout(1) = 'fields_pp/mesh.xyz'
  cfout(2) = 'fields_pp/'//cdump//'.xyz'


  write(*,'(2a)') 'Writing 3D dump ',cfout(2)
  open(20,file=cfout(1))
  open(21,file=cfout(2))

  ! write AVS field record
  open(60,file='fields3d.fld',position='append')
  write(60,'(3a)') 'time file='//cfout(2)//' filetype=ascii'
  write(60,'(3a)') 'variable 1 file='//cfout(2)//' filetype=ascii offset=0 stride=7'
  write(60,'(3a)') 'variable 2 file='//cfout(2)//' filetype=ascii offset=1 stride=7'
  write(60,'(3a)') 'variable 3 file='//cfout(2)//' filetype=ascii offset=2 stride=7'
  write(60,'(3a)') 'variable 4 file='//cfout(2)//' filetype=ascii offset=3 stride=7'
  write(60,'(3a)') 'variable 5 file='//cfout(2)//' filetype=ascii offset=4 stride=7'
  write(60,'(3a)') 'variable 6 file='//cfout(2)//' filetype=ascii offset=5 stride=7'
  write(60,'(3a)') 'variable 7 file='//cfout(2)//' filetype=ascii offset=6 stride=7'
  write(60,'(a3)') 'EOT'
  close(60)
  !  iskip3d = 2
  write(20,'(3i6)') ngx/iskip3d,ngy/iskip3d,ngz/iskip3d

  do k=1,ngz,iskip3d
     do j=1,ngy,iskip3d
        do i=1,ngx,iskip3d
           write(20,'(3f7.2)') i*dx+xmin, j*dy+ymin, k*dz+zmin
           write(21,'(7(1pe10.2))') & 
                abs(SUM(rho_ele(i:i+iskip3d-1,j,k))), &
                SUM(rho_ion(i:i+iskip3d-1,j,k))/iskip3d, &
                SUM(jx_ion(i:i+iskip3d-1,j,k))/iskip3d, &
                SUM(jy_ion(i:i+iskip3d-1,j,k))/iskip3d, &
                SUM(jz_ion(i:i+iskip3d-1,j,k))/iskip3d, & 
                SUM(Te(i:i+iskip3d-1,j,k))/iskip3d, & 
                SUM(Ti(i:i+iskip3d-1,j,k))/iskip3d 
        end do
     end do
  end do
  close(20)
  close(21)




  ! ====================================================================================
  !
  ! ELECTRON MOMENTUM SPACE
  !
  ! ====================================================================================
  !
  ! box limits 
  !  uxmin = -10.
  !  uxmax = 10.
  !  uymin = -5.
  !  uymax = 5.
  !  uzmin = -5.
  !  uzmax = 5.
  uxl = uxmax-uxmin
  uyl = uymax-uymin
  uzl = uzmax-uzmin

  dux = uxl/ngux
  duy = uyl/nguy
  duz = uzl/nguz
  rdx = 1./dux
  rdy = 1./duy
  rdz = 1./duz

  xl = xmax-xmin

  dx = xl/ngx

  du =  umevmax/.511/ngux  ! Bin size for energy spectra


  write(50,'(a8,a3,4(f12.3,a1))') 'UPA','"',uxmin,'/',uxmax,'/',uymin,'/',uymax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.3,a1))') 'UDR','"',uxmin+dux,'/',uxmax,'/',uymin+duy,'/',uymax,'"'  ! data region
  write(50,'(a8,a3,2(f12.3,a1))') 'UMESH','"',dux,'/',duy,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.3,a2))') 'UBOX','"X',xbox,'i/',xbox,'i"'   ! plot size in inches 
  write(50,'(a8,a3,2(f12.3,a1))') 'UAXES','"',uxtick,'/',uytick,'"'   ! tick intervals
  write(50,'(a8,a3,3(f12.3,a1))') 'FPMAP','"',0.0,'/',fpmax,'/',fpmax/50.,'"'  ! phase-space map
  write(50,'(a8,a3,f12.3,a1)') 'SSFP','"',fpmax/2,'"'  ! density scale ticks 

  write(50,'(a8,a3,4(f12.3,a1))') 'UXPA','"',xmin,'/',xmax,'/',uxmin,'/',uxmax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.3,a1))') 'UXPACONV','"',xmin*conv_length,'/',xmax*conv_length, &
       '/',uxmin,'/',uxmax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.3,a1))') 'UXDR','"',xmin+dx,'/',xmax,'/',uxmin+dux,'/',uxmax,'"'  ! data region
  write(50,'(a8,a3,2(f12.3,a1))') 'UXMESH','"',dx,'/',dux,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.3,a1))') 'UXAXES','"',xtick,'/',uxtick,'"'   ! tick intervals
  write(50,'(a8,a3,f12.3,a2)') 'UMEVMAX','"',umevmax,' "'  ! Max energy for fhot 

  fpxpy_ele(0:ngux+1,0:nguz+1) = 0.
  fpxpy_ele(0:ngux+1,0:nguy+1) = 0.
  fpxx_ele(0:ngx+1,0:ngux+1) = 0.
  fhote(0:ngux+1) = 0.

  do i=1,n

     xa=(ux(i)-uxmin)*rdx
     ya=(uy(i)-uymin)*rdy
     za=(uz(i)-uzmin)*rdz
     xxa = (x(i)-xmin)/dx

     gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)

     iu = (gamma-1.0)/du+1  ! Energy bin

     !  indices
     ii1 = xxa+1  ! space
     i1=xa+1      ! momenta
     j1=ya+1
     k1=za+1

     iu = min(max(0,iu),ngux+1)
     ii1 = min(max(0,ii1),ngx+1)

     i1 = min(max(0,i1),ngux+1)
     j1 = min(max(0,j1),nguy+1)
     k1 = min(max(0,k1),nguz+1)

     !  count charges at nearest point in phase space

     if (q(i)<0) then
        fpxx_ele(ii1,i1)=fpxx_ele(ii1,i1) + 1
        fpxpy_ele(i1,j1)=fpxpy_ele(i1,j1) + 1
        fpxpz_ele(i1,k1)=fpxpz_ele(i1,k1) + 1
        fhote(iu) = fhote(iu) + 1  ! Energy spectrum
     endif


  end do



  ng = (ngux+2)*(nguy+2)                         ! total # gridpoints

  ! output data files
  nout = 4
  cfout(1) = cdump//'/fpxpy_ele'
  cfout(2) = cdump//'/fpxpz_ele'
  cfout(3) = cdump//'/fhot_ele'
  cfout(4) = cdump//'/fpxx_ele'



  do i=0,nout-1
     write(*,'(2a)') 'Writing slice ',cfout(i+1)
     open(20+i,file=cfout(i+1))
  end do

  fpxpy_ele = max(fpxpy_ele,0.9)
  fpxpz_ele = max(fpxpz_ele,0.9)
  fpxx_ele = max(fpxx_ele,0.9)

  ! px-py, px-pz

  do j=1,nguy
     do i=1,ngux
        write(20,'(2f13.4,f16.6)') i*dux+uxmin, j*duy+uymin, log(fpxpy_ele(i,j))
        write(21,'(2f13.4,f16.6)') i*dux+uxmin, j*duy+uymin, log(fpxpz_ele(i,j))
     end do
  end do

  ! px-x

  do j=1,ngux
     do i=1,ngx
        write(23,'(2f13.4,f16.6)') i*dx+xmin, j*dux+uxmin, log(fpxx_ele(i,j))
     end do
  end do

  ! fhot

  fhote = max(fhote,0.1)
  write(22,'((f13.4,f16.6))') (i*du*0.511, fhote(i), i=1,ngux)

  do i=0,nout-1
     close(20+i)
  end do

  ! ====================================================================================
  !
  ! ION MOMENTUM SPACE
  !
  ! ====================================================================================

  ! box limits - ion momenta 
  uxl = uximax-uximin
  uyl = uyimax-uyimin
  uzl = uzimax-uzimin

  dux = uxl/ngux
  duy = uyl/nguy
  duz = uzl/nguz
  rdx = 1./dux
  rdy = 1./duy
  rdz = 1./duz
  xl = xmax-xmin

  dx = xl/ngx
  du =  uimevmax/.511/mass_ratio/ngux  ! Bin size for energy spectra
  !  aimax=90.
  !  aimin=-90.
  dalpha = 2*aimax/nalpha ! Bin size for angular distrib.

  write(50,'(a8,a3,4(f12.5,a1))') 'UIPA','"',uximin,'/',uximax,'/',uyimin,'/',uyimax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.5,a1))') 'UIDR','"',uximin+dux,'/',uximax,'/',uyimin+duy,'/',uyimax,'"'  ! data region
  write(50,'(a8,a3,2(f12.5,a1))') 'UIMESH','"',dux,'/',duy,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.5,a2))') 'UIBOX','"X',xbox,'i/',xbox,'i"'   ! plot size in inches 
  write(50,'(a8,a3,2(f12.5,a1))') 'UIAXES','"',uxitick,'/',uyitick,'"'   ! tick intervals
  write(50,'(a8,a3,3(f12.5,a1))') 'FPIMAP','"',0.0,'/',fpmax,'/',fpmax/50.,'"'  ! phase-space map
  write(50,'(a8,a3,f12.3,a1)') 'SSIFP','"',fpmax/2,'"'  ! density scale ticks 

  write(50,'(a8,a3,4(f12.5,a1))') 'UXIPA','"',xmin,'/',xmax,'/',uximin,'/',uximax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.5,a1))') 'UIPACONV','"',xmin*conv_length,'/',xmax*conv_length,&
       '/',uximin,'/',uximax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.5,a1))') 'UXIDR','"',xmin+dx,'/',xmax,'/',uximin+dux,'/',uximax,'"'  ! data region
  write(50,'(a8,a3,2(f12.5,a1))') 'UXIMESH','"',dx,'/',dux,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.5,a1))') 'UXIAXES','"',xtick,'/',uxitick,'"'   ! tick intervals
  write(50,'(a8,a3,4(f12.5,a1))') 'UZIPA','"',xmin,'/',xmax,'/',uzimin,'/',uzimax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.5,a1))') 'UZIPAC','"',xmin*conv_length,'/',xmax*conv_length, &
       '/',uzimin,'/',uzimax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.5,a1))') 'UZIDR','"',xmin+dx,'/',xmax,'/',uzimin+duz,'/',uzimax,'"'  ! data region
  write(50,'(a8,a3,2(f12.5,a1))') 'UZIMESH','"',dx,'/',duz,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.5,a1))') 'UZIAXES','"',xtick,'/',uzitick,'"'   ! tick intervals
  write(50,'(a8,a3,f12.3,a2)') 'UIMEVMAX','"',uimevmax,' "'  ! Max energy for fihot
  write(50,'(a8,a3,4(f12.5,a1))') 'AIPA','"',-aimax,'/',aimax,'/',-aimax,'/',aimax,'"'  ! plot area
  write(50,'(a8,a3,4(f12.5,a1))') 'AIDR','"',-aimax+dalpha,'/',aimax,'/',-aimax+dalpha,'/',aimax,'"'  ! data region
  write(50,'(a8,a3,2(f12.5,a1))') 'AIMESH','"',dalpha,'/',dalpha,'"'   ! mesh size 
  write(50,'(a8,a3,2(f12.5,a2))') 'AIBOX','"X',xbox,'i/',xbox,'i"'   ! plot size in inches 
  write(50,'(a8,a3,2(f12.5,a1))') 'AIAXES','"',aimax/6.,'/',aimax/6.,'"'   ! tick intervals
  write(50,'(a8,a3,f12.3,a2)') 'UIMEV1','"',uimev1,' "'  ! ion energy intervals 
  write(50,'(a8,a3,f12.3,a2)') 'UIMEV2','"',uimev2,' "'  ! for angular plots 

  fpxpy_ion(0:ngux+1,0:nguz+1) = 0.
  fpxpy_ion(0:ngux+1,0:nguy+1) = 0.
  fpxx_ion(0:ngx+1,0:ngux+1) = 0.
  fpzx_ion(0:ngx+1,0:nguz+1) = 0.
  fang_ion(0:nalpha+1,0:nalpha+1) = 0.
  uang_ion(0:nalpha+1,0:nalpha+1) = 0.
  fang_ion1(0:nalpha+1,0:nalpha+1) = 0.
  fang_ion2(0:nalpha+1,0:nalpha+1) = 0.
  fang_back(0:nalpha+1,0:nalpha+1) = 0.

  fhoti_for = 0.
  fhoti_back = 0.
  ntnsa=0
  utnsa=0.
  nfront = 0
  ufront = 0.
  nblow = 0
  ublow = 0.
  uf_max = 0.
  ur_max = 0.
  ub_max = 0.

  do i=1,n
     if (q(i)>0) then

        xa=(ux(i)-uximin)*rdx
        ya=(uy(i)-uyimin)*rdy
        za=(uz(i)-uzimin)*rdz
        xxa = (x(i)-xmin)/dx

        !  indices
        i1=xa+1
        j1=ya+1
        k1=za+1
        ii1 = xxa+1  ! space

        i1 = min(max(0,i1),ngux+1)
        j1 = min(max(0,j1),nguy+1)
        k1 = min(max(0,k1),nguz+1)
        ii1 = min(max(0,ii1),ngx+1)

        !  count charges at nearest point in phase space

        fpxx_ion(ii1,i1)=fpxx_ion(ii1,i1) + 1
        fpzx_ion(ii1,j1)=fpzx_ion(ii1,j1) + 1

        fpxpy_ion(i1,j1)=fpxpy_ion(i1,j1) + 1

        fpxpz_ion(i1,k1)=fpxpz_ion(i1,k1) + 1

        gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
        umev  = .511*mass_ratio*(gamma-1.0)

        iu = (gamma-1.0)/du+1  ! Energy bin
        iu = min(max(0,iu),ngux+1)  ! limits
        if (ux(i)>0 .and. umev > 1. ) then     ! Forward ions > 1MeV
           alphay=atan(uy(i)/ux(i))*180./pi  ! convert to degrees
           alphaz=atan(uz(i)/ux(i))*180./pi
           iaz = (alphaz-aimin)/dalpha+1  ! bin
           iaz = min(max(0,iaz),nalpha+1)
           iay = (alphay-aimin)/dalpha+1  ! bin
           iay = min(max(0,iay),nalpha+1)
           fang_ion(iay,iaz) = fang_ion(iay,iaz) + 1
           uang_ion(iay,iaz) = uang_ion(iay,iaz) + umev  ! Mean ion energy
        endif
        !        if (ux(i)>0 .and. umev > uimev1 .and. x(i)>40 &
        !		.and. z(i)>40 .and. z(i) < 140) then     ! Forward TNSA ions > uimev1 
        if (ux(i)>0 .and. umev > uimev1 .and. umev <=uimev2 ) then
           alphay=atan(uy(i)/ux(i))*180./pi  ! convert to degrees
           alphaz=atan(uz(i)/ux(i))*180./pi
           iaz = (alphaz-aimin)/dalpha+1  ! bin
           iaz = min(max(0,iaz),nalpha+1)
           iay = (alphay-aimin)/dalpha+1  ! bin
           iay = min(max(0,iay),nalpha+1)
           fang_ion1(iay,iaz) = fang_ion1(iay,iaz) + 1
	endif
        if (ux(i)>0 .and. umev > uimev1 .and. x(i) > rear_edge-10. )  then ! TNSA >Umev1 ions 
           !             .and. z(i)>40 .and. z(i) < 140 ) then     
           alphay=atan(uy(i)/ux(i))*180./pi  ! convert to degrees
           alphaz=atan(uz(i)/ux(i))*180./pi
           iaz = (alphaz-aimin)/dalpha+1  ! bin
           iaz = min(max(0,iaz),nalpha+1)
           iay = (alphay-aimin)/dalpha+1  ! bin
           iay = min(max(0,iay),nalpha+1)
           fang_ion2(iay,iaz) = fang_ion2(iay,iaz) + 1
	endif
        if (ux(i)<0 .and. umev > uimev1  ) then     ! Backward > uimev1 ions
           alphay=atan(-uy(i)/ux(i))*180./pi  ! convert to degrees
           alphaz=atan(-uz(i)/ux(i))*180./pi
           iaz = (alphaz-aimin)/dalpha+1  ! bin
           iaz = min(max(0,iaz),nalpha+1)
           iay = (alphay-aimin)/dalpha+1  ! bin
           iay = min(max(0,iay),nalpha+1)
           fang_back(iay,iaz) = fang_back(iay,iaz) + 1
        endif
        if (x(i) > rear_edge .and. umev > uimev1/2) then
           utnsa=utnsa + umev
           ntnsa=ntnsa + 1
	endif
	if (x(i) > rear_edge) then
           ur_max = max(ur_max,umev)
        endif

 	if ( x(i) < rear_edge-10. .and. umev>uimev1/2) then
           if (ux(i)>0) then  ! forward ions
              ufront = ufront + umev
              nfront = nfront + 1
              uf_max = max(uf_max,umev)
           else
              nblow = nblow+1
              ublow = ublow + umev
              ub_max = max(ub_max,umev)
           endif
        endif

        if (ux(i)>0) then
           fhoti_for(iu) = fhoti_for(iu) + 1
        else
           fhoti_back(iu) = fhoti_back(iu) + 1
        endif

     endif


  end do

  uang_ion = uang_ion/max(1.,fang_ion)   ! ave. KE angular distn 

  ng = (ngux+2)*(nguy+2)                         ! total # gridpoints

  ! output data files
  nout = 4
  cfout(1) = cdump//'/fpxpy_ion'
  cfout(2) = cdump//'/fpxpz_ion'
  cfout(3) = cdump//'/fpxx_ion'
  cfout(4) = cdump//'/fpzx_ion'

  cfout(5) = cdump//'/fhot_ion_forw'
  cfout(6) = cdump//'/fhot_ion_back'

  cfout(7) = cdump//'/fang_ion'
  cfout(8) = cdump//'/fang_ion1'
  cfout(9) = cdump//'/fang_ion2'
  cfout(10) = cdump//'/fang_back'


  do i=0,nout-1
     write(*,'(2a)') 'Writing slice ',cfout(i+1)
     open(20+i,file=cfout(i+1))
  end do

  fpxpy_ion = max(fpxpy_ion,0.9)
  fpxpz_ion = max(fpxpz_ion,0.9)
  fpxx_ion = max(fpxx_ion,0.9)
  fpzx_ion = max(fpzx_ion,0.9)

  do j=1,nguy
     do i=1,ngux
        write(20,'(2f13.4,f16.6)') i*dux+uximin, j*duy+uyimin, log(fpxpy_ion(i,j))
        write(21,'(2f13.4,f16.6)') i*dux+uximin, j*duy+uyimin, log(fpxpz_ion(i,j))
     end do
  end do


  ! px-x

  do j=1,ngux
     do i=1,ngx
        write(22,'(2f13.4,f16.6)') i*dx+xmin, j*dux+uximin, log(fpxx_ion(i,j))
     end do
  end do

  do j=1,nguz
     do i=1,ngx
        write(23,'(2f13.4,f16.6)') i*dx+xmin, j*duz+uzimin, log(fpzx_ion(i,j))
     end do
  end do


  do i=0,nout-1
     close(20+i)
  end do


  !  Close GMT grid_defs file
  close(50)
  !  Close GMT snapshot defs file
  close(80)
  ! Close protocol/info file
  close(10)
end program ppfields
