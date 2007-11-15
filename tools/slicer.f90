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
    ! 9/10 Added arrays for 3rd species (protons):
    !  rho_pro, Tp, g_pro, fpxx_pro, fhotp
    real, allocatable :: rho_ion(:,:,:), rho_ele(:,:,:), rho_pro(:,:,:) 
    real, allocatable :: jx_ele(:,:,:), jx_ion(:,:,:),&
        jy_ele(:,:,:), jy_ion(:,:,:), jz_ele(:,:,:), jz_ion(:,:,:)
    real, allocatable :: rho_hot(:,:,:), jx_hot(:,:,:), jy_hot(:,:,:), jz_hot(:,:,:)
    real, allocatable :: Egx(:,:,:), Egy(:,:,:), Egz(:,:,:), &
        phig(:,:,:), phi_whole(:,:,:)
    real, allocatable :: Te(:,:,:), Ti(:,:,:), Tp(:,:,:)
    real, allocatable :: g_ion(:,:,:), g_ele(:,:,:), g_pro(:,:,:), g_w(:,:,:)  ! particle weights for averages

    real, allocatable :: fpxpy_ion(:,:), fpxpy_ele(:,:)
    real, allocatable :: fpxpz_ion(:,:), fpxpz_ele(:,:)
    real, allocatable :: fpxx_ion(:,:), fpyy_ion(:,:), fpxx_pro(:,:), fpxx_ele(:,:), fpzx_ion(:,:), fpyy_pro(:,:)
    real, allocatable :: fang_ion(:,:), fang_ion1(:,:), fang_ion2(:,:), fang_back(:,:)
    real, allocatable :: uang_ion(:,:)

    real, allocatable :: x(:),y(:),z(:),ux(:),uy(:),uz(:),q(:),m(:),ex(:),ey(:),ez(:),bz(:),phi(:)
    integer, allocatable :: own(:), label(:), index(:)
    real, allocatable :: aveni(:),avene(:),avejxi(:),avejxe(:),aveex(:),aveti(:),avete(:),avephi(:)
    real, allocatable :: fhote(:), fhoti_for(:), fhoti_back(:), fhotp(:), avenhot(:)

    real, dimension(10) :: xslice, yslice, zslice

    real :: xmin, xmax, ymin, ymax, zmin, zmax
    real :: rdx, rdy, rdz, dx, dy, dz, cweight, jxweight, jyweight, jzweight
    real :: omega, tweight, exweight, eyweight, yya
    real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za, xl, yl, zl, fr1, fr2, xxa
    real :: uxmin, uxmax, uymin, uymax, uzmin, uzmax, uxl, uyl, uzl, dux, duy, duz, du
    real :: uximin, uximax, uyimin, uyimax, uzimin, uzimax, dalpha
    real :: temin, temax, timin, timax, jemax, jimax, emax, fpmax, rhomax, tcold
    real :: xbox, ybox, zbox, yshift, zshift, pshift, jevec, jivec, evec, epsilon
    real :: xtick, ytick, ztick, uxtick, uytick, uztick, uxitick, uyitick, uzitick
    real :: umevmax, uimevmax, umev, uimev1, uimev2, xrear

    integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, islice, jslice, kslice, nave, norm
    integer :: iu, ii1, jj1, kmax, nk
    character(30) :: cfile,cfile1, cfile2, cinfo, csnap
    character(30), dimension(25) :: cfout

    character(5) :: cme
    character(6) :: cdump, cvis
    logical :: found
    real :: xc1, rho_track, gamma, mass_ratio, phimin,phimax, logrho

    integer :: icm, icp, nout, nelecs, nions, nsel

    integer ::   start_step, npartr,  ner, nir, np_beamr, iconf, iens
    real :: xlr, ylr, zlr, boxr, epsr, thetar, tlaser, trun, lambda, Qs, sigma
    real :: aimax, aimin, alphaz, alphay
    real :: xpmin, xpmax, ypmin, ypmax, zpmin, zpmax, xptick, yptick, xd, yd, zd
    integer :: iaz, iay, iskip3d

    ! Conversion factors
    real :: cowp_micron  ! lengths c/wp -> microns
    real :: wpr_fs  ! time 1/wp -> fs
    real :: C_Ne ! # sim. particles -> real # particles
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
    read(20,'(10(9x,f12.5/))') uxmin, uxmax, uxtick, uymin, uymax, uytick, uzmin, uzmax, uztick, umevmax
    read(20,'(14(9x,f12.5/))') uximin, uximax, uxitick, uyimin, uyimax, uyitick, uzimin, uzimax, uzitick, uimevmax,aimin,aimax, uimev1, uimev2
    read(20,'(19(9x,f12.5/))') xslice(1),yslice(1), zslice(1), rear_edge, mass_ratio, rhomax, temin, temax, timin, timax, &
        jemax, jimax, emax, fpmax, jevec, jivec, evec, epsr, tcold
    read(20,'(6(9x,f12.5/))') xbox, ybox, zbox, yshift, zshift, pshift    ! plot dimensions and scale positions in inches
    read(20,'(9x,i6/)') iskip3d  ! skip stride for writing out 3D xyz plots
    close(20)
    write(*,*) n,ngx,ngy,ngz,ngux,nguy,nguz
    write(*,*) xmin,xmax,ymin,ymax,zmin,zmax
    write(*,*) uzmax, uzimax
    write(*,*) zslice(1), mass_ratio, sigma

    ! open GMT defs file 
    cfile = 'grid_defs.gmt'
    open(50,file=cfile)

    allocate ( x(n), y(n), z(n), ux(n), uy(n), uz(n), q(n), m(n), ex(n), ey(n), ez(n), bz(n), phi(n), own(n), label(n), index(n) )
    allocate ( rho_ion(0:ngx+1,0:ngy+1,0:ngz+1), rho_ele(0:ngx+1,0:ngy+1,0:ngz+1), &
        rho_pro(0:ngx+1,0:ngy+1,0:ngz+1), &
        jx_ele(0:ngx+1,0:ngy+1,0:ngz+1), jx_ion(0:ngx+1,0:ngy+1,0:ngz+1), &
        jz_ele(0:ngx+1,0:ngy+1,0:ngz+1), jz_ion(0:ngx+1,0:ngy+1,0:ngz+1), &
        jy_ele(0:ngx+1,0:ngy+1,0:ngz+1), jy_ion(0:ngx+1,0:ngy+1,0:ngz+1) )
    allocate ( Egx(0:ngx+1,0:ngy+1,0:ngz+1), Egy(0:ngx+1,0:ngy+1,0:ngz+1), &
        Egz(0:ngx+1,0:ngy+1,0:ngz+1), &
        Te(0:ngx+1,0:ngy+1,0:ngz+1), Ti(0:ngx+1,0:ngy+1,0:ngz+1),&
        Tp(0:ngx+1,0:ngy+1,0:ngz+1),&
        phig(0:ngx+1,0:ngy+1,0:ngz+1), phi_whole(0:ngx+1,0:ngy+1,0:ngz+1)) 
    allocate ( rho_hot(0:ngx+1,0:ngy+1,0:ngz+1), jx_hot(0:ngx+1,0:ngy+1,0:ngz+1), &
        jy_hot(0:ngx+1,0:ngy+1,0:ngz+1), jz_hot(0:ngx+1,0:ngy+1,0:ngz+1) )
    allocate ( g_ion(0:ngx+1,0:ngy+1,0:ngz+1), g_ele(0:ngx+1,0:ngy+1,0:ngz+1), &
        g_pro(0:ngx+1,0:ngy+1,0:ngz+1), g_w(0:ngx+1,0:ngy+1,0:ngz+1) )
    allocate ( fpxpy_ion(0:ngux+1,0:nguy+1), fpxpy_ele(0:ngux+1,0:nguy+1), &
        fpxpz_ion(0:ngux+1,0:nguz+1), fpxpz_ele(0:ngux+1,0:nguz+1), &
        fpxx_ion(0:ngx+1,0:ngux+1), fpxx_ele(0:ngx+1,0:ngux+1), &
        fpxx_pro(0:ngx+1,0:ngux+1), fpyy_pro(0:ngy+1,0:nguy+1), &
        fpzx_ion(0:ngx+1,0:nguz+1), fpyy_ion(0:ngy+1,0:nguy+1), &
        fang_ion(0:nalpha+1,0:nalpha+1),   fang_ion2(0:nalpha+1,0:nalpha+1),&
        fang_back(0:nalpha+1,0:nalpha+1), fang_ion1(0:nalpha+1,0:nalpha+1), &
        uang_ion(0:nalpha+1,0:nalpha+1), &
        fhote(0:ngux+1), fhoti_for(0:ngux+1), fhoti_back(0:ngux+1), fhotp(0:ngux+1) )
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

    read(80,'(7(9x,i8/),11(9x,f12.5/))')  &    ! info block - skip variable names
        start_step, npartr, &
        ner, nir, np_beamr, iconf, iens, &
        xlr, ylr, zlr, boxr, &
        epsilon, &  ! Coulomb smoothing radius
        thetar, &   ! Force clumping parameter
        tlaser, &            ! laser duration
        trun, &     ! Total simulation time
        omega, &    ! ratio laser frequency/ plasma frequency
        lambda, &   ! laser wavelength in microns
                                !       sigma, &    ! spotsize  TODO: put this back into info file!
        Qs          ! electron macrocharge

    close(80)

    ! Compute time, length, number conversion factors

    cowp_micron = lambda/2./pi*omega
    wpr_fs = 5./3./pi*lambda*omega
    C_Ne = 1.e9/(2.*pi)**3*lambda*omega*Qs
    C_pot = 0.511

    write(10,*) 'Conversion factors: '

    write(10,*) 'Length           ',cowp_micron
    write(10,*) 'Time             ',wpr_fs
    write(10,*) 'Particle number  ',C_Ne
    write(10,*) 'Potential  ',C_pot
    write(10,*)

    cfile = 'dumps/parts_dump.'//cdump
    write(*,*) 'Reading data from ',cfile  ! original ES format
    open(20,file=cfile)
    read(20,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),ex(i),ey(i),ez(i),phi(i),own(i),label(i),i=1,n)
    close(20)



    ! Do label check
    !  write (*,*) 'Doing label check ..'
    !  open(22,file="sorted.labels")
    !  call indsort(label,index,n)
    !  do i=1,n
    !     write(22,'(3i9)') i,label(i),label(index(i))
    !  end do
    !  close(22)
    !  write(*,*) '.. done'

    ! Dump for AVS-EXPRESS
    cfile = 'dumps/ions_dump.'//cdump
    cfile1 = 'dumps/elecs_dump.'//cdump
    cfile2 = 'dumps/protons_dump.'//cdump

    write(*,*) 'Writing ions to ',cfile
    write(*,*) 'Writing electrons to ',cfile1
    write(*,*) 'Writing protons to ',cfile2
    open(20,file=cfile)
    open(21,file=cfile1)
    open(22,file=cfile2)
    nsel=0
    do i=1,n
        gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
        umev  = .511*(gamma-1.0)
        if (q(i)>0 .and. m(i)/q(i)>1850. ) then
            nsel=nsel+1
            write(20,'(8(1pe12.4))') x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i)  ! ions
        else if (q(i)>0.) then
            write(22,'(8(1pe12.4))') x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i)  ! protons
        else
            write(21,'(8(1pe12.4))') x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i)  ! electrons
        endif
    end do

    close(20)
    close(21)
    close(22)
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

    ! grid defs for yz slices

    write(50,'(a8,a3,4(f12.3,a1))') 'YZPA','"',ymin,'/',ymax,'/',zmin,'/',zmax,'"'  ! plot area
    write(50,'(a8,a3,4(f12.3,a1))') 'YZDR','"',ymin+dy,'/',ymax,'/',zmin+dz,'/',zmax,'"'  ! data region 
    write(50,'(a8,a3,2(f12.3,a1))') 'YZMESH','"',dy,'/',dz,'"'   ! mesh size 
    write(50,'(a8,a3,2(f12.3,a2))') 'YZBOX','"X',ybox,'i/',zbox,'i"'   ! plot size in inches 
    write(50,'(a8,a3,2(f12.3,a1))') 'YZAXES','"',ytick,'/',ztick,'"'   ! tick intervals

    write(50,'(a8,a3,3(f12.3,a1))') 'RHOMAP','"',-rhomax/10.,'/',rhomax,'/',rhomax/50.,'"'  ! density map
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

    do k=0,ngz+1
        do j=0,ngy+1
            do i=0,ngx+1
                Te(i,j,k)=0.
                Ti(i,j,k)=0.
                Tp(i,j,k)=0.
                rho_ele(i,j,k) = 0.
                rho_ion(i,j,k) = 0.
                rho_pro(i,j,k) = 0.
                jx_ele(i,j,k) = 0.
                jx_ion(i,j,k) = 0.
                jy_ele(i,j,k) = 0.
                jy_ion(i,j,k) = 0.
                jz_ele(i,j,k) = 0.
                jz_ion(i,j,k) = 0.
                g_ele(i,j,k) = 0.
                g_ion(i,j,k) = 0.
                g_pro(i,j,k) = 0.
                Egx(i,j,k) = 0.
                Egy(i,j,k) = 0.
                Egz(i,j,k) = 0.
                phig(i,j,k) = 0.
                rho_hot(i,j,k) = 0.
                jx_hot(i,j,k) = 0.
                jy_hot(i,j,k) = 0.
                jz_hot(i,j,k) = 0.
            end do
        end do
    end do

    do i=1,n

        xa=(x(i)-xmin)*rdx
        ya=(y(i)-ymin)*rdy
        za=(z(i)-zmin)*rdz
        gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)

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
        fx2=i1-xa
        fx1=1.-fx2
        fy2=j1-ya
        fy1=1.-fy2
        fz2=k1-za
        fz1=1.-fz2

        !  gather charge at nearest grid points

        cweight = abs(q(i))*rdx*rdy*rdz/omega**2       ! charge weighting factor
        jxweight = cweight*ux(i)/gamma
        jyweight = cweight*uy(i)/gamma
        jzweight = cweight*uz(i)/gamma
        tweight = gamma-1.  ! K.E. of particle in keV
        fr1 = sqrt(fx1**2+fy1**2+fz1**2+epsilon**2)
        fr2 = sqrt(fx2**2+fy2**2+fz2**2+epsilon**2)


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

            if (511.*tweight <= tcold) then
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

                ! Hot electrons
                rho_hot(i1,j1,k1)=rho_hot(i1,j1,k1) + cweight*fx1*fy1*fz1
                rho_hot(i2,j1,k1)=rho_hot(i2,j1,k1) + cweight*fx2*fy1*fz1
                rho_hot(i1,j2,k1)=rho_hot(i1,j2,k1) + cweight*fx1*fy2*fz1
                rho_hot(i2,j2,k1)=rho_hot(i2,j2,k1) + cweight*fx2*fy2*fz1
                rho_hot(i1,j1,k2)=rho_hot(i1,j1,k2) + cweight*fx1*fy1*fz2
                rho_hot(i2,j1,k2)=rho_hot(i2,j1,k2) + cweight*fx2*fy1*fz2
                rho_hot(i1,j2,k2)=rho_hot(i1,j2,k2) + cweight*fx1*fy2*fz2
                rho_hot(i2,j2,k2)=rho_hot(i2,j2,k2) + cweight*fx2*fy2*fz2

                jx_hot(i1,j1,k1)=jx_hot(i1,j1,k1) + jxweight*fx1*fy1*fz1
                jx_hot(i2,j1,k1)=jx_hot(i2,j1,k1) + jxweight*fx2*fy1*fz1
                jx_hot(i1,j2,k1)=jx_hot(i1,j2,k1) + jxweight*fx1*fy2*fz1
                jx_hot(i2,j2,k1)=jx_hot(i2,j2,k1) + jxweight*fx2*fy2*fz1
                jx_hot(i1,j1,k2)=jx_hot(i1,j1,k2) + jxweight*fx1*fy1*fz2
                jx_hot(i2,j1,k2)=jx_hot(i2,j1,k2) + jxweight*fx2*fy1*fz2
                jx_hot(i1,j2,k2)=jx_hot(i1,j2,k2) + jxweight*fx1*fy2*fz2
                jx_hot(i2,j2,k2)=jx_hot(i2,j2,k2) + jxweight*fx2*fy2*fz2

                jy_hot(i1,j1,k1)=jy_hot(i1,j1,k1) + jyweight*fx1*fy1*fz1
                jy_hot(i2,j1,k1)=jy_hot(i2,j1,k1) + jyweight*fx2*fy1*fz1
                jy_hot(i1,j2,k1)=jy_hot(i1,j2,k1) + jyweight*fx1*fy2*fz1
                jy_hot(i2,j2,k1)=jy_hot(i2,j2,k1) + jyweight*fx2*fy2*fz1
                jy_hot(i1,j1,k2)=jy_hot(i1,j1,k2) + jyweight*fx1*fy1*fz2
                jy_hot(i2,j1,k2)=jy_hot(i2,j1,k2) + jyweight*fx2*fy1*fz2
                jy_hot(i1,j2,k2)=jy_hot(i1,j2,k2) + jyweight*fx1*fy2*fz2
                jy_hot(i2,j2,k2)=jy_hot(i2,j2,k2) + jyweight*fx2*fy2*fz2

                jz_hot(i1,j1,k1)=jz_hot(i1,j1,k1) + jzweight*fx1*fy1*fz1
                jz_hot(i2,j1,k1)=jz_hot(i2,j1,k1) + jzweight*fx2*fy1*fz1
                jz_hot(i1,j2,k1)=jz_hot(i1,j2,k1) + jzweight*fx1*fy2*fz1
                jz_hot(i2,j2,k1)=jz_hot(i2,j2,k1) + jzweight*fx2*fy2*fz1
                jz_hot(i1,j1,k2)=jz_hot(i1,j1,k2) + jzweight*fx1*fy1*fz2
                jz_hot(i2,j1,k2)=jz_hot(i2,j1,k2) + jzweight*fx2*fy1*fz2
                jz_hot(i1,j2,k2)=jz_hot(i1,j2,k2) + jzweight*fx1*fy2*fz2
                jz_hot(i2,j2,k2)=jz_hot(i2,j2,k2) + jzweight*fx2*fy2*fz2

            endif
            Te(i1,j1,k1)=Te(i1,j1,k1) + tweight*fx1*fy1*fz1
            Te(i2,j1,k1)=Te(i2,j1,k1) + tweight*fx2*fy1*fz1
            Te(i1,j2,k1)=Te(i1,j2,k1) + tweight*fx1*fy2*fz1
            Te(i2,j2,k1)=Te(i2,j2,k1) + tweight*fx2*fy2*fz1
            Te(i1,j1,k2)=Te(i1,j1,k2) + tweight*fx1*fy1*fz2
            Te(i2,j1,k2)=Te(i2,j1,k2) + tweight*fx2*fy1*fz2
            Te(i1,j2,k2)=Te(i1,j2,k2) + tweight*fx1*fy2*fz2
            Te(i2,j2,k2)=Te(i2,j2,k2) + tweight*fx2*fy2*fz2

        else if (q(i)>0 .and. m(i)/q(i)<1850.) then
            ! protons
            g_pro(i1,j1,k1)=g_pro(i1,j1,k1) + fx1*fy1*fz1  
            g_pro(i2,j1,k1)=g_pro(i2,j1,k1) + fx2*fy1*fz1
            g_pro(i1,j2,k1)=g_pro(i1,j2,k1) + fx1*fy2*fz1
            g_pro(i2,j2,k1)=g_pro(i2,j2,k1) + fx2*fy2*fz1
            g_pro(i1,j1,k2)=g_pro(i1,j1,k2) + fx1*fy1*fz2
            g_pro(i2,j1,k2)=g_pro(i2,j1,k2) + fx2*fy1*fz2
            g_pro(i1,j2,k2)=g_pro(i1,j2,k2) + fx1*fy2*fz2
            g_pro(i2,j2,k2)=g_pro(i2,j2,k2) + fx2*fy2*fz2

            rho_pro(i1,j1,k1)=rho_pro(i1,j1,k1) + cweight*fx1*fy1*fz1
            rho_pro(i2,j1,k1)=rho_pro(i2,j1,k1) + cweight*fx2*fy1*fz1
            rho_pro(i1,j2,k1)=rho_pro(i1,j2,k1) + cweight*fx1*fy2*fz1
            rho_pro(i2,j2,k1)=rho_pro(i2,j2,k1) + cweight*fx2*fy2*fz1
            rho_pro(i1,j1,k2)=rho_pro(i1,j1,k2) + cweight*fx1*fy1*fz2
            rho_pro(i2,j1,k2)=rho_pro(i2,j1,k2) + cweight*fx2*fy1*fz2
            rho_pro(i1,j2,k2)=rho_pro(i1,j2,k2) + cweight*fx1*fy2*fz2
            rho_pro(i2,j2,k2)=rho_pro(i2,j2,k2) + cweight*fx2*fy2*fz2

            ! proton temp
            Tp(i1,j1,k1)=Tp(i1,j1,k1) + tweight*fx1*fy1*fz1
            Tp(i2,j1,k1)=Tp(i2,j1,k1) + tweight*fx2*fy1*fz1
            Tp(i1,j2,k1)=Tp(i1,j2,k1) + tweight*fx1*fy2*fz1
            Tp(i2,j2,k1)=Tp(i2,j2,k1) + tweight*fx2*fy2*fz1
            Tp(i1,j1,k2)=Tp(i1,j1,k2) + tweight*fx1*fy1*fz2
            Tp(i2,j1,k2)=Tp(i2,j1,k2) + tweight*fx2*fy1*fz2
            Tp(i1,j2,k2)=Tp(i1,j2,k2) + tweight*fx1*fy2*fz2
            Tp(i2,j2,k2)=Tp(i2,j2,k2) + tweight*fx2*fy2*fz2

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

            jx_ion(i1,j1,k1)=jx_ion(i1,j1,k1) + jxweight*fx1*fy1*fz1
            jx_ion(i2,j1,k1)=jx_ion(i2,j1,k1) + jxweight*fx2*fy1*fz1
            jx_ion(i1,j2,k1)=jx_ion(i1,j2,k1) + jxweight*fx1*fy2*fz1
            jx_ion(i2,j2,k1)=jx_ion(i2,j2,k1) + jxweight*fx2*fy2*fz1
            jx_ion(i1,j1,k2)=jx_ion(i1,j1,k2) + jxweight*fx1*fy1*fz2
            jx_ion(i2,j1,k2)=jx_ion(i2,j1,k2) + jxweight*fx2*fy1*fz2
            jx_ion(i1,j2,k2)=jx_ion(i1,j2,k2) + jxweight*fx1*fy2*fz2
            jx_ion(i2,j2,k2)=jx_ion(i2,j2,k2) + jxweight*fx2*fy2*fz2

            jy_ion(i1,j1,k1)=jy_ion(i1,j1,k1) + jyweight*fx1*fy1*fz1
            jy_ion(i2,j1,k1)=jy_ion(i2,j1,k1) + jyweight*fx2*fy1*fz1
            jy_ion(i1,j2,k1)=jy_ion(i1,j2,k1) + jyweight*fx1*fy2*fz1
            jy_ion(i2,j2,k1)=jy_ion(i2,j2,k1) + jyweight*fx2*fy2*fz1
            jy_ion(i1,j1,k2)=jy_ion(i1,j1,k2) + jyweight*fx1*fy1*fz2
            jy_ion(i2,j1,k2)=jy_ion(i2,j1,k2) + jyweight*fx2*fy1*fz2
            jy_ion(i1,j2,k2)=jy_ion(i1,j2,k2) + jyweight*fx1*fy2*fz2
            jy_ion(i2,j2,k2)=jy_ion(i2,j2,k2) + jyweight*fx2*fy2*fz2

            jz_ion(i1,j1,k1)=jz_ion(i1,j1,k1) + jzweight*fx1*fy1*fz1
            jz_ion(i2,j1,k1)=jz_ion(i2,j1,k1) + jzweight*fx2*fy1*fz1
            jz_ion(i1,j2,k1)=jz_ion(i1,j2,k1) + jzweight*fx1*fy2*fz1
            jz_ion(i2,j2,k1)=jz_ion(i2,j2,k1) + jzweight*fx2*fy2*fz1
            jz_ion(i1,j1,k2)=jz_ion(i1,j1,k2) + jzweight*fx1*fy1*fz2
            jz_ion(i2,j1,k2)=jz_ion(i2,j1,k2) + jzweight*fx2*fy1*fz2
            jz_ion(i1,j2,k2)=jz_ion(i1,j2,k2) + jzweight*fx1*fy2*fz2
            jz_ion(i2,j2,k2)=jz_ion(i2,j2,k2) + jzweight*fx2*fy2*fz2

            ! ion temp
            Ti(i1,j1,k1)=Ti(i1,j1,k1) + tweight*fx1*fy1*fz1
            Ti(i2,j1,k1)=Ti(i2,j1,k1) + tweight*fx2*fy1*fz1
            Ti(i1,j2,k1)=Ti(i1,j2,k1) + tweight*fx1*fy2*fz1
            Ti(i2,j2,k1)=Ti(i2,j2,k1) + tweight*fx2*fy2*fz1
            Ti(i1,j1,k2)=Ti(i1,j1,k2) + tweight*fx1*fy1*fz2
            Ti(i2,j1,k2)=Ti(i2,j1,k2) + tweight*fx2*fy1*fz2
            Ti(i1,j2,k2)=Ti(i1,j2,k2) + tweight*fx1*fy2*fz2
            Ti(i2,j2,k2)=Ti(i2,j2,k2) + tweight*fx2*fy2*fz2


        endif
        ! electric fields - include all species
        Egx(i1,j1,k1)=Egx(i1,j1,k1) + ex(i)/fr1**2
        Egy(i1,j1,k1)=Egy(i1,j1,k1) + ey(i)/fr1**2
        Egz(i1,j1,k1)=Egz(i1,j1,k1) + ez(i)/fr1**2
        ! potential
        phig(i1,j1,k1)=phig(i1,j1,k1) + phi(i)/fr1
    end do

    ! normalise averaged quantities
    nelecs = SUM(g_ele(1:ngx,1:ngy,1:ngz))
    nions = SUM(g_ion(1:ngx,1:ngy,1:ngz))
    write(10,*) 'density integrals: ',nelecs, nions
    g_ele = max(1.,g_ele)
    g_ion = max(1.,g_ion)
    g_pro = max(1.,g_pro)

    Te = .511*Te/g_ele   ! Temperature in MeV (KE per particle)
    Ti = .511*mass_ratio*Ti/g_ion   ! Ion K.E. in MeV
    Tp = .511*1836*Tp/g_pro   ! Proton K.E. in MeV

    Egx = Egx/(g_ion)    ! Normalise fields/potential
    Egy = Egy/(g_ion)
    Egz = Egz/(g_ion)
    phig = phig/(g_ion)



    ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints

    ! output data file
    csnap = cdump//'/snap_defs.gmt'
    open(80,file=csnap)
    write(*,'(a)') 'Writing out snapshot defs to ',csnap
    write(80,'(a8,a3,f12.3,a2)') 'TLASER','"',tlaser,' "'  ! run time since laser switched on (timestamp)
    write(80,'(a8,a3,f12.3,a2)') 'TRUN','"',trun,' "'  ! total simulation time 



!  2D PLOTS: x-y plane

    nout = 17
    cfout(1) = cdump//'/xy_slice_eden'
    cfout(2) = cdump//'/xy_slice_iden'
    cfout(3) = cdump//'/xy_slice_jxe'
    cfout(4) = cdump//'/xy_slice_jye'
    cfout(5) = cdump//'/xy_slice_jxi'
    cfout(6) = cdump//'/xy_slice_jyi'
    cfout(7) = cdump//'/xy_slice_te'
    cfout(8) = cdump//'/xy_slice_ti'
    cfout(9) = cdump//'/xy_slice_ex'
    cfout(10) = cdump//'/xy_slice_ey'
    cfout(11) = cdump//'/xy_slice_phi'
    cfout(12) = cdump//'/xy_slice_nhot'
    cfout(13) = cdump//'/xy_slice_jxhot'
    cfout(14) = cdump//'/xy_slice_jyhot'
    cfout(15) = cdump//'/xy_slice_jzhot'
    cfout(16) = cdump//'/xy_slice_pden'
    cfout(17) = cdump//'/xy_slice_tp'

    do i=0,nout-1
        write(*,'(2a)') 'Writing slice ',cfout(i+1)
        open(20+i,file=cfout(i+1))
    end do



    kslice = (zslice(1)-zmin)*rdz
    !  sigma = 40.
    kmax = sigma*rdz/8    ! average over 1/8 laser spot\
    ! kmax=2
    nk = 2*kmax+1
    do j=1,ngy
        do i=1,ngx
            write(20,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(abs(rho_ele(i,j,kslice-kmax:kslice+kmax)))/nk
            !        logrho = log10(max(0.001,SUM(rho_ion(i,j,kslice-kmax:kslice+kmax))/nk))

            !        write(21,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, logrho
            write(21,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, SUM(rho_ion(i,j,kslice-kmax:kslice+kmax))/nk
            write(22,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(jx_ele(i,j,kslice-1:kslice+1))/3.
            write(23,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(jy_ele(i,j,kslice-1:kslice+1))/3.
            write(24,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(jx_ion(i,j,kslice-1:kslice+1))/3.
            write(25,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(jy_ion(i,j,kslice-1:kslice+1))/3.
            write(26,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(te(i,j,kslice-kmax:kslice+kmax))/nk
            write(27,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(ti(i,j,kslice-kmax:kslice+kmax))/nk
            write(28,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(Egx(i,j,kslice-1:kslice+1))/3.
            write(29,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(Egy(i,j,kslice-1:kslice+1))/3.
            write(30,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, phig(i,j,kslice)
            write(31,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(abs(rho_hot(i,j,kslice-kmax:kslice+kmax)))/nk
            write(32,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(jx_hot(i,j,kslice-1:kslice+1))/3.
            write(33,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(jy_hot(i,j,kslice-1:kslice+1))/3.
            write(34,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, jz_hot(i,j,kslice)
            write(35,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, SUM(rho_pro(i,j,kslice-kmax:kslice+kmax))/nk
            write(36,'(2f13.4,f16.6)') i*dx+xmin, j*dy+ymin, sum(tp(i,j,kslice-kmax:kslice+kmax))/nk

        end do
    end do

    do i=0,nout-1
        close(20+i)
    end do


!  2D PLOTS: x-z plane


    ! output data file
    nout = 15
    cfout(1) = cdump//'/xz_slice_eden'
    cfout(2) = cdump//'/xz_slice_iden'
    cfout(3) = cdump//'/xz_slice_jxe'
    cfout(4) = cdump//'/xz_slice_jze'
    cfout(5) = cdump//'/xz_slice_jxi'
    cfout(6) = cdump//'/xz_slice_jzi'
    cfout(7) = cdump//'/xz_slice_te'
    cfout(8) = cdump//'/xz_slice_ti'
    cfout(9) = cdump//'/xz_slice_ex'
    cfout(10) = cdump//'/xz_slice_ez'
    cfout(11) = cdump//'/xz_slice_phi'
    cfout(12) = cdump//'/xz_slice_nhot'
    cfout(13) = cdump//'/xz_slice_jxhot'
    cfout(14) = cdump//'/xz_slice_jyhot'
    cfout(15) = cdump//'/xz_slice_jzhot'

    do i=0,nout-1
        write(*,'(2a)') 'Writing slice ',cfout(i+1)
        open(20+i,file=cfout(i+1))
    end do



    jslice = (yslice(1)-ymin)*rdy
    do k=1,ngz
        do i=1,ngx
            write(20,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, abs(rho_ele(i,jslice,k))
            write(21,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, rho_ion(i,jslice,k)
            write(22,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jx_ele(i,jslice,k)
            write(23,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jz_ele(i,jslice,k)
            write(24,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jx_ion(i,jslice,k)
            write(25,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jz_ion(i,jslice,k)
            write(26,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, te(i,jslice,k)
            write(27,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, ti(i,jslice,k)
            write(28,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, Egx(i,jslice,k)
            write(29,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, Egz(i,jslice,k)
            write(30,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, phig(i,jslice,k)
            write(31,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, abs(rho_hot(i,jslice,k))
            write(32,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jx_hot(i,jslice,k)
            write(33,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jy_hot(i,jslice,k)
            write(34,'(2f13.4,f16.6)') i*dx+xmin, k*dz+zmin, jz_hot(i,jslice,k)
        end do
    end do

    do i=0,nout-1
        close(20+i)
    end do

!  2D PLOTS: y-z plane


    ! output data file
    nout = 15
    cfout(1) = cdump//'/yz_slice_eden'
    cfout(2) = cdump//'/yz_slice_iden'
    cfout(3) = cdump//'/yz_slice_jxe'
    cfout(4) = cdump//'/yz_slice_jze'
    cfout(5) = cdump//'/yz_slice_jxi'
    cfout(6) = cdump//'/yz_slice_jzi'
    cfout(7) = cdump//'/yz_slice_te'
    cfout(8) = cdump//'/yz_slice_ti'
    cfout(9) = cdump//'/yz_slice_ex'
    cfout(10) = cdump//'/yz_slice_ez'
    cfout(11) = cdump//'/yz_slice_phi'
    cfout(12) = cdump//'/yz_slice_nhot'
    cfout(13) = cdump//'/yz_slice_jxhot'
    cfout(14) = cdump//'/yz_slice_jyhot'
    cfout(15) = cdump//'/yz_slice_jzhot'

    do i=0,nout-1
        write(*,'(2a)') 'Writing slice ',cfout(i+1)
        open(20+i,file=cfout(i+1))
    end do



    islice = (xslice(1)-xmin)*rdx
    do k=1,ngz
        do j=1,ngy
            write(20,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, abs(rho_ele(islice,j,k))
            write(21,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, rho_ion(islice,j,k)
            write(22,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jx_ele(islice,j,k)
            write(23,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jz_ele(islice,j,k)
            write(24,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jx_ion(islice,j,k)
            write(25,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jz_ion(islice,j,k)
            write(26,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, te(islice,j,k)
            write(27,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, ti(islice,j,k)
            write(28,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, Egx(islice,j,k)
            write(29,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, Egz(islice,j,k)
            write(30,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, phig(islice,j,k)
            write(31,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, abs(rho_hot(islice,j,k))
            write(32,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jx_hot(islice,j,k)
            write(33,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jy_hot(islice,j,k)
            write(34,'(2f13.4,f16.6)') j*dy+ymin, k*dz+zmin, jz_hot(islice,j,k)
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
            avenhot(1:ngx) = avenhot(1:ngx) + abs(rho_hot(1:ngx,j,k)/norm)
            aveni(1:ngx) = aveni(1:ngx) + rho_ion(1:ngx,j,k)/norm
            avejxi(1:ngx) = avejxi(1:ngx) + jx_ion(1:ngx,j,k)/norm
            avejxe(1:ngx) = avejxe(1:ngx) + jx_ele(1:ngx,j,k)/norm
            aveex(1:ngx) = aveex(1:ngx) + egx(1:ngx,j,k)/norm
            avete(1:ngx) = avete(1:ngx) + 1.e3*te(1:ngx,j,k)/norm
            aveti(1:ngx) = aveti(1:ngx) + 1.e3*ti(1:ngx,j,k)/norm
            avephi(1:ngx) = avephi(1:ngx) + phig(1:ngx,j,k)/norm
        end do
    end do
    write(20,'(10(1pe16.4))') (i*dx+xmin, abs(avene(i)), aveni(i), avejxe(i), avejxi(i), aveex(i), &
        max(avete(i),1.e-5), max(aveti(i),1.e-5), avephi(i), avenhot(i), i=1,ngx )
    close(20)

    cfout(1) = 'fields_pp/mesh.xyz'
    cfout(2) = 'fields_pp/'//cdump//'.xyz'


    !  write(*,'(2a)') 'Writing 3D dump ',cfout(2)
    !  open(20,file=cfout(1))
    !  open(21,file=cfout(2))

    ! write AVS field record
    !  open(60,file='fields3d.fld',position='append')
    !  write(60,'(3a)') 'time file='//cfout(2)//' filetype=ascii'
    !  write(60,'(3a)') 'variable 1 file='//cfout(2)//' filetype=ascii offset=0 stride=7'
    !  write(60,'(3a)') 'variable 2 file='//cfout(2)//' filetype=ascii offset=1 stride=7'
    !  write(60,'(3a)') 'variable 3 file='//cfout(2)//' filetype=ascii offset=2 stride=7'
    !  write(60,'(3a)') 'variable 4 file='//cfout(2)//' filetype=ascii offset=3 stride=7'
    !  write(60,'(3a)') 'variable 5 file='//cfout(2)//' filetype=ascii offset=4 stride=7'
    !  write(60,'(3a)') 'variable 6 file='//cfout(2)//' filetype=ascii offset=5 stride=7'
    !  write(60,'(3a)') 'variable 7 file='//cfout(2)//' filetype=ascii offset=6 stride=7'
    !  write(60,'(a3)') 'EOT'
1   close(60)
    !  iskip3d = 2
    !  write(20,'(3i6)') ngx/iskip3d,ngy/iskip3d,ngz/iskip3d

    !  do k=1,ngz,iskip3d
    !     do j=1,ngy,iskip3d
    !       do i=1,ngx,iskip3d
    !           write(20,'(3f7.2)') i*dx+xmin, j*dy+ymin, k*dz+zmin
    !           write(21,'(7(1pe10.2))') & 
    !                abs(SUM(rho_ele(i:i+iskip3d-1,j,k))+SUM(rho_hot(i:i+iskip3d-1,j,k)))/iskip3d, &
    !                SUM(rho_ion(i:i+iskip3d-1,j,k))/iskip3d, &
    !                SUM(jx_ion(i:i+iskip3d-1,j,k))/iskip3d, &
    !                SUM(jy_ion(i:i+iskip3d-1,j,k))/iskip3d, &
    !                SUM(jz_ion(i:i+iskip3d-1,j,k))/iskip3d, & 
    !                SUM(Te(i:i+iskip3d-1,j,k))/iskip3d, & 
    !                SUM(Ti(i:i+iskip3d-1,j,k))/iskip3d 
    !        end do
    !     end do
    !  end do
    !  close(20)
    !  close(21)




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
    write(50,'(a8,a3,4(f12.3,a1))') 'UXPACONV','"',xmin*cowp_micron,'/',xmax*cowp_micron, &
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
    dy = yl/ngy
    du =  uimevmax/ngux  ! Bin size for energy spectra
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
    write(50,'(a8,a3,4(f12.5,a1))') 'UIPACONV','"',xmin*cowp_micron,'/',xmax*cowp_micron,&
        '/',uximin,'/',uximax,'"'  ! plot area
    write(50,'(a8,a3,4(f12.5,a1))') 'UXIDR','"',xmin+dx,'/',xmax,'/',uximin+dux,'/',uximax,'"'  ! data region
    write(50,'(a8,a3,2(f12.5,a1))') 'UXIMESH','"',dx,'/',dux,'"'   ! mesh size 
    write(50,'(a8,a3,2(f12.5,a1))') 'UXIAXES','"',xtick,'/',uxitick,'"'   ! tick intervals
    write(50,'(a8,a3,4(f12.5,a1))') 'UZIPA','"',zmin,'/',zmax,'/',uzimin,'/',uzimax,'"'  ! plot area
    write(50,'(a8,a3,4(f12.5,a1))') 'UZIPAC','"',xmin*cowp_micron,'/',xmax*cowp_micron, &
        '/',uzimin,'/',uzimax,'"'  ! plot area
    write(50,'(a8,a3,4(f12.5,a1))') 'UZIDR','"',zmin+dz,'/',zmax,'/',uzimin+duz,'/',uzimax,'"'  ! data region
    write(50,'(a8,a3,2(f12.5,a1))') 'UZIMESH','"',dz,'/',duz,'"'   ! mesh size 
    write(50,'(a8,a3,2(f12.5,a1))') 'UZIAXES','"',ztick,'/',uzitick,'"'   ! tick intervals
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
    fpxx_pro(0:ngx+1,0:ngux+1) = 0.
    fpyy_pro(0:ngy+1,0:nguy+1) = 0.
    fpzx_ion(0:ngx+1,0:nguz+1) = 0.
    fang_ion(0:nalpha+1,0:nalpha+1) = 0.
    uang_ion(0:nalpha+1,0:nalpha+1) = 0.
    fang_ion1(0:nalpha+1,0:nalpha+1) = 0.
    fang_ion2(0:nalpha+1,0:nalpha+1) = 0.
    fang_back(0:nalpha+1,0:nalpha+1) = 0.

    fhoti_for = 0.
    fhoti_back = 0.
    fhotp = 0.
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
        if (q(i)>0 ) then
            xa=(ux(i)-uximin)*rdx
            ya=(uy(i)-uyimin)*rdy
            za=(uz(i)-uzimin)*rdz
            xxa = (x(i)-xmin)/dx
            yya = (y(i)-ymin)/dy

            !  indices
            i1=xa+1
            j1=ya+1
            k1=za+1
            ii1 = xxa+1  ! space
            jj1 = yya+1  ! space

            i1 = min(max(0,i1),ngux+1)
            j1 = min(max(0,j1),nguy+1)
            k1 = min(max(0,k1),nguz+1)
            ii1 = min(max(0,ii1),ngx+1)
            jj1 = min(max(0,jj1),ngy+1)

            !  count charges at nearest point in phase space

            if (m(i)/q(i)>1850.) then
                ! heavy ions

                fpxx_ion(ii1,i1)=fpxx_ion(ii1,i1) + 1
                fpzx_ion(ii1,j1)=fpzx_ion(ii1,j1) + 1
                fpyy_ion(jj1,j1)=fpyy_ion(jj1,j1) + 1

                fpxpy_ion(i1,j1)=fpxpy_ion(i1,j1) + 1

                fpxpz_ion(i1,k1)=fpxpz_ion(i1,k1) + 1

                gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
                umev  = .511*mass_ratio*(gamma-1.0)

                iu = umev/du+1  ! Energy bin
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

            else
                ! forward protons
                fpxx_pro(ii1,i1)=fpxx_pro(ii1,i1) + 1
		! radial protons
                if (x(i) > rear_edge) fpyy_pro(jj1,j1)=fpyy_pro(jj1,j1) + 1
                gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
                umev  = .511*1836.*(gamma-1.0)
                iu = umev/du+1  ! Energy bin
                iu = min(max(0,iu),ngux+1)  ! limits
                if (ux(i)>0) fhotp(iu) = fhotp(iu) + 1
            endif

        endif


    end do

    uang_ion = uang_ion/max(1.,fang_ion)   ! ave. KE angular distn 
    write (10,'(//a20,2a10/a40/3(a20,2f10.3/))') &
        'Ion energies (MeV):','Average','Maximum',&
        '============================================================',&
        'Rear-side:', utnsa/max(ntnsa,1), ur_max, &
        'Front-side:', ufront/max(nfront,1), uf_max,&
        'Blowoff:',ublow/max(nblow,1), ub_max

    ng = (ngux+2)*(nguy+2)                         ! total # gridpoints

    ! output data files
    nout = 13
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

    cfout(11) = cdump//'/fpxx_pro'
    cfout(12) = cdump//'/fpyy_pro'
    cfout(13) = cdump//'/fhot_pro'


    do i=0,nout-1
        write(*,'(2a)') 'Writing slice ',cfout(i+1)
        open(20+i,file=cfout(i+1))
    end do

    fpxpy_ion = max(fpxpy_ion,0.9)
    fpxpz_ion = max(fpxpz_ion,0.9)
    fpxx_ion = max(fpxx_ion,0.9)
    fpzx_ion = max(fpzx_ion,0.9)
    fpxx_pro = max(fpxx_pro,0.9)
    fpyy_pro = max(fpyy_pro,0.9)
    fpyy_ion = max(fpyy_ion,0.9)

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
            write(30,'(2f13.4,f16.6)') i*dx+xmin, j*dux+uximin, log(fpxx_pro(i,j))
        end do
    end do

    do j=1,nguy
        do i=1,ngy
            write(23,'(2f13.4,f16.6)') i*dy+ymin, j*duy+uyimin, log(fpyy_ion(i,j))
            write(31,'(2f13.4,f16.6)') i*dy+ymin, j*duy+uyimin, log(fpyy_pro(i,j))
        end do
    end do


    ! energy spectra

    fhoti_for = max(fhoti_for,0.1)
    fhoti_back = max(fhoti_back,0.1)
    fhotp = max(fhotp,0.1)

    write(24,'((f13.4,f16.6))') (i*du, fhoti_for(i), i=1,ngux)
    write(25,'((f13.4,f16.6))') (i*du, fhoti_back(i), i=1,ngux)
    write(32,'((f13.4,f16.6))') (i*du, fhotp(i), i=1,ngux)

    ! far-field y-z angular distribution
    do j=1,nalpha
        do i=1,nalpha
            write(26,'(2f13.4,f16.6)') i*dalpha+aimin, j*dalpha+aimin, uang_ion(i,j)
            write(27,'(2f13.4,f16.6)') i*dalpha+aimin, j*dalpha+aimin, fang_ion1(i,j)
            write(28,'(2f13.4,f16.6)') i*dalpha+aimin, j*dalpha+aimin, fang_ion2(i,j)
            write(29,'(2f13.4,f16.6)') i*dalpha+aimin, j*dalpha+aimin, fang_back(i,j)
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



!  Index sort for 4-byte integer list

  subroutine indsort(iarr,list,n)
    implicit none

    integer, intent(in) :: n

    integer, dimension(n), intent(in) :: iarr
    integer, dimension(n), intent(inout) :: list
    integer :: i, indxt, ir, l, j
    integer :: q

    list(1:n) =  (/ (i,i=1,n) /)  

    if(n==1) return
    l= n/2 + 1
    ir = n

    do
       if (l>1) then
          l=l-1
          indxt = list(l)
          q = iarr(indxt)
       else
          indxt = list(ir)
          q = iarr(indxt)
          list(ir) = list(1)
          ir = ir - 1
          if (ir == 1) then
             list(1) = indxt
             return
          endif
       endif

       i = l
       j = l+l
       do while (j <= ir)
          if (j < ir) then
             if (iarr(list(j)) < iarr(list(j+1)) ) j=j+1
          endif
          if (q < iarr(list(j)) ) then
             list(i) = list(j)
             i=j
             j=j+j
          else
             j = ir+1
          endif
       end do
       list(i) = indxt
    end do

  end subroutine indsort



  subroutine swap_ab(p,q)
    integer*8 :: p,q, dum
    dum = p
    p=q
    q = dum
  end subroutine swap_ab


