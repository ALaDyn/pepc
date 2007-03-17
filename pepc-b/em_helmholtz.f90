
!>     ==================================
!>     
!>     Helmholtz solver for electromagnetic fields 
!>     
!>     Based on SPLIM project helmholtz.f90     
!>     ==================================


subroutine em_helmholtz(me,itime,n,dx,theta,a0,w0,rhoe,Az)

  implicit none
  integer, intent(in) :: n !< number of mesh points for 1D Helmholtz fields
  integer, intent(in) :: itime  !< current timestep 
  integer, intent(in) :: me  !< current rank 
  real, intent(in) :: theta  !< angle of incidence
  real, intent(in) ::  a0    !< laser amplitude vosc/c(t)
  real, intent(in) ::  w0    !< laser frequency (w/wp)
  real, intent(in) ::  dx    !< mesh spacing (c/wp)
  real*4, intent(in) :: rhoe(0:n+1)  !< cycle-averaged electron density
!  real*8, intent(out) :: Ezr(0:n+1), Byr(0:n+1), Azr(0:n+1), epond(0:n+1)
  complex, dimension(n) :: alpha,beta,gamma,y !< trisolver work arrays
  complex, dimension(0:n+1) :: Az,Ao	!< Vector potential
  complex, dimension(n) :: Ez, By, Bx !< Fields derived from Az
  real, dimension(n) :: eps  !< permittivity
  real :: rgam(0:n+1)  !<  relativistic gamma
  real :: err(0:n+1)   !<  error check
  complex :: yi,  aave, carg  !< complex args
  integer :: i,j,n_iter
  real :: pi, s2th, cth, errmax, pha !< constants
  real*8 :: lskin  !< skin depth
  real*8 :: ncrit  !< critical density

  real*8 :: phipon, epon_x, epon_y, epon_z !< pond. potential and fields

  real*8 :: xf, yf, zf, g0, nu_eff
  integer :: iplas, itav

  itav=1
  n_iter = 3
  ncrit = w0**2  ! critical density in norm. units
  nu_eff=0.2  ! Effective collision rate
  pi = asin(1.0)*2
  yi = (0.,1.)
  err=0.
  s2th=sin(pi/180*theta)**2
  cth = cos(pi/180*theta)


! Use gamma from previous step to speed up iteration
  ao=az
  rgam(0)=1.
  rgam(n+1) = 1.
  do i=1,n
     rgam(i) = sqrt(1 + 0.5*abs(ao(i))**2)
     az(i)=(0.,0.) 
  end do


  Az(0)=(0.,0.)
  Az(n+1)=(0.,0.)

  if (a0==0) return

  do j=1,n_iter
     do i=1,n
        !  coefficients as for s-pol light
        ! rhoe normalized to nc defined on own 1D grid along laser axis
        eps(i) = 1.-rhoe(i)/ncrit/rgam(i)/(1.+yi*nu_eff)
        y(i)=(0.,0.)
        alpha(i)=1
        beta(i)=-2 + w0**2*dx**2*(eps(i)-s2th)
        gamma(i)=1
     end do

     !  BCs - note additional factor w0=k0 for phase factors 
     y(1) = 2*yi*a0*sin(w0*dx*cth)
     carg = yi*w0*dx*cth
     beta(1) = beta(1) + cexp(carg)

     call trisolve(alpha,beta,gamma,y,Az(1:n),n-1,n)

     ! relativistic factor - average old and new potentials
     errmax=0.
     do i=1,n
        rgam(i) = sqrt(1 + 0.5*abs(az(i))**2) 
        err(i) = sqrt(abs(az(i)**2-ao(i)**2)/4/a0**2)
     end do

     ! BCs for next iterate (Laplacian of gamma)
     rgam(0) = 2*rgam(1) - rgam(2)
     rgam(n+1) = rgam(n)

     errmax = maxval(err(1:n))
     Ao = Az  ! Store previous iterate

  end do

iplas = n/2
if (me==0) write(*,'(i6,2f12.3)') itime,rhoe(iplas),abs(az(iplas))
if (itime .eq. itav .and. me==0) then
  write (*,'(a20,i2,a10,f12.5,a10,f12.3)') 'Iterate ',j,' error=',errmax,' amplitude',a0
  g0 = sqrt(1+a0**2/2)
  open (40,file='a_error.dat')
  write(40,'(a)') '! x, rho, eps, az/a0, gam/g0, err'
  write(40,'(6(1pe12.3))') (dx*i,rhoe(i),eps(i),abs(az(i))/a0,rgam(i)/g0,err(i),i=1,n)
  close(40)
  
endif

 
  ! Bcs 
  Az(0) = 2*Az(1) - Az(2)  ! gradient continuous
  Az(n+1) = Az(n) ! zero in solid


end subroutine em_helmholtz









