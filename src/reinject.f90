

!  ========================================
!
!     Thermal particle reinjection
!    - conserves  background temp. 
!    - pairwise injection for uy,uz to give zero transverse current
!
!  ========================================

 subroutine reinject(vt,idir,uxi,uyi,uzi)

 use utils

 implicit none
 real, parameter :: pi=3.14159265
 integer, parameter :: nsamp=30000
 integer :: idir
 real, save :: usamp(nsamp+1)
 real  :: uxi, uyi
 real, save :: theta
 real :: df,f0,A,u,du,g,uzi,g0,rs
 real :: ute,  pio2, umax
 real, intent(in) :: vt
 integer, save :: idum=-13
 integer, save :: isamp=1
 integer, save :: i2
 integer :: i,k, nsteps

!  conserve flux in x dirn: maintain const. Te

!  Tabulate momenta on first call
    ute = vt**2  ! Temperature normalised to mc^2
    pio2 = pi/2
   if (isamp.eq.1) then

        nsteps=100000
!  dimensionless electron/ion temp
        umax=6.*sqrt(ute + 4*ute**2)
        du=umax/nsteps
!  distn norm factor
        A = 1.0001*nsamp/2./pi/ute/(1+ute)
        f0=0.
        k=1

!  load loop
        do i=1,nsteps

!  cumulative integral of rel. distn function f(u)
	  u = i*du
          g = sqrt(u**2 + 1.0)
	  df = 2*pi*A*du*u*exp((1.-g)/ute)

	  f0 = f0+df

	  if(f0.ge.k) then
!  store u corresponding to integer values of f0
            usamp(k) = u
	    k=k+1
	  endif
        end do
!        write(40,*) (usamp(k),k=1,nsamp)
      endif

!  Flux in x-dirn: invert Int (vx f(ux) dux) directly
     
      rs=rano(idum)
      g0 = dmax1(1.d0,1.d0-ute*alog(rs))
      uzi = idir*sqrt(g0**2-1.d0)

!  pick random  u from sample:  
!  usamp contains integrated momentum flux Int(u f(u) du)
      if (mod(isamp,2).eq.1) then
       i2=min(1.*nsamp,nsamp*rano(idum)+1)
!  components for uy,uz
       theta=2*pi*rano(idum)
       uxi = usamp(i2)*cos(theta)
       uyi = usamp(i2)*sin(theta)
!  write(40,*) 'isamp=',isamp,' i2=',i2,' idum=',idum,'ux,uy=',uxi, uyi
      else
!  use -ve of previous sample
       uxi = -usamp(i2)*cos(theta)
       uyi = -usamp(i2)*sin(theta)
!  write(40,*) 'isamp=',isamp,' i2=',i2,' idum=',idum, 'ux,uy=',uxi, uyi
      endif
  
!  increment sample index
      isamp = isamp + 1
      if (isamp.gt.nsamp) isamp = 1

  end subroutine reinject
