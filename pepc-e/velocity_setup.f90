!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates functions for setting up particle velocities with diefferent models
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module velocity_setup
      use physvars
      implicit none
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public maxwell1
      public scramble_v
      public cold_start
      public rano

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   COLD_START
        !>   initialises velocity array slice to 0.
        !>   @param u array of velocities to be initialized
        !>   @param nmax total number of entries in u
        !>   @param i1 minimal index in u to be used
        !>   @param n maximum index in u to be used
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine cold_start(ux,uy,uz,nmax,i1,n)

		  implicit none
		  integer, intent(in) :: i1,n,nmax
          real*8 :: ux(nmax), uy(nmax), uz(nmax)

		  ux(i1:i1+n) = 0.
		  uy(i1:i1+n) = 0.
		  uz(i1:i1+n) = 0.

		end subroutine cold_start


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   MAXWELL1
        !>   initialises 1D Maxwellian velocity distribution
        !>   @param u array of velocities to be initialized
        !>   @param nmax total number of entries in u
        !>   @param i1 minimal index in u to be used
        !>   @param n maximum index in u to be used
        !>   @param vt desired average velocity of particles
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine maxwell1(u,nmax,i1,n,vt)
		  implicit none
		  integer, intent(in) ::  nmax, i1, n
		  real, intent(in) :: vt
		  real*8 :: u(nmax)
		  real, parameter :: pi=3.141592654

		  integer :: ip1, ip2, i, cntr, nv
		  real :: f0 ,df, v, dv, vmax, finf, vip, deltv

		  if (n.eq.0) return
		  nv = 30*n
		  vmax = 4.0
		  dv = vmax/nv
		  f0 = 0.5
		  finf = sqrt(pi/2.0)
		  cntr = 1
		  ip1 = n/2+i1-1
		  ip2 = n/2+i1

		  do  i=1,nv
		     v = vmax - (i-0.5)*dv
		     deltv = vmax/10.
		     df = exp( max(-30.0,-0.5*v**2) )*dv/finf*n/2.
		     f0 = f0 + df       ! integrate dist. fn.
		     if(f0.ge.cntr) then
		        vip = vt*(v-dv*(f0-cntr)/df)
		        u(ip1) = vip
		        u(ip2) = -vip
		        cntr = cntr + 1
		        ip1 = ip1-1
		        ip2 = ip2+1
		     endif
		  end do

		  !  odd one out
		  if (mod(n,2).eq.1) then
		     u(n)=0.
		  endif
		end subroutine maxwell1


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   SCRAMBLE_V
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine scramble_v(ux,uy,uz,nmax,i1,n)

		  implicit none

		  integer :: dum1, dum2, dum3
		  real*8 :: uxt, uyt, uzt
          real*8 :: ux(nmax), uy(nmax), uz(nmax)
		  integer :: i, j, k, kk, p, i1, n, n1,nmax

		  dum1 = -71 - 10*my_rank
		  dum2 = -113301 - 10*my_rank
		  dum3 = -8651 - 10*my_rank
		  !  exclude odd one out
		  if (mod(n,2).ne.0) then
		     n1=n-1
		  else
		    n1 = n
		  endif
		 write(15,*) 'scrambling ',n1,' particles'
		  !  scramble indices to remove correlation between ux,uy,uz
		  do i=1,n1
		     p=i+i1-1
		     j=int(n1*rano(dum1)+i1)
		     k=int(n1*rano(dum2)+i1)
		     kk=int(n1*rano(dum3)+i1)
		     uxt=ux(p)
		     uyt=uy(p)
		     uzt=uz(p)
		     ux(p)=ux(kk)
		     uy(p)=uy(j)
		     uz(p)=uz(k)
		     ux(kk)=uxt
		     uy(j)=uyt
		     uz(k)=uzt
		  end do

		end subroutine scramble_v


		  ! Random number scrambler
		  ! =======================
		  !
		  !  called with:
		  !               x=rano(iseed)
		  !
		  !  returns real number in the interval (0.0 - 1.0)
		  !  set iseed = -ve integer before using call
		  !
		  !  Example of calling routine:
		  !
		  !      subroutine coords
		  !      include 'common.h'
		  !
		  !      iseed1 = -11
		  !      iseed2 = -7
		  !
		  !
		  !      do i=1,n
		  ! x(i)=xlen*rano(iseed1)
		  !       y(i)=ylen*rano(iseed2)
		  !      end do
		  !
		  !
		  !      end


		  ! - routine taken from Numerical Recipies, p195
		  !
		  real*8 function rano(idum)
		    implicit none
		    integer :: idum
		    real*8, save :: dseed, dum
		    real*8, save :: v(97), y
		    integer, save :: iff, icall, i, j
		    data iff,icall/0,0/
		    if (idum.lt.0.or.iff.eq.0) then
		       !  first call generates new sequence
		       iff = 1
		       dseed=abs(idum)*1.0
		       idum=1
		       do  j=1,97
		          dum=genran(dseed)
		       end do
		       do j=1,97
		          v(j)=genran(dseed)
		       end do
		       y=genran(dseed)
		    endif

		    !  next index - make sure we don`t overstep array bounds if
		    !  generator returns a 0.0 or 1.0

		    j=max(mod(1+int(97.*y),98),1)
		    if(j.gt.97.or.j.lt.1) then
		       write (6,*) 'Call: ',icall
		       write (6,*) 'idum = ',idum,'j = ',j,' y= ',y
		       write (6,*) 'Random No. generator not initialised properly'
		       write (6,*) 'dummy =',dum,' dseed=',dseed
		       write (6,100) (i,v(i),i=1,97)
		100    format (i4,f10.6)
		       stop
		    endif
		    !  get next variate and generate new one to fill gap

		    y=v(j)
		    rano=y
		    v(j)=genran(dseed)
		    icall = icall + 1

		    return
		  end function rano


		  real*8 function genran (dseed)
		    real*8 ::  dseed
		    real*8 ::  d2p31m,d2p31
		    !                                  d2p31m=(2**31) - 1
		    !                                  d2p31 =(2**31)(or an adjusted value)
		    data               d2p31m/2147483647.0/
		    data               d2p31 /2147483648.0/

		    dseed = mod(16807.0*dseed,d2p31m)
		    genran = dseed / d2p31

		  end function genran



end module velocity_setup
