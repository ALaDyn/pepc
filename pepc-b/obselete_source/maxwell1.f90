!  ============
!
!   MAXWELL1 
!   initialises 1D Maxwellian velocity distribution
!
!  ============

subroutine maxwell1(u,nmax,i1,n,vt)
  implicit none
  integer, intent(in) ::  nmax, i1, n
  real, intent(in) :: vt
  real*8, intent(out) :: u(nmax)
  real, parameter :: pi=3.141592654

  integer :: ip1, ip2, i, cntr, nv
  real*8 :: f0 ,df, v, dv, vmax, finf, vip, deltv

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
     df = exp( max(-30.d0,-0.5*v**2) )*dv/finf*n/2.
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
