!  ==============================================
!
!     3v particle pusher
!
!  TM fields only (s-pol):  (Ex,0,Ez)  (0,By,0)
!  TE fields to follow (0,Ey,0) (Bx,0,Bz)
!  ==============================================


subroutine push_em(p_start,p_finish,dts)
  use physvars
  use treevars
  implicit none
  integer, intent(in) :: p_start, p_finish
  real, intent(in) :: dts
  integer :: p
  real :: beta, gam2,gam1,tt,sx, sy, sz, tz, ty, tx
  real :: uxd,uyp,uzp,uxp,uxm,uym,uzm,uzn1
  real :: exi, eyi, ezi, bxi, byi, bzi
  real :: phipon, ez_em, by_em, az_em, xd, yd, zd


  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  do p = p_start, p_finish
     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

     xd = x(p)-focus(1)
     yd = y(p)-focus(2)
     zd = z(p)-focus(3)

     ! evaluate external fields at particle positions

     if (beam_config.eq.4) then
! pond. standing wave on step-profile
        call empond(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,az_em,phipon)

     else if (beam_config.eq.6) then
! plane wave with Gaussian spot
        call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,az_em,phipon)
     endif

     !  Sum internal and external fields
     exi = ex(p)
     eyi = ey(p)
     ezi = ez(p)+ez_em
     bxi = 0.
     byi = by_em
     bzi = 0.

     ! transverse momentum from pz=az including thermal motion 
     !	uzi = -q(p)/m(p)*az_em + uz(p)


     !   first half-accn
     uxm = ux(p) + beta*exi
     uym = uy(p) + beta*eyi
     uzm = uz(p) + beta*ezi
     !     uzm = uz(p) - q(p)/m(p)*az_em
     !   rotation
     gam1=dsqrt(1.d0 + uxm**2 + uym**2 + uzm**2)
     ty = beta*byi/gam1
     tz = beta*bzi/gam1
     tt = 1.0 + ty**2+tz**2
     sy = 2.0*ty/tt
     sz = 2.0*tz/tt

     uxd = uxm + uym*tz - uzm*ty
     uyp = uym - uxd*sz
     uzp = uzm + uxd*sy
     uxp = uxd + uyp*tz - uzp*ty

     !   second half-accn
     ux(p) = uxp + beta*exi
     uy(p) = uyp + beta*eyi
     uz(p) = uzp + beta*ezi

  end do


end subroutine push_em








