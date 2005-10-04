!  ==============================================
!
!     Full 3v particle pusher (Boris rotation scheme)
!
!  TM fields only (s-pol):  (Ex,0,Ez)  (0,By,0)
!  TE fields  (0,Ey,0) (Bx,0,Bz)
!  ==============================================


subroutine push_full3v(p_start,p_finish,dts)
  use physvars
  use treevars
  implicit none
  integer, intent(in) :: p_start, p_finish
  real, intent(in) :: dts
  integer :: p
  real :: beta, gam2,gam1,tt, sx, sy, sz, tz, ty, tx
  real :: uxd, uyd, uzd, uyp, uzp, uxp, uxm, uym, uzm
  real :: exi, eyi, ezi, bxi, byi, bzi
  real :: phipon, ez_em, bx_em, by_em, az_em, xd, yd, zd


  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  do p = p_start, p_finish
     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

     xd = x(p)-focus(1)
     yd = y(p)-focus(2)
     zd = z(p)-focus(3)


     !  Sum internal and external fields 

!  TODO: need to include inductive dA/dt term in E-field (have to store prev. timestep of A(p))
 
     exi = ex(p)
     eyi = ey(p)
     ezi = ez(p)
     bxi = bx(p)
     byi = by(p)
     bzi = bz(p)



     !   first half-accn
     uxm = ux(p) + beta*exi
     uym = uy(p) + beta*eyi
     uzm = uz(p) + beta*ezi

     !   rotation
     gam1=dsqrt(1.d0 + uxm**2 + uym**2 + uzm**2)

     tx = beta*bxi/gam1
     ty = beta*byi/gam1
     tz = beta*bzi/gam1
     tt = 1.0 + tx**2 + ty**2 + tz**2

     sx = 2.0*tx/tt
     sy = 2.0*ty/tt
     sz = 2.0*tz/tt

     uxd = uxm + uym*tz - uzm*ty
     uyd = uym + uzm*tx - uxm*tz
     uzd = uzm + uxm*ty - uym*tx

     uxp = uxm + uyd*sz - uzd*sy
     uyp = uym + uzd*sx - uxd*sz
     uzp = uzm + uxd*sy - uyd*sx

     !   second half-accn
     ux(p) = uxp + beta*exi
     uy(p) = uyp + beta*eyi
     uz(p) = uzp + beta*ezi

  end do


end subroutine push_full3v

