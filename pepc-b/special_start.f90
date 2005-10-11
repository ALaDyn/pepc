! ==============================================
!
!                SPECIAL_START
!
!  Initialise set of particles with zero velocity
!
! ==============================================

subroutine special_start(iconf)

  use physvars
  use treevars
  use utils
  implicit none

  integer, intent(in) :: iconf
  integer :: i,iseed
  real :: rs,v0,gamma0,vt,xt,yt,zt,thetvel,phivel
  real :: dpx, vx_beam, vy_beam, vz_beam, st, ct, sp, cp, volb
  integer :: npb, p

  iseed = -17-me

  config: select case(iconf)
  case(1)  

     ! electron disc at x=0 with 2pi phase spread
     r_beam=sigma/2.
     gamma0=sqrt(1+vosc**2/2.)

     do while (i < nep)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = 2*pi*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do
     q(1:nep) = qe           ! plasma electrons
     m(1:nep) = mass_e            ! electron mass

  case(2) 

     ! electron disc at laser focus with 2pi phase spread
     r_beam=sigma/2.
     gamma0=sqrt(1+vosc**2/2.)


     x_crit = focus(1)
     do while (i < nep)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = focus(1)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do
     q(1:nep) = qe           ! plasma electrons
     m(1:nep) = mass_e            ! electron mass


  case(3)  ! random placement in box with vte in x,y plane for B test 
     do i=1,npp
        if (i<=nep) then
           vt=vte
           q(i) = qe        ! plasma electrons
           m(i) = mass_e 
        else
           vt = vti
           q(i) = qi        ! plasma ions (need Z* here)
           m(i) = mass_i    ! ion mass

        endif
        xt =  xl * rano(iseed) 
        yt =  yl * rano(iseed)
        zt =  zl * rano(iseed)
        x(i) = xt
        y(i) = yt
        z(i) = zt
        rs = rano(iseed)
        v0 = vt*sqrt(-2.0*log(rs))
        thetvel = 2*pi*rano(iseed)
        phivel = pi*rano(iseed)
        ux(i) = v0*cos(thetvel)
        uy(i) = v0*sin(thetvel)
        uz(i) = v0*cos(phivel)
     end do

  case(4)

     !  Special config for FMM comparison
     open(30,file='fmm_config.data')
     do i=1,npp
        read(30,*) x(i), y(i), z(i), q(i)
        m(i) = 1.
        ux(i)=0.
        uy(i)=0.
        uz(i)=0.
     end do
     force_const=1.
     dt=0.
     eps=0.
     close(30)

  case(5)

     write (*,*) 'Preparing particle beam'
     ! beam initialised along x-axis and rotated by theta, phi
     npb = max(ni,ne)  ! Take whichever species switched on
     np_beam=0 ! Discard normal beam mode (plasma+beam)
     npp = max(nep,nip)

     ct = cos(theta_beam)
     st = sin(theta_beam)
     cp = cos(phi_beam)
     sp = sin(phi_beam)
     vz_beam = u_beam*st
     vx_beam = u_beam*ct*cp
     vy_beam = u_beam*ct*sp
     Volb = pi*r_beam**2*x_beam
     qeb = Volb*rho_beam/npb    ! charge
     dpx=x_beam/npb           ! x-axis spacing
     i=0
     do while (i < npp)
        yt = r_beam*(2*rano(iseed)-1.)
        zt = r_beam*(2*rano(iseed)-1.)

        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i)= start_beam + x_beam*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = vx_beam
           uy(i) = vy_beam
           uz(i) = vz_beam

        endif
     end do
     q(1:npp) = qeb           ! plasma electrons
     m(1:npp) = abs(qeb*mass_beam)  ! electron mass

!     write (*,'(a30/(5f13.5))') 'Beam config:', &
!          (x(i),y(i),z(i),q(i),m(i),i=1,npp)
  end select config

  ex(1:npp) = 0.
  ey(1:npp) = 0.
  ez(1:npp) = 0.
  bx(1:npp) = 0.
  by(1:npp) = 0.
  bz(1:npp) = 0.
  ax(1:npp) = 0.
  ay(1:npp) = 0.
  az(1:npp) = 0.
  axo(1:npp) = 0.
  ayo(1:npp) = 0.
  azo(1:npp) = 0.
  pot(1:npp) = 0.


  pelabel(1:nep)       = me * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:npp) = ne + me * nip + (/(i, i = 1, nip)/) ! Ion labels
  pepid(1:npp) = me                ! processor ID
  work(1:npp) = 1.
  call push_full3v(1,npp,dt/2.)

end subroutine special_start



