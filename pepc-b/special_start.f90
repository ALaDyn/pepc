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
  integer :: i,iseed=-17
  real :: rs,v0,gamma0,vt,xt,yt,zt,thetvel,phivel

  r_beam=sigma/2.
  gamma0=sqrt(1+vosc**2/2.)

  config: select case(iconf)
    case(1)  

! electron disc at x=0 with 2pi phase spread

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

    case(2) 

  ! electron disc at laser focus with 2pi phase spread
 
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

    case(3)  ! random placement in box with vte in x,y plane for B test 
	do i=1,npp
          if (i<=nep) then
	    vt=vte
	  else
	    vt = vti
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
  q(1:nep) = qe           ! plasma electrons
  m(1:nep) = mass_e            ! electron mass
  q(nep + 1:npp)       = qi        ! plasma ions (need Z* here)
  m(nep + 1:npp)       = mass_i    ! ion mass
  pelabel(1:nep)       = me * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:npp) = ne + me * nip + (/(i, i = 1, nip)/) ! Ion labels
  pepid(1:npp) = me                ! processor ID
  work(1:npp) = 1.
  call push_full3v(1,npp,dt/2.)

end subroutine special_start



