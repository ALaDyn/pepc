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
  real :: gamma0,yt,zt

  r_beam=sigma/2.
  gamma0=sqrt(1+vosc**2/2.)

  config: select case(iconf)
  case(1)  ! electron disc at x=0 with 2pi phase spread

     do while (i < npp)
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
     x_crit = focus(1)
     ! electron disc at laser focus with 2pi phase spread
     do while (i < npp)
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

  end select config

  ex(1:nep) = 0.
  ey(1:nep) = 0.
  ez(1:nep) = 0.
  pot(1:nep) = 0.
  q(1:nep) = qe           ! plasma electrons
  m(1:nep) = mass_e            ! electron mass
  pepid(1:npp) = me                ! processor ID
  pelabel(1:nep) = me*nep + (/ (i,i=1,nep) /)     ! Electron labels
  work(1:npp) = 1.
  call push_full3v(1,npp,dt/2.)

end subroutine special_start



