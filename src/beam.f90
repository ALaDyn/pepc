! ==============================================
!
!                BEAM
!
!  Sets up particle beam
!
! ==============================================

subroutine beam
 
  use treevars
  use utils
  implicit none
  integer :: i, p, iseed1, iseed2
  real :: Volb, dpx, yt, zt
  real :: vx_beam, vy_beam, vz_beam, ct, cp, st, sp




  !  evaluate system constants  from inputs:
  !  beam cylinder volume:  r_beam is radius

  nb_pe = np_beam/num_pe  ! # beam particles to load per PE

  Volb = pi*r_beam**2*x_beam
  qeb = Volb*rho_beam/np_beam    ! charge
  dpx=x_beam/nb_pe           ! x-axis spacing
  iseed2 = -131 - np_beam -3*me ! seed
  iseed1 = -333 - np_beam -3*me

  ! proton beam: initialised along x-axis

  ! beam initialised along x-axis and rotated by theta, phi

  ct = cos(theta_beam)
  st = sin(theta_beam)
  cp = cos(phi_beam)
  sp = sin(phi_beam)
  vz_beam = u_beam*st
  vx_beam = u_beam*ct*cp
  vy_beam = u_beam*ct*sp

  i = 0

  do while (i < nb_pe)
     yt = r_beam*(2*rano(iseed2)-1.)          
     zt = r_beam*(2*rano(iseed1)-1.)  
     if (yt**2 + zt**2 <= r_beam**2 ) then
        i = i+1
        p = npp+i  ! put them after plasma ions
        x(p)= start_beam + dpx*i + dpx/num_pe*me
        y(p) = yt + yl/1.9
        z(p) = zt + zl/1.8
        q(p) = qeb
        m(p) = -qeb*mass_beam
        ux(p)=vx_beam
        uy(p)=vy_beam
        uz(p)=vz_beam
        pepid(p) = me                ! processor ID
        pelabel(p) = npart+me*nb_pe+i  ! labels
     endif
  end do

  npart = npart + np_beam  ! Augment particle numbers
  npp = npp + nb_pe

end subroutine beam









