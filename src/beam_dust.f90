! ==============================================
!
!                BEAM_DUST
!
!  Sets up spherical dust particle
!
! ==============================================

subroutine beam_dust

  use physvars
  use treevars
  use utils
  implicit none
  integer :: i, p, iseed1, iseed2
  real :: Volb, dpx, xt, yt, zt
  real :: vx_beam, vy_beam, vz_beam, ct, cp, st, sp




  !  evaluate system constants  from inputs:
  !  beam cylinder volume:  r_beam is radius


  Volb = 4*pi/3.*r_beam**3
  qeb = Volb*rho_beam/np_beam    ! charge

  if (me==0) then
     iseed1 = -333 - np_beam -3*me

     ! proton beam: initialised along x-axis

     ! beam initialised along x-axis and rotated by theta, phi
     write (*,*) 'Setting up dust particle'
     write (*,'(a20,f12.5)') 'Beam charge ',qeb
     write (*,'(a20,f12.5)') 'Ion charge ',qi
     ct = cos(theta_beam)
     st = sin(theta_beam)
     cp = cos(phi_beam)
     sp = sin(phi_beam)
     vz_beam = u_beam*st
     vy_beam = u_beam*ct*cp
     vx_beam = u_beam*ct*sp


     i = 0

     ! Put dust particle on root PE
     do while ( i < np_beam)  
        yt = r_beam*(2*rano(iseed1)-1.)          
        zt = r_beam*(2*rano(iseed1)-1.)  
        xt = r_beam*(2*rano(iseed1)-1.)  
        if (xt**2 + yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           p = npp+i  ! put them after plasma ions
           x(p)= xl/2. + xt        ! Assume disc or slab geometry and inject from side
           y(p) = yt + start_beam
           z(p) = zt + zl/2.
           q(p) = qeb
           m(p) = -qeb*mass_beam
           ux(p)=vx_beam
           uy(p)=vy_beam
           uz(p)=vz_beam
           pepid(p) = me                ! processor ID
           pelabel(p) = ne+ni+i  ! labels
        endif
     end do
     npp = npp + np_beam

  endif
  npart = npart + np_beam  ! Augment particle numbers

end subroutine beam_dust



