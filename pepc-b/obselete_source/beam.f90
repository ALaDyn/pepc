! ==============================================
!
!                BEAM
!
!  Sets up particle beam on CPU 0
!
! ==============================================

subroutine beam
 
  use physvars
  use treevars
  use utils
  implicit none
  integer :: i, p, iseed1, iseed2
  real :: Volb, dpx, yt, zt
  real :: vx_beam, vy_beam, vz_beam, ct, cp, st, sp




  !  evaluate system constants  from inputs:
  !  beam cylinder volume:  r_beam is radius
  

  if (my_rank.eq.0) then

     write(*,*) 'Setting up particle beam ',r_beam,' x ',x_beam
   Volb = pi*r_beam**2*x_beam
   qeb = Volb*rho_beam/np_beam    ! charge
   dpx=x_beam/np_beam           ! x-axis spacing
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

   do while (i < np_beam)
     yt = r_beam*(2*rano(iseed2)-1.)          
     zt = r_beam*(2*rano(iseed1)-1.)  
     if (yt**2 + zt**2 <= r_beam**2 ) then
        i = i+1
        p = npp+i  ! put them after plasma ions
        x(p)= start_beam + dpx*i
        y(p) = yt + yl/1.9
        z(p) = zt + zl/1.8
        q(p) = qeb
        m(p) = -qeb*mass_beam
        ux(p)=vx_beam
        uy(p)=vy_beam
        uz(p)=vz_beam
        pepid(p) = me                ! processor ID
        pelabel(p) = npart+np_beam+i  ! labels
     endif
   end do
  endif

  npp = npp + np_beam

  npart = npart + np_beam  ! Augment particle numbers for all CPUs
! zero fields
  Ex(1:npp) = 0
  Ey(1:npp) = 0
  Ez(1:npp) = 0
  Bx(1:npp) = 0
  By(1:npp) = 0
  Bz(1:npp) = 0
  Ax(1:npp) = 0
  Ay(1:npp) = 0
  Az(1:npp) = 0
  Axo(1:npp) = 0
  Ayo(1:npp) = 0
  Azo(1:npp) = 0
  pot(1:npp) = 0
  work(1:npp) = 1.   ! set work load balanced initially

end subroutine beam



