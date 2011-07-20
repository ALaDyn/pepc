!  ===============================================================
!
!                        ADD_ELECTRONS
!
!   $Revision: 608 $
!
!   Add neutralising electrons after ions-only const-temp eqm phase
!
!  ===============================================================

subroutine add_electrons

  use physvars
  use treevars
  use utils

  integer :: iseed1
  real :: xt, yt, zt




  ! Double # particles/PE
  ne = ni
  nip = npp
  nep = nip
  npart = ni+ne
  npp = nep+nip

  !  now have nep=nip;  ions (1:nip); electrons  (nip+1:npp) 
  !  reverse of usual order, but get mixed up by sorting later

  q(nip+1:npp) = qe        ! plasma electrons

  m(nip+1:npp) = mass_e      ! electron mass

  pepid(nip+1:npp) = me                ! processor ID


  pelabel(nip+1:npp) =  pelabel(1:nip)  ! Electron labels: 1->ne copied from ions
  pelabel(1:nip) = pelabel(1:nip) + ni  ! Augment ion labels: ne+1 -> npart

  ! zero accelerations - should really compute these for electrons
  Ex(nip+1:npp) = 0.
  Ey(nip+1:npp) = 0.
  Ez(nip+1:npp) = 0.
  pot(nip+1:npp) = 0.

  work(1:npp) = 1.   ! set work load balanced initially
  iseed1 = -7901-me

  !  Place electrons within 1/10 of ave. ion spacing in vicinity of ions
  xt = .1*a_ii*(2*rano(iseed1)-1.)
  yt = .1*a_ii*(2*rano(iseed1)-1.)          
  zt = .1*a_ii*(2*rano(iseed1)-1.)  
  x(nip+1:npp) = x(1:nip)+xt
  y(nip+1:npp) = y(1:nip)+yt
  z(nip+1:npp) = z(1:nip)+zt

  !  Set up thermal distribution
  if (vte > 0) then
     call maxwell1(ux,nppm,nip+1,nep,vte)
     call maxwell1(uy,nppm,nip+1,nep,vte)
     call maxwell1(uz,nppm,nip+1,nep,vte)
     call scramble_v(nip+1,nep)   ! remove x,y,z correlations
  else
     call cold_start(nip+1,nep)
  endif

  !  Set ion velocities to zero
  ux(1:nip) = 0.
  uy(1:nip) = 0.
  uz(1:nip) = 0.

end subroutine add_electrons
