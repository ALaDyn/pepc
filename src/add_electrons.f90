 
!  ===============================================================
!
!                        ADD_ELECTRONS
!
!   $Revision$
!
!   Add neutralising electrons after ions-only const-temp eqm phase
!
!  ===============================================================

subroutine add_electrons

  use treevars

  integer :: i,p
  real :: ratio_clamp  ! mass ratio used for NVT dynamics




  ! Double # particles/PE
  ne = ni
  nep = nip
  npart = ni+ne
  npp = nip + nep

  !  now have nep=nip;  ions (1:nip); electrons  (nip+1:npp) 
  !  reverse of usual order, but get mixed up by sorting later

  q(nip+1:npp) = qe        ! plasma electrons

  m(nip+1:npp) = mass_e      ! electron mass

  pepid(nip+1:npp) = me                ! processor ID


  pelabel(nip+1:npp) = me*nep + (/ (i,i=1,nep) /)  ! Electron labels: 1->ne
  pelabel(1:nip) = pelabel(1:nip) + ni             ! Augment ion labels: ne+1 -> npart

  ! zero accelerations - should really compute these for electrons
  ax(nip+1:npp) = 0.
  ay(nip+1:npp) = 0.
  az(nip+1:npp) = 0.
  pot(nip+1:npp) = 0.

  work(1:npp) = 1.   ! set work load balanced initially


!  Place electrons on top of ions
  x(nip+1:npp) = x(1:nip)
  y(nip+1:npp) = y(1:nip)
  z(nip+1:npp) = z(1:nip)

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

