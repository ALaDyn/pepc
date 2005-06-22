 
!  ===============================================================
!
!                       DOUBLE_TARGET 
!
!   $Revision$
!
!   Add secondary target shifted by displacement vector displace(1:3) 
!
!  ===============================================================


subroutine double_target

  use physvars
  use treevars
  use utils

  integer :: i,p, iseed1
  real :: xt, yt, zt



! Make copy of target displaced by vector displace(1:3)

  x(npp+1:2*npp) = x(1:npp) + displace(1)
  y(npp+1:2*npp) = y(1:npp) + displace(2)
  z(npp+1:2*npp) = z(1:npp) + displace(3)

  ux(npp+1:2*npp) = ux(1:npp) ! same velocities      
  uy(npp+1:2*npp) = uy(1:npp) ! 
  uz(npp+1:2*npp) = uz(1:npp) !
 
  q(npp+1:2*npp) = q(1:npp) ! same charge      
  m(npp+1:2*npp) = m(1:npp)  ! same mass

  pepid(npp+1:2*npp) = me                ! processor ID


  pelabel(npp+1:2*npp) =  pelabel(1:npp) + npart  ! labels: electrons_1, ions_1, electrons_2, ions_2 

  ! zero accelerations - should really compute these for electrons
  Ex(npp+1:2*npp) = 0.
  Ey(npp+1:2*npp) = 0.
  Ez(npp+1:2*npp) = 0.
  pot(npp+1:2*npp) = 0.
  Ax(1:2*npp)  = 0.
  Ay(1:2*npp)  = 0.
  Az(1:2*npp)  = 0.

  work(1:2*npp) = 1.   ! set work load balanced initially
  iseed1 = -7901-me



  ! Double # particles/PE
  nip = 2*nip
  nep = 2*nep
  ni = 2*ni
  ne = 2*ne
  npart = ni+ne
  npp = 2*npp
end subroutine double_target

