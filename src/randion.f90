! ==============================================
!
!                RANDION
!
!  Sets up 3D random distribution of particles
!
! ==============================================

subroutine randion

  use treevars
  use utils
  implicit none
  integer :: i,j
  integer :: idum, iseed1, iseed2, iseed3, i1, n1,p, k, nsphere, nx, ny
  real :: r(npp), phi(npp), the(npp), xt, yt, zt, radius, dpx, s

  iseed1 = -11-me      ! Select seed depending on PE
  iseed2 = -10011 - me
  iseed3 = -30013 - me
  write (ipefile,'(a,3i8)') 'Seeds: ',iseed1,iseed2, iseed3

  !  Initialise particles according to initial_config
  !       0 = random slab
  !       1 = random sphere
  !       2 = random disc
  !       3 = random wire
  !       4 = ion crystal slab

  nsphere = npp

  if (initial_config == 4) then
     s = (Vplas/ni)**(1./3.)
     nx = x_plasma/s
     ny = y_plasma/s
     if (nx*ny**2 /= ni .and. me==0) then
        write(*,*) '# ions no good for crystal config - check ni or dimensions'
     endif

  else
     dpx = x_plasma/npp
  endif

  p = 0
  do while (p < nsphere)

     ! sphere
     if (initial_config ==1 ) then
        xt = r_sphere*(2*rano(iseed1)-1.)
        yt = r_sphere*(2*rano(iseed2)-1.)          
        zt = r_sphere*(2*rano(iseed3)-1.)  
        if (xt**2 + yt**2 + zt**2 <= r_sphere**2 ) then
           p = p+1
           x(p) = xt + plasma_centre(1)
           y(p) = yt + plasma_centre(2)
           z(p) = zt + plasma_centre(3)

        endif

        ! disc
     else if (initial_config == 2) then
        yt = r_sphere*(2*rano(iseed2)-1.)          
        zt = r_sphere*(2*rano(iseed1)-1.)  


        if (yt**2 + zt**2 <= r_sphere**2 ) then
           p = p+1
           x(p)= plasma_centre(1) + dpx*p + dpx/num_pe*me-x_plasma/2.
           y(p) = yt + plasma_centre(2)
           z(p) = zt + plasma_centre(3)
        endif

        ! wire
     else if (initial_config == 3) then
        xt = r_sphere*(2*rano(iseed2)-1.)          
        yt = r_sphere*(2*rano(iseed1)-1.)  


        if (xt**2 + yt**2 <= r_sphere**2 ) then
           p = p+1
           z(p)= plasma_centre(3) + dpx*p + dpx/num_pe*me-x_plasma/2.
           y(p) = yt + plasma_centre(2)
           x(p) = xt + plasma_centre(1)
        endif

!     else if (initial_config ==4) then

!        if (p <= nep) then 
           ! random electrons
!           p=p+1
!           x(p) = x_plasma*rano(iseed1) + plasma_centre(1) -x_plasma/2.
!           y(p) = y_plasma*rano(iseed2) + plasma_centre(2) -y_plasma/2.         
!           z(p) = y_plasma*rano(iseed3) + plasma_centre(3) -y_plasma/2. 
!        else
!           p = p+1
!           x(p) = dpx*p + plasma_centre(1) -x_plasma/2.

        ! slab by default
     else
        p=p+1
        x(p) = x_plasma*rano(iseed1) + plasma_centre(1) -x_plasma/2.
        y(p) = y_plasma*rano(iseed2) + plasma_centre(2) -y_plasma/2.         
        z(p) = z_plasma*rano(iseed3) + plasma_centre(3) -z_plasma/2. 

     endif

  end do

  ! scramble to remove correlations
  iseed3=-17-4*me
  n1=nsphere
  !  exclude odd one out
  if (mod(nsphere,2).ne.0) then
     n1=n1-1
  endif

  do  i=1,n1
     k=n1*rano(iseed3)+1
     xt=x(i)
     yt=y(i)
     zt=z(i)
     x(i)=x(k)
     y(i)=y(k)
     z(i)=z(k)
     x(k)=xt
     y(k)=yt
     z(k)=zt
  end do



  q(1:nep) = qe           ! plasma electrons
  q(nep+1:npp) = qi        ! plasma ions (need Z* here)
  m(1:nep) = mass_e            ! electron mass
  m(nep+1:npp) = mass_i      ! ion mass
  pepid(1:npp) = me                ! processor ID
  pelabel(1:nep) = me*nep + (/ (i,i=1,nep) /)     ! Electron labels
  pelabel(nep+1:npp) = ne + me*nip + (/ (i,i=1,nip) /)       ! Ion labels

! zero accelerations
  ax(1:npp) = 0
  ay(1:npp) = 0
  az(1:npp) = 0
  work(1:npp) = 1.   ! set work load balanced initially

end subroutine randion

