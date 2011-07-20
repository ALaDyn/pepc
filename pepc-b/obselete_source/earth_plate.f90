 
!  ===============================================================
!
!                          EARTH_PLATE 
!
!    Emulate grounded earth plate at z=0 
!
!  ===============================================================

subroutine earth_plate

  use physvars
  use treevars
  use utils
  implicit none

  integer :: p
  integer, save :: iseed1=5
  integer, save :: iseed2=7
  real :: x_limit, y_limit, r_limit, xt, yt, zt, xt2, yt2

  do p=1,npp
     xt = x(p) - plasma_centre(1)  ! shift origin to 1st target centre
     yt = y(p) - plasma_centre(2)
     zt = z(p) - plasma_centre(3)

     xt2 = xt - displace(1)
     yt2 = yt - displace(2)

     r_limit = r_sphere
     x_limit = x_plasma/2.
     y_limit = y_plasma/2.

     if (target_geometry==1) then
        ! sphere
        !	      constrained = (xt**2 + yt**2 + zt**2 <= r_limit**2) 


     else if (target_geometry==2) then

     else if (target_geometry==3) then
        ! wire

        if ( q(p)>0 .and. zt < -x_limit ) then
 ! Reflect ions that hit earth plate
           z(p) = z(p) + 2*(-x_limit - zt)
           uz(p) = -uz(p)

        else if ( q(p)<0 .and. zt <  -x_limit ) then 

 ! Absorb electrons and reinject with thermal velocity
 	   call reinject(vte,1,ux(p),uy(p),uz(p))
!write(*,*) 'new velocities:', ux(p),uy(p),uz(p)
           z(p)= plasma_centre(3) -x_plasma/2 + uz(p)*dt/2. 

           if (pelabel(p) <= npart/2 .and. xt**2+yt**2 > r_sphere**2) then
 ! Wire 1 electrons falling outside wire1 get reinjected inside wire1
             do while (xt**2+yt**2 > r_sphere**2)
                xt = r_sphere*(2*rano(iseed2)-1.)
                yt = r_sphere*(2*rano(iseed1)-1.)
             end do

             y(p) = yt + plasma_centre(2)
             x(p) = xt + plasma_centre(1)

           else if (pelabel(p) <= npart/2 .and. xt**2+yt**2 <= r_sphere**2) then
  ! Wire 1 electrons falling inside get put back in  wire 2
	     y(p) = y(p) + displace(2)
   	     x(p) = x(p) + displace(1)
!             pelabel(p) = pelabel(p) + npart/2  ! Relabel as wire2 electron
! write(*,*) 'wire 1 reflect:',x(p),y(p),z(p),ux(p),uy(p),uz(p),pelabel(p)
 
           else if (pelabel(p) > npart/2 .and. xt2**2+yt2**2 > r_sphere**2) then
 ! Wire 2 electrons falling outside wire2 get reinjected inside wire2
             do while (xt2**2+yt2**2 > r_sphere**2)
                xt2 = r_sphere*(2*rano(iseed2)-1.)
                yt2 = r_sphere*(2*rano(iseed1)-1.)
             end do
             y(p) = yt2 + plasma_centre(2)+displace(2)
             x(p) = xt2 + plasma_centre(1)+displace(1)

           else if (pelabel(p) > npart/2 .and. xt2**2+yt2**2 <= r_sphere**2) then
  ! Wire 2 electrons falling inside get put back in  wire 1
	     y(p) = y(p) - displace(2)
   	     x(p) = x(p) - displace(1)
 !            pelabel(p) = pelabel(p) - npart/2  ! Relabel as wire1 electron
 ! write(*,*) 'wire 2 reflect:',x(p),y(p),z(p),ux(p),uy(p),uz(p),pelabel(p)
           endif

        endif


     else
 
    endif
 
end do

end subroutine earth_plate
