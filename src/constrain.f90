 
!  ===============================================================
!
!                           CONSTRAIN
!
!   Constrain particle movement to pre-defined geometry
!
!  ===============================================================

subroutine constrain

  use treevars
  use utils
  implicit none

  integer :: i,p
  real :: gamma, hx, hy, hz, phi, delta_r
  real :: x_limit, y_limit, z_limit, r_limit, xt, yt, zt

  do p=1,npp
     xt = x(p) - plasma_centre(1)  ! shift origin to target centre
     yt = y(p) - plasma_centre(2)
     zt = z(p) - plasma_centre(3)


     r_limit = r_sphere
     x_limit = x_plasma/2.
     y_limit = y_plasma/2.

     if (initial_config==1) then
        ! sphere
        !	      constrained = (xt**2 + yt**2 + zt**2 <= r_limit**2) 


     else if (initial_config==2) then
        ! disc
        if ( yt**2 + zt**2 > r_limit**2) then
           phi = phase(yt,zt)
           hy = cos(phi) ! direction cosine
           hz = sin(phi)
           delta_r = sqrt(yt**2+zt**2)-r_limit
           y(p) = y(p) - 2*hy*delta_r  ! reflect particle back into circle
           z(p) = z(p) - 2*hz*delta_r
        endif
        if ( xt < -x_limit ) x(p) = x(p) + 2*(-x_limit - xt) 
        if ( xt >  x_limit ) x(p) = x(p) - 2*(xt - x_limit)


     else if (initial_config==3) then
        ! wire

        if ( yt**2 + xt**2 > r_limit**2) then
           phi = phase(xt,yt)
           hx = cos(phi) ! direction cosine
           hy = sin(phi)
           delta_r = sqrt(xt**2+yt**2)-r_limit
           x(p) = x(p) - 2*hx*delta_r  ! reflect particle back into circle
           y(p) = y(p) - 2*hy*delta_r
           ux(p) = -ux(p)
           uy(p) = -uy(p)
        endif
        if ( zt < -x_limit ) then
           z(p) = z(p) + 2*(-x_limit - zt)! reflect back into column
           uz(p) = -uz(p)
        else if ( zt >  x_limit ) then 
           z(p) = z(p) - 2*(zt - x_limit)
           uz(p) = -uz(p)

        endif


     else
        ! slab: periodic
        if ( xt < -x_limit ) then
           x(p) = x(p) + 2*(-x_limit-xt)
           ux(p) = -ux(p)
        else if ( xt >  x_limit ) then
           x(p) = x(p) - 2*(xt-x_limit)
           ux(p) = -ux(p)
        endif

        if ( yt < -y_limit ) then
           y(p) = y(p) + 2*(-y_limit-yt) 
           uy(p) = -uy(p)
        else if ( yt >  y_limit ) then
           y(p) = y(p) - 2*(yt-y_limit)
           uy(p) = -uy(p)
        endif
        if ( zt < -y_limit )then
           z(p) = z(p) + 2*(-y_limit-zt)
           uz(p) = -uz(p)
        else if ( zt >  y_limit ) then
           z(p) = z(p) - 2*(zt-y_limit)
           uz(p) = -uz(p)
        endif
     endif

  end do

end subroutine constrain
