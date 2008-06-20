
!> MAC_CHOOSE
!>
!>  Picks one of several multipole acceptance criteria (MACs)
!>
!> mac=	0  Barnes Hut s/d
!>     	1  Minimum distance - replace d with shortest distance to cell box
!>	2  Bmax - replace s with longest sep between coq and box corner
!>	3  Upper bound: Use E-field from previous step to estimate acceptable error and equiv dmin 
!>	10 Periodic variant of s/d - use nearest min. image and return its quadrant for force sum.
!>         - returns the index offset of the nearest image's cell as an integer(3) array in neighbour

subroutine mac_choose(p,p_ex_p,p_ey_p,p_ez_p,np_local,walk_node,walk_key,walk_abs_charge,boxl2,theta2, mac, mac_ok, neighbour)
  use treevars
  implicit none
 


  real, intent(in) :: boxl2 !< square of the boxlength of the cell in question
  real, intent(in) :: theta2 !< square of the multipole acceptance parameter
  logical, intent(out) :: mac_ok !< on exit: if cell needs not be resolved: .true. else .false.
  integer, intent(out) :: neighbour(3) !< only for periodic mode: contains offsets to nearest  image cell (i,j,k)
  integer, intent(in) :: mac
  integer, intent(in) :: p
  real*8, intent(in) :: walk_abs_charge
  integer*8,intent(in) :: walk_key 
  integer, intent(in) :: walk_node
  integer, intent(in) :: np_local
  real*8, intent(in) :: p_ex_p, p_ey_p, p_ez_p  
  real*8 :: dist2,s2
  real*8 :: deltax, deltay, deltaz
  real*8 :: Ldiv2
  integer :: ix,iy,iz,nbits,i,j,k
!  integer, intent(in) :: npshort
  real :: xt, yt,zt
  real :: boxl
  real :: x_wn,y_wn,z_wn,b_max2,b_temp,e_dx,e_dy,e_dz
  real :: x_min,y_min,z_min,x_max,y_max,z_max,x_len,y_len,z_len,px,py,pz,md2
  real,dimension(3) :: x_val,y_val,z_val
!  real :: delta_x,delta_y,delta_z 
  real :: alpha,field_old,dist
  real*8 :: d,dx,dy,dz
  real*8 :: quad_mp,eps2
  real*8 :: B2,rc
  
  
  eps2 = 10**(-8)
  
  if (mac.ne.0) then
  !  get levels of twigs
  nbits = log( 1.*walk_key )/log(8.)
! should use:
!  nbits = node_level(walk_node)
  neighbour = 0

! TODO: need to prestore boxl and corner coord as tree node property
! Do this at end of tree_fill and then for new nodes added during tree walk
! - is there a faster bit-op to extract these?

!  Ldiv2 = periodic_L*0.5
  boxl =  sqrt(boxl2)

  ix = SUM( (/ (2**i*ibits( walk_key,3*i,1 ), i=0,nbits-1) /) )
  iy = SUM( (/ (2**i*ibits( walk_key,3*i+1,1 ), i=0,nbits-1) /) )
  iz = SUM( (/ (2**i*ibits( walk_key,3*i+2,1 ), i=0,nbits-1) /) )

 endif

 
  select case(mac) 
  case(0)                               ! BH-MAC
     ! write(*,*) "MAC0"
     dx = x(p) - xcoc( walk_node )      ! Separations
     dy = y(p) - ycoc( walk_node )
     dz = z(p) - zcoc( walk_node )
     
       
     dist2 = theta2*(dx*dx+dy*dy+dz*dz)

     mac_ok = (dist2 > boxl2)
!     mac_ok = .false.
  case(1)                                  !MD-MAC                
    ! write(*,*) "MAC1"
     px = x(p)
     py = y(p) 
     pz = z(p) 
  
   
     xt=ix*boxl + xmin
     yt=iy*boxl + ymin
     zt=iz*boxl + zmin
     
     
     x_min = xt
     y_min = yt
     z_min = zt
     x_max = xt + boxl
     y_max = yt + boxl
     
     !Nur eins von beiden!!!!

     z_max = zt + boxl !3D
     !z_max = zt       !2D
     
     x_val = (/x_min-px,px-x_max,0.0/)
     y_val = (/y_min-py,py-y_max,0.0/)
     z_val = (/z_min-pz,pz-z_max,0.0/)
     

     x_len = maxval(x_val)
     y_len = maxval(y_val)
     z_len = maxval(z_val)

     md2=x_len*x_len+y_len*y_len+z_len*z_len
   
     dist2= theta2*md2
     
     mac_ok = (dist2 > boxl2)
     
       
  case(2)                               !b_max-MAC
    
     !write(*,*)"MAC2"
     x_wn = xcoc( walk_node )
     y_wn = ycoc( walk_node )
     z_wn = zcoc( walk_node )

     dx = x(p) - x_wn      
     dy = y(p) - y_wn
     dz = z(p) - z_wn

     
     xt=ix*boxl + xmin
     yt=iy*boxl + ymin
     zt=iz*boxl + zmin
 !    write(*,*) xmin,ymin,zmin
     
!     b_max2 = 0
!     b_temp = 0

!     do i=1,2
!        do j =1,2
!           do k =1,2
!              e_dx = xt+((i-1)*boxl) - x_wn      
!              e_dy = yt+((j-1)*boxl) - y_wn
!              e_dz = zt+((k-1)*boxl) - z_wn  !!3 D
 !             e_dz = zt - z_wn  !!2 D
!             b_temp = e_dx*e_dx+e_dy*e_dy+e_dz*e_dz    
!              if(b_max2 .le. b_temp ) then
!                 b_max2 = b_temp
!              endif
!           enddo
!        enddo
!     enddo

     
     !write(*,*)"WN",x_wn,y_wn,z_wn,"Box",xt,yt,zt,boxl
 
     if ((x_wn - (xt+(boxl/2))) < 0) then
       
        e_dx = xt + boxl - x_wn
     else
        e_dx = xt - x_wn
     end if

     if ((y_wn - (yt+(boxl/2))) < 0) then
        e_dy = yt + boxl - y_wn
     else
        e_dy = yt - y_wn
     end if
    
     if ((z_wn - (zt+(boxl/2))) < 0) then
        e_dz = zt + boxl - z_wn
     else
        e_dz = zt - z_wn
     end if

     b_max2 = e_dx*e_dx + e_dy*e_dy + e_dz*e_dz
     
  !   write(*,*) "theta2",theta2,"me",me
     
     dist2 = theta2*( dx*dx + dy*dy + dz*dz)
    ! write(*,*) dist2,b_max2,boxl
   
     mac_ok = (dist2 > b_max2)
!     if(me == 2 .and. mac_ok ) write(*,*) mac_ok
      
!     write(*,*) mac_ok, theta2
!     write(*,*) x_wn,y_wn,z_wn,"teilchen",x(p),y(p),z(p)
!     write(*,*)"ddd",dx,dy,dz
!     write(*,*) dist2, b_max2
  
  case(3) ! MAC E-Feld
!     write(*,*) "MAC3"
     alpha = theta2

     field_old = sqrt(p_ex_p**2 + p_ey_p**2 + p_ez_p**2)

     dx = x(p) - xcoc( walk_node )      
     dy = y(p) - ycoc( walk_node )
     dz = z(p) - zcoc( walk_node )
     
     
     dist = sqrt((dx*dx+dy*dy+dz*dz)) !Distance between COC and particle
     
     mac_ok = (dist**5 * field_old * alpha >  walk_abs_charge * boxl**3) 


  case(4)     !erw b_max !! doesn't work
     alpha = 1000


     dx = x(p) - xcoc( walk_node )      
     dy = y(p) - ycoc( walk_node )
     dz = z(p) - zcoc( walk_node )

 
     xt=ix*boxl + xmin
     yt=iy*boxl + ymin
     zt=iz*boxl + zmin
     
    
     x_wn = xcoc( walk_node )
     y_wn = ycoc( walk_node )
     z_wn = zcoc( walk_node )

     b_max2 = 0
     b_temp = 0

     do i=1,2
        do j =1,2
           do k =1,2
              e_dx = xt+((i-1)*boxl) - x_wn      
              e_dy = yt+((j-1)*boxl) - y_wn
              e_dz = zt+((k-1)*boxl) - z_wn  !!3 D
             ! e_dz = zt - z_wn  !!2 D
              b_temp = e_dx*e_dx+e_dy*e_dy+e_dz*e_dz    
              if(b_max2 .le. b_temp ) then
                 b_max2 = b_temp
              endif
           enddo
        enddo
     enddo
  
     B2 = b_max2*abs_charge(walk_node)
     rc = sqrt(b_max2)/2 + sqrt(sqrt(3*B2/alpha)+(b_max2/4))
     
     mac_ok = (dx**2+dy**2+dz**2 > rc**2)

     
!  case(10)                              ! BH-MAC periodic
!     dx = x(p) - xcoc( walk_node )      ! Separations
!     dy = y(p)  - ycoc( walk_node )
!     dz = z(p) - zcoc( walk_node )
     
!     if (dx.gt.Ldiv2) then
!        deltax = dx-periodic_L
!        neighbour(1) = -1
!     else if (dx.le.(-Ldiv2)) then
!        deltax = dx+periodic_L
!        neighbour(1) = +1
!     else
!        deltax = dx
!        neighbour(1) = 0
!     end if
     
!     if (dy.gt.Ldiv2) then
!        deltay = dy-periodic_L
!        neighbour(2) = -1
!     else if (dy.le.(-Ldiv2)) then
!        deltay = dy+periodic_L
!        neighbour(2) = +1
!     else
!        deltay = dy
!        neighbour(2) = 0
!     end if
     
!     if (dz.gt.Ldiv2) then
!        deltaz = dz-periodic_L
!        neighbour(3) = -1
!     else if (dz.le.(-Ldiv2)) then
!        deltaz = dz+periodic_L
!        neighbour(3) = +1
!     else
!        deltaz = dz
!        neighbour(3) = 0
!     end if
!     dist2 = theta2*(deltax*deltax+deltay*deltay+deltaz*deltaz)
!     mac_ok = (dist2 > boxl2) 


  end select

end subroutine mac_choose
