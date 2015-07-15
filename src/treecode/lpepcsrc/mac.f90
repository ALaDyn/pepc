! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2015 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!


!> MAC_CHOOSE
!>
!>  Picks one of several multipole acceptance criteria (MACs)
!>
!> mac=	0  Barnes Hut s/d
!>     	1  Minimum distance - replace d with shortest distance to cell box
!>	2  Bmax - replace s with longest sep between coq and box corner
!>	3  Upper bound: Use E-field from previous step to estimate acceptable error and equiv dmin 
!>	10 Periodic variant of s/d - use nearest min. image and return its quadrant for force sum.

subroutine mac_choose(p,p_ex_p,p_ey_p,p_ez_p,walk_node,walk_key,walk_abs_charge,boxl2,theta2, mac, mac_ok, vbox)
  use treevars
  use module_spacefilling
  implicit none
 


  real, intent(in) :: boxl2 !< square of the boxlength of the cell in question
  real, intent(in) :: theta2 !< square of the multipole acceptance parameter
  logical, intent(out) :: mac_ok !< on exit: if cell needs not be resolved: .true. else .false.
  integer, intent(in) :: vbox(3) !< only for periodic mode: contains vector to the actually processed near-field box
  integer, intent(in) :: mac
  integer, intent(in) :: p
  real*8, intent(in) :: walk_abs_charge
  integer*8,intent(in) :: walk_key 
  integer, intent(in) :: walk_node
  real*8, intent(in) :: p_ex_p, p_ey_p, p_ez_p  
  real*8 :: dist2
  integer :: ix,iy,iz,nbits,i,j,k
  real*8 :: xt, yt,zt
  real :: boxl
  real*8 :: x_wn,y_wn,z_wn,dist,field_old,e_dx,e_dy,e_dz
  real*8 :: b_max2,b_temp
  real*8 :: x_min,y_min,z_min,x_max,y_max,z_max,x_len,y_len,z_len,px,py,pz,md2
  real*8,dimension(3) :: x_val,y_val,z_val
  real :: alpha
  real*8 :: dx,dy,dz
  real*8 :: eps2
  real*8 :: B2,rc
  
  
  eps2 = 10**(-8)
  
  if (mac.ne.0) then
  !  get levels of twigs
  nbits = level_from_key(walk_key)
! should use:
!  nbits = node_level(walk_node)

! TODO: need to prestore boxl and corner coord as tree node property
! Do this at end of tree_fill and then for new nodes added during tree walk
! - is there a faster bit-op to extract these?

!  Ldiv2 = periodic_L*0.5
  boxl =  sqrt(boxl2)

  ix = int(SUM( (/ (2**i*ibits( walk_key,3*i  ,1 ), i=0,nbits-1) /) ))
  iy = int(SUM( (/ (2**i*ibits( walk_key,3*i+1,1 ), i=0,nbits-1) /) ))
  iz = int(SUM( (/ (2**i*ibits( walk_key,3*i+2,1 ), i=0,nbits-1) /) ))

 endif

 
  select case(mac) 
  case(0)                               ! BH-MAC
     ! write(*,*) "MAC0"
     dx = x(p) - ( tree_nodes( walk_node )%coc(1) + vbox(1) )     ! Separations
     dy = y(p) - ( tree_nodes( walk_node )%coc(2) + vbox(2) )
     dz = z(p) - ( tree_nodes( walk_node )%coc(3) + vbox(3) )
       
     dist2 = theta2*(dx*dx+dy*dy+dz*dz)

     mac_ok = (dist2 > boxl2)
!     mac_ok = .false.
  case(1)                                  !MD-MAC                
    ! write(*,*) "MAC1"
     px = x(p)
     py = y(p) 
     pz = z(p) 
  
   
     xt= ( ix*boxl + vbox(1) ) + boxmin(1)
     yt= ( iy*boxl + vbox(2) ) + boxmin(2)
     zt= ( iz*boxl + vbox(3) ) + boxmin(3)
     
     
     x_min = xt
     y_min = yt
     z_min = zt
     x_max = xt + boxl
     y_max = yt + boxl
     
     !Nur eins von beiden!!!!

     z_max = zt + boxl !3D
     !z_max = zt       !2D
     
     x_val = (/x_min-px,px-x_max,0.0_8/)
     y_val = (/y_min-py,py-y_max,0.0_8/)
     z_val = (/z_min-pz,pz-z_max,0.0_8/)
     

     x_len = maxval(x_val)
     y_len = maxval(y_val)
     z_len = maxval(z_val)

     md2=x_len*x_len+y_len*y_len+z_len*z_len
   
     dist2= theta2*md2
     
     mac_ok = (dist2 > boxl2)
     
       
  case(2)                               !b_max-MAC
    
     !write(*,*)"MAC2"
     x_wn = tree_nodes( walk_node )%coc(1) + vbox(1)
     y_wn = tree_nodes( walk_node )%coc(2) + vbox(2)
     z_wn = tree_nodes( walk_node )%coc(3) + vbox(3)

     dx = x(p) - x_wn      
     dy = y(p) - y_wn
     dz = z(p) - z_wn

     
     xt= ( ix*boxl + vbox(1) ) + boxmin(1)
     yt= ( iy*boxl + vbox(2) ) + boxmin(2)
     zt= ( iz*boxl + vbox(3) ) + boxmin(3)
     !    write(*,*) boxmin
     
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

     dx = x(p) - ( tree_nodes( walk_node )%coc(1) + vbox(1) )
     dy = y(p) - ( tree_nodes( walk_node )%coc(2) + vbox(2) )
     dz = z(p) - ( tree_nodes( walk_node )%coc(3) + vbox(3) )
     
     dist = sqrt((dx*dx+dy*dy+dz*dz)) !Distance between COC and particle
     
     mac_ok = (dist**5 * field_old * alpha >  walk_abs_charge * boxl**3) 


  case(4)     !erw b_max !! does not work
     alpha = 1000


     dx = x(p) - ( tree_nodes( walk_node )%coc(1) + vbox(1) )
     dy = y(p) - ( tree_nodes( walk_node )%coc(2) + vbox(2) )
     dz = z(p) - ( tree_nodes( walk_node )%coc(3) + vbox(3) )


     xt= ( ix*boxl + vbox(1) ) + boxmin(1)
     yt= ( iy*boxl + vbox(2) ) + boxmin(2)
     zt= ( iz*boxl + vbox(3) ) + boxmin(3)

    
     x_wn = tree_nodes( walk_node )%coc(1) + vbox(1)
     y_wn = tree_nodes( walk_node )%coc(2) + vbox(2)
     z_wn = tree_nodes( walk_node )%coc(3) + vbox(3)

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
  
     B2 = b_max2*tree_nodes( walk_node )%abs_charge
     rc = sqrt(b_max2)/2 + sqrt(sqrt(3*B2/alpha)+(b_max2/4))
     
     mac_ok = (dx**2+dy**2+dz**2 > rc**2)

  end select

end subroutine mac_choose
