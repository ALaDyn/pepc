!  ===================================================================
!
!                              SUM_FORCE
!
!   Calculate forces of particle from interaction list
!
!  ===================================================================

subroutine sum_force( p, n, inode, sumfx, sumfy, sumfz, sumphi )
  use treevars
  implicit none
  integer, intent(in) :: p  ! particle label 
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  integer :: jnode, i

  real :: dx,dy,dz,d,dx2,dy2,dz2 
 real :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
 real :: fsx,fsy,fsz,phi

  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: eps2

  eps2=eps**2
  sumfx = 0
  sumfy = 0
  sumfz = 0
  sumphi = 0

  do i=1,n
     
  !  preprocess distances
     jnode = inode(i)
     dx = x(p) - xcoc( jnode )
     dy = y(p) - ycoc( jnode )
     dz = z(p) - zcoc( jnode ) 




     d = sqrt(dx**2+dy**2+dz**2+eps2)
     rd3 = 1./d**3
     rd5 = 1./d**5
     rd7 = 1./d**7

     dx2 = dx**2
     dy2 = dy**2
     dz2 = dz**2
     dx3 = dx**3
     dy3 = dy**3
     dz3 = dz**3

     fd1 = -3.*dx2*rd5 + rd3
     fd2 = -3.*dy2*rd5 + rd3
     fd3 = -3.*dz2*rd5 + rd3
     fd4 = -3.*dx*dy*rd5
     fd5 = -3.*dy*dz*rd5
     fd6 = -3.*dx*dz*rd5

     !  sum force, potential to quadrupole
     !       f(p)=f(p)+Fpm+Fpd+Fpq

     ! potential
     sumphi = sumphi + charge( jnode )/d    &                           !  monopole term
          ! 
     + (dx*xdip( jnode ) + dy*ydip( jnode ) + dz*zdip( jnode ))*rd3  &    !  dipole 
          -0.5*fd1*xxquad( jnode ) - 0.5*fd2*yyquad( jnode ) - 0.5*fd3*zzquad( jnode )  &  !  quadrupole
          - fd4*xyquad( jnode ) + fd5*yzquad( jnode ) + fd6*zxquad( jnode )      


     !  forces

     sumfx = sumfx + charge( jnode )*dx*rd3 &      ! monopole term
          !
     - ( fd1*xdip( jnode ) + fd4*ydip( jnode ) + fd6*zdip( jnode ) )  &   !  dipole term
          !
     + ( -9.*dx*rd5 + 15.*dx3*rd7 )*0.5*xxquad( jnode ) &
          + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*xyquad( jnode ) &
          + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*zxquad( jnode ) &   !   quadrupole term
          + ( 15*dx*dy*dz*rd7 )*yzquad( jnode ) &
          + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*0.5*yyquad( jnode ) &
          + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*0.5*zzquad( jnode ) 

     sumfy = sumfy + charge( jnode )*dy*rd3 &
          - ( fd2*ydip( jnode ) + fd4*xdip( jnode ) + fd5*zdip( jnode ) ) &
          + ( -9.*dy*rd5 + 15.*dy3*rd7 )*0.5*yyquad( jnode ) &
          + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*xyquad( jnode ) &
          + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*yzquad( jnode ) &
          + ( 15.*dx*dy*dz*rd7 )*zxquad( jnode ) &
          + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*0.5*xxquad( jnode ) &
          + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*0.5*zzquad( jnode ) 

     sumfz = sumfz + charge( jnode )*dz*rd3 &
          - ( fd3*zdip( jnode ) + fd5*ydip( jnode ) + fd6*xdip( jnode ) ) &
          + ( -9.*dz*rd5 + 15.*dz3*rd7 )*0.5*zzquad( jnode ) &
          + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*zxquad( jnode ) &
          + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*yzquad( jnode ) &
          + ( 15.*dx*dy*dz*rd7 )*xyquad( jnode ) &
          + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*0.5*yyquad( jnode ) &
          + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*0.5*xxquad( jnode ) 

  end do


end subroutine sum_force
