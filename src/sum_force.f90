!  ===================================================================
!
!                              SUM_FORCE
!
!   Calculate forces of particle from interaction list
!
!  ===================================================================

subroutine sum_force( p, n, inode, eps, sumfx, sumfy, sumfz, sumphi )
  use treevars
  implicit none
  integer, intent(in) :: p  ! particle label 
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  real, intent(in) :: eps ! smoothing parameter
  integer :: jnode, i

  real :: rd,dx,dy,dz,d,dx2,dy2,dz2 
 real :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
 real :: fsx,fsy,fsz,phi
 real, dimension(n) :: q_loc, xdl, ydl, zdl, xxql, yyql, zzql, xyql, yzql, zxql 
  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: eps2

  eps2=eps**2
  sumfx = 0
  sumfy = 0
  sumfz = 0
  sumphi = 0

! copy multipole moments into stride 1 array
  do i=1,n
    jnode=inode(i)
    q_loc(i) = charge(jnode)
    xdl(i) = xdip(jnode)
    ydl(i) = ydip(jnode)
    zdl(i) = zdip(jnode)
    xxql(i) = xxquad(jnode)
    yyql(i) = yyquad(jnode)
    zzql(i) = zzquad(jnode)
    xyql(i) = xyquad(jnode)
    yzql(i) = yzquad(jnode)
    zxql(i) = zxquad(jnode)
  end do

  do i=1,n
     
  !  preprocess distances
     jnode = inode(i)
     dx = x(p) - xcoc( jnode )
     dy = y(p) - ycoc( jnode )
     dz = z(p) - zcoc( jnode ) 

     d = sqrt(dx**2+dy**2+dz**2+eps2)
     rd = 1./d
     rd3 = rd**3
     rd5 = rd**5
     rd7 = rd**7

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

     ! potential

     sumphi = sumphi + q_loc(i)*rd    &                           !  monopole term
          ! 
     + (dx*xdl(i) + dy*ydl(i) + dz*zdl(i))*rd3  &    !  dipole 
          -0.5*fd1*xxql(i) - 0.5*fd2*yyql(i) - 0.5*fd3*zzql(i)  &  !  quadrupole
          + fd4*xyql(i) + fd5*yzql(i) + fd6*zxql(i)      

     !  forces

     sumfx = sumfx + q_loc(i)*dx*rd3 &      ! monopole term
          !
     - ( fd1*xdl(i) + fd4*ydl(i) + fd6*zdl(i) )  &   !  dipole term
          !
     + ( -9.*dx*rd5 + 15.*dx3*rd7 )*0.5*xxql(i) &
          + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*xyql(i) &
          + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*zxql(i) &   !   quadrupole term
          + ( 15*dx*dy*dz*rd7 )*yzql(i) &
          + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*0.5*yyql(i) &
          + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*0.5*zzql(i) 

     sumfy = sumfy + q_loc(i)*dy*rd3 &
          - ( fd2*ydl(i) + fd4*xdl(i) + fd5*zdl(i) ) &
          + ( -9.*dy*rd5 + 15.*dy3*rd7 )*0.5*yyql(i) &
          + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*xyql(i) &
          + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*yzql(i) &
          + ( 15.*dx*dy*dz*rd7 )*zxql(i) &
          + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*0.5*xxql(i) &
          + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*0.5*zzql(i) 

     sumfz = sumfz + q_loc(i)*dz*rd3 &
          - ( fd3*zdl(i) + fd5*ydl(i) + fd6*xdl(i) ) &
          + ( -9.*dz*rd5 + 15.*dz3*rd7 )*0.5*zzql(i) &
          + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*zxql(i) &
          + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*yzql(i) &
          + ( 15.*dx*dy*dz*rd7 )*xyql(i) &
          + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*0.5*yyql(i) &
          + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*0.5*xxql(i) 

  end do



end subroutine sum_force
