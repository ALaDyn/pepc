!  ===================================================================
!
!                              SUM_FORCE
!
!   Calculate forces of particle from interaction list
!
!  ===================================================================

subroutine sum_force( p, n, inode, sumfx, sumfy, sumfz, sumphi )
  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  integer, intent(in) :: p  ! particle label 
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode

  real, dimension(1:n) :: dx,dy,dz,d,dx2,dy2,dz2 
  real, dimension(1:n) :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
  real, dimension(1:n) :: fsx,fsy,fsz,phi

  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: eps2

!VAMPINST subroutine_start
       CALL VTENTER(IF_sum_force,VTNOSCL,VTIERR)
!      write(*,*) 'VT: sum_force S>',VTIERR,
!     *    IF_sum_force,ICLASSH
!
  eps2=eps**2

  !  preprocess distances - loops implicit from 1 -> n

  dx = x(p) - xcoc( inode )
  dy = y(p) - ycoc( inode )
  dz = z(p) - zcoc( inode ) 


  !  sum force to quadrupole
  !       f(p)=f(p)+Fpm+Fpd+Fpq

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

  ! potential
  phi =   charge( inode )/d    &                           !  monopole term
       ! 
  + (dx*xdip( inode ) + dy*ydip( inode ) + dz*zdip( inode ))*rd3  &    !  dipole 
       -0.5*fd1*xxquad( inode ) - 0.5*fd2*yyquad( inode ) - 0.5*fd3*zzquad( inode )  &  !  quadrupole
       - fd4*xyquad( inode ) + fd5*yzquad( inode ) + fd6*zxquad( inode )      


  !  forces

  fsx =   charge( inode )*dx*rd3 &      ! monopole term
       !
  - ( fd1*xdip( inode ) + fd4*ydip( inode ) + fd6*zdip( inode ) )  &   !  dipole term
       !
  + ( -9.*dx*rd5 + 15.*dx3*rd7 )*0.5*xxquad( inode ) &
       + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*xyquad( inode ) &
       + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*zxquad( inode ) &   !   quadrupole term
       + ( 15*dx*dy*dz*rd7 )*yzquad( inode ) &
       + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*0.5*yyquad( inode ) &
       + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*0.5*zzquad( inode ) 

  fsy =    charge( inode )*dy*rd3 &
       - ( fd2*ydip( inode ) + fd4*xdip( inode ) + fd5*zdip( inode ) ) &
       + ( -9.*dy*rd5 + 15.*dy3*rd7 )*0.5*yyquad( inode ) &
       + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*xyquad( inode ) &
       + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*yzquad( inode ) &
       + ( 15.*dx*dy*dz*rd7 )*zxquad( inode ) &
       + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*0.5*xxquad( inode ) &
       + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*0.5*zzquad( inode ) 

  fsz =     charge( inode )*dz*rd3 &
       - ( fd3*zdip( inode ) + fd5*ydip( inode ) + fd6*xdip( inode ) ) &
       + ( -9.*dz*rd5 + 15.*dz3*rd7 )*0.5*zzquad( inode ) &
       + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*zxquad( inode ) &
       + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*yzquad( inode ) &
       + ( 15.*dx*dy*dz*rd7 )*xyquad( inode ) &
       + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*0.5*yyquad( inode ) &
       + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*0.5*xxquad( inode ) 



  !  sum partial forces over interaction list

  sumfx = SUM(fsx)
  sumfy = SUM(fsy)
  sumfz = SUM(fsz)
  sumphi = SUM(phi)

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: sum_force S<',VTIERR,ICLASSH
!
end subroutine sum_force
