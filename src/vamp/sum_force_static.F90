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
  integer, dimension(1:nintmax) ::  inode

  real, dimension(1:nintmax) :: dx,dy,dz,d,dx2,dy2,dz2 
  real, dimension(1:nintmax) :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
  real, dimension(1:nintmax) :: fsx,fsy,fsz,phi

  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: eps2

!VAMPINST subroutine_start
       CALL VTENTER(IF_sum_force,VTNOSCL,VTIERR)
!      write(*,*) 'VT: sum_force S>',VTIERR,
!     *    IF_sum_force,ICLASSH
!
  eps2=eps**2

  !  preprocess distances - loops implicit from 1 -> n

  dx(1:n) = x(p) - xcoc( inode(1:n) )
  dy(1:n) = y(p) - ycoc( inode(1:n) )
  dz(1:n) = z(p) - zcoc( inode(1:n) ) 


  !  sum force to quadrupole
  !       f(p)=f(p)+Fpm+Fpd+Fpq

  d(1:n) = sqrt(dx(1:n)**2+dy(1:n)**2+dz(1:n)**2+eps2)
  rd3(1:n) = 1./d(1:n)**3
  rd5(1:n) = 1./d(1:n)**5
  rd7(1:n) = 1./d(1:n)**7

  dx2(1:n) = dx(1:n)**2
  dy2(1:n) = dy(1:n)**2
  dz2(1:n) = dz(1:n)**2
  dx3(1:n) = dx(1:n)**3
  dy3(1:n) = dy(1:n)**3
  dz3(1:n) = dz(1:n)**3

  fd1(1:n) = -3.*dx2(1:n)*rd5(1:n) + rd3(1:n)
  fd2(1:n) = -3.*dy2(1:n)*rd5(1:n) + rd3(1:n)
  fd3(1:n) = -3.*dz2(1:n)*rd5(1:n) + rd3(1:n)
  fd4(1:n) = -3.*dx(1:n)*dy(1:n)*rd5(1:n)
  fd5(1:n) = -3.*dy(1:n)*dz(1:n)*rd5(1:n)
  fd6(1:n) = -3.*dx(1:n)*dz(1:n)*rd5(1:n)

  ! potential
  phi(1:n) =   charge( inode(1:n) )/d(1:n)    &                !  monopole term
       ! 
  + (dx(1:n)*xdip( inode(1:n) ) + dy(1:n)*ydip( inode(1:n) ) + dz(1:n)*zdip( inode(1:n) ))*rd3(1:n)  &    !  dipole 
   -0.5*fd1(1:n)*xxquad( inode(1:n) ) - 0.5*fd2(1:n)*yyquad( inode(1:n) ) - 0.5*fd3(1:n)*zzquad( inode(1:n) )  &  !  quadrupole
  - fd4(1:n)*xyquad( inode(1:n) ) + fd5(1:n)*yzquad( inode(1:n) ) + fd6(1:n)*zxquad( inode(1:n) )      


  !  forces

  fsx(1:n) =   charge( inode(1:n) )*dx(1:n)*rd3(1:n) &      ! monopole term
       !
  - ( fd1(1:n)*xdip( inode(1:n) ) + fd4(1:n)*ydip( inode(1:n) ) + fd6(1:n)*zdip( inode(1:n) ) )  &   !  dipole term
       !
  + ( -9.*dx(1:n)*rd5(1:n) + 15.*dx3(1:n)*rd7(1:n) )*0.5*xxquad( inode(1:n) ) &
       + ( -3.*dy(1:n)*rd5(1:n) + 15.*dy(1:n)*dx2(1:n)*rd7(1:n) )*xyquad( inode(1:n) ) &
       + ( -3.*dz(1:n)*rd5(1:n) + 15.*dz(1:n)*dx2(1:n)*rd7(1:n) )*zxquad( inode(1:n) ) &   !   quadrupole term
       + ( 15*dx(1:n)*dy(1:n)*dz(1:n)*rd7(1:n) )*yzquad( inode(1:n) ) &
       + ( -3.*dx(1:n)*rd5(1:n) + 15.*dx(1:n)*dy2(1:n)*rd7(1:n) )*0.5*yyquad( inode(1:n) ) &
       + ( -3.*dx(1:n)*rd5(1:n) + 15.*dx(1:n)*dz2(1:n)*rd7(1:n) )*0.5*zzquad( inode(1:n) ) 

  fsy(1:n) =    charge( inode(1:n) )*dy(1:n)*rd3(1:n) &
       - ( fd2(1:n)*ydip( inode(1:n) ) + fd4(1:n)*xdip( inode(1:n) ) + fd5(1:n)*zdip( inode(1:n) ) ) &
       + ( -9.*dy(1:n)*rd5(1:n) + 15.*dy3(1:n)*rd7(1:n) )*0.5*yyquad( inode(1:n) ) &
       + ( -3.*dx(1:n)*rd5(1:n) + 15.*dx(1:n)*dy2(1:n)*rd7(1:n) )*xyquad( inode(1:n) ) &
       + ( -3.*dz(1:n)*rd5(1:n) + 15.*dz(1:n)*dy2(1:n)*rd7(1:n) )*yzquad( inode(1:n) ) &
       + ( 15.*dx(1:n)*dy(1:n)*dz(1:n)*rd7(1:n) )*zxquad( inode(1:n) ) &
       + ( -3.*dy(1:n)*rd5(1:n) + 15.*dy(1:n)*dx2(1:n)*rd7(1:n) )*0.5*xxquad( inode(1:n) ) &
       + ( -3.*dy(1:n)*rd5(1:n) + 15.*dy(1:n)*dz2(1:n)*rd7(1:n) )*0.5*zzquad( inode(1:n) ) 

  fsz(1:n) =     charge( inode(1:n) )*dz(1:n)*rd3(1:n) &
       - ( fd3(1:n)*zdip( inode(1:n) ) + fd5(1:n)*ydip( inode(1:n) ) + fd6(1:n)*xdip( inode(1:n) ) ) &
       + ( -9.*dz(1:n)*rd5(1:n) + 15.*dz3(1:n)*rd7(1:n) )*0.5*zzquad( inode(1:n) ) &
       + ( -3.*dx(1:n)*rd5(1:n) + 15.*dx(1:n)*dz2(1:n)*rd7(1:n) )*zxquad( inode(1:n) ) &
       + ( -3.*dy(1:n)*rd5(1:n) + 15.*dy(1:n)*dz2(1:n)*rd7(1:n) )*yzquad( inode(1:n) ) &
       + ( 15.*dx(1:n)*dy(1:n)*dz(1:n)*rd7(1:n) )*xyquad( inode(1:n) ) &
       + ( -3.*dz(1:n)*rd5(1:n) + 15.*dz(1:n)*dy2(1:n)*rd7(1:n) )*0.5*yyquad( inode(1:n) ) &
       + ( -3.*dz(1:n)*rd5(1:n) + 15.*dz(1:n)*dx2(1:n)*rd7(1:n) )*0.5*xxquad( inode(1:n) ) 



  !  sum partial forces over interaction list

  sumfx = SUM(fsx(1:n))
  sumfy = SUM(fsy(1:n))
  sumfz = SUM(fsz(1:n))
  sumphi = SUM(phi(1:n))

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: sum_force S<',VTIERR,ICLASSH
!
end subroutine sum_force
