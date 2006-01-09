!  ===================================================================
!
!                              SUM_FORCE
!
!   Calculate forces of particle from interaction list
!
!  ===================================================================

subroutine sum_force( p, n, inode, eps, sumfx, sumfy, sumfz, sumphi, load )
  use treevars
  implicit none
  integer, intent(in) :: p  ! particle label 
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  real, intent(in) :: eps ! smoothing parameter
  real, intent(out) :: load ! work load for particle p
  integer :: jnode, i,j,k 

  real*8 :: rd,dx,dy,dz,d,dx2,dy2,dz2 
 real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
 real*8 :: fsx,fsy,fsz,phi
 real*8, dimension(n*10) :: mult 
 real*8, dimension(n*3) :: coc
  real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: eps2

  eps2=eps**2
  sumfx = 0.
  sumfy = 0.
  sumfz = 0.
  sumphi = 0.

! copy multipole moments into stride 1 array
  do j=1,n
    jnode=inode(j)
    i = (j-1)*10+1
    k = (j-1)*3+1
    coc(k) = xcoc(jnode)
    coc(k+1) = ycoc(jnode)
    coc(k+2) = zcoc(jnode)
    mult(i) = charge(jnode)
    mult(i+1) = xdip(jnode)
    mult(i+2) = ydip(jnode)
    mult(i+3) = zdip(jnode)
    mult(i+4) = xxquad(jnode)
    mult(i+5) = yyquad(jnode)
    mult(i+6) = zzquad(jnode)
    mult(i+7) = xyquad(jnode)
    mult(i+8) = yzquad(jnode)
    mult(i+9) = zxquad(jnode)
  end do

  call delay(me,0)
!  write(*,*) p,' x_p=',x(p)
  do j=1,n
     
  !  preprocess distances
     i = 10*(j-1) + 1  ! multipole index
     k = 3*(j-1) + 1  ! coc index
!     jnode=inode(j)
!     dx = x(p) - xcoc(jnode)
!     dy = y(p) - ycoc(jnode)
!     dz = z(p) - zcoc(jnode) 
     dx = x(p) - coc(k)
     dy = y(p) - coc(k+1)
     dz = z(p) - coc(k+2) 

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

     sumphi = sumphi + mult(i)*rd    &                           !  monopole term
          ! 
     + (dx*mult(i+1) + dy*mult(i+2) + dz*mult(i+3))*rd3  &    !  dipole 
     !        Dx             Dy            Dz
          -0.5*fd1*mult(i+4) - 0.5*fd2*mult(i+5) - 0.5*fd3*mult(i+6)  &  !  quadrupole
          !           Qxx                 Qyy                 Qzz
          + fd4*mult(i+7) + fd5*mult(i+8) + fd6*mult(i+9)      
     !             Qxy            Qyz             Qzx

     !  forces

     sumfx = sumfx + mult(i)*dx*rd3 &      ! monopole term
          !
     - ( fd1*mult(i+1) + fd4*mult(i+2) + fd6*mult(i+3) )  &   !  dipole term
          !
     + ( -9.*dx*rd5 + 15.*dx3*rd7 )*0.5*mult(i+4) &
          + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*mult(i+7) &
          + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*mult(i+9) &   !   quadrupole term
          + ( 15*dx*dy*dz*rd7 )*mult(i+8) &
          + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*0.5*mult(i+5) &
          + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*0.5*mult(i+6) 

     sumfy = sumfy + mult(i)*dy*rd3 &
          - ( fd2*mult(i+2) + fd4*mult(i+1) + fd5*mult(i+3) ) &
          + ( -9.*dy*rd5 + 15.*dy3*rd7 )*0.5*mult(i+5) &
          + ( -3.*dx*rd5 + 15.*dx*dy2*rd7 )*mult(i+7) &
          + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*mult(i+8) &
          + ( 15.*dx*dy*dz*rd7 )*mult(i+9) &
          + ( -3.*dy*rd5 + 15.*dy*dx2*rd7 )*0.5*mult(i+4) &
          + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*0.5*mult(i+6) 

     sumfz = sumfz + mult(i)*dz*rd3 &
          - ( fd3*mult(i+3) + fd5*mult(i+2) + fd6*mult(i+1) ) &
          + ( -9.*dz*rd5 + 15.*dz3*rd7 )*0.5*mult(i+6) &
          + ( -3.*dx*rd5 + 15.*dx*dz2*rd7 )*mult(i+9) &
          + ( -3.*dy*rd5 + 15.*dy*dz2*rd7 )*mult(i+8) &
          + ( 15.*dx*dy*dz*rd7 )*mult(i+7) &
          + ( -3.*dz*rd5 + 15.*dz*dy2*rd7 )*0.5*mult(i+5) &
          + ( -3.*dz*rd5 + 15.*dz*dx2*rd7 )*0.5*mult(i+4) 

!     write(*,'(i5,a5,f12.3,a5,f12.3,a5,f12.3)') jnode,' q_j=',charge(jnode),' x_j=',xcoc(jnode),' d=',d

  end do

  ! work load - can make this for accurate by sorting leaf/twig terms

  work(p) = n    ! store for next domain decomp
  load = work(p)  ! return to calling routine
end subroutine sum_force







