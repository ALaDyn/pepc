module treevars

  implicit none

  integer, parameter :: N = 500000
  real*8 :: x(N), y(N), z(N)
  real*8 :: charge(N), xcoc(N), ycoc(N), zcoc(N)
  real*8 :: xdip(N), ydip(N), zdip(N)
  real*8 :: xxquad(N), xyquad(N), yyquad(N)
  real*8 :: yzquad(N), zxquad(N), zzquad(N)
  integer :: access1(N), access2(N)


  contains
    subroutine my_rand(val)
      
      implicit none
      real*8, intent(inout) :: val
      real :: tmp
      
      call random_number(tmp)
      val = tmp

    end subroutine my_rand

end module treevars


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is directly involved in force calculation
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_calc_force
  implicit none
  save

contains


  subroutine calc_force_coulomb(p, inode, vbox, eps, sumfx, sumfy, sumfz, sumphi)
    use treevars
    implicit none

    include 'mpif.h'

    integer, intent(in) :: p  !< particle label
    integer, intent(in) :: inode !< index of particle to interact with
    real*8, intent(in) :: vbox(3) !< vector to neighbour box that is currently processed
    real, intent(in) :: eps ! smoothing parameter
    real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

    real*8 :: rd,dx,dy,dz,d,dx2,dy2,dz2
    real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
    real :: eps2

    eps2=eps**2
    sumfx = 0.
    sumfy = 0.
    sumfz = 0.
    sumphi = 0.

    !  preprocess distances
    dx = x(p) - ( xcoc(inode) + vbox(1) )
    dy = y(p) - ( ycoc(inode) + vbox(2) )
    dz = z(p) - ( zcoc(inode) + vbox(3) )


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

    fd1 = 3.*dx2*rd5 - rd3
    fd2 = 3.*dy2*rd5 - rd3
    fd3 = 3.*dz2*rd5 - rd3
    fd4 = 3.*dx*dy*rd5
    fd5 = 3.*dy*dz*rd5
    fd6 = 3.*dx*dz*rd5

    ! potential

    sumphi = sumphi + charge(inode)*rd    &                           !  monopole term
                                !
         + (dx*xdip(inode) + dy*ydip(inode) + dz*zdip(inode))*rd3  &    !  dipole
                                !     Dx             Dy            Dz
         + 0.5*fd1*xxquad(inode) + 0.5*fd2*yyquad(inode) + 0.5*fd3*zzquad(inode)  &  !  quadrupole
                                !           Qxx                 Qyy                 Qzz
         + fd4*xyquad(inode) + fd5*yzquad(inode) + fd6*zxquad(inode)
    !   Qxy            Qyz             Qzx

    !  forces

    sumfx = sumfx + charge(inode)*dx*rd3 &      ! monopole term
                                !
         + fd1*xdip(inode) + fd4*ydip(inode) + fd6*zdip(inode)   &   !  dipole term
                                !
         + (15.*dx3*rd7 - 9.*dx*rd5 )*0.5*xxquad(inode) &     !
         + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*xyquad(inode) &     !
         + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*zxquad(inode) &     !   quadrupole term
         + ( 15*dx*dy*dz*rd7 )*yzquad(inode) &                !
         + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*0.5*yyquad(inode) & !
         + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*0.5*zzquad(inode)   !

    sumfy = sumfy + charge(inode)*dy*rd3 &
         + fd2*ydip(inode) + fd4*xdip(inode) + fd5*zdip(inode)  &
         + ( 15.*dy3*rd7 - 9.*dy*rd5 )*0.5*yyquad(inode) &
         + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*xyquad(inode) &
         + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*yzquad(inode) &
         + ( 15.*dx*dy*dz*rd7 )*zxquad(inode) &
         + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*0.5*xxquad(inode) &
         + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*0.5*zzquad(inode)

    sumfz = sumfz + charge(inode)*dz*rd3 &
         + fd3*zdip(inode) + fd5*ydip(inode) + fd6*xdip(inode)  &
         + ( 15.*dz3*rd7 - 9.*dz*rd5 )*0.5*zzquad(inode) &
         + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*zxquad(inode) &
         + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*yzquad(inode) &
         + ( 15.*dx*dy*dz*rd7 )*xyquad(inode) &
         + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*0.5*yyquad(inode) &
         + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*0.5*xxquad(inode)

  end subroutine calc_force_coulomb

end module module_calc_force



program force_test

  use treevars

  use module_calc_force

  implicit none

  include 'mpif.h'

  integer :: cnt, cnt2, p1, p2

  real*8 :: t1, t2

  real*8 :: tmp

  real*8 :: ex(N), ey(N), ez(N), phi(N)
  real*8 :: integral

  integer :: ierr

  do cnt=1, N

     call my_rand(x(cnt))
     call my_rand(y(cnt))
     call my_rand(z(cnt))

     call my_rand(charge(cnt))
     call my_rand(xcoc(cnt))
     call my_rand(ycoc(cnt))
     call my_rand(zcoc(cnt))

     call my_rand(xdip(cnt))
     call my_rand(ydip(cnt))
     call my_rand(zdip(cnt))

     call my_rand(xxquad(cnt))
     call my_rand(xyquad(cnt))
     call my_rand(yyquad(cnt))

     call my_rand(yzquad(cnt))
     call my_rand(zxquad(cnt))
     call my_rand(zzquad(cnt))

     call my_rand(tmp)
!     access1(cnt) = int(tmp*N)
     access1(cnt) = cnt

     call my_rand(tmp)
!     access2(cnt) = int(tmp*N)
     access2(cnt) = cnt
     

     ex(cnt) = 0.0
     ey(cnt) = 0.0
     ez(cnt) = 0.0
     phi(cnt) = 0.0

  end do

  call MPI_Init(ierr)

  t1 = MPI_Wtime()

  do cnt=1, N

     p1 = access1(cnt)
     p2 = access2(cnt)

     call calc_force_coulomb(p1, p2, [0._8,0._8,0._8], 1e-6, ex(p1), ey(p1), ez(p1), phi(p1))
     
     if(modulo(cnt, N/10) .eq. 0) write(*,*) "progress [%]: ", (100.0*cnt)/N, p1, p2, ex(cnt), phi(cnt)

  end do

  t2 = MPI_Wtime()

  call MPI_Finalize(ierr)

  integral = 0.0

  do cnt=1, N
     integral = integral + ex(cnt) + ey(cnt) + ez(cnt) + phi(cnt)
  end do

  write(*,*) "fertig", integral
  write(*,*) "time[s] ", t2-t1
  write(*,*) "FP ops [M] ", 179*N/1000000.0
  write(*,*) "perf [MF] ", (179.0*N) / (t2-t1) / 1000000.0
  write(*,*) "perf [%]  ", (179.0*N) / (t2-t1) / 1000000.0 /3400. * 100.

end program force_test
