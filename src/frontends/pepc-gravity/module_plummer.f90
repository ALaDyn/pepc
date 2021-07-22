! Generate Plummer model initial conditions for particles,
! scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
! See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.

module plummer
  use module_pepc_kinds
  use module_pepc_types

  implicit none

  real*8,  parameter :: PI = 3.14159265358979323846_8
  real*8,  parameter :: MFRAC = 0.999_8 ! mass cut off at MFRAC of total
  real*8,  parameter :: TWOTO31 = 2147483648.0_8
  integer, parameter :: MULT = 1103515245
  integer, parameter :: MASK = int(Z'7FFFFFFF')
  integer, parameter :: ADD = 12345
  integer, parameter :: seed = 123
  
  integer:: A, B, randx, lastrand

  contains

  subroutine pranset(seed)
    implicit none
    integer, intent(in):: seed
    A = 1
    B = 0
    randx = IAND(A*seed+B,  MASK)
    A = IAND(MULT*A, MASK)
    B = IAND(MULT*B+ADD, MASK)
  end subroutine

  !< \returns A random double in [0, 1.0)
  real*8 function prand()
    implicit none
    lastrand = randx
    randx = IAND(A*randx+B, MASK)
    prand = DBLE(lastrand) / TWOTO31
  end function 

  !< \returns A random double in [xl, xh)
  real*8 function xrand(xl, xh)
    implicit none
    real*8, intent(in):: xl, xh
    xrand = xl + (xh - xl) * prand()
  end function

  !< Pick a random point on a sphere of specified radius.
  !! \param p Coordinate vector chosen
  !! \param rad Radius of chosen point 
  subroutine pickshell(p, rad)
    implicit none
    real*8, intent(in) :: rad
    real*8, intent(inout) :: p(3)
    real*8 :: rsq, rsc

    do 
      p(1) = xrand(-1._8, 1._8)
      p(2) = xrand(-1._8, 1._8)
      p(3) = xrand(-1._8, 1._8)
      rsq = DOT_PRODUCT(p, p)
      if (rsq <= 1.0) exit
    end do
    
    rsc = rad / DSQRT(rsq)
    p = p * rsc
  end subroutine

  !< Generate particles and put them in an array
  !! \param p An already allocated array of length \p np 
  !! \param tnp Total number of particles
  !! \param np Number of particles on this rank
  !! \param lb Lower bound id of particles on this rank
  !! particle ids (labels) run from 0 to tnp-1. On this rank, they are in range of [lb, lb+np)
  subroutine generate_particles(p, tnp, np, lb)
    implicit none
    integer(kind_particle), intent(in) :: tnp, np, lb 
    type(t_particle), intent(inout) :: p(np)

    type(t_particle) :: tmp1, tmp2
    integer(kind_particle) :: ip, id1, id2, ub, halftnp
    real*8 :: r, v, x, y
    real*8 :: cmr(3), cmv(3)
   
    real*8, parameter :: offset = 4._8 
    real*8, parameter :: rsc = 9*PI / 16
    real*8, parameter :: vsc = DSQRT(1.0/rsc)

    ub = lb + np ! This rank only needs to store particles in [lb, ub)

    halftnp = tnp/2 + MOD(tnp, 2_kind_particle)

    call pranset(seed)

    cmr = 0._8
    cmv = 0._8

    do ip = 0, halftnp-1
      ! Generate a particle in the first half particles 
      id1 = ip
      tmp1%label = id1
      tmp1%data%q = 1._8/tnp
      tmp1%work = 1._8
            
      do
        r = 1 / DSQRT(xrand(0._8, MFRAC)**(-2._8/3._8) - 1._8)
        if (r <= 9.0) exit
      end do

      call pickshell(tmp1%x, rsc*r)

      !if (ip .eq. 0) write(*, '(3f10.6)') tmp1%x 

      do
        x = xrand(0._8, 1._8)
        y = xrand(0._8, 0.1_8)
        if (y <= x*x*((1-x*x)**3.5_8)) exit
      end do
        
      v = DSQRT(2._8)*x / (1+r*r)**0.25_8

      call pickshell(tmp1%data%v, vsc*v)
    
      cmr = cmr + tmp1%x
      cmv = cmv + tmp1%data%v

      if (lb <= id1 .and. id1 < ub) p(id1-lb+1) = tmp1
    
      ! Generate a particle in the second half particles 
      id2 = id1 + halftnp
      if (id2 < tnp) then
        tmp2%label = id2
        tmp2%data%q = 1._8/tnp
        tmp2%work = 1._8
        tmp2%x = tmp1%x + offset
        tmp2%data%v = tmp1%data%v
    
        cmr = cmr + tmp2%x
        cmv = cmv + tmp2%data%v
        if (lb <= id2 .and. id2 < ub) p(id2-lb+1) = tmp2
      end if
    end do

    cmr = cmr / tnp
    cmv = cmv / tnp 

    do ip = 1, np
      p(ip)%x = p(ip)%x - cmr
      p(ip)%data%v = p(ip)%data%v - cmv
    end do
    
  end subroutine

end module 
