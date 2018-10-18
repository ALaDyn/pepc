! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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

!>
!> PRNG module
!> Do not trust the random properties! The PRNG used is limited to 4 bytes of randomness
!> and is not made to be called from multiple threads/ranks. We store its state per
!> thread. This will probably not produce random streams across threads but will be
!> callable from different threads w/o locking.
!>
!> Call each RNG with random_state_t%idum < 0 to seed them individually per thread.
!> If this is not done, all random states will be in sync. Also, setting idum = 0 is
!> equivalent to -1
!>
module module_random
  implicit none

  private

  integer, parameter :: IM1  = 2147483563
  integer, parameter :: IM2  = 2147483399
  real,    parameter :: AM   = 1.0/IM1
  integer, parameter :: IMM1 = IM1-1
  integer, parameter :: IA1  = 40014
  integer, parameter :: IA2  = 40692
  integer, parameter :: IQ1  = 53668
  integer, parameter :: IQ2  = 52774
  integer, parameter :: IR1  = 12211
  integer, parameter :: IR2  = 3791
  integer, parameter :: NTAB = 32
  integer, parameter :: NDIV = 1+IMM1/NTAB
  real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
  real,    parameter :: RNMX = 1.0 - eps_

  type, public :: random_state_t
    integer :: idum  = -1
    integer :: idum2 =  123456789
    integer :: iy    =  0
    integer :: iv(NTAB)
  end type

  interface random
    module procedure random4, random8, random16, random4_v, random8_v, random16_v
  end interface
  ! We have one common interface to 'produce' pseudo-random numbers for arrays and scalars
  ! of 4, 8, and 16 bytes length.
  ! !!! NOTE !!!
  ! The PRNG only uses 4 bytes to generate the random numbers, so anything longer than
  ! 4 bytes will work, but not be random to the same extent!
  
  public :: random

  contains

  subroutine random4(rnd, rnd_stat)
    implicit none
    real*4, intent(out) :: rnd
    type(random_state_t), intent(inout) :: rnd_stat

    rnd  = par_rand(rnd_stat)
  end subroutine random4

  subroutine random8(rnd, rnd_stat)
    implicit none
    real*8, intent(out) :: rnd
    type(random_state_t), intent(inout) :: rnd_stat

    rnd  = par_rand(rnd_stat)
  end subroutine random8

  subroutine random16(rnd, rnd_stat)
    implicit none
    real*16, intent(out) :: rnd
    type(random_state_t), intent(inout) :: rnd_stat

    rnd = par_rand(rnd_stat)
  end subroutine random16

  subroutine random4_v(array, rnd_stat)
    implicit none
    real*4, intent(out) :: array(:)
    type(random_state_t), intent(inout) :: rnd_stat
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand(rnd_stat)
    end do
  end subroutine random4_v

  subroutine random8_v(array, rnd_stat)
    implicit none
    real*8, intent(out) :: array(:)
    type(random_state_t), intent(inout) :: rnd_stat
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand(rnd_stat)
    end do
  end subroutine random8_v

  subroutine random16_v(array, rnd_stat)
    implicit none
    real*16, intent(out) :: array(:)
    type(random_state_t), intent(inout) :: rnd_stat
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand(rnd_stat)
    end do
  end subroutine random16_v

  !>
  !> Portable Random Number Generator, see Numerical Recipes
  !>
  !> The parameter iseed contains the state of the RNG and is intent(inout).
  !> To enable being called from multiple threads and MPI ranks, each 
  !> using a 'distinct' seed. If the seed is not initialised per thread/state,
  !> all will have the same sequence.
  !>
  function par_rand(iseed)
    implicit none
    real :: par_rand
    type(random_state_t), intent(inout) :: iseed

    integer :: j, k

    associate(idum => iseed%idum, idum2 => iseed%idum2, iy => iseed%iy, iv => iseed%iv)

       if (idum <=0) then
          if (-idum < 1) then
             idum = 1
          else
             idum = -idum
          endif

          idum2 = idum

          do j = NTAB+7,0,-1
             k = idum/IQ1
             idum = IA1 * (idum-k*IQ1) - k*IR1
             if (idum < 0 ) idum = idum + IM1

             if (j<NTAB) iv(j+1) = idum

          end do
          iy = iv(1)
       end if

       k = idum/IQ1
       idum = IA1 * (idum-k*IQ1) - k*IR1
       if (idum < 0) idum = idum + IM1

       k = idum2/IQ2
       idum2 = IA2 * (idum2-k*IQ2) - k*IR2
       if (idum2 < 0) idum2 = idum2 + IM2

       j = iy/NDIV + 1
       iy = iv(j)-idum2
       iv(j) = idum

       if (iy < 1) iy = iy + IMM1
       par_rand = AM*iy

    end associate

    if (par_rand > RNMX) par_rand = RNMX

  end function par_rand
end module
