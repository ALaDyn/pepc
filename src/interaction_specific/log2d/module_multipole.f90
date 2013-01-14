! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

#include "multipole.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Provides functions and operators used in the multipole / Taylor 
!> expansion of the 2D logarithmic potential following [GR1987]:
!>
!> A configuration of line charges of strength q_i perpendicular to the 
!> x-y plane at positions z_i = x_i + j y_i with j^2 = -1 creates a
!> potential:
!>
!> Phi(z) = - Sum_i q_i log(z - z_i)
!>
!> (the physical potential corresponds only to the real part of Phi) and a
!> force field:
!>
!> E(z) = - grad Phi(z) = - (Re { Phi'(z) }, - Im { Phi'(z) })
!>
!> For |z| > max_i |z_i|, the potential can be expanded as follows:
!>
!> Phi(z)  ~= - [ Q log(z) + Sum_k=1^p omega_k M_k(z) ]
!> Phi'(z) ~= - [ Q M_1(z) + Sum_k=1^p omega_k M'_k(z) ]
!>
!> with the total charge Q = Sum_i q_i and the multipole moments:
!>
!> omega_k = - Sum_i q_i O_k(z) / k
!>
!> O and M are the monomials:
!>
!> O_k(z) = z^k    M_k(z) = z^(-k)
!>
!> [GR1987] L. Greengard and V. Rokhlin, A Fast Algorithm for Particle
!>          Simulations, Journal of Computational Physics 73, 325-348 
!>          (1987)
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_multipole
  implicit none

  contains

  !
  ! The following four monomials are actually implemented as CPP macros in
  ! multipole.h to ensure inlining.
  !

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> Monomial: O_k(z) = z^k
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !complex(kind = 8) function OMultipole(k, z)
  !  implicit none
  !  integer, intent(in) :: k
  !  complex(kind = 8), intent(in) :: z
  !
  !  OMultipole = z**k
  !end function OMultipole

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> Monomial: O'_k(z) = k z^(k - 1)
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !complex(kind = 8) function OMultipolePrime(k, z)
  !  implicit none
  !  integer, intent(in) :: k
  !  complex(kind = 8), intent(in) :: z
  !
  !  OMultipolePrime = k * z**(k - 1)
  !end function OMultipolePrime

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> Monomial: M_k(z) = z^(-k)
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !complex(kind = 8) function MTaylor(k, z)
  !  implicit none
  !  integer, intent(in) :: k
  !  complex(kind = 8), intent(in) :: z
  !
  !  MTaylor = z**(-k)
  !end function MTaylor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> Monomial: M'_k(z) = -k z^(-k - 1)
  !>
  !> Used to compute the force field.
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !complex(kind = 8) function MTaylorPrime(k, z)
  !  implicit none
  !  integer, intent(in) :: k
  !  complex(kind = 8), intent(in) :: z
  !
  !  MTaylorPrime = -k * z**(-k-1)
  !end function MTaylorPrime

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> A is the higher order part of the M2M operator that shifts higher order
  !> multipole moments:
  !>
  !> Q' = Q
  !>
  !> omega_k' = A_Q,k(z0) Q + Sum_l=1^k A_k,l(z0) omega_l
  !>                                    ^^^^^^^^^
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind = 8) function ATranslate(k, l, z0)
    implicit none
    integer, intent(in) :: k, l
    complex(kind = 8), intent(in) :: z0

    ATranslate = OMultipole(k - l, z0) * BinomialCoefficient(k - 1, l - 1)
  end function ATranslate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> B is the part of the M2L operator that converts higher order multipole
  !> moments to Taylor coefficients:
  !>
  !> Phi(z) = Sum_l=0^q mu_l O_l(z)
  !>
  !> mu_l ~= B_Q,l(z0) Q + Sum_k=1,p B_l,k(z0) omega_k
  !>                                 ^^^^^^^^^
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind = 8) function BConvert(l, k, z0)
    implicit none
    integer, intent(in) :: k, l
    complex(kind = 8), intent(in) :: z0

    BConvert = (-1)**k * BinomialCoefficient(l + k - 1, k - 1) * MTaylor(l + k, z0)
  end function BConvert

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> Binomial coefficients: BinomialCoefficient(n, k) = n! / [ k! (n - k)! ]
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8 function BinomialCoefficient(n, k)
    implicit none
    integer, intent(in) :: n, k

    BinomialCoefficient = factorial(n) / (factorial(k) * factorial(n - k))
  end function BinomialCoefficient

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> Factorial function: factorial(n) = n! = 1 x 2 x ... x n
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8 function factorial(n)
    implicit none
    integer, intent(in) :: n

    select case (n)
      case ( 0)
        factorial =                            1._8
      case ( 1)
        factorial =                            1._8
      case ( 2)
        factorial =                            2._8
      case ( 3)
        factorial =                            6._8
      case ( 4)
        factorial =                           24._8
      case ( 5)
        factorial =                          120._8
      case ( 6)
        factorial =                          720._8
      case ( 7)
        factorial =                         5040._8
      case ( 8)
        factorial =                        40320._8
      case ( 9)
        factorial =                       362880._8
      case (10)
        factorial =                      3628800._8
      case (11)
        factorial =                     39916800._8
      case (12)
        factorial =                    479001600._8
      case (13)
        factorial =                   6227020800._8
      case (14)
        factorial =                  87178291200._8
      case (15)
        factorial =                1307674368000._8
      case (16)
        factorial =               20922789888000._8
      case (17)
        factorial =              355687428096000._8
      case (18)
        factorial =             6402373705728000._8
      case (19)
        factorial =           121645100408832000._8
      case (20)
        factorial =          2432902008176640000._8
      case (21)
        factorial =         51090942171709440000._8
      case (22)
        factorial =       1124000727777607680000._8
      case (23)
        factorial =      25852016738884976640000._8
      case (24)
        factorial =     620448401733239439360000._8
      case default
        factorial = gamma(real(n+1, kind = 8))
    end select
  end function factorial

end module module_multipole
