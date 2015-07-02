! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

!!!!!!!!!!!  SHORTCUT


module module_shortcut
  use module_pepc_kinds
  use module_pepc_types
  implicit none
  save

    real(kind_physics), parameter :: dotnine          =  0.9_kind_physics
    real(kind_physics), parameter :: dotone           =  0.1_kind_physics
    real(kind_physics), parameter :: half             =  0.5_kind_physics
    real(kind_physics), parameter :: quarter          =  0.25_kind_physics
    real(kind_physics), parameter :: tentominusfive   =  1.E-5_kind_physics
    real(kind_physics), parameter :: tentominusseven  =  1.E-7_kind_physics
    real(kind_physics), parameter :: tentominuseight  =  1.E-8_kind_physics
    real(kind_physics), parameter :: zero             =  0._kind_physics
    real(kind_physics), parameter :: one              =  1._kind_physics
    real(kind_physics), parameter :: two              =  2._kind_physics
    real(kind_physics), parameter :: three            =  3._kind_physics
    real(kind_physics), parameter :: four             =  4._kind_physics
    real(kind_physics), parameter :: five             =  5._kind_physics
    real(kind_physics), parameter :: six              =  6._kind_physics
    real(kind_physics), parameter :: seven            =  7._kind_physics
    real(kind_physics), parameter :: eight            =  8._kind_physics
    real(kind_physics), parameter :: nine             =  9._kind_physics
    real(kind_physics), parameter :: ten              =  10._kind_physics
    real(kind_physics), parameter :: twelve           =  12._kind_physics
    real(kind_physics), parameter :: fifteen          =  15._kind_physics
    real(kind_physics), parameter :: twenty           =  20._kind_physics
    real(kind_physics), parameter :: twentyfour       =  24._kind_physics
    real(kind_physics), parameter :: thirtyfive       =  35._kind_physics
    real(kind_physics), parameter :: hundredfive      =  105._kind_physics
    real(kind_physics), parameter :: pi               =  two*acos(zero)
    real(kind_physics), parameter :: oneoverpi        =  one/pi
    real(kind_physics), parameter :: c                =  one
    real(kind_physics), parameter :: prec             = 1.E-16_kind_physics

end module module_shortcut
