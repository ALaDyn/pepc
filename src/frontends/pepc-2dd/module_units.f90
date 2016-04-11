! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2016 Juelich Supercomputing Centre, 
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
!>  global listing of physical units
!>
module pepca_units
      implicit none


      ! basic units (here: length=microns, time=ns)
      real*8, public, parameter :: pi                 = acos(-1.0)                  !< ludolfs number :-)
      real*8, public, parameter :: unit_epsilon0      =       1./(4.*pi)            !< vacuum permittivity
      real*8, public, parameter :: unit_hbar          =       1.                    !< reduced planck constant
      real*8, public, parameter :: unit_me            =       1./2.                 !< electron mass
      real*8, public, parameter :: unit_qe            =      -sqrt(2.)              !< electron charge
      real*8, public, parameter :: unit_kB            =       1.
      ! natural constants and conversion factors (from http://physics.nist.gov/cuu/Constants/index.html)
      real*8, public, parameter :: unit_mp_over_me    =   1836.15267245             !< mass ratio between electron and proton
      ! derived units
      real*8, public, parameter :: unit_4piepsilon0   =       4.*pi*unit_epsilon0   !< 4*pi*epsilon0
      real*8, public, parameter :: unit_alpha         =       1./137.0359895        !< fine structure constant
      real*8, public, parameter :: unit_c             =       unit_qe*unit_qe/(unit_hbar*unit_alpha) !< speed of light
      real*8, public, parameter :: unit_c2            =           unit_c*unit_c     !< (speed of light)^2
      real*8, public, parameter :: unit_qp            = -unit_qe                    !< proton charge
      real*8, public, parameter :: unit_mp            =  unit_mp_over_me*unit_me    !< proton mass
      ! conversion factors
      real*8, public, parameter :: unit_length_m_per_simunit      = 5.29177249e-11                    !< simulation unit length in m
      real*8, public, parameter :: unit_length_micron_per_simunit = unit_length_m_per_simunit/1.e-6   !< simulation unit length in microns
      real*8, public, parameter :: unit_time_s_per_simunit        = 4.8377687E-17                     !< simulation unit time in seconds
      real*8, public, parameter :: unit_time_fs_per_simunit       = unit_time_s_per_simunit/1.e-15    !< simulation unit time in femotseconds
      real*8, public, parameter :: unit_time_ns_per_simunit       = unit_time_s_per_simunit/1.e-9     !< simulation unit time in nanoseconds

!      contains


end module
