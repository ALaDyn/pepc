! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  global listing of physical units
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_units
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! basic units (here: atomic rydberg system)
      real, public, parameter :: pi                 = acos(-1.0)                  !< ludolfs number :-)
      real, public, parameter :: unit_epsilon0      =       1./(4.*pi)            !< vacuum permittivity
      real, public, parameter :: unit_hbar          =       1.                    !< reduced planck constant
      real, public, parameter :: unit_me            =       1./2.                 !< electron mass
      real, public, parameter :: unit_qe            = -sqrt(2.)                   !< electron charge
      real, public, parameter :: unit_kB            =       1.
      ! natural constants and conversion factors (from http://physics.nist.gov/cuu/Constants/index.html)
      real, public, parameter :: unit_mp_over_me    =    1836.15267247            !< mass ratio between electron and proton
      real, public, parameter :: unit_alpha         =       1./137.035999679      !< fine-structure constant
      real, public, parameter :: unit_NA            =       6.02214179E+23        !< Avogadro constant
      real, public, parameter :: unit_abohr_in_nm   =       5.2917720859E-2       !< Bohr Radius in nm
      real, public, parameter :: unit_abohr_in_m    =      unit_abohr_in_nm*1.E-9 !< Bohr Radius in m
      real, public, parameter :: unit_Ryd_in_eV     =      13.60569193            !< Rydberg energy in eV
      real, public, parameter :: unit_t0_in_s       =       4.8377687E-17         !< fundamental time unit of Rydberg atomic system in s
      real, public, parameter :: unit_P0_in_W       =      45.059494E-3           !< fundamental power unit of Rydberg atomic system in Watts
      real, public, parameter :: unit_E0_in_Vpercm  =       3.6360903E9           !< fundamental unit of the electrical field strrength of Rydberg atomic system in V/cm
      real, public, parameter :: unit_kB_in_eVperK  =       8.617343E-5           !< Boltzman constant in eV/K
      real, public, parameter :: unit_kB_in_JperK   =       1.3806504E-23         !< Boltzman constant in J/K
      real, public, parameter :: unit_kB_in_RydperK = unit_kB_in_eVperK / unit_Ryd_in_eV !< Boltzman constant in Ryd/K
      real, public, parameter :: unit_t0_in_fs      = unit_t0_in_s / 1.E-15       !< conversion factor from wp^-1 to fs
      ! derived units
      real, public, parameter :: unit_4piepsilon0   =       4.*pi*unit_epsilon0   !< 4*pi*epsilon0
      real, public, parameter :: unit_c             =          (unit_qe*unit_qe)/(unit_hbar*unit_alpha) !< speed of light
      real, public, parameter :: unit_c2            =           unit_c*unit_c     !< (speed of light)^2
      real, public, parameter :: unit_hbarc         =           unit_hbar*unit_c  !< hbar*c
      real, public, parameter :: unit_mu0           =       1./(unit_epsilon0*unit_c*unit_c) !< vacuum permeability
      real, public, parameter :: unit_qp            = -unit_qe                    !< proton charge
      real, public, parameter :: unit_mp            =  unit_mp_over_me*unit_me    !< proton mass
      real, public, parameter :: unit_abohr         = (unit_4piepsilon0*unit_hbar*unit_hbar)/(unit_me*unit_qe*unit_qe) !< bohr radius
      real, public, parameter :: unit_ERydberg      = (unit_qe*unit_qe)/(unit_4piepsilon0*2.*unit_abohr) !< rydberg energy
      real, public, parameter :: unit_t0            =  unit_hbar/unit_ERydberg    !< basic time unit
      real, public, parameter :: unit_Z0            = sqrt(unit_mu0/unit_epsilon0)!< characteristic impedance of vacuum


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      contains


end module module_units
