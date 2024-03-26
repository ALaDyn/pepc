! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2023-2023 Juelich Supercomputing Centre,
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
!> A module to extend timing with user specified abbreviations.
!>
!> You can use own (frontend-defined) timer constants in the range
!> `t_userdefined_first` .. `t_userdefined_last`, e.g.
!>
!>      call timer_start(t_userdefined + 0)
!>      call timer_start(t_userdefined + 7)
!>
!> etc.
!>
!> This module only defines named constants for better readability.

! The following is a diagram of the hierarchy of a subset of the timers
! provided by this module.
!
!
!  |
!  +-> t_remesh
!  |   +
!  |   |
!  |   +-> t_remesh_interpol
!  |   |
!  |   +-> t_remesh_sort
!  |
!  +-> t_particle_count
!  |
!  +-> t_io
!  |   +
!  |   |
!  |   +-> t_checkpoint
!  |   |
!  |   +-> t_dump
!  |   |
!  |   +-> t_dump_results

module module_user_timings
   use module_timings
   implicit none

   !integer, parameter :: t_userdefined_first = 60
   !integer, parameter :: t_userdefined_last = 90
   integer, parameter :: t_remesh          = t_userdefined_first + 1 !&
   integer, parameter :: t_remesh_interpol = t_userdefined_first + 2 !&
   integer, parameter :: t_remesh_sort     = t_userdefined_first + 3 !&
   integer, parameter :: t_particle_count  = t_userdefined_first + 4 !& Not an acutal timer, rather used to track the number of particles from the timings
   integer, parameter :: t_io              = t_userdefined_first + 5 !&
   integer, parameter :: t_checkpoint      = t_userdefined_first + 6 !&
   integer, parameter :: t_dump            = t_userdefined_first + 7 !&
   integer, parameter :: t_dump_results    = t_userdefined_first + 8 !&

end module module_user_timings
