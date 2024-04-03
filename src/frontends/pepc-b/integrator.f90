! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

subroutine integrator

   use module_physvars
   use module_integration_scheme
   use module_particle_boundaries
   implicit none

   ! Integration schemes (selected via scheme)
   !      1 = NVE - total energy conserved
   !      2 = NVT - global Te, Ti cons
   !      3 = global NVT electron Te conserved; ions frozen
   !      4 = local NVT: each PE keeps Te clamped; ions frozen
   !      5 = local NVT, ions only; electrons added at end of run
   !      6 = Full EM pusher (all E, B components), total energy conserved
   !      7 = ES, Nonrelativistic push, total energy conserved
   !      8 = 2v EM push, TE fields (Ex, Ey, Bz)

   if (my_rank .eq. 0) then
      write (6, *) 'Push with scheme', scheme
      write (24, *) 'Push with scheme', scheme
   end if

   ! Perform 1/2 step backwards to centre velocities at t=0
   if (itime .eq. 0) then

      centre: select case (scheme)

      case (1, 2, 3, 4, 6)
         call velocities(1, np_local, -dt / 2.)  ! pure ES, NVT ensembles

      case (7)
         call push_full3v(1, np_local, -dt / 2.)  ! full EM pusher (all E, B components)

      case (8)
         call push_TE(1, np_local, -dt / 2.)  ! full EM pusher (all E, B components)

      case default
         ! do nothing!

      end select centre
      return
   end if

   pusher: select case (scheme)

   case (1, 2, 3, 4)
      call velocities(1, np_local, dt)  ! pure ES, NVT ensembles
      call push_x(1, np_local, dt)  ! update positions

   case (5)
      call velocities(1, np_local, dt)  ! ion quiet start NVT ensemble
      call push_restrict(1, np_local, dt, a_ii / 3.)  ! restricted update

   case (6)
      call velocities(1, np_local, dt)  ! nonrelativistic push
      call push_nonrel(1, np_local, dt)

   case (7)
      call push_full3v(1, np_local, dt)  ! full EM pusher (all E, B components)
      call push_x(1, np_local, dt)  ! update positions

   case (8)
      call push_TE(1, np_local, dt)  ! full EM pusher (all E, B components)
      call push_nonrel(1, np_local, dt)
   case default
      ! do nothing!

   end select pusher

end subroutine integrator
