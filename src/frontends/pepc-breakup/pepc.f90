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

program pepc
   ! pepc modules
   use module_pepc
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use module_debug
   use interactions_integrator

   ! frontend helper routines
   use helper
   implicit none

   ! control variable
   logical :: doDiag

   ! initialize pepc library and MPI
   call pepc_initialize("pepc-essential", my_rank, n_ranks, .true.)

   root = my_rank .eq. 0

   call timer_start(t_user_total)
   call timer_start(t_user_init)
   call set_parameter()
   call init_particles(particles)
   call timer_stop(t_user_init)

   if (root) write (*, '(a,es12.4)') " === init time [s]: ", timer_read(t_user_init)

   ! shift velocity to half a step back, before initial condition
   ! In order for that to happen, initial field config has to be known first
   call pepc_particleresults_clear(particles)
   call pepc_grow_tree(particles)
   call pepc_traverse_tree(particles)
   call boris_scheme(particles, -dt/2.0_8, external_e, buffer, electron_num)

   ! free tree specific allocations
   call pepc_timber_tree()

   do step = 0, nt - 1
      if (root) then
         write (*, *) " "
         write (*, '(a,i12)') " ====== computing step  :", step
         write (*, '(a,es12.4)') " ====== simulation time :", step*dt
      end if

      call timer_start(t_user_step)

      doDiag = MOD(step, diag_interval) .eq. 0

      call allocate_ll_buffer(electron_num, buffer)
      call boris_scheme(particles, dt, external_e, buffer, electron_num)
      call deallocate_ll_buffer(buffer)

      call pepc_particleresults_clear(particles)
      call pepc_grow_tree(particles)
      np = size(particles, kind=kind(np))
      if (root) write (*, '(a,es12.4)') " ====== tree grow time  :", timer_read(t_fields_tree)
      call pepc_traverse_tree(particles)
      if (root) write (*, '(a,es12.4)') " ====== tree walk time  :", timer_read(t_fields_passes)

      if (doDiag .and. domain_output) call write_domain(particles)

      if (dbg(DBG_STATS)) call pepc_statistics(step)
      call pepc_timber_tree()

      if (doDiag .and. particle_test) call test_particles()
      if (doDiag .and. particle_output) call write_particles(particles)

      if (reflecting_walls) call filter_particles(particles)

      call timer_stop(t_user_step)
      if (root) write (*, '(a,es12.4)') " == time in step [s]                              : ", timer_read(t_user_step)

      call timings_GatherAndOutput(step, 0, 0 == step)

   end do

   deallocate (particles)

   call timer_stop(t_user_total)

   if (root) then
      write (*, *) " "
      write (*, '(a)') " ===== finished pepc simulation"
      write (*, '(a,es12.4)') " ===== total run time [s]: ", timer_read(t_user_total)
   end if

   ! cleanup pepc and MPI
   call pepc_finalize()

end program pepc