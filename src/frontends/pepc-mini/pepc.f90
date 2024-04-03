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

program pepc

   ! pepc modules
   use module_pepc
   use module_pepc_types
   use module_timings
   use module_debug

   ! frontend helper routines
   use helper
   implicit none

   ! timing variables
   real*8 :: timer(5)
   real*8 :: t1, t2

   ! control variable
   logical :: doDiag

  !!! initialize pepc library and MPI
   call pepc_initialize("pepc-mini", my_rank, n_ranks, .true.)

   root = my_rank .eq. 0

   timer(1) = get_time()

   call set_parameter()

   call init_particles(particles)

   timer(2) = get_time()

   if (root) write (*, '(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)

   do step = 0, nt
      if (root) then
         write (*, *) " "
         write (*, '(a,i12)') " ====== computing step  :", step
         write (*, '(a,es12.4)') " ====== simulation time :", step * dt
      end if

      timer(3) = get_time()

      doDiag = MOD(step, diag_interval) .eq. 0

      call pepc_particleresults_clear(particles)
      t1 = get_time()
      call pepc_grow_tree(particles)
      np = size(particles, kind=kind(np))
      t2 = get_time()
      if (root) write (*, '(a,es12.4)') " ====== tree grow time  :", t2 - t1
      t1 = get_time()
      call pepc_traverse_tree(particles)
      t2 = get_time()
      if (root) write (*, '(a,es12.4)') " ====== tree walk time  :", t2 - t1

      call apply_external_field()

      if (doDiag .and. domain_output) call write_domain(particles)

      if (doDiag .and. particle_probe) call compute_field()

      if (dbg(DBG_STATS)) call pepc_statistics(step)
      call pepc_timber_tree()
      !call pepc_restore_particles(np, particles)

      if (doDiag .and. particle_test) call test_particles()

      if (doDiag .and. particle_output) call write_particles(particles)

      call push_particles(particles)

      if (particle_filter) call filter_particles(particles)

      timer(4) = get_time()
      if (root) write (*, '(a,es12.4)') " == time in step [s]                              : ", timer(4) - timer(3)

      call timings_GatherAndOutput(step, 0)

   end do

   deallocate (particles)

   timer(5) = get_time()

   if (root) then
      write (*, *) " "
      write (*, '(a)') " ===== finished pepc simulation"
      write (*, '(a,es12.4)') " ===== total run time [s]: ", timer(5) - timer(1)
   end if

  !!! cleanup pepc and MPI
   call pepc_finalize()

end program pepc

