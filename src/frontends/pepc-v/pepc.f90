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

! ==============================================================
!
!
!                  PEPC-V
!
!    Parallel Efficient Parallel Coulomb-solver: Vortex particles
!
!  ==============================================================

#ifndef TEST_INTERACTION
program pepcv

   use physvars
   use manipulate_particles
   use module_pepc
   use module_timings
   use files
   use diagnostics
   use module_pepc_kinds
   use module_pepc_types
   use module_interaction_specific, only: theta2
   use iso_fortran_env
   implicit none

   integer(kind_particle) :: i
   real :: trun                     ! total run time including restarts and offset
   integer :: itime, stage, t_flag
   integer, parameter :: t_remesh = t_userdefined_first + 1

   ! Allocate array space for tree
   call pepc_initialize("pepc-v", my_rank, n_cpu, .true.)
   call pepc_read_parameters_from_first_argument()

   ! Set up O/P files
   call openfiles

   ! Time stamp
   if (my_rank .eq. 0) call stamp(output_unit, 1)
   if (my_rank .eq. 0) call stamp(run_unit, 1)

   ! Each CPU gets copy of initial data
   call pepc_setup(itime, trun)

   ! Set up particles
   call special_start()

   call pepc_prepare(3_kind_dim)

   ! Loop over all timesteps
   do while (itime .lt. nt)

      itime = itime + 1
      if (my_rank .eq. 0) write (*, '(a5,i8,a3,i8,a7,i8,a)') 'Step', itime, ' of', nt, ', using', n, ' particles -------------'

      ! Runge-Kutta time integration
      do stage = 1, rk_stages

         if (my_rank .eq. 0) write (*, *) 'Time:', trun, '/', te

         call timer_start(t_tot)

         call pepc_particleresults_clear(vortex_particles)

         if (theta2 .gt. 0.0) then
            call pepc_grow_and_traverse(vortex_particles, itime, no_dealloc=.false., no_restore=.true.)
            np = size(vortex_particles, kind=kind(np))
         else
            call direct_sum(np, vortex_particles, vortex_particles%results, my_rank, n_cpu)
         end if

         do i = 1, np
            vortex_particles(i)%results%u(1:3) = vortex_particles(i)%results%u(1:3) * force_const
            vortex_particles(i)%results%af(1:3) = vortex_particles(i)%results%af(1:3) * force_const
            vortex_particles(i)%results%div = vortex_particles(i)%results%div * force_const
         end do

         !call verify_direct()

         !if (stage == rk_stages)  call dump(itime, trun)

         call push_rk2(stage)

         if (stage .lt. rk_stages) then
            call timer_stop(t_tot)   ! total loop time without diags
            call timings_LocalOutput(itime, stage)
            call timings_GatherAndOutput(itime, stage)
         end if

         flush (output_unit)

         trun = trun + dt / rk_stages

      end do

      ! dump, if needed (need to do this before remeshing or we will loose velocity information)
      call dump(itime, trun)

      ! if remeshing is requested at all and if it is time right now, do it!
      if ((rem_freq .gt. 0) .and. (mod(itime, rem_freq) .eq. 0)) then

         if (my_rank .eq. 0) write (*, '("PEPC-V | ", a)') 'Starting remeshing...'
         call timer_start(t_remesh)

         call remeshing()

         call timer_stop(t_remesh)
         if (my_rank .eq. 0) write (*, '("PEPC-V | ", a,f12.8,a)') 'Finished remeshing after ', timer_read(t_remesh), ' seconds'
         t_flag = -rk_stages

      else

         t_flag = rk_stages

      end if

      call timer_stop(t_tot)   ! total loop time incl. remeshing if requested
      call timings_LocalOutput(itime, t_flag)
      call timings_GatherAndOutput(itime, t_flag)

      ! Some linear diagnostics
      call linear_diagnostics(itime, trun)
      call divergence_diag()

   end do

   call dump_results()

   ! deallocate array space for particles
   call cleanup(my_rank, n_cpu)

   ! Time stamp
   if (my_rank .eq. 0) call stamp(output_unit, 2)
   if (my_rank .eq. 0) call stamp(run_unit, 2)

   ! Tidy up O/P files
   call closefiles

   ! cleanup of lpepc static data
   call pepc_finalize()

end program pepcv
#endif
