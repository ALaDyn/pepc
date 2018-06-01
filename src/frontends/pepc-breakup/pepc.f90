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
   use module_checkpoint
   use interactions_integrator
   use rng_wrapper

   ! frontend helper routines
   use helper
   implicit none
   include 'mpif.h'
   ! control variable
   logical :: doDiag

   ! initialize pepc library and MPI
   call pepc_initialize("pepc-breakup", my_rank, n_ranks, .true.)

  !  seed = (/ 0, 0 /)
  !  call frand123Init( state, my_rank, 0, seed )

   root = my_rank .eq. 0

   call timer_start(t_user_total)
   call timer_start(t_user_init)
   call set_parameter()

   if (resume == 1) then
     np = tnp/n_ranks
     if (my_rank < MOD(tnp, 1_kind_particle*n_ranks)) np = np + 1
     call read_particles_mpiio(itime_in, MPI_COMM_WORLD, checkin_step, tnp, particles, checkpoint_file, &
                               int(np))
     call write_particles(particles)
   else
     itime_in = 0
     call init_particles(particles)
    !  call torus_diagnostic_grid(major_radius, minor_radius, 7, particles)
   end if

   !========================read cross section data======================
   allocate(CS_tables)
   CS_guide => CS_tables
   ! IMPORTANT NOTE: the order of set_cross_section_table must correspond to the case orders in collision_update()
  !  call getcwd(file_path)
   file_path = "../src/frontends/pepc-breakup/cross_sections/"
   call set_cross_section_table(trim(file_path)//"total_scattering_H2.txt", CS_guide, 11, 1)

   ! NOTE: Future prospect: add proper function to maximize the collision freq. over energy (assuming initial density is highest, hence constant)
   !       look at some external function DDFSA or DFSA.
   call determine_absolute_max_CS(CS_tables, abs_max_CS)
   call deallocate_CS_buffer(CS_tables)

   allocate(CS_tables)
   CS_guide => CS_tables
   call set_cross_section_table(trim(file_path)//"elastic_scattering_H2.txt", CS_guide, 11, 0)
   call set_cross_section_table(trim(file_path)//"rotational_excitation_J_0_2.txt", CS_guide, 12, 0)
   call set_cross_section_table(trim(file_path)//"vibrational_excitation_v_0_1.txt", CS_guide, 13, 0)
   call set_cross_section_table(trim(file_path)//"nondissociative_ionization_H2+.txt", CS_guide, 14, 0)
   call set_cross_section_table(trim(file_path)//"dissociative_ionization_H+.txt", CS_guide, 15, 1)
   allocate(cross_sections_vector(5))
   !=====================================================================
   call timer_stop(t_user_init)

   if (root) write (*, '(a,es12.4)') " === init time [s]: ", timer_read(t_user_init)

   ! shift velocity to half a step back, before initial condition
   ! In order for that to happen, initial field config has to be known first
   call pepc_particleresults_clear(particles)
   call pepc_grow_tree(particles)
   call pepc_traverse_tree(particles)

   neutral_density = calculate_neutral_density(pressure, init_temperature)
   if (root) write (*, '(a,es12.4)') " === number density of neutrals: ", neutral_density
   E_q_dt_m = (e*(1.0e12))/(4.0*pi*eps_0*e_mass*c)

   ! number of injected electrons per time step
   electron_num = 20

   do i = 1, size(particles)
      call particle_EB_field(particles(i), external_e)
      call boris_velocity_update(particles(i), -dt*0.5_8)
      V_loop = -1.*V_loop
      call particle_EB_field(particles(i), -external_e)
      V_loop = -1.*V_loop
   end do

   ! free tree specific allocations
   call pepc_timber_tree()

   do step = 0, nt - 1
      if (root) then
         write (*, *) " "
         write (*, '(a,i12)') " ====== computing step  :", step
         write (*, '(a,es12.4)') " ====== simulation time :", step*dt
      end if

      call timer_start(t_user_step)

      doDiag = MOD(step+1, diag_interval) .eq. 0

      call allocate_ll_buffer(electron_num, buffer)
      particle_guide => buffer
      call timer_start(t_boris)

      new_particle_cnt = 0
      swapped_num = 0
      charge_count = 0.0
      total_charge_count = 0.0
      do i = 1, size(particles)
         if (i > (size(particles)-swapped_num)) then
           EXIT
         end if
         call filter_and_swap(particles, 2, i, swapped_num)

         call particle_EB_field(particles(i), external_e)
         call boris_velocity_update(particles(i), dt)
         call particle_pusher(particles(i), dt)

         if (particles(i)%data%species == 0) then
           call collision_update(particles(i), particle_guide, new_particle_cnt, electron_num, cross_sections_vector)
         end if
        !  call test_ionization(particles(i), particle_guide, new_particle_cnt, electron_num)
      end do

      call MPI_REDUCE(charge_count, total_charge_count, 3, MPI_KIND_PHYSICS, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) then
        print *, "SUMMED CHARGE COUNT: ", total_charge_count(1), total_charge_count(2), total_charge_count(3)
        ! total_charge_count(1) = -1.0*electron_num
        ! call write_text_output(total_charge_count(1), total_charge_count(2), step)
      end if

      if (root) then
        if ((new_particle_cnt > 0) .or. (swapped_num /= 0 .or. electron_num /= 0)) then
           call extend_particles_list_swap_inject(particles, buffer, new_particle_cnt, electron_num, swapped_num)
          !  print *, "extending particle list successful! New size: ", my_rank, size(particles)
        end if
      else
        if ((new_particle_cnt > 0) .or. (swapped_num /= 0)) then
          call extend_particles_list_swap(particles, buffer, new_particle_cnt, swapped_num)
        end if
      end if

      np = size(particles)
      tnp = 0
      call MPI_REDUCE(np, tnp, 1, MPI_KIND_PARTICLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) print *, "Total particles: ", tnp

      call timer_stop(t_boris)
      if (root) write (*, '(a,es12.4)') " ====== boris_scheme [s]:", timer_read(t_boris)
      call deallocate_ll_buffer(buffer)

      ! if (doDiag .and. particle_output) call write_particles(particles)

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

      if (doDiag .and. particle_mpi_output) then
        call MPI_BCAST(tnp, 1, MPI_KIND_PARTICLE, 0, MPI_COMM_WORLD, ierr)
        call write_particles_mpiio(MPI_COMM_WORLD, step+itime_in+1, tnp, particles, checkpoint_file)
      end if

      if (reflecting_walls) call filter_particles(particles)

      call timer_stop(t_user_step)
      if (root) write (*, '(a,es12.4)') " == time in step [s]                              : ", timer_read(t_user_step)

      call timings_GatherAndOutput(step, 0, 0 == step)

   end do

   call deallocate_CS_buffer(CS_tables)
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
