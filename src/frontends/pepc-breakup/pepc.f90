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
   use treevars, only: maxlevel
   use module_pepc
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use module_debug
   use module_checkpoint
   use module_spacefilling
   use diagnostics
   use interactions_integrator
   use iso_fortran_env
   use rng_wrapper

   ! frontend helper routines
   use helper
   use particles_resize
   use mpi
   use omp_lib
   implicit none
   ! control variable
   logical :: doDiag

   ! initialize pepc library and MPI
   call pepc_initialize("pepc-breakup", my_rank, n_ranks, .true.)
   call register_density_diag_type()

   root = my_rank .eq. 0

   call timer_start(t_user_total)
   call timer_start(t_user_init)
   call set_parameter()

#ifdef _OPENMP
  call OMP_set_num_threads(init_omp_threads)
#endif

   !=====================prepare array for density diagnostics==================
   if (density_output) then
     call init_diagnostic_verts(density_verts, minimum_x, minimum_y, minimum_z, &
                                x_length, y_length, z_length)
     if (root) then
       allocate(final_density(size(density_verts)*n_ranks))
     else
       allocate(final_density(1))
     end if
   end if

   ! NOTE: Random123 seeding
   ! Seeding procedure for Random123 (any expression that generates integer unique to the process works)
   ! Only 1 instance of seed is required for the whole of pepc-breakup. Counter will be incremented
   ! at every call of "gen_norm_double_rng()", as well as assigning a random key.
   ! Excerpt from Random123 documentation (https://www.thesalmons.org/john/random123/releases/1.11.2pre/docs/):
   !   The result is highly sensitive to small changes in the inputs, so that the sequence of values produced by simply
   !   incrementing the counter (or key) is effectively indistinguishable from a sequence of samples of a uniformly distributed
   !   random variable.

   if (resume == 1) then
     np = tnp/n_ranks
     if (my_rank < MOD(tnp, 1_kind_particle*n_ranks)) np = np + 1
     call read_particles_mpiio(itime_in, MPI_COMM_WORLD, checkin_step, tnp, particles, checkpoint_file, &
                               int(np))
     call write_particles(particles)

     ! ctr_s(1) = (my_rank + 1)*np
     ! ctr_s(2:4) = CEILING(particles(1)%x*1e5, kind=int32)
     ! key_s(1) = (my_rank + 1)*n_ranks
     ! key_s(2:4) = CEILING(particles(1)%data%v*1e17, kind=int32)

   else
     itime_in = 0

     ctr_s(1) = (my_rank + 1)*np
     ctr_s(2) = MOD(CEILING(sqrt(2.0)*10**(my_rank+11), kind=int32),100**(my_rank+1))
     ctr_s(3) = MOD(CEILING(0.5*(1+sqrt(5.0))*10**(my_rank+11), kind=int32),100**(my_rank+1))
     ctr_s(4) = MOD(CEILING(sqrt(3.0)*10**(my_rank+11), kind=int32),100**(my_rank+1))
     key_s(1) = (my_rank + 1)*n_ranks
     key_s(2) = MOD(CEILING(sqrt(2.0)*10**(my_rank+7), kind=int32),100**(my_rank+1))
     key_s(3) = MOD(CEILING(0.5*(1+sqrt(5.0))*10**(my_rank+7), kind=int32),100**(my_rank+1))
     key_s(4) = MOD(CEILING(sqrt(3.0)*10**(my_rank+7), kind=int32),100**(my_rank+1))

     call init_particles(particles, sim_type)
     ! call torus_diagnostic_xz_grid(major_radius, minor_radius, 8, particles)
     ! call torus_diagnostic_xz_breakdown(major_radius, minor_radius, 9, particles)
     ! call torus_diagnostic_xy_grid(major_radius, minor_radius, 10, particles, 1.0_8)
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
   call set_cross_section_table(trim(file_path)//"Ext_elastic_scattering_H2.txt", CS_guide, 11, 0)
   call set_cross_section_table(trim(file_path)//"rotational_excitation_J_0_2.txt", CS_guide, 12, 0)
   call set_cross_section_table(trim(file_path)//"vibrational_excitation_v_0_1.txt", CS_guide, 13, 0)
   call set_cross_section_table(trim(file_path)//"Ext_nondissociative_ionization_H2+.txt", CS_guide, 14, 0)
   call set_cross_section_table(trim(file_path)//"Ext_dissociative_ionization_H+.txt", CS_guide, 15, 1)
   total_cross_sections = 5
  !  allocate(flow_count(3))
  !  allocate(total_flow_count(3))
   call set_Xi_table(trim(file_path)//"Ohkri_Xi_H2.txt", 101, Xi_table)
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

   !$OMP PARALLEL DO if(np/omp_threads > 10) default(private) shared(particles, external_e, dt, V_loop)
   do i = 1, size(particles)
      call particle_EB_field(particles(i), external_e)
      call boris_velocity_update(particles(i), -dt*0.5_8)
      V_loop = -1.*V_loop
      call particle_EB_field(particles(i), -external_e)
      V_loop = -1.*V_loop
   end do
   !$OMP END PARALLEL DO

   ! free tree specific allocations
   call pepc_timber_tree()

   ! NOTE: Possible error in reported charge values! due to positive charges generated
   !       below the anode, which is then counted! Causes underestimation of
   !       reported electron values at anode, notable at high E/p ranges.

   do step = 0, nt - 1
      if (root) then
         write (*, *) " "
         write (*, '(a,i12)') " ====== computing step  :", step
         write (*, '(a,es12.4)') " ====== simulation time :", step*dt
      end if

      call timer_start(t_user_step)

      doDiag = MOD(step+itime_in+1, diag_interval) .eq. 0
      np = size(particles)

      call timer_start(t_boris)
      allocate(new_particles_offset(init_omp_threads, 4))
      allocate(thread_charge_count(init_omp_threads, 3))
      allocate(generic_array(init_omp_threads))
      rank_charge_count = 0.0
      total_charge_count = 0.0
      thread_charge_count = 0.0
      new_particles_offset = 0
      generic_array = 0

      !$OMP PARALLEL if(np/init_omp_threads > 10) default(private) &
      !$OMP shared(init_omp_threads, particles, new_particles_offset, np) &
      !$OMP shared(gathered_new_buffer, external_e, dt, V_loop, d, E_q_dt_m) &
      !$OMP shared(rank_charge_count, thread_charge_count, flt_geom, my_rank) &
      !$OMP shared(abs_max_CS, neutral_density, CS_tables, B0, B_p, major_radius) &
      !$OMP shared(minor_radius, plasma_dimensions, generic_array) &
      !$OMP shared(total_cross_sections, step, omp_threads)

      ! NOTE: counter and key for Random123 is redefined on thread basis.
#ifdef _OPENMP
      omp_threads = init_omp_threads
      thread_id = OMP_GET_THREAD_NUM()
      if (np/init_omp_threads .le. 10) then
        thread_id = 0
        omp_threads = 1
      end if
#else
      thread_id = 0
      omp_threads = 1
#endif
      local_size = np/omp_threads
      IStart = thread_id*local_size + 1
      ctr_s(1) = (thread_id + 1)*np
      ctr_s(2:4) = CEILING(particles(IStart)%x*1e5, kind=int32)
      key_s(1) = (thread_id + 1)*omp_threads
      key_s(2:4) = CEILING(particles(IStart)%data%v*1e17, kind=int32)

      if (thread_id .eq. (omp_threads-1)) then
        local_size = np - local_size*(omp_threads-1)
      end if
      IStop = IStart + local_size-1

      call allocate_ll_buffer(50, buffer)
      particle_guide => buffer

      new_particle_cnt = 0
      swapped_num = 0
      charge_count = 0.0
      break = 0
      ! flow_count = 0.0
      ! total_flow_count = 0.0
      do i = IStart, IStop
         call filter_and_swap(particles, flt_geom, i, IStart, IStop, swapped_num, charge_count, break)

!====================== save the resolved 'e' from tree traverse================
         traversed_e = particles(i)%results%e

         call particle_EB_field(particles(i), external_e)
         call boris_velocity_update(particles(i), dt)
         call particle_pusher(particles(i), dt)

         if (particles(i)%data%species == 0) then
            collision_checks = floor(abs(particles(i)%data%q))
            j = 0
            do while (j < collision_checks)
              call collision_update(particles(i), particle_guide, new_particle_cnt, electron_num, total_cross_sections, ctr_s, key_s)
              j = j + 1
            end do
         end if
         ! call test_ionization(particles(i), particle_guide, new_particle_cnt, electron_num)

!=============== revert 'e' to state before addition of external field==========
         particles(i)%results%e = traversed_e

         if (break .eq. 1) then
           EXIT
         end if
      end do

      new_particles_offset((thread_id + 1),1) = new_particle_cnt
      new_particles_offset((thread_id + 1),2) = swapped_num
      new_particles_offset((thread_id + 1),3) = IStart
      new_particles_offset((thread_id + 1),4) = IStop - swapped_num - IStart
      thread_charge_count((thread_id + 1),1) = charge_count(1)
      thread_charge_count((thread_id + 1),2) = charge_count(2)
      thread_charge_count((thread_id + 1),3) = charge_count(3)
      !$OMP BARRIER

      !================Gathering totals of filtered and new particles===========
      if (thread_id .eq. 0) then
        do i = 2, omp_threads
          new_particles_offset(i,1) = new_particles_offset(i,1) + new_particles_offset((i-1),1)
          new_particles_offset(i,2) = new_particles_offset(i,2) + new_particles_offset((i-1),2)
          thread_charge_count(i,1) = thread_charge_count(i,1) + thread_charge_count((i-1),1)
          thread_charge_count(i,2) = thread_charge_count(i,2) + thread_charge_count((i-1),2)
          thread_charge_count(i,3) = thread_charge_count(i,3) + thread_charge_count((i-1),3)
          generic_array(i) = new_particles_offset(i,1)
        end do

      !============Overwrite original particles(:), discarding swapped ones=====
        do i = 1,omp_threads-1
          CStart = new_particles_offset((i+1),3) - new_particles_offset(i,2)
          CStop = CStart + new_particles_offset((i+1),4)
          IStart = new_particles_offset((i+1),3)
          IStop = IStart + new_particles_offset((i+1),4)
          particles(CStart:CStop) = particles(IStart:IStop)
        end do

        generic_array(1) = new_particles_offset(1,1)
        rank_charge_count = thread_charge_count(omp_threads,:)

        if (new_particles_offset(omp_threads,1) .ne. 0) then
          allocate(gathered_new_buffer(new_particles_offset(omp_threads,1)))
        else
          allocate(gathered_new_buffer(1))
        end if
      end if
      !$OMP BARRIER

      if (new_particles_offset(omp_threads,1) .ne. 0) then
        call gather_ll_buffers_omp(buffer, generic_array, gathered_new_buffer, thread_id, omp_threads)
      end if
      call deallocate_ll_buffer(buffer)
      !$OMP END PARALLEL
      deallocate(generic_array)

      new_particle_cnt = new_particles_offset(omp_threads,1)
      swapped_num = new_particles_offset(omp_threads,2)
      deallocate(new_particles_offset)
      deallocate(thread_charge_count)

      call MPI_REDUCE(rank_charge_count, total_charge_count, 3, MPI_KIND_PHYSICS, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) then
        print *, "SUMMED CHARGE COUNT: ", total_charge_count(1), total_charge_count(2), total_charge_count(3)
      end if

      if (root) then
        if ((new_particle_cnt > 0) .or. (swapped_num /= 0 .or. electron_num /= 0)) then
           ! call extend_particles_list_swap_inject(particles, buffer, new_particle_cnt, electron_num, swapped_num)
           call extend_particles_swap_inject_omp(particles, gathered_new_buffer, new_particle_cnt, swapped_num, electron_num)
          !  print *, "extending particle list successful! New size: ", my_rank, size(particles)
        end if
      else
        if ((new_particle_cnt > 0) .or. (swapped_num /= 0)) then
          ! call extend_particles_list_swap(particles, buffer, new_particle_cnt, swapped_num)
          call extend_particles_swap_omp(particles, gathered_new_buffer, new_particle_cnt, swapped_num)
        end if
      end if

      ! np = size(particles)
      ! tnp = 0
      ! call MPI_REDUCE(np, tnp, 1, MPI_KIND_PARTICLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ! if (root) print *, "Total particles: ", tnp

      deallocate(gathered_new_buffer)
      call timer_stop(t_boris)
      if (root) write (*, '(a,es12.4)') " ====== boris_scheme [s]:", timer_read(t_boris)

      ! if (doDiag .and. particle_output) call write_particles(particles)

      !=====================Particle Merging====================================
      if (tnp > 400000) then
        sibling_upper_limit = 500
        merge_ratio = 0.75
        call compute_particle_keys(global_tree%bounding_box, particles)
        call sort_particles_by_key(particles) !Counter act jumbling by filter_and_swap(), as well as new particles.
        ! call determine_siblings_at_level(particles, sibling_cnt, unique_parents, 4_kind_level)

        !NOTE: sibling_cnt is allocated here.
        !      sibling_upper_limit is also updated to the max no. of actual siblings across all parents.
        call defined_siblings_number_grouping(particles, sibling_upper_limit, sibling_cnt, unique_parents)
        do i = 1, unique_parents
          call sort_sibling_species(particles, sibling_cnt, i, unique_parents)
        end do

        call allocate_ll_buffer(sibling_upper_limit, buffer)
        particle_guide => buffer
        new_particle_cnt = 0
        do i = 1, unique_parents
          !NOTE: actual merging. Include check, if particles(i)%data%mp_int1 == -1, don't merge!
          ! print *, "Merging ", i, "of ", unique_parents, " unique parents."
          call momentum_partition_merging(particles, sibling_cnt, sibling_upper_limit, &
                                               i, particle_guide, new_particle_cnt)
        end do
        call merge_replace_particles_list(particles, buffer, new_particle_cnt)
        call deallocate_ll_buffer(buffer)
        deallocate(sibling_cnt)
      end if
      !=========================================================================

      np = size(particles)
      tnp = 0
      call MPI_REDUCE(np, tnp, 1, MPI_KIND_PARTICLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) print *, "Total particles: ", tnp

      call pepc_particleresults_clear(particles)

      ! if (mod(step,5) .eq. 0) then
        call pepc_grow_tree(particles)
        np = size(particles, kind=kind(np))
        if (root) write (*, '(a,es12.4)') " ====== tree grow time  :", timer_read(t_fields_tree)
        call pepc_traverse_tree(particles)
        if (root) write (*, '(a,es12.4)') " ====== tree walk time  :", timer_read(t_fields_passes)

        if (doDiag .and. domain_output) call write_domain(particles)
        if (dbg(DBG_STATS)) call pepc_statistics(step)
        call pepc_timber_tree()
      ! end if

!=============================Writing output files==============================
      if (doDiag .and. particle_output) call write_particles(particles)

      ! NOTE: if density diagnostic is on, do these
      if (doDiag .and. density_output) then
         call timer_start(t_interpolate)
         do i = 1, size(particles)
           call density_interpolation(particles(i), density_verts)
         end do

         call MPI_GATHER(density_verts, size(density_verts), MPI_TYPE_density, &
                         final_density, size(density_verts), MPI_TYPE_density, 0, &
                         MPI_COMM_WORLD, ierr)
         if (root) then
           call clear_density_results(density_verts)

           do ir = 0, n_ranks-1
             do iv = 1, size(density_verts)
               cnt = ir*size(density_verts) + iv
               density_verts(iv)%q_density = final_density(cnt)%q_density + density_verts(iv)%q_density
               density_verts(iv)%J_density = final_density(cnt)%J_density + density_verts(iv)%J_density
               final_density(cnt)%q_density = 0.0
               final_density(cnt)%J_density = 0.0
             end do
           end do

           call write_densities(density_verts)
         end if
         call clear_density_results(density_verts)
         call timer_stop(t_interpolate)
         if (root) write (*, '(a,es12.4)') " ====== density interpolation [s]:", timer_read(t_interpolate)
      end if

      if (doDiag .and. particle_mpi_output) then
        call MPI_BCAST(tnp, 1, MPI_KIND_PARTICLE, 0, MPI_COMM_WORLD, ierr)
        call write_particles_mpiio(MPI_COMM_WORLD, step+itime_in+1, tnp, particles, checkpoint_file)
        call write_updated_resume_variables(step+itime_in+1)
      end if

      call timer_stop(t_user_step)
      if (root) write (*, '(a,es12.4)') " == time in step [s]                              : ", timer_read(t_user_step)
      prev_t_user_step = timer_read(t_user_step)

      call timings_GatherAndOutput(step, 0, 0 == step)

      call timer_stop(t_user_total)
      current_wall_time = timer_read(t_user_total)
      ! if (root) print *, timer_read(t_user_total)
      call timer_resume(t_user_total)
!============================preempt_checkpointing here=========================
      call preempt_checkpointing(current_wall_time, prev_t_user_step, doDiag, particles, step, break_loop)
      if (break_loop) EXIT

   end do

   call deallocate_CS_buffer(CS_tables)

   if (density_output) then
     deallocate(final_density)
     deallocate(density_verts)
   end if
   deallocate(particles)

   call timer_stop(t_user_total)

   if (root) then
      write (*, *) " "
      write (*, '(a)') " ===== finished pepc simulation"
      write (*, '(a,es12.4)') " ===== total run time [s]: ", timer_read(t_user_total)
   end if

   ! cleanup pepc and MPI
   call free_density_diag_type()
   call pepc_finalize()

end program pepc
