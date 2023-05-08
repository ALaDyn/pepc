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
     if (mesh_mode == 0) then
       call init_diagnostic_verts(density_verts, minimum_x, minimum_y, minimum_z, &
                                  x_length, y_length, z_length)
     elseif(mesh_mode == 1) then
       call read_msh_file(mesh_name, density_verts, connectivity_tets, my_rank)
       N_element = size(connectivity_tets(:,1))
     end if

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
     ! call write_particles(particles)

     ! ctr_s(1) = (my_rank + 1)*np
     ! ctr_s(2:4) = CEILING(particles(1)%x*1e5, kind=int32)
     ! key_s(1) = (my_rank + 1)*n_ranks
     ! key_s(2:4) = CEILING(particles(1)%data%v*1e17, kind=int32)

   else
     itime_in = 0
     ! call init_particles(particles, sim_type)
     ! call torus_diagnostic_xz_grid(major_radius, minor_radius, 8, particles)
     ! call torus_diagnostic_xz_breakdown(major_radius, minor_radius, 43, particles)
     call circle_points(particles, 5.8_kind_physics, 0.0_kind_physics, 1.0_kind_physics, 256)
     ! call torus_diagnostic_xy_grid(major_radius, minor_radius, 10, particles, 1.0_8)
   end if

   if (root) print *, "Particle initialisation done."

   E_q_dt_m = (e*(1.0e12))/(4.0*pi*eps_0*e_mass*c)
   call poloidal_B_grid(B_pol_grid, 200, 200, 4.05_kind_physics, 1.75_kind_physics, &
                        3.5_kind_physics, 3.5_kind_physics)
   !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
   do i = 1, size(particles)
     call particle_EB_field(particles(i), external_e, B_pol_grid)
   end do

   call write_particles(particles)
   
   np = size(particles)
   tnp = 0
   call MPI_ALLREDUCE(np, tnp, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, ierr)
   if (root) print *, "tnp: ", tnp, np
   if (root) print*, 'init_x    ', 'init_y   ', 'init_z   ', 'length   ', 'init_Bmag   ', 'final_Bmag   '
!   call write_particles_mpiio(MPI_COMM_WORLD, 1, tnp, particles, checkpoint_file)
   allocate(global_table2(6, tnp))

   call system_clock(time_start, rate)
   !$OMP PARALLEL  default(none) &
   !$OMP shared(init_omp_threads, np) &
   !$OMP shared(V_loop, d, E_q_dt_m) &
   !$OMP shared(my_rank, global_table2) &
   !$OMP shared(major_radius) &
   !$OMP shared(minor_radius, plasma_dimensions) &
   !$OMP shared(omp_threads, allowed_wall_time) &
   !$OMP private(dummy, start_i, neutral_density, local_min_x) &
   !$OMP firstprivate(B_pol_grid, particles, time_start, rate)
   dummy = OMP_GET_NUM_THREADS()
   start_i = OMP_GET_THREAD_NUM()
   if (start_i == 0) print *, "Total threads: ", dummy

   !$OMP DO SCHEDULE(DYNAMIC,1)
   do i = 1, np
     ! call particle_EB_field(particles(i), external_e, B_pol_grid)
     ! Using local_table2(:,:) as a dummy carrier. Data structure as follows:
     ! local_table2(1:3,:) = particle's initial coordinate.
     ! local_table2(4,:)   = connection length
     ! local_table2(5,:)   = particle's init B_mag
     ! local_table2(6,:)        = particle's final B_mag
     ! Using global_table2(:,:) as a global collector.
     ! neutral_density as a dummy for B_mag.
     ! local_min_x as a dummy for origin of torus.
     ! print *, start_i, i
     local_min_x = 0.0_kind_physics
     neutral_density = sqrt(dot_product(particles(i)%data%b,particles(i)%data%b))
     global_table2(1,i) = particles(i)%x(1)
     global_table2(2,i) = particles(i)%x(2)
     global_table2(3,i) = particles(i)%x(3)
     global_table2(5,i) = neutral_density
     call streakline_integral_timing(particles(i), major_radius, minor_radius, local_min_x, global_table2(4,i), time_start, rate, allowed_wall_time)
     neutral_density = sqrt(dot_product(particles(i)%data%b,particles(i)%data%b))
     global_table2(6,i) = neutral_density

     
     print*, start_i,  global_table2(1,i), global_table2(2,i), global_table2(3,i), global_table2(4,i), global_table2(5,i), global_table2(6,i) 
   end do
   !$OMP END PARALLEL
   step = 0
   call write_particles_mpiio(MPI_COMM_WORLD, step+itime_in+1, tnp, particles, checkpoint_file)

   if (root) print *, "Done Write!"

   ! call MPI_GATHER(local_table2, size(local_table2), MPI_KIND_PHYSICS, &
   !                 global_table2, size(local_table2), MPI_KIND_PHYSICS, 0, &
   !                 MPI_COMM_WORLD, ierr)
   if (root) print *, "after GATHER"
   call connection_length_output(global_table2, 23)
  
   if (root) print *, "done output"
   deallocate(global_table2)
!=============================Writing output files==============================
!      if (doDiag .and. particle_output) then
!   call charge_poloidal_distribution(particles, local_table2, global_table2, tnp, itime_in + step + 1)
!        call write_particles(particles)
!        write(file_name, '(A6,I10.10,A4)') 'table_', itime_in + step + 1, '.txt'
!        file_name = trim(file_name)
        ! call toroidal_weight_distribution(particles, local_table2, 1000)
        ! call gather_weights_tables(local_table2, global_table2, 55, file_name)
        ! call toroidal_max_weights(particles, local_table1D, local_table1D_1, 1000)
        ! call gather_minmaxWeights_tables(local_table1D, local_table1D_1, global_table1D, global_table1D_1, 55, file_name)
        ! call unit_vector_distribution(particles, local_table2, 50, 100)
        ! call gather_spherical_angle_tables(local_table2, global_table2, 55, file_name)
        ! call V_mean_phi_distribution(particles, local_table2, 1000)
        ! call V_mean_phi_distribution_gather(local_table2, global_table2, 55, file_name)
!        call V_par_perp_calculation(particles, local_table2)
!        call V_par_perp_histogram(local_table2, 100, tnp, 55, file_name)
!      end if

   deallocate(B_pol_grid)
   deallocate(particles)

   ! cleanup pepc and MPI
   call free_density_diag_type()
   call pepc_finalize()

end program pepc
