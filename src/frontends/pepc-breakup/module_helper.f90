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

!>
!> helper module
!>
module helper
   use module_box
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use iso_fortran_env
   use rng_wrapper
   use, intrinsic :: iso_c_binding, only: c_double

   implicit none

   ! timing variables
   integer, parameter :: t_user_total = t_userdefined_first
   integer, parameter :: t_user_init = t_userdefined_first + 1
   integer, parameter :: t_user_step = t_userdefined_first + 2
   integer, parameter :: t_user_directsum = t_userdefined_first + 3
   integer, parameter :: t_user_particleio = t_userdefined_first + 4
   integer, parameter :: t_boris = t_userdefined_first + 5
   integer, parameter :: t_interpolate = t_userdefined_first + 6

   ! MPI variables
   integer(kind_pe) :: my_rank, n_ranks, ierr
   logical :: root

   ! time variables
   real*8 :: dt
   integer :: step

   ! control variables
   integer :: nt ! number of timesteps
   integer(kind_particle) :: tnp ! total number of particles
   integer(kind_particle) :: np ! local number of particles
   logical :: particle_output ! turn vtk output on/off
   logical :: domain_output ! turn vtk output on/off
   logical :: particle_mpi_output !turn mpi IO on/off
   logical :: particle_test ! check tree code results against direct summation
   logical :: reflecting_walls ! reflect particles at walls
   integer :: diag_interval ! number of timesteps between all diagnostics and IO
   real(kind_physics) :: plasma_dimensions(3) ! size of the simulation box
   real(kind_physics) :: external_e(3), traversed_e(3)

   integer, parameter :: particle_direct = -1 ! number of particle for direct summation

   type, public :: linked_list_elem
      type(linked_list_elem), pointer :: next
      type(t_particle), allocatable :: tmp_particles(:)
   end type linked_list_elem

   type, public :: linked_list_CS
      type(linked_list_CS), pointer :: next_CS
      real(kind_physics), allocatable :: CS(:,:)
   end type linked_list_CS

   type, public :: diag_vertex
      real(kind_physics) :: x(1:3) = 0.0
      real(kind_physics) :: q_density(1:2) = 0.0
      real(kind_physics) :: J_density(1:3) = 0.0
   end type diag_vertex

   ! particle data (position, velocity, mass, charge)
   type(t_particle), allocatable   :: particles(:)
   real(kind_physics), allocatable :: direct_L2(:)

   ! particle merging variables
   type(t_particle), allocatable   :: merged_particles(:)
   type(linked_list_elem), pointer :: merge_buffer_0, merge_buffer_1
   integer, allocatable :: sibling_cnt(:)
   integer :: unique_parents, merged_cnt, actual_parts_cnt(3), total_actual_parts(3)
   integer :: sibling_upper_limit, collision_checks
   real(kind_physics) :: merge_ratio,  local_min_x(3), local_max_x(3), min_x(3), max_x(3)
   type(t_box) :: bounding_box
   integer(kind_particle) :: last_merge_tnp

   ! buffer to record newly generated particles & related counters
   type(linked_list_elem), pointer :: buffer, particle_guide, slab_particles
   type(t_particle), allocatable   :: gathered_new_buffer(:)
   integer :: electron_num, i, j, new_particle_cnt, local_electron_num, &
              swapped_num, break, virtual_particle_cnt(2), H1s, H2s
   real(kind_physics) :: rank_charge_count(5), charge_count(5), &
                         total_charge_count(5), seed_dl(3)
   integer, allocatable :: new_particles_offset(:,:), generic_array(:)
   real(kind_physics), allocatable :: thread_charge_count(:,:)

   ! variables related to OMP
   integer :: init_omp_threads, omp_threads, thread_id, local_size, IStart, IStop, CStart, CStop

   integer, public :: MPI_TYPE_density

   ! density diagnostics related variables
   logical :: density_output
   type(diag_vertex), allocatable  :: density_verts(:), final_density(:)
   real(kind_physics) :: minimum_x, minimum_y, minimum_z, x_length, y_length, z_length
   real(kind_physics) :: dx, dy ,dz, s_min_x, s_min_y, s_min_z
   integer :: x_cell, y_cell, z_cell, iv, ir, cnt, mesh_mode,  N_element
   integer(kind_particle), allocatable :: connectivity_array(:), connectivity_tets(:,:)
   character(255) :: mesh_name, file_name, coil_data_file

   ! general diagnostics
   logical :: slice_parts
   real(kind_physics), allocatable :: local_table1D(:), global_table1D(:), local_ion_instance(:)
   real(kind_physics), allocatable :: local_table1D_1(:), global_table1D_1(:)
   real(kind_physics), allocatable :: local_table2(:,:), global_table2(:,:)
   real(kind_physics), allocatable :: glob_ion_instance(:), thread_ion_instance(:,:), rank_ion_instance(:)

   ! variables for random number generation
   integer :: dummy
   integer, allocatable :: r_seeds(:)
   integer(kind=int32) :: ctr_s(4), key_s(4)
   real(kind_physics):: rand_num(8), dummy_val(8)

   ! variables related to cross sections and probabilistic collisions
   real(kind_physics), allocatable :: Xi_table(:,:)
   integer :: total_cross_sections, eirene_cross_sections, coll_type, start_i
   real(kind_physics) :: abs_max_CS, neutral_density, init_temperature, pressure
   real(kind_physics), allocatable :: eirene_coeffs1(:), eirene_coeffs2(:)
   integer, allocatable :: outcomes(:)

   ! lookup tables for cross section data
   character(255) :: file_path
   type(linked_list_CS), pointer :: CS_tables, CS_guide, CS_total_scatter
   real(kind_physics), allocatable :: slopes(:)

   ! variables describing external fields
   real(kind_physics) :: d, major_radius, minor_radius, B0, B_p, V_loop, Itf
   real(kind_physics), allocatable :: B_pol_grid(:,:)

   ! checkpoint related variables
   character(255) :: checkpoint_file, i_wall_time
   integer :: checkin_step, resume, itime_in, sim_type, mode, flt_geom
   real(kind_physics) :: allowed_wall_time
   real(kind_physics) :: current_wall_time, prev_t_user_step = 0.0
   logical :: break_loop

   ! particle merging related variables
   integer(kind_key) :: dummy_key
   integer(kind_key), allocatable :: key_array(:)
   real(kind_physics), allocatable :: energy_group_levels(:)
   real(kind_physics) :: last_v(3), stored_vel(3)
   integer :: new_mass, old_part_cnt, tmp_buff_pos, stored_i, steps_since_last

   ! constants & scaling factors
   real(kind_physics), parameter :: c = 299792458.0_kind_physics ! m/s
   real(kind_physics), parameter :: e_mass = 510998.9461_kind_physics ! eV/c^2
   real(kind_physics), parameter :: e = 1.6021766208e-19_kind_physics ! Coulomb
   real(kind_physics), parameter :: eps_0 = 8.85418781762e-12_kind_physics ! Coulomb/Vm
   real(kind_physics), parameter :: mu_0 = 12.566370614e-7_kind_physics ! Newton/Ampere^2
   real(kind_physics), parameter :: pi = 3.141592653589793238462643383279502884197_kind_physics
   real(kind_physics), parameter :: kboltzmann = 1.38064852e-23_kind_physics ! J K^(-1)
   real(kind_physics) :: E_q_dt_m

   interface random
      module procedure random8, random16
   end interface

contains

   subroutine set_parameter()

      use module_pepc
      use module_interaction_specific, only: theta2, eps2, force_law, include_far_field_if_periodic
      implicit none

      integer, parameter :: fid = 12
      character(255)     :: para_file
      logical            :: read_para_file

      namelist /pepcbreakup/ resume, itime_in, init_omp_threads, i_wall_time, slice_parts, density_output, &
         mesh_mode, mesh_name, x_cell, y_cell, z_cell, minimum_x, minimum_y, minimum_z, x_length, &
         y_length, z_length, sim_type, mode, coil_data_file, Itf, d, electron_num, tnp, H1s, H2s, nt, dt, &
         particle_output, domain_output, particle_mpi_output, reflecting_walls, &
         particle_test, diag_interval, plasma_dimensions, init_temperature, &
         pressure, external_e, major_radius, minor_radius, B0, B_p, V_loop

      ! set default parameter values
      resume = 0
      itime_in = 0
      init_omp_threads = 1
      i_wall_time = '00:00:00'
      slice_parts = .false.
      density_output = .false.
      mesh_mode = 0
      mesh_name = './'
      x_cell = 0
      y_cell = 0
      z_cell = 0
      minimum_x = 0
      minimum_y = 0
      minimum_z = 0
      x_length = 0
      y_length = 0
      z_length = 0
      sim_type = 0
      mode = 0
      coil_data_file = './'
      Itf = 0.0
      d = 0.0015
      electron_num = 0
      tnp = 10000
      H1s = 0
      H2s = 0
      nt = 25
      dt = 1e-2
      particle_test = .false.
      particle_output = .false.
      domain_output = .false.
      particle_mpi_output = .true.
      reflecting_walls = .false.
      diag_interval = 1
      plasma_dimensions = (/1.0_8, 1.0_8, 1.0_8/)
      external_e = (/0.0_8, 0.0_8, 0.0_8/)
      electron_num = 0
      omp_threads = 1

      ! read in namelist file
      call pepc_read_parameters_from_first_argument(read_para_file, para_file)

      if (read_para_file) then
         if (root) write (*, '(a)') " == reading parameter file, section pepcbreakup: ", para_file
         open (fid, file=para_file)
         read (fid, NML=pepcbreakup)
         close (fid)
      else
         if (root) write (*, *) " == no param file, using default parameter "
      end if

      if (root) then
         write (*, '(a,i12)') " == resume from previous runs?          : ", resume
         write (*, '(a,i12)') " == resume from time step               : ", itime_in
         write (*, '(a,a8)')   " == allowed wall time                   : ", i_wall_time
         write (*, '(a,i12)') " == total number of particles           : ", tnp
         write (*, '(a,i12)') " == number of time steps                : ", nt
         write (*, '(a,es12.4)') " == time step                           : ", dt
         write (*, '(a,i12)') " == diag & IO interval                  : ", diag_interval
         write (*, '(a,l12)') " == particle test                       : ", particle_test
         write (*, '(a,l12)') " == particle output                     : ", particle_output
         write (*, '(a,l12)') " == density interpolation               : ", density_output
         write (*, '(a,i12)') " == mesh mode                           : ", mesh_mode
         write (*, '(a,a100)')   " == .msh filepath                       : ", mesh_name
         write (*, '(a,a100)')   " == coil filepath                       : ", coil_data_file
         write (*, '(a,l12)') " == domain output                       : ", domain_output
         write (*, '(a,l12)') " == particle mpi output                 : ", particle_mpi_output
         write (*, '(a,l12)') " == reflecting walls                    : ", reflecting_walls
         write (*, '(a,3(es12.4))') " == plasma dimensions                   : ", plasma_dimensions
         write (*, '(a,es12.4)') " == initial electron temperature(K)     : ", init_temperature
         write (*, '(a,es12.4)') " == pressure(Pa)                        : ", pressure
         write (*, '(a,3(es12.4))') " == external electric field(V/m)        : ", external_e
         write (*, '(a,es12.4)') " == major radius(m)                     : ", major_radius
         write (*, '(a,es12.4)') " == minor radius(m)                     : ", minor_radius
         write (*, '(a,es12.4)') " == toroidal magnetic field strength (T): ", B0
         write (*, '(a,es12.4)') " == poloidal magnetic field strength (T): ", B_p
         write (*, '(a,es12.4)') " == toroidal loop voltage(V)            : ", V_loop
      end if

      ! NOTE: Scale the appropriate read-in variables!
      !       1. multiply 4*pi*eps_0*(1e-12 sec*light speed)^2/(elementary charge) to electric field (V/m)
      !       2. multiply 4*pi*eps_0*(1e-12 sec*light speed)/(elementary charge) to voltage
      !       3. multiply ((1.e-12)*(c**2))/e_mass to B field (Tesla)
      !       4. divide length variables with (1.e-12*c)
      external_e = external_e*4.0*pi*eps_0*((c*1e-12)**2)/e
      V_loop = V_loop*4.0*pi*eps_0*((c*1e-12))/e
      B0 = B0 * ((1.e-12)*(c**2))/e_mass
      B_p = B_p * ((1.e-12)*(c**2))/e_mass
      plasma_dimensions = plasma_dimensions/(c*1e-12)
      major_radius = major_radius/(c*1e-12)
      minor_radius = minor_radius/(c*1e-12)

      read(i_wall_time(1:2),*)  i
      allowed_wall_time = i*3600.0_8
      read(i_wall_time(4:5),*)  i
      allowed_wall_time = allowed_wall_time + i*60.0_8
      read(i_wall_time(7:8),*)  i
      allowed_wall_time = allowed_wall_time + i*1.0_8

      if (sim_type == 0) then
        V_loop = 0.0
        B0 = 0.0
        B_p = 0.0
        Itf = 0.0
        flt_geom = 2
      else if (sim_type == 1) then
        external_e = 0.0
        electron_num = 0
        flt_geom = 3
      end if

      virtual_particle_cnt(1) = H1s
      virtual_particle_cnt(2) = H2s

      call pepc_prepare(3_kind_dim)
   end subroutine set_parameter

   function cross_product(vector1, vector2) result(vector_ans)
      implicit none
      real*8, intent(in) :: vector1(3), vector2(3)
      real*8 :: vector_ans(3)

      vector_ans(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
      vector_ans(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3)
      vector_ans(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)
   end function cross_product

   function Rodriguez_rotation(angle, rot_axis, vec_in) result(vec_out)
     implicit none
     real(kind_physics), intent(in), dimension(:) :: vec_in, rot_axis
     real(kind_physics), intent(in) :: angle
     real(kind_physics) :: vec_out(3), k_cross_v(3)
     real(kind_physics) :: k_dot_v, cos_angle
     ! rot_axis == k, vec_in  == v
     cos_angle = cos(angle)
     k_cross_v = cross_product(rot_axis, vec_in)
     k_dot_v = dot_product(rot_axis, vec_in) * (1. - cos_angle)

     vec_out = vec_in*cos_angle + k_cross_v * sin(angle) + rot_axis * k_dot_v
   end function Rodriguez_rotation

   ! Given the velocity vector of particle, its azimuthal and inclination angle is calculated.
   ! can be 'adapted' to calculate azimuthal position of particle in torus.
   ! Phi angle is aligned 0 deg from east, counterclock-wise. east coincides with +x-axis
   subroutine angles_calc(vec_in, vel_mag, theta, phi)
     implicit none
     real(kind_physics), dimension(:), intent(in):: vec_in
     real(kind_physics), intent(in) :: vel_mag
     real(kind_physics), intent(out) :: theta, phi

     if ((vec_in(1) /= 0.0) .and. (vec_in(2) /= 0.0)) then
       if ((vec_in(2) > 0.0) .and. (vec_in(1) < 0.0)) then
         phi = pi + atan(vec_in(2)/vec_in(1))
       else if ((vec_in(2) < 0.0) .and. (vec_in(1)> 0.0)) then
         phi = 2. * pi + atan(vec_in(2)/vec_in(1))
       else if ((vec_in(2) > 0.0) .and. (vec_in(1) > 0.0)) then
         phi = atan(vec_in(2)/vec_in(1))
       else
         phi = pi + atan(vec_in(2)/vec_in(1))
       end if
     else if (vec_in(1) == 0.0) then
       if (vec_in(2) > 0.0) then
         phi = 0.5*pi
       else
         phi = 1.5*pi
       end if
     else if (vec_in(2) == 0.0) then
       if (vec_in(1) > 0.0) then
         phi = 0.0
       else
         phi = pi
       end if
     end if

     theta = acos(vec_in(3)/vel_mag)
   end subroutine angles_calc

   function torus_geometry(mode) result(pos)
     implicit none
     real(kind_physics) :: pos(3), ran(3)
     integer, intent(in) :: mode
     real(kind_physics) :: theta, phi, l, r, dist, inplane_dist, x, y, z, ran_dist, dist_outer, dist_inner
     integer :: rep

     select case(mode)
     case(0) ! random torus distribution
       rep = 1

       do while (rep .eq. 1)
         dist_outer = major_radius + minor_radius
         dist_inner = major_radius - minor_radius

         ! Generating Random Number between [0,1]
         dummy = gen_norm_double_rng(ctr_s, key_s, rand_num)

         x = (2.*rand_num(1) - 1.) * dist_outer
         y = (2.*rand_num(2) - 1.) * dist_outer
         ran_dist = sqrt(x**2 + y**2)

         if ((ran_dist .LE. dist_outer) .and. (ran_dist .GE. dist_inner)) then
           z = (2.*rand_num(3) - 1.) * minor_radius
           inplane_dist = sqrt((ran_dist - major_radius)**2 + z**2)

           if (inplane_dist .LE. minor_radius) then
             pos(1) = x
             pos(2) = y
             pos(3) = z
             rep = 0
           else
             rep = 1
           end if
         else
           rep = 1
         end if
       end do

     case(1) ! plane distribution at phi = pi/2

       ! Generating Random Number between [0,1]
       dummy = gen_norm_double_rng(ctr_s, key_s, rand_num)

       r = rand_num(1) * minor_radius
       theta  = rand_num(3) * 2. * pi
       l = major_radius + r*sin(theta)
       pos(3) = r * cos(theta)

       phi = rand_num(2) * pi
       pos(1) = l
       pos(2) = 0.0

     case(2) ! at the major axis of torus
       l = 1.0

       ! Generating Random Number between [0,1]
       dummy = gen_norm_double_rng(ctr_s, key_s, rand_num)

       if (rand_num(2) > 0.5) then
          l = -1.0
       end if

       pos(1) = (2*rand_num(1) - 1.0) * major_radius
       pos(2) = l * sqrt(major_radius*major_radius - pos(1)*pos(1))
       pos(3) = 0.0

     end select

   end function torus_geometry

   function thermal_velocity_mag(mass, temp) result(velocity)
     implicit none
     real(kind_physics), intent(in) :: mass, temp
     real(kind_physics) :: velocity, kb_ev

     kb_ev = 8.61733035e-5_kind_physics !boltzmann constant in eV/K
     velocity = sqrt(8.*temp*kb_ev/(mass*e_mass*pi)) !dimensionless velocity
    !  print *, velocity
   end function thermal_velocity_mag

   subroutine init_particles(p, geom)
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      integer :: rc, geom
      real*8 :: dummy
      real(kind_physics) :: rand_scale, magnitude, ez(3), toroidal_vec(3)

      if (root) write (*, '(a)') " == [init] init particles "

      ! set initially number of local particles
      np = tnp/n_ranks
      if (my_rank < MOD(tnp, 1_kind_particle*n_ranks)) np = np + 1
      !1_kind_particle means a value of 1, enforced into specified precision defined by kind_particle

      allocate (p(np), stat=rc)
      if (rc .ne. 0) write (*, *) " === particle allocation error!"

      allocate (direct_L2(np), stat=rc)
      if (rc .ne. 0) write (*, *) " === direct_L2 allocation error!"
      direct_L2 = -1.0_8

      ! set random seed
      dummy = par_rand(1*my_rank)

      select case(geom)
      case(0)
        do ip = 1, np
           p(ip)%label = 0 !my_rank*(tnp/n_ranks) + ip - 1
           p(ip)%data%q = -1.0_8
           p(ip)%data%m = 1.0_8
           p(ip)%data%b = 0.0_8
           p(ip)%data%species = 0

           p(ip)%data%age = 0.0_8

          ! NOTE: square plane distribution
           call random(p(ip)%x)
           p(ip)%x(1) = p(ip)%x(1)*plasma_dimensions(1) ! p(ip)%x(1)*0.8*plasma_dimensions(1) - plasma_dimensions(1)*0.4
           p(ip)%x(2) = p(ip)%x(2)*plasma_dimensions(2) ! p(ip)%x(2)*0.8*plasma_dimensions(2) - plasma_dimensions(2)*0.4
           p(ip)%x(3) = p(ip)%x(3)*plasma_dimensions(3)
           ! p(ip)%x(3) = -0.01

           ! magnitude = MOD(real(ip), 100.0) !sqrt((2*100.0/e_mass))
           magnitude = sqrt((2*17.0/e_mass)) ! 0.0_kind_physics
           ! call random_number(rand_scale)
           rand_scale = 1.0_8
           p(ip)%data%v = 0.0_kind_physics
           p(ip)%data%v(1) = 1.0*magnitude*rand_scale
           p(ip)%data%f_b = 0.0_kind_physics
           p(ip)%data%f_e = 0.0_kind_physics
           p(ip)%data%mp_int1 = 0
           p(ip)%label = 0
           p(ip)%work = 1.0_8
        end do

      case(1)
        do ip = 1, np
           p(ip)%label = 0 !my_rank*(tnp/n_ranks) + ip - 1
           if (MOD(ip,2) .eq. 0) then
               p(ip)%data%q = 1.0_8
               p(ip)%data%m = 3673.43889456_8
               p(ip)%data%species = 2
           else
               p(ip)%data%q = -1.0_8
               p(ip)%data%m = 1.0_8
               p(ip)%data%species = 0
           end if
           p(ip)%data%age = 0.0_8

          ! NOTE: torus distribution
           p(ip)%x = torus_geometry(mode)

           ez = 0.0
           ez(3) = 1.0
           toroidal_vec(1) = p(ip)%x(2)*ez(3) - p(ip)%x(3)*ez(2)
           toroidal_vec(2) = p(ip)%x(3)*ez(1) - p(ip)%x(1)*ez(3)
           toroidal_vec(3) = p(ip)%x(1)*ez(2) - p(ip)%x(2)*ez(1)
           toroidal_vec = toroidal_vec/sqrt(dot_product(toroidal_vec, toroidal_vec))

           magnitude = sqrt((2*0.02/e_mass))
           call random_number(rand_scale)
           rand_scale = 1.0_8
           p(ip)%data%v = -1.0*magnitude*toroidal_vec*rand_scale
           ! p(ip)%data%v = 0.0_8
           p(ip)%data%b = 0.0_kind_physics
           p(ip)%data%f_b = 0.0_kind_physics
           p(ip)%data%f_e = 0.0_kind_physics
           p(ip)%data%mp_int1 = 0
           p(ip)%label = 0

           p(ip)%work = 1.0_8
        end do
      end select
   end subroutine init_particles

   subroutine push_particles(p)
      use module_mirror_boxes
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      real*8  :: fact

      if (root) write (*, '(a)') " == [pusher] push particles "

      fact = dt

      do ip = 1, np
         p(ip)%data%v = p(ip)%data%v + fact*p(ip)%data%q/p(ip)%data%m*p(ip)%results%e
         p(ip)%x = p(ip)%x + dt*p(ip)%data%v
      end do
   end subroutine push_particles

   subroutine filter_particles(p)
      use mpi
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      integer :: id, ncoll, ncoll_total, ierr

      ncoll = 0

      do ip = 1, np
         do id = 1, 3
            if (p(ip)%x(id) < 0.0_8) then
               p(ip)%x(id) = -p(ip)%x(id)
               p(ip)%data%v(id) = -p(ip)%data%v(id)
               ncoll = ncoll + 1
            else if (p(ip)%x(id) > plasma_dimensions(id)) then
               p(ip)%x(id) = 2*plasma_dimensions(id) - p(ip)%x(id)
               p(ip)%data%v(id) = -p(ip)%data%v(id)
               ncoll = ncoll + 1
            end if
         end do
      end do

      call mpi_reduce(ncoll, ncoll_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) write (*, '(a,i12)') " == [filter] total number of wall collisions      : ", ncoll_total
   end subroutine filter_particles

   subroutine test_particles()
      use module_pepc_types
      use module_directsum
      use mpi
      implicit none

      integer(kind_particle), allocatable   :: tindx(:)
      real*8, allocatable                   :: trnd(:)
      type(t_particle_results), allocatable :: trslt(:)
      integer(kind_particle)                :: tn, tn_global, ti
      integer                               :: rc
      real(kind_physics)                    :: L2sum_local, L2sum_global, L2

      call timer_start(t_user_directsum)

      if (allocated(direct_L2)) then
         deallocate (direct_L2)
      end if
      allocate (direct_L2(np))
      direct_L2 = -1.0_8

      if (particle_direct .eq. -1) then
         tn = np
      else
         tn = particle_direct/n_ranks
         if (my_rank .eq. (n_ranks - 1)) tn = tn + MOD(particle_direct, n_ranks)
      endif

      allocate (tindx(tn), trnd(tn), trslt(tn))

      if (particle_direct .eq. -1) then
         do ti = 1, tn
            tindx(ti) = ti
         enddo
      else
         call random(trnd(1:tn))

         tindx(1:tn) = int(trnd(1:tn)*(np - 1)) + 1
      endif

      call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

      L2sum_local = 0.0
      L2sum_global = 0.0
      do ti = 1, tn
         L2 = &
            (particles(tindx(ti))%results%e(1) - trslt(ti)%e(1))**2 + &
            (particles(tindx(ti))%results%e(2) - trslt(ti)%e(2))**2 + &
            (particles(tindx(ti))%results%e(3) - trslt(ti)%e(3))**2
         L2sum_local = L2sum_local + L2
         direct_L2(tindx(ti)) = L2
      end do

      call MPI_ALLREDUCE(tn, tn_global, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(L2sum_local, L2sum_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      L2sum_global = sqrt(L2sum_global)/tn_global

      call timer_stop(t_user_directsum)
      if (root) then
         write (*, '(a,i12)') " == [direct test] number tested particles         : ", tn_global
         write (*, '(a,es12.4)') " == [direct test] L2 error in probed particles    : ", L2sum_global
         write (*, '(a,es12.4)') " == [direct test] time in test [s]                : ", timer_read(t_user_directsum)
      end if

      deallocate (tindx)
      deallocate (trnd)
      deallocate (trslt)
   end subroutine test_particles

   integer function vtk_step_of_step(step) result(vtk_step)
      use module_vtk
      implicit none

      integer, intent(in) :: step

      if (step .eq. 0) then
         vtk_step = VTK_STEP_FIRST
      else if (nt - 1 - step < diag_interval) then
         vtk_step = VTK_STEP_LAST
      else
         vtk_step = VTK_STEP_NORMAL
      endif
   end function vtk_step_of_step

   subroutine write_text_output(anode_charge_count, cathode_charge_count, step)
     implicit none
     real(kind_physics), intent(in) :: anode_charge_count, cathode_charge_count
     integer, intent(in) :: step

    !  open(12, file = fname, action='WRITE')
     if (step == 0) then
       write(12, *) 'step ', 'anode_count ', 'cathode_count'
     end if
     write(12, *) step, anode_charge_count, cathode_charge_count
    !  close(12)
   end subroutine write_text_output

   subroutine vtk_write_densities(fname, step, time, vtk_step, vertices, &
                                  ncell, npart, offset, vtk_type, helper_func, coord_scale)
     use module_vtk
     use module_pepc_types
     use module_interaction_specific_types
     use mpi
     implicit none

     character(*), intent(in) :: fname
     integer, intent(in) :: step, vtk_step, vtk_type, offset
     real*8, intent(in) :: time
     type(diag_vertex), intent(in) :: vertices(:)
     integer(kind_particle), intent(in) :: ncell, npart

     interface
       subroutine helper_func(d, r, vtkf)
         use module_interaction_specific_types, only: t_particle_data, t_particle_results
         use module_vtk, only: vtkfile_unstructured_grid
         implicit none
         type(t_particle_data), intent(in) :: d(:)
         type(t_particle_results), intent(in) :: r(:)
         type(vtkfile_unstructured_grid), intent(inout) :: vtkf
       end subroutine helper_func
     end interface

     optional :: helper_func
     real(kind_physics), intent(in), optional :: coord_scale

     integer(kind_particle) :: i
     integer(kind_pe) :: mpi_rank, mpi_size
     integer(kind_default) :: ierr
     type(vtkfile_unstructured_grid) :: vtk

     call vtk%create(fname, step, time, vtk_step)
       call vtk%write_headers(npart, ncell)
       call vtk%startpoints()
         call vtk%write_data_array("xyz", vertices(:)%x(1), vertices(:)%x(2), vertices(:)%x(3), coord_scale)
       call vtk%finishpoints()
       call vtk%startpointdata()
         call vtk%write_data_array("electron density", vertices(:)%q_density(1))
         call vtk%write_data_array("ion density", vertices(:)%q_density(2))
         call vtk%write_data_array("current density", vertices(:)%J_density(1),&
                                    vertices(:)%J_density(2), vertices(:)%J_density(3))
       call vtk%finishpointdata()
       call vtk%startcells()
         call vtk%write_data_array("connectivity", connectivity_array)
         call vtk%write_data_array("offsets", [((i*offset),i=1,ncell)])
         call vtk%write_data_array("types", [(vtk_type,i=1,ncell)])
       call vtk%finishcells()
       print *, "Done writing Cell"
       call vtk%startcelldata()
        ! no cell data here as cells correspond to points anyway, in case of problems use PointDataToCellData Filter in Paraview
       call vtk%finishcelldata()
       call vtk%write_final()
     call vtk%close()
   end subroutine vtk_write_densities

   subroutine write_densities(vertices, mode)
     use module_vtk
     use mpi
     implicit none

     type(diag_vertex), intent(in) :: vertices(:)
     integer, intent(in) :: mode
     integer :: vtk_step, offset
     integer(kind_particle) :: ncell, npart

     vtk_step = vtk_step_of_step(step)
     if (mode == 0) then
       ncell = x_cell*y_cell*z_cell
       npart = (x_cell + 1)*(y_cell + 1)*(z_cell + 1)
       offset = 8
       call vtk_write_densities("densities", step, dt*step, vtk_step, vertices, &
                                ncell, npart, offset, VTK_VOXEL)
     elseif (mode == 1) then
       ncell = N_element
       npart = size(vertices)
       offset = 4
       call vtk_write_densities("densities", step, dt*step, vtk_step, vertices, &
                                ncell, npart, offset, VTK_TETRA)
     end if
     ! call vtk_write_densities("densities", step, dt*step, vtk_step, vertices)
   end subroutine write_densities

   subroutine write_updated_resume_variables(resume_step)
     implicit none
     integer, intent(in) :: resume_step

     open(23, file = 'update_variable.txt', action='WRITE')
     write(23, '(i1)') 1
     write(23, '(i10)') resume_step
     write(23, '(i10)') tnp
     write(23, '(i10)') virtual_particle_cnt(1)
     write(23, '(i10)') virtual_particle_cnt(2)
     close(23)
   end subroutine write_updated_resume_variables

   subroutine preempt_checkpointing(current_wall_time, previous_step_duration, doDiag, particles, step, break)
     use module_checkpoint
     use mpi
     implicit none

     real(kind_physics), intent(in) :: current_wall_time, previous_step_duration
     logical, intent(in) :: doDiag
     type(t_particle), allocatable, intent(in) :: particles(:)
     integer, intent(in) ::step
     real(kind_physics) :: difference_to_wall_time, buffer_multiplier, time_check
     logical, intent(out) :: break

     break = .false.
     buffer_multiplier = 20.0

     difference_to_wall_time = allowed_wall_time - current_wall_time
     time_check = buffer_multiplier*previous_step_duration
     if (my_rank == 0) print *, 'Time left till wall time (seconds): ', difference_to_wall_time

     if (difference_to_wall_time < time_check) then
       if (my_rank == 0) then
         print *, "Pre-empted Checkpointing triggered!"
       end if
       if(.not. doDiag) then
         call write_particles(particles)
         if (my_rank == 0) call write_updated_resume_variables(step+itime_in+1)
         call MPI_BCAST(tnp, 1, MPI_KIND_PARTICLE, 0, MPI_COMM_WORLD, ierr)
         call write_particles_mpiio(MPI_COMM_WORLD, step+itime_in+1, tnp, particles, checkpoint_file)
       end if
       break = .true.
     end if
   end subroutine preempt_checkpointing

   subroutine write_particles(p)
      use module_vtk_helpers
      use mpi
      implicit none

      type(t_particle), intent(in) :: p(:)

      integer :: vtk_step, temp_step

      temp_step = (step + 1)/1000

      call timer_start(t_user_particleio)
      vtk_step = vtk_step_of_step(temp_step)
      call vtk_write_particles("particles", MPI_COMM_WORLD, temp_step, dt*step, vtk_step, p, coulomb_and_l2)
      call timer_stop(t_user_particleio)
      if (root) write (*, '(a,es12.4)') " == [write particles] time in vtk output [s]      : ", timer_read(t_user_particleio)

   contains

      subroutine coulomb_and_l2(d, r, vtkf)
         use module_vtk
         use module_interaction_specific_types
         implicit none

         type(t_particle_data), intent(in) :: d(:)
         type(t_particle_results), intent(in) :: r(:)
         type(vtkfile_unstructured_grid), intent(inout) :: vtkf

         call vtk_write_particle_data_results(d, r, vtkf)
         if (particle_test) call vtkf%write_data_array("L2 error", direct_L2(:))
      end subroutine
   end subroutine write_particles

   subroutine write_domain(p)
      use module_vtk
      use module_vtk_helpers
      use module_pepc, only: global_tree
      implicit none

      type(t_particle), allocatable, intent(in) :: p(:)

      integer :: vtk_step

      ! output of tree diagnostics
      vtk_step = vtk_step_of_step(step)
      call vtk_write_branches(step, dt*step, vtk_step, global_tree)
      call vtk_write_leaves(step, dt*step, vtk_step, global_tree)
      call vtk_write_spacecurve(step, dt*step, vtk_step, p)
   end subroutine write_domain

   subroutine random_gauss(list)
      implicit none

      real(kind_physics), intent(inout) :: list(:)

      real(kind_physics) :: v(2), pi, r, p
      integer :: n, i

      pi = 2.0_kind_physics*acos(0.0_kind_physics)
      n = size(list)

      do i = 1, n, 2

         call random(v)

         r = sqrt(-2.0_8*log(v(1)))
         p = 2.0_8*pi*v(2)

         list(i) = r*sin(p)
         if ((i + 1) <= n) list(i + 1) = r*cos(p)

      end do
   end subroutine

   subroutine random8(array)
      implicit none
      real*8 :: array(:)
      integer :: i

      do i = 1, size(array)
         array(i) = par_rand()
      end do
   end subroutine random8

   subroutine random16(array)
      implicit none
      real*16 :: array(:)
      integer :: i

      do i = 1, size(array)
         array(i) = par_rand()
      end do
   end subroutine random16

   !>
   !> portable random number generator, see numerical recipes
   !> check for the random numbers:
   !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
   !> the parameter iseed is optional
   !>
   function par_rand(iseed)
      implicit none
      real :: par_rand
      integer, intent(in), optional :: iseed

      integer, parameter :: IM1 = 2147483563
      integer, parameter :: IM2 = 2147483399
      real, parameter :: AM = 1.0/IM1
      integer, parameter :: IMM1 = IM1 - 1
      integer, parameter :: IA1 = 40014
      integer, parameter :: IA2 = 40692
      integer, parameter :: IQ1 = 53668
      integer, parameter :: IQ2 = 52774
      integer, parameter :: IR1 = 12211
      integer, parameter :: IR2 = 3791
      integer, parameter :: NTAB = 32
      integer, parameter :: NDIV = 1 + IMM1/NTAB
      real, parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
      real, parameter :: RNMX = 1.0 - eps_

      integer :: j, k
      integer, volatile, save :: idum = -1
      integer, volatile, save :: idum2 = 123456789
      integer, volatile, save :: iy = 0
      integer, volatile, save :: iv(NTAB)

      if (idum <= 0 .or. present(iseed)) then
         if (present(iseed)) then
            idum = iseed
         else
            if (-idum < 1) then
               idum = 1
            else
               idum = -idum
            endif
         endif

         idum2 = idum

         do j = NTAB + 7, 0, -1
            k = idum/IQ1
            idum = IA1*(idum - k*IQ1) - k*IR1
            if (idum < 0) idum = idum + IM1

            if (j < NTAB) iv(j + 1) = idum

         end do
         iy = iv(1)
      end if

      k = idum/IQ1
      idum = IA1*(idum - k*IQ1) - k*IR1
      if (idum < 0) idum = idum + IM1

      k = idum2/IQ2
      idum2 = IA2*(idum2 - k*IQ2) - k*IR2
      if (idum2 < 0) idum2 = idum2 + IM2

      j = iy/NDIV + 1
      iy = iv(j) - idum2
      iv(j) = idum

      if (iy < 1) iy = iy + IMM1
      par_rand = AM*iy
      if (par_rand > RNMX) par_rand = RNMX
   end function par_rand
end module
