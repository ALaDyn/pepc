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
!> interactions module
!>

module interactions_integrator
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use helper
   use particles_resize
   use rng_wrapper
   use elliptic_integrals
   implicit none

contains
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

   function matrix_vector_multiplication(mat, vec) result(ans)
      implicit none
      real*8, allocatable, intent(in) :: mat(:, :), vec(:)
      real*8 :: ans(size(mat, 1))
      integer :: i, j

      ans = 0.0_8
      do i = 1, size(mat, 1)
         do j = 1, size(mat, 2)
            ans(i) = ans(i) + mat(i, j)*vec(j)
         end do
      end do
   end function matrix_vector_multiplication

   function current_loop_B(particle_pos, r_pos, coil_z_origin, coil_radius, I) &
     result(B_out)
     implicit none
     real(kind_physics), intent(in) :: particle_pos(3), r_pos, coil_z_origin
     real(kind_physics), intent(in) :: coil_radius, I
     real(kind_physics) :: r_vec(2), Br, Bz, Bphi, B_out(3), z
     real(kind_physics) :: k, Kk, Ek, r_pos2, R2, z2, denom1, denom2, prefac

     z = particle_pos(3) - coil_z_origin
     z2 = z*z
     r_pos2 = r_pos*r_pos
     R2 = coil_radius*coil_radius

     denom1 = sqrt(z2 + (coil_radius + r_pos)**2)
     denom2 = z2 + (r_pos - coil_radius)**2
     prefac = mu_0*I/(2*pi*denom1)

     k = sqrt(4*r_pos*coil_radius/(denom1*denom1))

     Kk = elliptic_fk(k)
     Ek = elliptic_ek(k)

     Br = ((z2 + r_pos2 + R2)*Ek/denom2 - Kk)*prefac*z/r_pos
     Bz = prefac*((R2 - z2 - r_pos2)*Ek/denom2 + Kk)

     r_vec = particle_pos(1:2)
     r_vec = r_vec/r_pos

     B_out(1:2) = Br*r_vec
     B_out(3) = Bz

     B_out = B_out*((1.e-12)*(c**2))/e_mass ! non-dimensionalise
   end function current_loop_B

   subroutine particle_EB_field(particle, E_field)
      implicit none
      type(t_particle), intent(inout) :: particle
      real(kind_physics), intent(in) :: E_field(3)
      real(kind_physics) :: Coord(3), B_field(3), ez(3), Pol_B_field(3), field_vector(3)
      real(kind_physics) :: R, PF_final(3), PF_temp(3), Itf, Ipf, PF_unit_vector(3)

      ez = 0.0
      ez(3) = 1.0
      B_field = 0.0
      field_vector = 0.0
      Coord = particle%x*(c*1e-12)

      field_vector = cross_product(particle%x,ez)
      field_vector = field_vector/sqrt(dot_product(field_vector, field_vector))

      Itf = 7.524e7_8 !3.62e7_8
      R = sqrt(Coord(1)**2 + Coord(2)**2)
      B_field = field_vector * mu_0 * Itf/(2.0_8 * pi *R) ! field_vector*B0*major_radius/R
      B_field = B_field*((1.e-12)*(c**2))/e_mass ! non-dimensionalise

      PF_final = 0.0
      Ipf = 3.0e6_8

      ! 4 correction coils surrounding breakdown region, centred around major_radius
      ! Coils are directly solving Biot-Savart Law
      PF_temp = current_loop_B(Coord, R, 3.0_8, 5.8_8, 2357.3662872100125_8) !(JET-like: 35943.1422902196_8) ! 4.0e4_8) !3.0e6_8) !
      PF_final = PF_final + PF_temp

      PF_temp = current_loop_B(Coord, R, -3.0_8, 5.8_8, 2357.3662872100125_8) !(JET-like: 35943.1422902196_8) !4.0e4_8) ! 3.0e6_8) !
      PF_final = PF_final + PF_temp

      PF_temp = current_loop_B(Coord, R, 0.0_8, 3.0_8, 1.0e5_8)!(JET-like: 1.77e5_8) !1.5e5_8) ! 7.0e6_8) !
      PF_final = PF_final + PF_temp

      PF_temp = current_loop_B(Coord, R, 0.0_8, 8.3_8, 1.5e4_8)!(JET-like: 1.77e4_8)!11467.320098924624_8) ! 115453.38348883369_8) !
      PF_final = PF_final + PF_temp

      ! PF_unit_vector = PF_final/sqrt(dot_product(PF_final, PF_final))
      ! PF_temp = 0.0005 * ((1.e-12)*(c**2))/e_mass * PF_unit_vector
      Pol_B_field = PF_final! + PF_temp

      R = R/(c*1e-12)
      particle%results%e = particle%results%e + field_vector*V_loop/(2.*pi*R) + E_field
      particle%data%b = 0.0_kind_physics !Pol_B_field + B_field
   end subroutine particle_EB_field

   subroutine test_ionization(particle, guide, new_particle, electron_count)
      implicit none
      type(t_particle), intent(inout) :: particle
      type(linked_list_elem), pointer, intent(inout) :: guide
      integer, intent(inout) :: new_particle, electron_count
      integer :: buffer_pos, ll_elem_gen

      if (particle%data%age > 6000) then
         ll_elem_gen = MOD(new_particle, size(guide%tmp_particles))
         buffer_pos = ll_elem_gen + 1
         ! print *, "Buffer_pos: ", buffer_pos, size(guide%tmp_particles)
         if (ll_elem_gen == 0 .and. new_particle > 0) then
            ! print *, "Extending Temp_array!!!!!"
            call allocate_ll_buffer(size(guide%tmp_particles), guide%next)
            guide => guide%next
         end if

         guide%tmp_particles(buffer_pos)%x(1) = particle%x(1)
         guide%tmp_particles(buffer_pos)%x(2) = particle%x(2)
         guide%tmp_particles(buffer_pos)%x(3) = particle%x(3)
         guide%tmp_particles(buffer_pos)%work = 1.0_8 !particle%work
         guide%tmp_particles(buffer_pos)%data%q = 1.0

         guide%tmp_particles(buffer_pos)%data%v = particle%data%v*0.005

         guide%tmp_particles(buffer_pos)%data%m = 100.0
         guide%tmp_particles(buffer_pos)%data%age = 0.0
         guide%tmp_particles(buffer_pos)%results%e = 0.0
         guide%tmp_particles(buffer_pos)%results%pot = 0.0

         particle%data%v = particle%data%v*0.5
         particle%data%age = 0.0

         new_particle = new_particle + 1
      end if
   end subroutine test_ionization

   subroutine boris_velocity_update(particle, dt)
      implicit none
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: dt
      real*8 :: Vm(3), Vd(3), Vp(3), tan_w(3), sin_w(3), V_cross(3), V_init(3)
      real*8 :: half_dt, q_m_ratio
      V_init = particle%data%v
      half_dt = dt*0.5_8

      q_m_ratio = particle%data%q*half_dt/particle%data%m

      Vm = particle%data%v + E_q_dt_m*particle%results%e*q_m_ratio
      tan_w = particle%data%b*q_m_ratio
      sin_w = 2.0_8*tan_w/(1 + dot_product(tan_w, tan_w))

      V_cross = cross_product(Vm, tan_w)
      Vd = Vm + V_cross
      V_cross = cross_product(Vd, sin_w)
      Vp = Vm + V_cross

      particle%data%v = Vp + E_q_dt_m*particle%results%e*q_m_ratio

      particle%data%f_e = 2.0*E_q_dt_m*particle%results%e*q_m_ratio*particle%data%m*e_mass/(c*dt*1e-12)
      !particle%data%q*particle%results%e/particle%data%m
      !particle%data%f_b = particle%data%m*e_mass*(particle%data%v - V_init)/(c*dt*1e-12) - particle%data%f_e

      ! To isolate the Delta V due to E field and B field, propose to compute
      ! D_V = V_final - V_init
      ! D_Ve = 2.0*E_q_dt_m*particle%results%e*q_m_ratio
      ! D_Vb = D_V - D_Ve
   end subroutine boris_velocity_update

   subroutine particle_pusher(particle, dt)
      implicit none
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: dt
      real*8 :: dist_vec(3)

      dist_vec = particle%data%v*dt
      particle%x = particle%x + dist_vec
      particle%data%f_b(1) = particle%data%f_b(1) + sqrt(dot_product(dist_vec,dist_vec))
      particle%data%f_b(2) = 0.0_kind_physics
      particle%data%f_b(3) = 0.0_kind_physics

      particle%data%age = particle%data%age + dt
   end subroutine particle_pusher

   function calculate_neutral_density(pressure, temperature) result(density)
     implicit none
     real(kind_physics), intent(in) :: pressure, temperature
     real(kind_physics) :: density

     density = pressure/(temperature * kboltzmann)

     ! Dimensionless scaling
     density = density * (1e-12 * c)**3

   end function calculate_neutral_density

   subroutine set_cross_section_table(fname, guide_CS, file_id, is_last)
     implicit none
     type(linked_list_CS), pointer, intent(inout) :: guide_CS
     character(len = *), intent(in) :: fname
     integer, intent(in) :: file_id, is_last
     integer :: entries, i

     open(file_id,file=fname,action='READ')
     read(file_id,*) ! Skipping first line
     read(file_id,*) entries ! Number of entries in data

     allocate(guide_CS%CS(entries,2))

     do i = 1,entries
       read(file_id,*) guide_CS%CS(i,1), guide_CS%CS(i,2)
     end do

     ! scale table(i,2) to become dimensionless, scaling is 1m = non_dim_l * (1e-12 * c)
     guide_CS%CS(:,2) = guide_CS%CS(:,2)*10000.0/(c**2)

     if (is_last == 0) then
       allocate(guide_CS%next_CS)
       guide_CS => guide_CS%next_CS
      !  nullify(guide_CS%next_CS)
     end if

     nullify(guide_CS%next_CS)
     close(file_id)
   end subroutine set_cross_section_table

   subroutine set_eirene_coeffs(fname, file_id, coeffs)
     implicit none
     real(kind_physics), allocatable, intent(inout) :: coeffs(:)
     real(kind_physics) :: emin, emax
     character(len = *), intent(in) :: fname
     integer, intent(in) :: file_id
     integer :: entries, i

     open(file_id,file=fname,action='READ')
     read(file_id,*) ! Skipping first line
     read(file_id,*) entries, emin, emax ! Number of entries in data

     allocate(coeffs(entries+2))
     coeffs(1) = emin
     coeffs(2) = emax

     do i = 3, entries+2
       read(file_id,*) coeffs(i)
     end do
     close(file_id)
   end subroutine set_eirene_coeffs

   subroutine Eirene_fit(coeffs, E, cross_sec, CS_index)
     real(kind_physics), allocatable, intent(in) :: coeffs(:)
     real(kind_physics), dimension(:), intent(inout) :: cross_sec
     real(kind_physics), intent(in) :: E
     real(kind_physics) :: ln_E, ln_sigma, Emin, Emax
     integer, intent(inout) :: CS_index
     integer :: i, n, starting_index
     ! E in eV, Emin & Emax describe the valid region. Any E outside of the range
     ! gives 0.0 cross section value. Note that the computed sigma by Eirene
     ! is at the unit of cm^2. Remember to non-dimensionalise!
     ln_sigma = 0.0
     ln_E = log(E)
     n = size(coeffs)
     Emin = coeffs(1)
     Emax = coeffs(2)
     starting_index = 3

     if ((E < Emin) .or. (E > Emax)) then
       cross_sec(CS_index) = 0.0
     else
       do i = starting_index, n
         ln_sigma = ln_sigma + coeffs(i)*ln_E**(i-starting_index)
       end do
       ln_sigma = ln_sigma - log(1e-16)
       cross_sec(CS_index) = cross_sec(CS_index-1) + exp(ln_sigma)*10000.0/(c**2) ! non-dimensionalise.
     end if
     CS_index = CS_index + 1
   end subroutine Eirene_fit

   subroutine determine_absolute_max_CS(guide_CS, max_CS)
     implicit none
     type(linked_list_CS), pointer, intent(in) :: guide_CS
     real(kind_physics), intent(inout) :: max_CS
     type(linked_list_CS), pointer :: temp_guide

     temp_guide => guide_CS
     max_CS = maxval(temp_guide%CS(:,2))
    !  max_CS = 0.0
    !  temp_guide => guide_CS
    !  do while (associated(temp_guide))
    !    max_CS = max_CS + maxval(temp_guide%CS(:,2))
    !    temp_guide => temp_guide%next_CS
    !  end do
   end subroutine determine_absolute_max_CS

   subroutine determine_cross_sections(particle, sigma_vec, guide_CS)
     implicit none
     type(t_particle), intent(in) :: particle
     type(linked_list_CS), pointer, intent(in) :: guide_CS
     real(kind_physics), dimension(:), intent(inout) :: sigma_vec
     integer :: i, table_rows, CS_table_no
     real(kind_physics) :: energy, scaled_mass
     type(linked_list_CS), pointer :: temp_guide

     scaled_mass = particle%data%m/abs(particle%data%q)

     CS_table_no = 1
     temp_guide => guide_CS
     do while (associated(temp_guide))
       table_rows = size(temp_guide%CS,1)
       energy = 0.5 * scaled_mass * e_mass * dot_product(particle%data%v,particle%data%v)
      !  print *, "energy: ", energy, dot_product(particle%data%v,particle%data%v)

       i = 1
       ! if energy is less/more than range of temp_guide%CS, return 0.0
       if (energy < temp_guide%CS(i,1) .or. energy > temp_guide%CS(table_rows,1)) then
         if (CS_table_no .ne. 1) then
           sigma_vec(CS_table_no) = sigma_vec(CS_table_no - 1)
         else
           sigma_vec(CS_table_no) = 0.0_kind_physics
         end if
       else
         do while (energy > temp_guide%CS(i,1))
           i = i + 1 ! this will return the 'i' that is immediately higher than 'energy'
         end do

         if (CS_table_no .ne. 1) then
           sigma_vec(CS_table_no) = (temp_guide%CS(i,2) - temp_guide%CS(i-1,2))/(temp_guide%CS(i,1) - temp_guide%CS(i-1,1)) &
                                * (energy-temp_guide%CS(i-1,1)) + temp_guide%CS(i-1,2) + sigma_vec(CS_table_no - 1)
         else
           sigma_vec(CS_table_no) = (temp_guide%CS(i,2) - temp_guide%CS(i-1,2))/(temp_guide%CS(i,1) - temp_guide%CS(i-1,1)) &
                                * (energy-temp_guide%CS(i-1,1)) + temp_guide%CS(i-1,2)
         end if
         ! Subsequent values in Sigma vector is a cumulation of previous cross sections
         ! multiplying each element in the sigma vector with local_density & velocity magnitude
         ! will give collision frequency.
       end if
       temp_guide => temp_guide%next_CS
       CS_table_no = CS_table_no + 1
     end do
   end subroutine determine_cross_sections

   subroutine set_Xi_table(fname, file_id, Xi_t)
     implicit none
     character(len = *), intent(in) :: fname
     real(kind_physics), intent(inout), allocatable :: Xi_t(:,:)
     integer, intent(in) :: file_id
     integer :: entries, i

     open(file_id,file=fname,action='READ')
     read(file_id,*) ! Skipping first line
     read(file_id,*) entries ! Number of entries in data

     allocate(Xi_t(entries,2))

     do i = 1,entries
       read(file_id,*) Xi_t(i,1), Xi_t(i,2)
     end do
     close(file_id)
   end subroutine set_Xi_table

   subroutine determine_xi(KE, xi)
     implicit none
     real(kind_physics), intent(in) :: KE
     real(kind_physics), intent(out) :: xi
     integer :: i

     i = 1
     do while (KE .gt. Xi_table(i,1))
       i = i + 1
     end do

     xi = (Xi_table(i,2) - Xi_table(i-1,2))/(Xi_table(i,1) - Xi_table(i-1,1)) &
          * (KE - Xi_table(i-1,1)) + Xi_table(i-1,2)
   end subroutine determine_xi

   subroutine momentum_partition_merging(particles, sibling_cnt, sibling_upper_limit, &
                                        parent_no, merged_guide, merged_cnt)
     ! NOTE: make sure the particle list is sorted by species within the sibling extent before using this subroutine!
     ! Quadrants  : 1 - +x +y +z
     !            : 2 - +x +y -z
     !            : 3 - +x -y +z
     !            : 4 - +x -y -z
     !            : 5 - -x +y +z
     !            : 6 - -x +y -z
     !            : 7 - -x -y +z
     !            : 8 - -x -y -z
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, allocatable, intent(in) :: sibling_cnt(:)
     integer, intent(in) :: sibling_upper_limit
     integer, intent(in) :: parent_no
     type(linked_list_elem), pointer, intent(inout) :: merged_guide
     integer, intent(inout) :: merged_cnt
     integer :: ll_elem_gen, istart, istop, i, j, species, direction, direction_cnt(8), buffer_pos
     type(t_particle), allocatable :: directional_buffer(:,:)

     allocate(directional_buffer(8,sibling_upper_limit))

     istart = 1
     do i = 2, parent_no
       istart = istart + sibling_cnt(i-1)
     end do
     if (parent_no .eq. 1) istart = 1
     istop = istart + sibling_cnt(parent_no) - 1

     ! count the number of particles of the same species, in each of the 8 directions.
     ! merging subroutine is called once for every species.
     species = 0
     direction_cnt = 0
     do i = istart, istop
       if (particles(i)%data%species .ne. species) then
         ! print *, "Different species detected: ", particles(i)%data%species, direction_cnt
         do j = 1, 8
           if (direction_cnt(j) .ne. 0) then
             call energy_elastic_merging(j, directional_buffer, direction_cnt(j), merged_guide, &
                                  merged_cnt, species, particles(i)%data%mp_int1)
           end if
         end do
         species = particles(i)%data%species
         direction_cnt = 0
       end if

       ! NOTE: momentum direction filtering
       if (particles(i)%data%v(1) > 0.0_kind_physics) then
         if (particles(i)%data%v(2) > 0.0_kind_physics) then
           if (particles(i)%data%v(3) > 0.0_kind_physics) then
             direction = 1
           else
             direction = 2
           end if
         else
           if (particles(i)%data%v(3) > 0.0_kind_physics) then
             direction = 3
           else
             direction = 4
           end if
         end if
       else
         if (particles(i)%data%v(2) > 0.0_kind_physics) then
           if (particles(i)%data%v(3) > 0.0_kind_physics) then
             direction = 5
           else
             direction = 6
           end if
         else
           if (particles(i)%data%v(3) > 0.0_kind_physics) then
             direction = 7
           else
             direction = 8
           end if
         end if
       end if

       direction_cnt(direction) = direction_cnt(direction) + 1
       directional_buffer(direction,direction_cnt(direction)) = particles(i)
     end do
     ! print *, "Final species: ", particles(i)%data%species, direction_cnt
     ! Final merge for the last species
     do i = 1, 8
       if (direction_cnt(i) .ne. 0) then
         call energy_elastic_merging(i, directional_buffer, direction_cnt(i), merged_guide, &
                              merged_cnt, particles(istop)%data%species, particles(istart)%data%mp_int1)
       end if
     end do
     ! print *, "Parent no: ", parent_no, istart, istop, "Finished merging.", size(particles)

     deallocate(directional_buffer)
   end subroutine momentum_partition_merging

   subroutine age_elastic_merging(direction, directional_buffer, direction_cnt, merged_guide, merged_cnt, species, parent_key)
     !NOTE: in order to not merge particles that are .le. 2 particles lying in the same momentum partition,
     !       most direct way of doing it is to copy the particles into respective directional buffer.
     !       don't sum them up yet! information will be lost if done so.
     implicit none
     integer, intent(in) :: direction, direction_cnt, species, parent_key
     type(t_particle), allocatable, intent(in) :: directional_buffer(:,:)
     type(linked_list_elem), pointer, intent(inout) :: merged_guide
     integer, intent(inout) :: merged_cnt
     integer :: i, j, buffer_pos, ll_elem_gen, filled, filled_ion, m_i, remainder, &
                extent
     type(t_particle), allocatable :: temp_particle(:), ion_temp_particle(:), merged_buffer(:)
     type(t_particle) :: sum_particle

     allocate(temp_particle(direction_cnt))
     allocate(ion_temp_particle(direction_cnt))
     allocate(merged_buffer(4))
     m_i = 0
     ! print *, "ll_elem count: ", CEILING(REAL(merged_cnt/size(merged_guide%tmp_particles))), buffer_pos

     !NOTE: use f_e as a dummy variable to store v^2, for energy balance during elastic merging.
     !      summing up the particles to one. used later for elastic merging.
     if (direction_cnt > 2) then
       ! variables used to count number of filled buffer. Sort to old or new particles.
       !NOTE: any new particles generated from ionisation events is marked with label = -1.0
       !      don't merge this with older particles.
       filled = 0
       filled_ion = 0
       do i = 1, direction_cnt
         if (directional_buffer(direction,i)%label .ge. 0) then
           filled = filled + 1
           temp_particle(filled) = directional_buffer(direction,i)
         else if(directional_buffer(direction,i)%label .lt. 0) then
           filled_ion = filled_ion + 1
           ion_temp_particle(filled_ion) = directional_buffer(direction,i)
         end if
       end do

       if (filled > 2) then
         call resolve_elastic_merge_momentum(temp_particle, filled, merged_buffer, m_i, 0)
       else
         do i = 1, filled
           m_i = m_i + 1
           merged_buffer(m_i) = temp_particle(i)
           merged_buffer(m_i)%label = 0
         end do
       end if

       if (filled_ion > 2) then
         call resolve_elastic_merge_momentum(ion_temp_particle, filled_ion, merged_buffer, m_i, 0)
       else
         do i = 1, filled_ion
           m_i = m_i + 1
           merged_buffer(m_i) = ion_temp_particle(i)
           merged_buffer(m_i)%label = 0
         end do
       end if
       ! after merging, 'label' variable in particle will have a new functionality.
       ! set label = 1 if 1 additional check during 'collision_update' is required.

     else
       do i = 1, direction_cnt
         m_i = m_i + 1
         merged_buffer(m_i) = directional_buffer(direction,i)
         merged_buffer(m_i)%label = 0
       end do
     end if

     !NOTE: now that merged_buffer is filled, copy to merged_guide ll_element.
     if (m_i .ne. 0) then
       ll_elem_gen = MOD(merged_cnt, size(merged_guide%tmp_particles))
       buffer_pos = ll_elem_gen + 1
       ! print *, "buffer_pos: ", buffer_pos
       if (ll_elem_gen == 0 .and. merged_cnt > 0) then
          call allocate_ll_buffer(size(merged_guide%tmp_particles), merged_guide%next)
          merged_guide => merged_guide%next
          nullify (merged_guide%next)
       end if

       ll_elem_gen = ll_elem_gen + m_i
       if (ll_elem_gen <= size(merged_guide%tmp_particles)) then
         remainder = buffer_pos + m_i - 1
         merged_guide%tmp_particles(buffer_pos:remainder) = merged_buffer(1:m_i)
         ! print *, "ll_elem_gen", ll_elem_gen, buffer_pos, remainder, m_i, parent_key
         ! print *, "Compare x: ", m_i,  merged_guide%tmp_particles(buffer_pos)%x, merged_buffer(1)%x
       else
         extent = size(merged_guide%tmp_particles) - buffer_pos + 1
         remainder = ll_elem_gen - size(merged_guide%tmp_particles)
         merged_guide%tmp_particles(buffer_pos:size(merged_guide%tmp_particles)) = &
         merged_buffer(1:extent)

         call allocate_ll_buffer(size(merged_guide%tmp_particles), merged_guide%next)
         merged_guide => merged_guide%next
         ! print *, "extension", ll_elem_gen, buffer_pos, remainder, (extent + 1), m_i, parent_key
         merged_guide%tmp_particles(1:remainder) = merged_buffer((extent + 1):m_i)
         nullify (merged_guide%next)
       end if
       merged_cnt = merged_cnt + m_i
       ! print *, "ll_elem_gen", ll_elem_gen, size(merged_guide%tmp_particles), merged_cnt, parent_key
     end if

     deallocate(temp_particle)
     deallocate(ion_temp_particle)
     deallocate(merged_buffer)
   end subroutine age_elastic_merging

   subroutine energy_elastic_merging(direction, directional_buffer, direction_cnt, merged_guide, merged_cnt, species, parent_key)
     !NOTE: in order to not merge particles that are .le. 2 particles lying in the same momentum partition,
     !       most direct way of doing it is to copy the particles into respective directional buffer.
     !       don't sum them up yet! information will be lost if done so.
     implicit none
     integer, intent(in) :: direction, direction_cnt, species, parent_key
     type(t_particle), allocatable, intent(in) :: directional_buffer(:,:)
     type(linked_list_elem), pointer, intent(inout) :: merged_guide
     integer, intent(inout) :: merged_cnt
     integer :: i, j, buffer_pos, ll_elem_gen, m_i, remainder, extent, filtered_instance, f_i, &
                merge_instance, merge_collector_size, IStart, IStop, j_start
     integer, allocatable :: grouped_count(:)
     real(kind_physics) :: kin_e, weight, vel(3)
     real(kind_physics), allocatable :: energy_threshold(:), max_weight(:)
     type(t_particle), allocatable :: energy_collector(:,:), merged_buffer(:), pass_buffer(:)
     type(t_particle) :: sum_particle

     filtered_instance = 6
     merge_collector_size = 4
     ! print *, "Direction ", direction, " of 8, starting merged_buffer allocation .", direction_cnt
     allocate(energy_threshold(filtered_instance))
     allocate(max_weight(filtered_instance))
     allocate(merged_buffer(direction_cnt))

     max_weight = 1.0_kind_physics

     m_i = 0
     ! print *, "ll_elem count: ", CEILING(REAL(merged_cnt/size(merged_guide%tmp_particles))), buffer_pos
     energy_threshold(1) = 0.0
     energy_threshold(2) = 10.0
     energy_threshold(3) = 20.0
     energy_threshold(4) = 30.0
     energy_threshold(5) = 40.0
     energy_threshold(6) = 60.0
     f_i = filtered_instance

     if (direction_cnt > 2) then
       allocate(grouped_count(filtered_instance))
       allocate(energy_collector(filtered_instance, direction_cnt))

       ! Sorting particle based on energy
       grouped_count = 0
       do i = 1, direction_cnt
         weight = abs(directional_buffer(direction,i)%data%q)
         vel = directional_buffer(direction,i)%data%v
         kin_e = 0.5*(directional_buffer(direction,i)%data%m/weight)*e_mass*dot_product(vel,vel)

         if (kin_e >= energy_threshold(filtered_instance)) then
           f_i = filtered_instance
         else
           do while (kin_e < energy_threshold(f_i))
             f_i = f_i - 1
           end do
         end if
         if(f_i < 1 .or. f_i > filtered_instance) print *, "Checking f_i: ", f_i
         grouped_count(f_i) = grouped_count(f_i) + 1
         energy_collector(f_i, grouped_count(f_i)) = directional_buffer(direction,i)
         if (weight > max_weight(f_i)) then
           max_weight(f_i) = weight
         end if

         ! Don't merge the particles with less than 10.0eV
         if (species .eq. 0) then
           if (f_i .eq. 1) then
             m_i = m_i + 1
             merged_buffer(m_i) = directional_buffer(direction,i)
             merged_buffer(m_i)%label = 0
           end if
         end if
       end do

       ! print *, "parent_key: ", parent_key, direction, grouped_count
       if (species .eq. 0) then
         j_start = 2
       else
         j_start = 1
       end if

       ! Merge the rest according to the grouped energy levels.
       do j = j_start, filtered_instance
         if(grouped_count(j) .ne. 0) then
           if (grouped_count(j) > 2) then
             if (max_weight(j) > 1.0_kind_physics) then
               allocate(pass_buffer(grouped_count(j)))
               pass_buffer = energy_collector(j,1:grouped_count(j))
               call sort_particles_by_weight(pass_buffer)
               energy_collector(j,1:grouped_count(j)) = pass_buffer
               deallocate(pass_buffer)
             end if

             merge_instance = FLOOR(0.5*grouped_count(j)*(1.0 - merge_ratio))
             if (merge_ratio < 0.5) then
               merge_instance = FLOOR(merge_ratio*grouped_count(j)*0.5)
               merge_collector_size = FLOOR(1.0*grouped_count(j)/merge_instance)
             end if
             remainder = grouped_count(j) - merge_collector_size*merge_instance

             do i = 1, merge_instance
               IStart = (i - 1)*merge_collector_size + 1
               IStop = IStart + merge_collector_size - 1
               allocate(pass_buffer(merge_collector_size))
               pass_buffer = energy_collector(j,IStart:IStop)
               call resolve_elastic_merge_momentum(pass_buffer, merge_collector_size, merged_buffer, m_i, j)
               deallocate(pass_buffer)
             end do

             IStart = grouped_count(j) - remainder + 1
             IStop = grouped_count(j)
             if (merge_ratio >= 0.5 .or. remainder <= 2) then
               do i = IStart, IStop
                 m_i = m_i + 1
                 merged_buffer(m_i) = energy_collector(j,i)
                 merged_buffer(m_i)%label = 0
               end do
             else if (merge_ratio < 0.5 .and. remainder > 2) then
               allocate(pass_buffer(remainder))
               pass_buffer = energy_collector(j,IStart:IStop)
               call resolve_elastic_merge_momentum(pass_buffer, remainder, merged_buffer, m_i, j)
               deallocate(pass_buffer)
             end if
             ! print *, "n_t: ", grouped_count(j), " target: ", merge_ratio*grouped_count(j), &
             !          " remainder: ", remainder, " merged p: ", merge_instance*2, ". Diff: ", &
             !          merge_ratio*grouped_count(j) - (merge_instance*2 + remainder)
           else
             do i = 1, grouped_count(j)
               m_i = m_i + 1
               merged_buffer(m_i) = energy_collector(j,i)
               merged_buffer(m_i)%label = 0
             end do
           end if
         end if
       end do

       deallocate(grouped_count)
       deallocate(energy_collector)
     else
       do i = 1, direction_cnt
         m_i = m_i + 1
         merged_buffer(m_i) = directional_buffer(direction,i)
         merged_buffer(m_i)%label = 0
       end do
     end if
     deallocate(energy_threshold)
     deallocate(max_weight)

     !NOTE: now that merged_buffer is filled, copy to merged_guide ll_element.
     if (m_i .ne. 0) then
       ll_elem_gen = MOD(merged_cnt, size(merged_guide%tmp_particles))
       buffer_pos = ll_elem_gen + 1
       ! print *, "buffer_pos: ", buffer_pos
       if (ll_elem_gen == 0 .and. merged_cnt > 0) then
          call allocate_ll_buffer(size(merged_guide%tmp_particles), merged_guide%next)
          merged_guide => merged_guide%next
          nullify (merged_guide%next)
       end if

       ll_elem_gen = ll_elem_gen + m_i
       if (ll_elem_gen <= size(merged_guide%tmp_particles)) then
         remainder = buffer_pos + m_i - 1
         merged_guide%tmp_particles(buffer_pos:remainder) = merged_buffer(1:m_i)
         ! print *, "ll_elem_gen", ll_elem_gen, buffer_pos, remainder, m_i, parent_key
       else
         extent = size(merged_guide%tmp_particles) - buffer_pos + 1
         remainder = ll_elem_gen - size(merged_guide%tmp_particles)
         merged_guide%tmp_particles(buffer_pos:size(merged_guide%tmp_particles)) = &
         merged_buffer(1:extent)

         call allocate_ll_buffer(size(merged_guide%tmp_particles), merged_guide%next)
         merged_guide => merged_guide%next
         ! print *, "extension", ll_elem_gen, buffer_pos, remainder, (extent + 1), m_i, parent_key
         merged_guide%tmp_particles(1:remainder) = merged_buffer((extent + 1):m_i)
         nullify (merged_guide%next)
       end if
       merged_cnt = merged_cnt + m_i
       ! print *, "ll_elem_gen", ll_elem_gen, size(merged_guide%tmp_particles), merged_cnt, parent_key
     end if
     deallocate(merged_buffer)
     ! print *, "merged_buffer Deallocation done"
   end subroutine energy_elastic_merging

   subroutine resolve_elastic_merge_momentum(input_buffer, input_cnt, merged_buffer, m_i, group_index)
     implicit none
     type(t_particle), allocatable, intent(in) :: input_buffer(:)
     integer, intent(in) :: input_cnt, group_index
     type(t_particle), allocatable, intent(inout) :: merged_buffer(:)
     integer, intent(inout) :: m_i
     type(t_particle) :: sum_particle
     real(kind_physics) :: d(3), max_vec(3), v_squared, vel_mag2, vel_mag, unit_vec(3), &
                           Chi, rot_axis(3), w1, w2, total_weight, ave_E, mass, charge, &
                           v_diff, e_diff
     integer :: i, j, species, parent_key, correct_sign

     species = input_buffer(1)%data%species
     parent_key = input_buffer(1)%data%mp_int1

     if (species == 0) then
       mass = 1.0_kind_physics
       charge = -1.0_kind_physics
     else if (species == 1) then
       mass = 1836.21957489_kind_physics
       charge = 1.0_kind_physics
     else if (species == 2) then
       mass = 3673.43889456_kind_physics
       charge = 1.0_kind_physics
     end if
     !initialise sum_particle
     sum_particle%data%q = 0.0_kind_physics
     sum_particle%data%v = 0.0_kind_physics
     sum_particle%data%m = 0.0_kind_physics
     sum_particle%data%f_e = 0.0_kind_physics
     sum_particle%data%f_b = 0.0_kind_physics
     sum_particle%x = 0.0_kind_physics

     ! initialise d(3) vector, essentially the unit vector of the max extent in x, y, z, direction
     d = 0.0
     max_vec = 0.0
     total_weight = 0.0
     do i = 1, input_cnt
       v_squared = abs(input_buffer(i)%data%q)*dot_product(input_buffer(i)%data%v, input_buffer(i)%data%v)
       sum_particle%data%q = sum_particle%data%q + input_buffer(i)%data%q
       sum_particle%data%v = sum_particle%data%v + abs(input_buffer(i)%data%q)*input_buffer(i)%data%v
       sum_particle%data%m = sum_particle%data%m + input_buffer(i)%data%m
       sum_particle%x = sum_particle%x + input_buffer(i)%x
       sum_particle%data%f_e(1) = sum_particle%data%f_e(1) + v_squared
       total_weight = total_weight + abs(input_buffer(i)%data%q)

       do j = 1, 3
         if (abs(input_buffer(i)%data%v(j)) > max_vec(j)) then
           max_vec(j) = abs(input_buffer(i)%data%v(j))
           d(j) = input_buffer(i)%data%v(j)
         end if
       end do
     end do
     d = d/sqrt(dot_product(d,d))
     sum_particle%x = sum_particle%x/abs(sum_particle%data%q)
     ave_E = 0.5*sum_particle%data%m*e_mass*sum_particle%data%f_e(1)/(total_weight*total_weight)

     max_vec = 0.0_kind_physics
     vel_mag2 = dot_product(sum_particle%data%v, sum_particle%data%v)
     vel_mag = sqrt(vel_mag2)

     if (vel_mag2 .ne. 0.0_kind_physics) then
       w1 = CEILING(1.0*total_weight/2.0)
       if (w1 < (total_weight - vel_mag2/sum_particle%data%f_e(1))) then
         w1 = CEILING(total_weight - vel_mag2/sum_particle%data%f_e(1))
       end if

       unit_vec = sum_particle%data%v/vel_mag
       Chi = acos(sqrt(total_weight/w1 - (total_weight*total_weight/w1 - total_weight)*sum_particle%data%f_e(1)/vel_mag2))
       ! if (isnan(Chi)) print *, "Chi calculation failed!"

       rot_axis = cross_product(sum_particle%data%v, d)
       rot_axis = rot_axis/sqrt(dot_product(rot_axis, rot_axis))
       max_vec = Rodriguez_rotation(Chi, rot_axis, unit_vec)

       ! Check if the resulting velocity vector is in the same momentum quadrant.
       ! if not, flip the rotation direction.
       correct_sign = 1
       do j = 1, 3
         if (max_vec(j)*sum_particle%data%v(j) .lt. 0) then
           correct_sign = 0
         end if
       end do
       if (correct_sign == 0) then
         max_vec = Rodriguez_rotation(-1.0*Chi, rot_axis, unit_vec)
       end if

       ! Fill in values for merged particle 1
       m_i = m_i + 1
       merged_buffer(m_i)%x = input_buffer(1)%x
       merged_buffer(m_i)%work = 1.0_8
       merged_buffer(m_i)%data%v = dot_product(sum_particle%data%v, max_vec)/total_weight * max_vec
       merged_buffer(m_i)%data%q = charge*w1
       merged_buffer(m_i)%data%m = mass*w1
       merged_buffer(m_i)%data%b = 0.0
       merged_buffer(m_i)%data%f_e = 0.0
       merged_buffer(m_i)%data%f_b = 0.0
       merged_buffer(m_i)%data%f_b(2) = ave_E
       merged_buffer(m_i)%data%age = 0.0
       merged_buffer(m_i)%label = group_index
       merged_buffer(m_i)%data%species = species
       merged_buffer(m_i)%data%mp_int1 = parent_key

       ! Fill in values for merged particle 2
       w2 = total_weight - w1
       m_i = m_i + 1
       merged_buffer(m_i)%x = input_buffer(2)%x
       merged_buffer(m_i)%work = 1.0_8
       merged_buffer(m_i)%data%v = (sum_particle%data%v - w1 * merged_buffer(m_i - 1)%data%v)/w2
       merged_buffer(m_i)%data%q = charge*w2
       merged_buffer(m_i)%data%m = mass*w2
       merged_buffer(m_i)%data%b = 0.0
       merged_buffer(m_i)%data%f_e = 0.0
       merged_buffer(m_i)%data%f_b = 0.0
       merged_buffer(m_i)%data%f_b(2) = ave_E
       merged_buffer(m_i)%data%age = 0.0
       merged_buffer(m_i)%label = group_index
       merged_buffer(m_i)%data%species = species
       merged_buffer(m_i)%data%mp_int1 = parent_key

       ! merged_buffer(m_i)%data%b = sum_particle%data%v - (w1*merged_buffer(m_i - 1)%data%v + w2*merged_buffer(m_i)%data%v)
       ! v_diff = sqrt(dot_product(merged_buffer(m_i)%data%b, merged_buffer(m_i)%data%b))
       ! if (v_diff > 1e-16_kind_physics) print *, "Noticeable Velocity Error!", v_diff
       !
       ! e_diff = sum_particle%data%f_e(1) - &
       !          (w1*dot_product(merged_buffer(m_i - 1)%data%v, merged_buffer(m_i - 1)%data%v) + &
       !           w2*dot_product(merged_buffer(m_i)%data%v, merged_buffer(m_i)%data%v))
       ! if (e_diff > 1e-16_kind_physics) print *, "Noticeable Energy Error!", e_diff
       !
       ! merged_buffer(m_i)%data%b = 0.0
     else
       w1 = CEILING(1.0*total_weight/2.0)
       ! Fill in values for merged particle 1
       m_i = m_i + 1
       merged_buffer(m_i)%x = input_buffer(1)%x
       merged_buffer(m_i)%work = 1.0_8
       merged_buffer(m_i)%data%v = 0.0_kind_physics
       merged_buffer(m_i)%data%q = charge*w1
       merged_buffer(m_i)%data%m = mass*w1
       merged_buffer(m_i)%data%b = 0.0
       merged_buffer(m_i)%data%f_e = 0.0
       merged_buffer(m_i)%data%f_b = 0.0
       merged_buffer(m_i)%data%age = 0.0
       merged_buffer(m_i)%label = group_index
       merged_buffer(m_i)%data%species = species
       merged_buffer(m_i)%data%mp_int1 = parent_key

       ! Fill in values for merged particle 2
       w2 = total_weight - w1
       m_i = m_i + 1
       merged_buffer(m_i)%x = input_buffer(2)%x
       merged_buffer(m_i)%work = 1.0_8
       merged_buffer(m_i)%data%v = 0.0_kind_physics
       merged_buffer(m_i)%data%q = charge*w2
       merged_buffer(m_i)%data%m = mass*w2
       merged_buffer(m_i)%data%b = 0.0
       merged_buffer(m_i)%data%f_e = 0.0
       merged_buffer(m_i)%data%f_b = 0.0
       merged_buffer(m_i)%data%age = 0.0
       merged_buffer(m_i)%label = group_index
       merged_buffer(m_i)%data%species = species
       merged_buffer(m_i)%data%mp_int1 = parent_key
     end if
   end subroutine resolve_elastic_merge_momentum

   subroutine add_particle(guide, particle, new_particle, buffer_pos, ran, type)
     ! NOTE: A big assumption is made here: Without considering the rovibrational states
     !       of the reaction products, generated particle doesn`t store their own internal
     !       energy value. Instead, the reaction energy is entirely converted into kinetic E.
     !       More traditionally, this redistribution of energy is termed as Kinetic Energy Release
     !       in mass spetroscopy
     implicit none
     type(t_particle), intent(in) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     integer, intent(inout) :: new_particle, buffer_pos
     real(kind_physics), intent(in) :: ran(8)
     integer, intent(in) :: type
     integer :: ll_elem_gen

     ll_elem_gen = MOD(new_particle, size(guide%tmp_particles))
     buffer_pos = ll_elem_gen + 1
     if (ll_elem_gen == 0 .and. new_particle > 0) then
        ! print *, "Extending Temp_array!!!!!"
        call allocate_ll_buffer(size(guide%tmp_particles), guide%next)
        guide => guide%next
        nullify (guide%next)
     end if

     select case(type)
     case(0) ! Electron

       guide%tmp_particles(buffer_pos)%data%q = -1.0
       guide%tmp_particles(buffer_pos)%data%m = 1.0_8
       guide%tmp_particles(buffer_pos)%data%species = 0

     case(1) ! H+

       guide%tmp_particles(buffer_pos)%data%q = 1.0
       guide%tmp_particles(buffer_pos)%data%m = 1836.21957489_8 ! times greater than mass of electron
       guide%tmp_particles(buffer_pos)%data%species = 1

     case(2) ! H2+

       guide%tmp_particles(buffer_pos)%data%q = 1.0
       guide%tmp_particles(buffer_pos)%data%m = 3673.43889456_8 ! times greater than mass of electron
       guide%tmp_particles(buffer_pos)%data%species = 2

     case(3) ! H atom

       guide%tmp_particles(buffer_pos)%data%q = 0.0
       guide%tmp_particles(buffer_pos)%data%m = 1837.21957489_8 ! times greater than mass of electron
       guide%tmp_particles(buffer_pos)%data%species = 3

     end select

     if (type .eq. 0) then
       guide%tmp_particles(buffer_pos)%x(1) = particle%x(1) + (2.0*ran(6) - 1.0)*seed_dl(1)
       guide%tmp_particles(buffer_pos)%x(2) = particle%x(2) + (2.0*ran(7) - 1.0)*seed_dl(2)
       guide%tmp_particles(buffer_pos)%x(3) = particle%x(3) + (2.0*ran(8) - 1.0)*seed_dl(3)
     else
       guide%tmp_particles(buffer_pos)%x = particle%x
     end if
     guide%tmp_particles(buffer_pos)%work = 1.0_8
     guide%tmp_particles(buffer_pos)%data%age = 0.0_8
     guide%tmp_particles(buffer_pos)%data%b = 0.0_kind_physics
     guide%tmp_particles(buffer_pos)%data%f_e = 0.0_kind_physics
     guide%tmp_particles(buffer_pos)%data%f_b = 0.0_kind_physics
     guide%tmp_particles(buffer_pos)%results%pot = 0.0_kind_physics
     guide%tmp_particles(buffer_pos)%data%mp_int1 = 0

     guide%tmp_particles(buffer_pos)%label = -1
     guide%tmp_particles(buffer_pos)%data%v = 0.0_8

     new_particle = new_particle + 1
   end subroutine add_particle

   subroutine collision_update(particle, guide, new_particle, electron_count, CS_numbers, r123_ctr, r123_key, charge_count)
     implicit none
     type(t_particle), intent(inout) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     integer, intent(inout) :: new_particle, electron_count
     integer(kind=int32), intent(inout) :: r123_ctr(4), r123_key(4)
     integer, intent(in) :: CS_numbers
     real(kind_physics), intent(inout) :: charge_count(5)
     real(kind_physics), dimension(:), allocatable :: CS_vector
     real(kind_physics) :: vel_mag, nu_prime, reduced_vel_mag, R_J02, V_V01, IE_H2_ion, AE_H_ion, theta, phi, polar_theta, polar_phi
     real(kind_physics) :: rot_axis(3), temp_vel(3), temp_vel1(3), reduced_incident(3), cos_theta, temp_vel_mag, temp_vel1_mag
     real(kind_physics) :: H2_mass, K_E, cos_Chi, Chi, unit_inc_vel(3), chi_rotation_axis(3), prefac, scatter_loss, xi, rand_num(8)
     real(kind_physics) :: weight, scaled_mass, Disso_H1, Disso_H2
     integer :: buff_pos, i, CS_index

     ! TODO: not a huge fan of the next 4 lines. Rather make those parameters which might help optimisation
     R_J02 = 0.0441 ! eV, transition energy associated with rotational excitation from J = 0 -> 2
     V_V01 = 0.516 ! eV, transition energy associated with vibrational excitation from V = 0 -> 1
     IE_H2_ion = 15.426 ! 15.283 eV, Ionization energy of H2+ [Source: T.E.Sharp Atomic Data 2, 119-169 (1971)]. 15.426 by Yoon 2008
     AE_H_ion = 18.075 ! eV, Appearance energy of H+, 18.1 by Yoon 2008
     Disso_H1 = 4.47787 !eV, H dissociation energy H(1s) + H(1s)
     Disso_H2 = 14.676 !eV, H dissociation energy H(1s) + H(2s)
     H2_mass = 3674.43889456_8 ! H2/e mass ratio
     ! [definition of Appearance energy vs Ionization energy from Mass Spectroscopy (2011) by Gross J.H.]

     ! NOTE: Currently, calculation of nu_prime involves obtaining abs_max_CS,
     !       which is the max of total cross section. nu_prime is then obtained
     !       by multiplying abs_max_CS with the particle`s velocity magnitude and
     !       constant local density. nu_prime will thus be always larger than actual nu!
     ! TODO: not at all critical, but there are many ways of doing the next line, would be interesting to see if there is a fastest
     ! way of computing the length and weigh that against readability, perhaps a separate inlinable function if there is a faster
     ! way?
     weight = abs(particle%data%q)
     scaled_mass = particle%data%m/weight
     vel_mag = sqrt(dot_product(particle%data%v,particle%data%v))
     K_E = 0.5*scaled_mass*e_mass*vel_mag**2
     nu_prime = abs_max_CS * vel_mag * neutral_density

     ! if ((abs(particle%data%q) .gt. 1.0_kind_physics)) then
     !   K_E = particle%data%f_b(2)
     ! end if
    !  call determine_cross_sections(particle, CS_vector, CS_tables)
    !  CS_vector = CS_vector * vel_mag * neutral_density
    !  print *, nu_prime, (1 - exp(-1*nu_prime*dt)), CS_vector/nu_prime

     ! Generating Random Number between [0,1]
     dummy = gen_norm_double_rng(r123_ctr, r123_key, rand_num)

    !  print *, "rand: ", rand_num(1), " expression: ", (1 - exp(-1*nu_prime*dt))!(1 - exp(-1*CS_vector(size(CS_vector))*dt))
     if (rand_num(1) < (1 - exp(-1*nu_prime*dt))) then ! type of collision determined if satisfied
       allocate(CS_vector(total_cross_sections))
       call determine_cross_sections(particle, CS_vector, CS_tables)
       call determine_xi(K_E, xi)
       CS_index = total_cross_sections - eirene_cross_sections + 1
       call Eirene_fit(eirene_coeffs1, K_E, CS_vector, CS_index)
       call Eirene_fit(eirene_coeffs2, K_E, CS_vector, CS_index)

       !======================== For use in random scatter ================
       polar_theta = 2.0*pi*rand_num(4)
       polar_phi = acos(2.0*rand_num(3) - 1.)

       !=====================For use in energy dependent scatter =========
       ! cos_Chi = (2. + K_E - 2.*(1. + K_E)**rand_num(5))/K_E ! [Vahedi & Surendra]
       cos_Chi = 1 - (2.*rand_num(5)*(1. - xi))/(1. + xi*(1. - 2.*rand_num(5))) ! [Ohkrimovsky 2002]
       Chi = acos(cos_Chi)
       unit_inc_vel = particle%data%v/vel_mag
       scatter_loss = 2.*(1. - cos_Chi)/H2_mass ! lost energy calculation

       ! Convert CS_vector (cross section of all reactions) to Collision freq., nu.
       CS_vector = CS_vector * vel_mag * neutral_density !6545520.13889 test value of constant local_number_density (at 0.001Pa)

       ! TODO: this to me looks like it could be a simple do-loop
       ! TODO: also looks like it can be combined with the select case below via a rather long/convoluted if-elseif construct that
       ! may read more easily
       i = 1
       do while (rand_num(2) > (CS_vector(i)/nu_prime))
         i = i + 1
         if (i > size(CS_vector)) then
           i = 0
           EXIT
         end if
       end do
       deallocate(CS_vector)
     else ! otherwise, no collision happened
       i = 0
     end if

     select case(i)
     case(0) ! null collision, no update performed
      !  print *, "null coll!"
     case(1) ! elastic scattering (no additional electron, no byproducts)
       ! update velocity to indicate scattering into random angle
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*scatter_loss/e_mass)
       ! particle%data%v(1) = vel_mag * sin(polar_phi) * cos(polar_theta)
       ! particle%data%v(2) = vel_mag * sin(polar_phi) * sin(polar_theta)
       ! particle%data%v(3) = vel_mag * cos(polar_phi)

       ! =================Vahedi Surendra 1995 / Ohkrimovsky 2002 ==============
       rot_axis = 0.0
       rot_axis(1) = 1.0
       polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
       temp_vel = cross_product(unit_inc_vel, rot_axis)
       temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
       prefac = sin(Chi)/sin(polar_phi)

       particle%data%v = reduced_vel_mag*(cos_Chi*unit_inc_vel + &
                         sin(polar_theta)*prefac*temp_vel + &
                         cos(polar_theta)*prefac*temp_vel1)

       ! =======Anti Parallel Scattering ============================
       ! particle%data%v = -1.0*reduced_vel_mag*unit_inc_vel

       particle%data%age = 0.0_8
      !  print *, "elastic coll.!"

     case(2) ! rotational excitation of H2 molecule, electron will lose the transition energy
       ! update velocity to indicate scattering into random angle
       ! TODO: will be less costly to keep 2/e_mass as parameter and multiplying with it - not sure the compile will pick up on this
       ! and do the short-cut for you. Better still: add it to the energies straight away.
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*(R_J02 + scatter_loss)/e_mass)

       ! particle%data%v(1) = reduced_vel_mag * sin(polar_phi) * cos(polar_theta)
       ! particle%data%v(2) = reduced_vel_mag * sin(polar_phi) * sin(polar_theta)
       ! particle%data%v(3) = reduced_vel_mag * cos(polar_phi)

       ! =================Vahedi Surendra 1995 / Ohkrimovsky 2002 ==============
       rot_axis = 0.0
       rot_axis(1) = 1.0
       polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
       temp_vel = cross_product(unit_inc_vel, rot_axis)
       temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
       prefac = sin(Chi)/sin(polar_phi)

       particle%data%v = reduced_vel_mag*(cos_Chi*unit_inc_vel + &
                        sin(polar_theta)*prefac*temp_vel + &
                        cos(polar_theta)*prefac*temp_vel1)

       ! ==========================Anti Parallel Scatter =======================
       ! particle%data%v = -1.0*reduced_vel_mag*unit_inc_vel

       particle%data%age = 0.0_8
      !  print *, "rotational exci.!"

     case(3) ! vibrational excitation of H2 molecule, electron will lose the transition energy
       ! update velocity to indicate scattering into random angle
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*(V_V01 + scatter_loss)/e_mass)

       ! particle%data%v(1) = reduced_vel_mag * sin(polar_phi) * cos(polar_theta)
       ! particle%data%v(2) = reduced_vel_mag * sin(polar_phi) * sin(polar_theta)
       ! particle%data%v(3) = reduced_vel_mag * cos(polar_phi)

       ! =================Vahedi Surendra 1995 / Ohkrimovsky 2002 ==============
       rot_axis = 0.0
       rot_axis(1) = 1.0
       polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
       temp_vel = cross_product(unit_inc_vel, rot_axis)
       temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
       prefac = sin(Chi)/sin(polar_phi)

       particle%data%v = reduced_vel_mag*(cos_Chi*unit_inc_vel + &
                         sin(polar_theta)*prefac*temp_vel + &
                         cos(polar_theta)*prefac*temp_vel1)

       ! ==========================Anti Parallel Scatter =======================
       ! particle%data%v = -1.0*reduced_vel_mag*unit_inc_vel

       particle%data%age = 0.0_8
      !  print *, "vibrational exci.!"

     case(4) ! nondissociative ionization (1 additional electron, 1 byproduct)
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*(IE_H2_ion + scatter_loss)/e_mass)
       reduced_incident = particle%data%v * (reduced_vel_mag/vel_mag)
       call angles_calc(reduced_incident, reduced_vel_mag, theta, phi)
      !  print *, "incident momentum: ", scaled_mass* reduced_incident

       call add_particle(guide, particle, new_particle, buff_pos, rand_num, 0)
       ! ===========================Random Scatter==============================
       ! temp_vel(1) = sin(polar_phi*0.5) * cos(polar_theta)
       ! temp_vel(2) = sin(polar_phi*0.5) * sin(polar_theta)
       ! temp_vel(3) = cos(polar_phi*0.5)

       ! rot_axis(1) = 0.0
       ! rot_axis(2) = 1.0
       ! rot_axis(3) = 0.0
       ! temp_vel = Rodriguez_rotation(theta, rot_axis, temp_vel)
       ! rot_axis(2) = 0.0
       ! rot_axis(3) = 1.0
       ! temp_vel = Rodriguez_rotation(phi, rot_axis, temp_vel)

       ! =================Vahedi Surendra 1995 / Ohkrimovsky 2002 ==============
       rot_axis = 0.0
       rot_axis(1) = 1.0
       polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
       temp_vel = cross_product(unit_inc_vel, rot_axis)
       temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
       prefac = sin(Chi)/sin(polar_phi)

       temp_vel = (cos_Chi*unit_inc_vel + sin(polar_theta)*prefac*temp_vel + &
                  cos(polar_theta)*prefac*temp_vel1)

       ! NOTE: Scaling proper velocity magnitudes w.r.t K_E & momentum
       cos_theta = dot_product(temp_vel, reduced_incident)/reduced_vel_mag
       temp_vel_mag = (2.*scaled_mass/(guide%tmp_particles(buff_pos)%data%m + scaled_mass)) &
                      * reduced_vel_mag * cos_theta
       temp_vel1_mag = reduced_vel_mag*sqrt(scaled_mass**2 + 2. * scaled_mass * guide%tmp_particles(buff_pos)%data%m &
                       * (1. - 2. * (cos_theta**2)) + guide%tmp_particles(buff_pos)%data%m**2)/(guide%tmp_particles(buff_pos)%data%m + scaled_mass)
       guide%tmp_particles(buff_pos)%data%v = temp_vel * temp_vel_mag

       temp_vel1 = reduced_incident - (guide%tmp_particles(buff_pos)%data%m/scaled_mass) * &
                   guide%tmp_particles(buff_pos)%data%v

       particle%data%v = temp_vel1
       particle%data%age = 0.0_8
      !  print *, "outgoing momentum: ", scaled_mass*temp_vel1 + guide%tmp_particles(buff_pos)%data%v * guide%tmp_particles(buff_pos)%data%m
      !  print *, "incident kinetic energy: ", 0.5*reduced_vel_mag**2
      !  print *, "outgoing kinetic energy: ", 0.5*(temp_vel1_mag**2 + temp_vel_mag**2)

       call add_particle(guide, particle, new_particle, buff_pos, rand_num, 2)
      !  print *, "nondissoc.!"

     case(5) ! dissociative ionization (1 additional electron, 2 byproducts), Hydrogen atom is ignored!
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*(AE_H_ion+scatter_loss)/e_mass)
       reduced_incident = particle%data%v * (reduced_vel_mag/vel_mag)
       call angles_calc(reduced_incident, reduced_vel_mag, theta, phi)
      !  print *, "incident momentum: ", scaled_mass * reduced_incident

       call add_particle(guide, particle, new_particle, buff_pos, rand_num, 0)
       ! temp_vel(1) = sin(polar_phi*0.5) * cos(polar_theta)
       ! temp_vel(2) = sin(polar_phi*0.5) * sin(polar_theta)
       ! temp_vel(3) = cos(polar_phi*0.5)

       ! rot_axis(1) = 0.0
       ! rot_axis(2) = 1.0
       ! rot_axis(3) = 0.0
       ! temp_vel = Rodriguez_rotation(theta, rot_axis, temp_vel)
       ! rot_axis(2) = 0.0
       ! rot_axis(3) = 1.0
       ! temp_vel = Rodriguez_rotation(phi, rot_axis, temp_vel)

       rot_axis = 0.0
       rot_axis(1) = 1.0
       polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
       temp_vel = cross_product(unit_inc_vel, rot_axis)
       temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
       prefac = sin(Chi)/sin(polar_phi)

       temp_vel = (cos_Chi*unit_inc_vel + sin(polar_theta)*prefac*temp_vel + &
                  cos(polar_theta)*prefac*temp_vel1)

       ! NOTE: Scaling proper velocity magnitudes w.r.t K_E & momentum
       cos_theta = dot_product(temp_vel, reduced_incident)/reduced_vel_mag
       temp_vel_mag = (2.*scaled_mass/(guide%tmp_particles(buff_pos)%data%m + scaled_mass)) &
                      * reduced_vel_mag * cos_theta
       temp_vel1_mag = reduced_vel_mag*sqrt(scaled_mass**2 + 2. * scaled_mass * guide%tmp_particles(buff_pos)%data%m &
                       * (1. - 2. * (cos_theta**2)) + guide%tmp_particles(buff_pos)%data%m**2)/(guide%tmp_particles(buff_pos)%data%m + scaled_mass)
       guide%tmp_particles(buff_pos)%data%v = temp_vel * temp_vel_mag

       temp_vel1 = reduced_incident - (guide%tmp_particles(buff_pos)%data%m/scaled_mass) * &
                   guide%tmp_particles(buff_pos)%data%v

       particle%data%v = temp_vel1
       particle%data%age = 0.0_8
      !  print *, "outgoing momentum: ", scaled_mass*temp_vel1 + guide%tmp_particles(buff_pos)%data%v * guide%tmp_particles(buff_pos)%data%m
      !  print *, "incident kinetic energy: ", 0.5*reduced_vel_mag**2
      !  print *, "outgoing kinetic energy: ", 0.5*(temp_vel1_mag**2 + temp_vel_mag**2)

       call add_particle(guide, particle, new_particle, buff_pos, rand_num, 1)
      !  print *, "dissoc.!"

    case(6) ! H2 molecule dissociation into H(1s) atoms (no additional electron, no byproducts)
       ! update velocity to indicate scattering into random angle
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*(Disso_H1 + scatter_loss)/e_mass)
       ! particle%data%v(1) = vel_mag * sin(polar_phi) * cos(polar_theta)
       ! particle%data%v(2) = vel_mag * sin(polar_phi) * sin(polar_theta)
       ! particle%data%v(3) = vel_mag * cos(polar_phi)

       ! =================Vahedi Surendra 1995 / Ohkrimovsky 2002 ==============
       rot_axis = 0.0
       rot_axis(1) = 1.0
       polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
       temp_vel = cross_product(unit_inc_vel, rot_axis)
       temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
       prefac = sin(Chi)/sin(polar_phi)

       particle%data%v = reduced_vel_mag*(cos_Chi*unit_inc_vel + &
                         sin(polar_theta)*prefac*temp_vel + &
                         cos(polar_theta)*prefac*temp_vel1)

       ! =======Anti Parallel Scattering ============================
       ! particle%data%v = -1.0*reduced_vel_mag*unit_inc_vel

       particle%data%age = 0.0_8
      !  print *, "elastic coll.!"
      !NOTE: remember to add a H atom counter, note that it is required to consider the
      !      enveloping Openmp operation, as well as the MPI implementation!
      !      Consider extending 'thread_charge_count' by 1 field.
      charge_count(4) = charge_count(4) + 2

    case(7) ! H2 molecule dissociation into H(1s) & H(2s) atoms (no additional electron, no byproducts)
      ! update velocity to indicate scattering into random angle
      reduced_vel_mag = sqrt(vel_mag**2 - 2.*(Disso_H2 + scatter_loss)/e_mass)
      ! particle%data%v(1) = vel_mag * sin(polar_phi) * cos(polar_theta)
      ! particle%data%v(2) = vel_mag * sin(polar_phi) * sin(polar_theta)
      ! particle%data%v(3) = vel_mag * cos(polar_phi)

      ! =================Vahedi Surendra 1995 / Ohkrimovsky 2002 ==============
      rot_axis = 0.0
      rot_axis(1) = 1.0
      polar_phi = acos(dot_product(unit_inc_vel,rot_axis))
      temp_vel = cross_product(unit_inc_vel, rot_axis)
      temp_vel1 = cross_product(unit_inc_vel, -1.*temp_vel)
      prefac = sin(Chi)/sin(polar_phi)

      particle%data%v = reduced_vel_mag*(cos_Chi*unit_inc_vel + &
                        sin(polar_theta)*prefac*temp_vel + &
                        cos(polar_theta)*prefac*temp_vel1)

      ! =======Anti Parallel Scattering ============================
      ! particle%data%v = -1.0*reduced_vel_mag*unit_inc_vel

      particle%data%age = 0.0_8
      !   print *, "elastic coll.!"
      !NOTE: remember to add a H atom counter, note that it is required to consider the
      !      enveloping Openmp operation, as well as the MPI implementation!
      !      Consider extending 'thread_charge_count' by 1 field.
      charge_count(4) = charge_count(4) + 1
      charge_count(5) = charge_count(5) + 1

     end select
     particle%label = 0
   end subroutine collision_update

   recursive subroutine filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, charge_count, break_check)
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, intent(inout) :: swapped_cnt, break_check
     real(kind_physics), intent(inout) :: charge_count(5)
     integer, intent(in) :: geometry, current_index, head_i, tail_i
     integer :: target_swap, init_swap_cnt!, buffer_size
     real(kind_physics) :: center_pos(3)
     real(kind_physics) :: radius, cyl_radius, plate_radius, cyl_length, box_dim, &
                           x, y, z, box_l, m_radius, buffer_zone

     ! buffer_size = size(particles)
     ! target_swap = buffer_size - swapped_cnt

     center_pos = 0.0
     center_pos(3) = 0.0

     x = particles(current_index)%x(1) - center_pos(1)
     y = particles(current_index)%x(2) - center_pos(2)

     init_swap_cnt = swapped_cnt
     target_swap = tail_i - swapped_cnt

     select case(geometry)
     case(1) ! box boundary
       box_dim = 1.0/(c*1e-12)
       box_l = 0.5/(c*1e-12)

       if (current_index .lt. target_swap) then
         if ((abs(x) > box_dim*0.5) .or. (abs(y) > box_dim*0.5)) then
           particles(current_index) = particles(target_swap)
           swapped_cnt = swapped_cnt + 1
         else
           if (particles(current_index)%x(3) > 0.0) then
             charge_count(1) = charge_count(1) + particles(current_index)%data%q
             particles(current_index) = particles(target_swap)
             swapped_cnt = swapped_cnt + 1
           else if (particles(current_index)%x(3) < -box_l) then
             charge_count(2) = charge_count(2) + particles(current_index)%data%q
             particles(current_index) = particles(target_swap)
             swapped_cnt = swapped_cnt + 1
           end if
         end if

         if (init_swap_cnt /= swapped_cnt) then
           call filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, charge_count, break_check)
         end if
       else
         if ((abs(x) > box_dim*0.5) .or. (abs(y) > box_dim*0.5)) then
           swapped_cnt = swapped_cnt + 1
         else
           if (particles(current_index)%x(3) > 0.0) then
             charge_count(1) = charge_count(1) + particles(current_index)%data%q
             swapped_cnt = swapped_cnt + 1
           else if (particles(current_index)%x(3) < -box_l) then
             charge_count(2) = charge_count(2) + particles(current_index)%data%q
             swapped_cnt = swapped_cnt + 1
           end if
         end if
         break_check = 1
       end if

     case(2) ! cylinder boundary
       cyl_radius = 0.05/(c*1e-12)
       plate_radius = 0.025/(c*1e-12)
       cyl_length = d/(c*1e-12)
       radius = sqrt(x**2 + y**2)

       if (current_index .lt. target_swap) then
         if (radius < cyl_radius) then
           if (particles(current_index)%x(3) < -cyl_length)then
             if (radius < plate_radius) then
               charge_count(2) = charge_count(2) + particles(current_index)%data%q
               particles(current_index) = particles(target_swap)
               swapped_cnt = swapped_cnt + 1
             else
               particles(current_index) = particles(target_swap)
               swapped_cnt = swapped_cnt + 1
             end if
           else if (particles(current_index)%x(3) > 0.0) then
             if (radius < plate_radius) then
               charge_count(1) = charge_count(1) + particles(current_index)%data%q
               if (particles(current_index)%data%q < 0.0) then
                 charge_count(3) = charge_count(3) + 1
               end if
               particles(current_index) = particles(target_swap)
               swapped_cnt = swapped_cnt + 1
             else
               particles(current_index) = particles(target_swap)
               swapped_cnt = swapped_cnt + 1
             end if
           end if
         else
           particles(current_index) = particles(target_swap)
           swapped_cnt = swapped_cnt + 1
         end if

         if (init_swap_cnt /= swapped_cnt) then
           call filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, charge_count, break_check)
         end if
       else
         if (radius < cyl_radius) then
           if (particles(current_index)%x(3) < -cyl_length)then
             if (radius < plate_radius) then
               charge_count(2) = charge_count(2) + particles(current_index)%data%q
               swapped_cnt = swapped_cnt + 1
             else
               swapped_cnt = swapped_cnt + 1
             end if
           else if (particles(current_index)%x(3) > 0.0) then
             if (radius < plate_radius) then
               charge_count(1) = charge_count(1) + particles(current_index)%data%q
               if (particles(current_index)%data%q < 0.0) then
                 charge_count(3) = charge_count(3) + 1
               end if
               swapped_cnt = swapped_cnt + 1
             else
               swapped_cnt = swapped_cnt + 1
             end if
           end if
         else
           swapped_cnt = swapped_cnt + 1
         end if
         break_check = 1
       end if

     case (3) ! torus shape
       x = particles(current_index)%x(1)
       y = particles(current_index)%x(2)
       z = particles(current_index)%x(3)
       radius = major_radius - sqrt(x**2 + y**2)
       m_radius = sqrt(radius**2 + z**2)

       buffer_zone = 0.75/(c*1e-12)

       if (current_index .lt. target_swap) then
         if (m_radius > (minor_radius+buffer_zone)) then
           ! if (particles(current_index)%data%q < 0.0_8) then
           !   charge_count(1) = charge_count(1) + particles(current_index)%data%q
           ! else
           !   charge_count(2) = charge_count(2) + particles(current_index)%data%q
           ! end if
           particles(current_index) = particles(target_swap)
           swapped_cnt = swapped_cnt + 1
         end if

         if (init_swap_cnt /= swapped_cnt) then
           call filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, charge_count, break_check)
         end if
       else
         if (m_radius > (minor_radius+buffer_zone)) then
           ! if (particles(current_index)%data%q < 0.0_8) then
           !   charge_count(1) = charge_count(1) + particles(current_index)%data%q
           ! else
           !   charge_count(2) = charge_count(2) + particles(current_index)%data%q
           ! end if
           swapped_cnt = swapped_cnt + 1
         end if
         break_check = 1
       end if

     end select
   end subroutine
end module
