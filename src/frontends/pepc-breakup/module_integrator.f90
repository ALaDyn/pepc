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

   ! Given the velocity vector of particle, it's azimuthal and inclination angle is calculated.
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
      particle%data%b = Pol_B_field + B_field
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
     real(kind_physics) :: energy
     type(linked_list_CS), pointer :: temp_guide

     CS_table_no = 1
     temp_guide => guide_CS
     do while (associated(temp_guide))
       table_rows = size(temp_guide%CS,1)
       energy = 0.5 * particle%data%m * e_mass * dot_product(particle%data%v,particle%data%v)
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

   subroutine add_particle(guide, particle, new_particle, buffer_pos, type)
     ! NOTE: A big assumption is made here: Without considering the rovibrational states
     !       of the reaction products, generated particle doesn't store their own internal
     !       energy value. Instead, the reaction energy is entirely converted into kinetic E.
     !       More traditionally, this redistribution of energy is termed as Kinetic Energy Release
     !       in mass spetroscopy
     implicit none
     type(t_particle), intent(in) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     integer, intent(inout) :: new_particle, buffer_pos
     integer, intent(in) :: type
     integer :: ll_elem_gen

     ll_elem_gen = MOD(new_particle, size(guide%tmp_particles))
     buffer_pos = ll_elem_gen + 1
     if (ll_elem_gen == 0 .and. new_particle > 0) then
        ! print *, "Extending Temp_array!!!!!"
        call allocate_ll_buffer(size(guide%tmp_particles), guide%next)
        guide => guide%next
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

     guide%tmp_particles(buffer_pos)%x = particle%x
     guide%tmp_particles(buffer_pos)%work = 1.0_8
     guide%tmp_particles(buffer_pos)%data%age = 0.0_8

     guide%tmp_particles(buffer_pos)%label = 0
     guide%tmp_particles(buffer_pos)%data%v = 0.0_8

     new_particle = new_particle + 1
   end subroutine add_particle

   subroutine collision_update(particle, guide, new_particle, electron_count, CS_numbers)
     implicit none
     type(t_particle), intent(inout) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     integer, intent(inout) :: new_particle, electron_count
     integer, intent(in) :: CS_numbers
     real(kind_physics), dimension(:), allocatable :: CS_vector
     real(kind_physics) :: vel_mag, nu_prime, reduced_vel_mag, R_J02, V_V01, IE_H2_ion, AE_H_ion, theta, phi, polar_theta, polar_phi
     real(kind_physics) :: rot_axis(3), temp_vel(3), temp_vel1(3), reduced_incident(3), cos_theta, temp_vel_mag, temp_vel1_mag
     real(kind_physics) :: H2_mass, K_E, cos_Chi, Chi, unit_inc_vel(3), chi_rotation_axis(3), prefac, scatter_loss, xi
     integer :: buff_pos, i

     ! TODO: not a huge fan of the next 4 lines. Rather make those parameters which might help optimisation
     R_J02 = 0.0441 ! eV, transition energy associated with rotational excitation from J = 0 -> 2
     V_V01 = 0.516 ! eV, transition energy associated with vibrational excitation from V = 0 -> 1
     IE_H2_ion = 15.426 ! 15.283 eV, Ionization energy of H2+ [Source: T.E.Sharp Atomic Data 2, 119-169 (1971)]. 15.426 by Yoon 2008
     AE_H_ion = 18.075 ! eV, Appearance energy of H+, 18.1 by Yoon 2008
     H2_mass = 3674.43889456_8 ! H2/e mass ratio
     ! [definition of Appearance energy vs Ionization energy from Mass Spectroscopy (2011) by Gross J.H.]

     ! NOTE: Currently, calculation of nu_prime involves obtaining abs_max_CS,
     !       which is the max of total cross section. nu_prime is then obtained
     !       by multiplying abs_max_CS with the particle's velocity magnitude and
     !       constant local density. nu_prime will thus be always larger than actual nu!
     ! TODO: not at all critical, but there are many ways of doing the next line, would be interesting to see if there is a fastest
     ! way of computing the length and weigh that against readability, perhaps a separate inlinable function if there is a faster
     ! way?
     vel_mag = sqrt(dot_product(particle%data%v,particle%data%v))
     K_E = 0.5*particle%data%m*vel_mag**2
     nu_prime = abs_max_CS * vel_mag * neutral_density
    !  call determine_cross_sections(particle, CS_vector, CS_tables)
    !  CS_vector = CS_vector * vel_mag * neutral_density
    !  print *, nu_prime, (1 - exp(-1*nu_prime*dt)), CS_vector/nu_prime

     ! Generating Random Number between [0,1]
     dummy = gen_norm_double_rng(ctr_s, key_s, rand_num)

    !  print *, "rand: ", rand_num(1), " expression: ", (1 - exp(-1*nu_prime*dt))!(1 - exp(-1*CS_vector(size(CS_vector))*dt))
     if (rand_num(1) < (1 - exp(-1*nu_prime*dt))) then ! type of collision determined if satisfied
       allocate(CS_vector(CS_numbers))
       call determine_cross_sections(particle, CS_vector, CS_tables)
       call determine_xi(K_E, xi)

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
      !  print *, "incident momentum: ", particle%data%m * reduced_incident

       call add_particle(guide, particle, new_particle, buff_pos, 0)
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
       temp_vel_mag = (2.*particle%data%m/(guide%tmp_particles(buff_pos)%data%m + particle%data%m)) &
                      * reduced_vel_mag * cos_theta
       temp_vel1_mag = reduced_vel_mag*sqrt(particle%data%m**2 + 2. * particle%data%m * guide%tmp_particles(buff_pos)%data%m &
                       * (1. - 2. * (cos_theta**2)) + guide%tmp_particles(buff_pos)%data%m**2)/(guide%tmp_particles(buff_pos)%data%m + particle%data%m)
       guide%tmp_particles(buff_pos)%data%v = temp_vel * temp_vel_mag

       temp_vel1 = reduced_incident - (guide%tmp_particles(buff_pos)%data%m/particle%data%m) * &
                   guide%tmp_particles(buff_pos)%data%v

       particle%data%v = temp_vel1
       particle%data%age = 0.0_8
      !  print *, "outgoing momentum: ", particle%data%m*temp_vel1 + guide%tmp_particles(buff_pos)%data%v * guide%tmp_particles(buff_pos)%data%m
      !  print *, "incident kinetic energy: ", 0.5*reduced_vel_mag**2
      !  print *, "outgoing kinetic energy: ", 0.5*(temp_vel1_mag**2 + temp_vel_mag**2)

       call add_particle(guide, particle, new_particle, buff_pos, 2)
      !  print *, "nondissoc.!"

     case(5) ! dissociative ionization (1 additional electron, 2 byproducts), Hydrogen atom is ignored!
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*(AE_H_ion+scatter_loss)/e_mass)
       reduced_incident = particle%data%v * (reduced_vel_mag/vel_mag)
       call angles_calc(reduced_incident, reduced_vel_mag, theta, phi)
      !  print *, "incident momentum: ", particle%data%m * reduced_incident

       call add_particle(guide, particle, new_particle, buff_pos, 0)
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
       temp_vel_mag = (2.*particle%data%m/(guide%tmp_particles(buff_pos)%data%m + particle%data%m)) &
                      * reduced_vel_mag * cos_theta
       temp_vel1_mag = reduced_vel_mag*sqrt(particle%data%m**2 + 2. * particle%data%m * guide%tmp_particles(buff_pos)%data%m &
                       * (1. - 2. * (cos_theta**2)) + guide%tmp_particles(buff_pos)%data%m**2)/(guide%tmp_particles(buff_pos)%data%m + particle%data%m)
       guide%tmp_particles(buff_pos)%data%v = temp_vel * temp_vel_mag

       temp_vel1 = reduced_incident - (guide%tmp_particles(buff_pos)%data%m/particle%data%m) * &
                   guide%tmp_particles(buff_pos)%data%v

       particle%data%v = temp_vel1
       particle%data%age = 0.0_8
      !  print *, "outgoing momentum: ", particle%data%m*temp_vel1 + guide%tmp_particles(buff_pos)%data%v * guide%tmp_particles(buff_pos)%data%m
      !  print *, "incident kinetic energy: ", 0.5*reduced_vel_mag**2
      !  print *, "outgoing kinetic energy: ", 0.5*(temp_vel1_mag**2 + temp_vel_mag**2)

       call add_particle(guide, particle, new_particle, buff_pos, 1)
      !  print *, "dissoc.!"

     end select
   end subroutine collision_update

   recursive subroutine filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, break_check)
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, intent(inout) :: swapped_cnt, break_check
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
           call filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, break_check)
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
           call filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, break_check)
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
           call filter_and_swap(particles, geometry, current_index, head_i, tail_i, swapped_cnt, break_check)
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
