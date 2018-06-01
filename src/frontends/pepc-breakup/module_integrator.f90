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

   function circular_poloidal_field(particle_pos, vector_parallel_maj_radius, distance_xy) result(Poloidal_field)
     implicit none
     real(kind_physics), intent(in) :: particle_pos(3), vector_parallel_maj_radius(3)
     real(kind_physics), intent(in) :: distance_xy
     real(kind_physics) :: major_r_vec(3), a, Poloidal_field(3)

     ! first calculate the vector from center of cylinder (2D slice of torus)
     ! to the particle (in plane), defined as major_r_vec
     major_r_vec(1) = particle_pos(1)
     major_r_vec(2) = particle_pos(2)
     major_r_vec(3) = 0.0

     major_r_vec = major_r_vec*major_radius/distance_xy
     major_r_vec = particle_pos - major_r_vec
     a = sqrt(particle_pos(3)**2 + (major_radius - distance_xy)**2)

     Poloidal_field = cross_product(vector_parallel_maj_radius,major_r_vec)
     Poloidal_field = Poloidal_field/sqrt(dot_product(Poloidal_field, Poloidal_field))
     Poloidal_field = Poloidal_field*(0.0001*((1.e-12)*(c**2))/e_mass + a*B_p/minor_radius)
   end function circular_poloidal_field

   function Amperes_law(particle_pos, vector_parallel_maj_radius, distance_xy, coil_ver, coil_hor, I_pf) &
     result(Poloidal_field)
     implicit none
     real(kind_physics), intent(in) :: particle_pos(3), vector_parallel_maj_radius(3)
     real(kind_physics), intent(in) :: distance_xy, coil_ver, coil_hor, I_pf
     real(kind_physics) :: major_r_vec(3), a, Poloidal_field(3), mu_0, B_mag

     Poloidal_field = 0.0
     mu_0 = ((1./c)**2)/eps_0
     ! first calculate vector point from PF coil to particle position in 2D plane
     ! defined as major_r_vec
     major_r_vec(1) = particle_pos(1)
     major_r_vec(2) = particle_pos(2)
     major_r_vec(3) = 0.0

     ! coil position is calculated, remember to scale coil_hor & coil_ver
     major_r_vec = major_r_vec*coil_hor/(distance_xy * (c*1e-12))
     major_r_vec(3) = coil_ver/(c*1e-12)

     ! vector pointing from coil to particle position and distance between them is calculated
     major_r_vec = particle_pos - major_r_vec
     a = sqrt(dot_product(major_r_vec,major_r_vec)) ! NOTE: 'a' is dimensionless

     ! unit vector direction of resulting poloidal field is calculated
     Poloidal_field = cross_product(vector_parallel_maj_radius,major_r_vec)
     Poloidal_field = Poloidal_field/sqrt(dot_product(Poloidal_field, Poloidal_field))

     ! Calculate strength of magnetic field by PF coil
     a = a * (c*1e-12)
     B_mag = mu_0*I_pf/(2.*pi*a) ! NOTE: Calculation done with Tesla unit
     B_mag = B_mag * ((1.e-12)*(c**2))/e_mass ! NOTE: scale to dimensionless

     ! final vector of PF B field
     Poloidal_field = Poloidal_field * B_mag
   end function Amperes_law

   subroutine particle_EB_field(particle, E_field)
      implicit none
      type(t_particle), intent(inout) :: particle
      real(kind_physics), intent(in) :: E_field(3)
      real(kind_physics) :: B_field(3), ez(3), Pol_B_field(3), field_vector(3)
      real(kind_physics) :: R, PF_final(3), PF_temp(3), disp, Ipf, PF_unit_vector(3)

      ez = 0.0
      ez(3) = 1.0
      B_field = 0.0
      field_vector = 0.0

      field_vector = cross_product(particle%x,ez)
      field_vector = field_vector/sqrt(dot_product(field_vector, field_vector))

      R = sqrt(particle%x(1)**2 + particle%x(2)**2)
      B_field = field_vector*B0*major_radius/R

      ! Option 2: a circular poloidal magnetic field, centred around major_radius
      ! Pol_B_field = circular_poloidal_field(particle%x, field_vector, R)

      PF_final = 0.0
      disp = 10.0_8
      Ipf = 3.4e8_8
      ! Option 3: 4 correction coils surrounding breakdown region, centred around major_radius
      PF_temp = Amperes_law(particle%x, field_vector, R, disp, 2.9_8, Ipf)
      PF_final = PF_final + PF_temp
      PF_temp = Amperes_law(particle%x, field_vector, R, 0.0_8, 2.9_8 - disp, Ipf)
      PF_final = PF_final + PF_temp
      PF_temp = Amperes_law(particle%x, field_vector, R, -disp, 2.9_8, Ipf)
      PF_final = PF_final + PF_temp
      PF_temp = Amperes_law(particle%x, field_vector, R, 0.0_8, 2.9_8 + disp, Ipf)
      PF_final = PF_final + PF_temp
      PF_unit_vector = PF_final/sqrt(dot_product(PF_final, PF_final))
      PF_temp = 0.0005 * ((1.e-12)*(c**2))/e_mass * PF_unit_vector
      Pol_B_field = PF_final + PF_temp

      particle%results%e = particle%results%e + field_vector*V_loop/(2.*pi*R) + E_field
      ! print *, particle%results%e*0.0160217662080007054395368083795655167047391940703667
      particle%data%b = Pol_B_field !B_field +
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
      real*8 :: Vm(3), Vd(3), Vp(3), tan_w(3), sin_w(3), V_cross(3)
      real*8 :: half_dt, q_m_ratio

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
   end subroutine boris_velocity_update

   subroutine particle_pusher(particle, dt)
      implicit none
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: dt
      real*8 :: dist_vec(3)

      dist_vec = particle%data%v*dt
      particle%x = particle%x + dist_vec
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

     guide%tmp_particles(buffer_pos)%data%v = 0.0_8

     new_particle = new_particle + 1
   end subroutine add_particle

   subroutine collision_update(particle, guide, new_particle, electron_count, CS_vector)
     implicit none
     type(t_particle), intent(inout) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     integer, intent(inout) :: new_particle, electron_count
     real(kind_physics), dimension(:), intent(inout) :: CS_vector
     real(kind_physics) :: vel_mag, nu_prime, reduced_vel_mag, R_J02, V_V01, IE_H2_ion, AE_H_ion, theta, phi, polar_theta, polar_phi
     real(kind_physics) :: rot_axis(3), temp_vel(3), temp_vel1(3), reduced_incident(3), cos_theta, temp_vel_mag, temp_vel1_mag
     integer :: buff_pos, i

     ! TODO: not a huge fan of the next 4 lines. Rather make those parameters which might help optimisation
     R_J02 = 0.0441 ! eV, transition energy associated with rotational excitation from J = 0 -> 2
     V_V01 = 0.516 ! eV, transition energy associated with vibrational excitation from V = 0 -> 1
     IE_H2_ion = 15.426 ! 15.283 eV, Ionization energy of H2+ [Source: T.E.Sharp Atomic Data 2, 119-169 (1971)]. 15.426 by Yoon 2008
     AE_H_ion = 18.075 ! eV, Appearance energy of H+, 18.1 by Yoon 2008
     ! [definition of Appearance energy vs Ionization energy from Mass Spectroscopy (2011) by Gross J.H.]

     ! NOTE: Currently, calculation of nu_prime involves obtaining abs_max_CS,
     !       which is the sum of all max value of cross section data of all considered
     !       reactions, disregarding the associated eV. nu_prime is then obtained
     !       by multiplying abs_max_CS with the particle's velocity magnitude and
     !       constant local density. nu_prime will thus be always larger than actual nu!
     ! TODO: not at all critical, but there are many ways of doing the next line, would be interesting to see if there is a fastest
     ! way of computing the length and weigh that against readability, perhaps a separate inlinable function if there is a faster
     ! way?
     vel_mag = sqrt(dot_product(particle%data%v,particle%data%v))
     nu_prime = abs_max_CS * vel_mag * neutral_density
    !  call determine_cross_sections(particle, CS_vector, CS_tables)
    !  CS_vector = CS_vector * vel_mag * neutral_density
    !  print *, nu_prime, (1 - exp(-1*nu_prime*dt)), CS_vector/nu_prime

     ! Seeding procedure for RNG (any expression that generates integer unique to the process works)
     ctr_s(1) = (my_rank + 1)*(nt - step)
     ctr_s(2:4) = CEILING(particle%x*1e5)
     key_s(1) = (my_rank + 1)*(step + 1)
     key_s(2:4) = CEILING(particle%data%v*1e17)

     ! Generating Random Number between [0,1]
     dummy = gen_norm_double_rng(ctr_s, key_s, rand_num)

    !  print *, "rand: ", rand_num(1), " expression: ", (1 - exp(-1*nu_prime*dt))!(1 - exp(-1*CS_vector(size(CS_vector))*dt))
     if (rand_num(1) < (1 - exp(-1*nu_prime*dt))) then ! type of collision determined if satisfied
       call determine_cross_sections(particle, CS_vector, CS_tables)

       ! Convert CS_vector (cross section of all reactions) to Collision freq., nu.
       CS_vector = CS_vector * vel_mag * neutral_density !6545520.13889 test value of constant local_number_density (at 0.001Pa)

       ! TODO: this to me looks like it could be a simple do-loop
       ! TODO: also looks like it can be combined with the select case below via a rather long/convoluted if-elseif construct that
       ! may read more easily
       i = 1
       ! TODO: CHECK: why is this rand_num(2) and not (1) again?
       do while (rand_num(2) > (CS_vector(i)/nu_prime))
         i = i + 1
         if (i > size(CS_vector)) then
           i = 0
           EXIT
         end if
       end do

     else ! otherwise, no collision happened
       i = 0
     end if

     polar_theta = 2.0*pi*rand_num(4)
     polar_phi = acos(2.0*rand_num(3) - 1.)

     select case(i)
     case(0) ! null collision, no update performed
      !  print *, "null coll!"
     case(1) ! elastic scattering (no additional electron, no byproducts)
       ! update velocity to indicate scattering into random angle
       ! TODO: ERROR: this distribution will not be uniform on the unit sphere, instead it will cause a bias towards front/back scattering
       ! TODO: c.f. http://mathworld.wolfram.com/SpherePointPicking.html
       ! TODO: c.f. https://corysimon.github.io/articles/uniformdistn-on-sphere/
       particle%data%v(1) = vel_mag * sin(polar_phi) * cos(polar_theta)
       particle%data%v(2) = vel_mag * sin(polar_phi) * sin(polar_theta)
       particle%data%v(3) = vel_mag * cos(polar_phi)
       particle%data%age = 0.0_8
      !  print *, "elastic coll.!"

     case(2) ! rotational excitation of H2 molecule, electron will lose the transition energy
       ! update velocity to indicate scattering into random angle
       ! TODO: will be less costly to keep 2/e_mass as parameter and multiplying with it - not sure the compile will pick up on this
       ! and do the short-cut for you. Better still: add it to the energies straight away.
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*R_J02/e_mass)
       particle%data%v(1) = reduced_vel_mag * sin(polar_phi) * cos(polar_theta)
       particle%data%v(2) = reduced_vel_mag * sin(polar_phi) * sin(polar_theta)
       particle%data%v(3) = reduced_vel_mag * cos(polar_phi)
       particle%data%age = 0.0_8
      !  print *, "rotational exci.!"

     case(3) ! vibrational excitation of H2 molecule, electron will lose the transition energy
       ! update velocity to indicate scattering into random angle
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*V_V01/e_mass)
       particle%data%v(1) = reduced_vel_mag * sin(polar_phi) * cos(polar_theta)
       particle%data%v(2) = reduced_vel_mag * sin(polar_phi) * sin(polar_theta)
       particle%data%v(3) = reduced_vel_mag * cos(polar_phi)
       particle%data%age = 0.0_8
      !  print *, "vibrational exci.!"

     case(4) ! nondissociative ionization (1 additional electron, 1 byproduct)
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*IE_H2_ion/e_mass)
       reduced_incident = particle%data%v * (reduced_vel_mag/vel_mag)
       call angles_calc(reduced_incident, reduced_vel_mag, theta, phi)
      !  print *, "incident momentum: ", particle%data%m * reduced_incident

       call add_particle(guide, particle, new_particle, buff_pos, 0)
       temp_vel(1) = sin(polar_phi*0.5) * cos(polar_theta)
       temp_vel(2) = sin(polar_phi*0.5) * sin(polar_theta)
       temp_vel(3) = cos(polar_phi*0.5)

       rot_axis(1) = 0.0
       rot_axis(2) = 1.0
       rot_axis(3) = 0.0
       temp_vel = Rodriguez_rotation(theta, rot_axis, temp_vel)
       rot_axis(2) = 0.0
       rot_axis(3) = 1.0
       temp_vel = Rodriguez_rotation(phi, rot_axis, temp_vel)

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
       reduced_vel_mag = sqrt(vel_mag**2 - 2.*AE_H_ion/e_mass)
       reduced_incident = particle%data%v * (reduced_vel_mag/vel_mag)
       call angles_calc(reduced_incident, reduced_vel_mag, theta, phi)
      !  print *, "incident momentum: ", particle%data%m * reduced_incident

       call add_particle(guide, particle, new_particle, buff_pos, 0)

       temp_vel(1) = sin(polar_phi*0.5) * cos(polar_theta)
       temp_vel(2) = sin(polar_phi*0.5) * sin(polar_theta)
       temp_vel(3) = cos(polar_phi*0.5)

       rot_axis(1) = 0.0
       rot_axis(2) = 1.0
       rot_axis(3) = 0.0
       temp_vel = Rodriguez_rotation(theta, rot_axis, temp_vel)
       rot_axis(2) = 0.0
       rot_axis(3) = 1.0
       temp_vel = Rodriguez_rotation(phi, rot_axis, temp_vel)

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

   subroutine extend_particles_list_add_e(particles, buffer, new_particles_size, electron_number)
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:)
      type(linked_list_elem), pointer, intent(inout) :: buffer
      type(linked_list_elem), pointer :: temp_guide
      integer, intent(in) :: new_particles_size, electron_number
      integer :: new_size, head, tail, ll_buffer_size, i, remainder, ll_elem_cnt, &
                 temp_size
      real(kind_physics) :: center_pos(3), plane_orient(3)

      ll_buffer_size = size(buffer%tmp_particles)
      ll_elem_cnt = CEILING(real(new_particles_size)/real(ll_buffer_size))
      remainder = MOD(new_particles_size, ll_buffer_size)

      temp_size = size(particles) + new_particles_size
      new_size = temp_size + electron_number
      i = 1
      head = size(particles) + 1
      tail = size(particles) + ll_buffer_size

      print *, "Linked List elements: ", ll_elem_cnt, size(particles), new_particles_size + electron_number, remainder

      call resize_array(particles, new_size)

      ! NOTE: if there is no new particle being generated, the following code section will cause fatal error!
      if (new_particles_size /= 0) then
        temp_guide => buffer
        do while (associated(temp_guide))
           if (i /= ll_elem_cnt) then
              particles(head:tail) = temp_guide%tmp_particles
           else if ((i == ll_elem_cnt) .and. (remainder == 0)) then
              particles(head:temp_size) = temp_guide%tmp_particles
           else
              particles(head:temp_size) = temp_guide%tmp_particles(1:(remainder))
           end if

           head = tail + 1
           tail = tail + ll_buffer_size
           i = i + 1

           temp_guide => temp_guide%next
        end do
        nullify (temp_guide)
      end if

      center_pos = 0.5/(c*1e-12)
      center_pos(3) = 0.0
      plane_orient = 0.0
      plane_orient(3) = -1.0
      call injected_electrons(center_pos, plane_orient, 2, 0.15/(c*1e-12), &
                                particles, temp_size + 1)
   end subroutine extend_particles_list_add_e

   subroutine extend_particles_list_swap_inject(particles, buffer, new_particles_size, electron_number, swapped_cnt)
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:)
      type(linked_list_elem), pointer, intent(inout) :: buffer
      type(linked_list_elem), pointer :: temp_guide
      integer, intent(in) :: new_particles_size, electron_number, swapped_cnt
      integer :: new_size, head, tail, ll_buffer_size, i, remainder, ll_elem_cnt, &
                 temp_size
      real(kind_physics) :: center_pos(3), plane_orient(3)

      ll_buffer_size = size(buffer%tmp_particles)
      ll_elem_cnt = CEILING(real(new_particles_size)/real(ll_buffer_size))
      remainder = MOD(new_particles_size, ll_buffer_size)

      temp_size = size(particles) + new_particles_size - swapped_cnt
      new_size = temp_size + electron_number
      i = 1
      head = size(particles) - swapped_cnt + 1
      tail = size(particles) - swapped_cnt + ll_buffer_size

      ! print *, "Linked List elements: ", ll_elem_cnt, size(particles), new_particles_size + electron_number, remainder

      call resize_array(particles, new_size)

      ! NOTE: if there is no new particle being generated, the following code section will cause fatal error!
      if (new_particles_size /= 0) then
        temp_guide => buffer
        do while (associated(temp_guide))
           if (i /= ll_elem_cnt) then
              particles(head:tail) = temp_guide%tmp_particles
           else if ((i == ll_elem_cnt) .and. (remainder == 0)) then
              particles(head:temp_size) = temp_guide%tmp_particles
           else
              particles(head:temp_size) = temp_guide%tmp_particles(1:(remainder))
           end if

           head = tail + 1
           tail = tail + ll_buffer_size
           i = i + 1

           temp_guide => temp_guide%next
        end do
        nullify (temp_guide)
      end if

      if (electron_number /= 0) then
        center_pos = 0.0
        center_pos(3) = -0.01
        plane_orient = 0.0
        plane_orient(3) = -1.0
        call injected_electrons(center_pos, plane_orient, 2, 0.025/(c*1e-12), &
                                  particles, temp_size + 1)
      end if
   end subroutine extend_particles_list_swap_inject

   subroutine injected_electrons(center_pos, plane, geometry, inlet_size, &
                                   particles_list, starting_index)
     ! NOTE: supports only planes in principal directions
     implicit none
     real(kind_physics), intent(in) :: plane(3), center_pos(3), inlet_size
     integer, intent(in) :: geometry, starting_index
     type(t_particle), allocatable, intent(inout) :: particles_list(:)
     real(kind_physics) :: ran(3), magnitude, theta, u !, mu, sigma2, res_mag, a
     integer :: i, j, index

    !  magnitude = thermal_velocity_mag(1.0_8, 1773.15_kind_physics)
     magnitude = sqrt((2*1.0)/e_mass)
     ! NOTE: below are the expressions for mean and variance of Maxwellian distribution.
    !  mu = magnitude
    !  a = sqrt(1773.15_kind_physics*8.61733035e-5_kind_physics/e_mass)
    !  sigma2 = a**2 * (3.*pi - 8.)/pi

    !  print *, "mean & variance: ", mu, sigma2
     select case(geometry)
     case(1) ! square plane
       do index = starting_index, size(particles_list)
         call random_number(ran)
        !  call frand123NormDouble( state, mu, sigma2, res_mag )
         particles_list(index)%data%q = -1.0_8
         particles_list(index)%data%m = 1.0_8
         particles_list(index)%data%species = 0
         particles_list(index)%data%age = 0.0_8
         particles_list(index)%work = 1.0_8

         particles_list(index)%x = center_pos
         particles_list(index)%data%v = 0.0

         do i = 1, 3
           if (plane(i) == 0.0) then
             particles_list(index)%x(i) = ran(i) * inlet_size + center_pos(i) - inlet_size*0.5
           else
             particles_list(index)%data%v(i) = plane(i) * magnitude
           end if
         end do
       end do

     case(2) ! circle plane
       do index = starting_index, size(particles_list)
         call random_number(ran)
        !  call frand123NormDouble( state, mu, sigma2, res_mag )
         particles_list(index)%data%q = -1.0_8
         particles_list(index)%data%m = 1.0_8
         particles_list(index)%data%species = 0
         particles_list(index)%data%age = 0.0_8
         particles_list(index)%work = 1.0_8

         particles_list(index)%x = center_pos
         particles_list(index)%data%v = plane * magnitude

         theta = ran(2)*2.*pi
         u = ran(1) + ran(3)
         if (u > 1.) then
           u = 2. - u
         end if

         particles_list(index)%x(1) = u * sin(theta) * inlet_size + center_pos(1)
         particles_list(index)%x(2) = u * cos(theta) * inlet_size + center_pos(2)
        !  print *, "res_mag = ", res_mag
        !  j = 0
        !  do i = 1, 3
        !    if ((plane(i) <= 0.00001) .and. (j == 0)) then
        !      particles_list(index)%x(i) = u * sin(theta) * inlet_size + center_pos(i)
        !      j = 1
        !    else if ((plane(i) <= 0.00001) .and. (j == 1)) then
        !      particles_list(index)%x(i) = u * cos(theta) * inlet_size + center_pos(i)
        !    end if
        !  end do
       end do
     end select
   end subroutine injected_electrons

   recursive subroutine filter_and_swap(particles, geometry, current_index, swapped_cnt)
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, intent(inout) :: swapped_cnt
     integer, intent(in) :: geometry, current_index
     integer :: buffer_size, target_swap, init_swap_cnt
     real(kind_physics) :: center_pos(3)
     real(kind_physics) :: radius, cyl_radius, plate_radius, cyl_length, box_dim, x, y, box_l

     buffer_size = size(particles)
     target_swap = buffer_size - swapped_cnt

     center_pos = 0.0
     center_pos(3) = 0.0

     x = particles(current_index)%x(1) - center_pos(1)
     y = particles(current_index)%x(2) - center_pos(2)

     init_swap_cnt = swapped_cnt

     select case(geometry)
     case(1) ! box boundary
       box_dim = 1.0/(c*1e-12)
       box_l = 0.5/(c*1e-12)

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
         call filter_and_swap(particles, geometry, current_index, swapped_cnt)
       end if

     case(2) ! cylinder boundary
       cyl_radius = 0.05/(c*1e-12)
       plate_radius = 0.025/(c*1e-12)
       cyl_length = 0.0011870845/(c*1e-12)
       radius = sqrt(x**2 + y**2)

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
             particles(current_index) = particles(target_swap)
             swapped_cnt = swapped_cnt + 1
             if (particles(current_index)%data%q < 0.0) then
               charge_count(3) = charge_count(3) + 1
             end if
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
         call filter_and_swap(particles, geometry, current_index, swapped_cnt)
       end if

     end select
   end subroutine
end module
