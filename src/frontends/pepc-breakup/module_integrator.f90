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

   subroutine particle_EB_field(particle, E_field)
      implicit none
      type(t_particle), intent(inout) :: particle
      real(kind_physics) :: B_field(3), E_field(3)
      ! integer :: i

      B_field = 0.0

      !TODO add details to calculate the external E_field & B_field experienced by particles
      particle%results%e = particle%results%e + E_field
      particle%data%b = B_field

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

         ! call rotation(particle%data%v, 25.0, 25.0, 25.0)
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

   subroutine set_cross_section_table(fname, guide_CS, file_id, is_last)
     implicit none
     type(linked_list_CS), pointer, intent(inout) :: guide_CS
     character(len = *), intent(in) :: fname
     integer, intent(in) :: file_id, is_last
     integer :: entries, i
     type(linked_list_CS), pointer :: temp_guide

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
       nullify(guide_CS%next_CS)
     end if

     close(file_id)
   end subroutine set_cross_section_table

   subroutine determine_absolute_max_CS(guide_CS, max_CS)
     implicit none
     type(linked_list_CS), pointer, intent(in) :: guide_CS
     real(kind_physics), intent(inout) :: max_CS
     type(linked_list_CS), pointer :: temp_guide

     max_CS = 0.0
     temp_guide => guide_CS
     do while (associated(temp_guide))
       max_CS = max_CS + maxval(temp_guide%CS(:,2))
       temp_guide => temp_guide%next_CS
     end do
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

   subroutine add_particle(guide, particle, vel_mag, new_particle, rand1, rand2, buffer_pos, type)
     implicit none
     type(t_particle), intent(in) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     real(kind_physics), intent(in) :: vel_mag, rand1, rand2
     integer, intent(inout) :: new_particle, buffer_pos
     integer, intent(in) :: type
     real(kind_physics) :: vel_update, theta, phi
     integer :: ll_elem_gen

     ll_elem_gen = MOD(new_particle, size(guide%tmp_particles))
     buffer_pos = ll_elem_gen + 1
     if (ll_elem_gen == 0 .and. new_particle > 0) then
        ! print *, "Extending Temp_array!!!!!"
        call allocate_ll_buffer(size(guide%tmp_particles), guide%next)
        guide => guide%next
     end if
     ! theta defines inclination, phi defines azimuth in spherical coord. system
     theta = rand1*pi
     phi = rand2*pi*2.0

     select case(type)
     case(0) ! Electron
       vel_update = vel_mag

       guide%tmp_particles(buffer_pos)%data%q = -1.0
       guide%tmp_particles(buffer_pos)%data%m = 1.0_8
       guide%tmp_particles(buffer_pos)%data%species = 0

     case(1) ! H+
       vel_update = vel_mag/1836.21957489_8

       guide%tmp_particles(buffer_pos)%data%q = 1.0
       guide%tmp_particles(buffer_pos)%data%m = 1836.21957489_8 ! times greater than mass of electron
       guide%tmp_particles(buffer_pos)%data%species = 1

     case(2) ! H2+
       vel_update = vel_mag/3673.43889456

       guide%tmp_particles(buffer_pos)%data%q = 1.0
       guide%tmp_particles(buffer_pos)%data%m = 3673.43889456_8 ! times greater than mass of electron
       guide%tmp_particles(buffer_pos)%data%species = 2

     case(3) ! H atom
       vel_update = vel_mag/1837.21957489_8

       guide%tmp_particles(buffer_pos)%data%q = 0.0
       guide%tmp_particles(buffer_pos)%data%m = 1837.21957489_8 ! times greater than mass of electron
       guide%tmp_particles(buffer_pos)%data%species = 3

     end select

     guide%tmp_particles(buffer_pos)%x = particle%x
     guide%tmp_particles(buffer_pos)%work = 1.0_8
     guide%tmp_particles(buffer_pos)%data%age = 0.0_8

     guide%tmp_particles(buffer_pos)%data%v(1) = vel_update*sin(theta)*cos(phi)
     guide%tmp_particles(buffer_pos)%data%v(2) = vel_update*sin(theta)*sin(phi)
     guide%tmp_particles(buffer_pos)%data%v(3) = vel_update*cos(theta)

     new_particle = new_particle + 1
   end subroutine add_particle

   subroutine collision_update(particle, guide, new_particle, electron_count, CS_vector)
     implicit none
     type(t_particle), intent(inout) :: particle
     type(linked_list_elem), pointer, intent(inout) :: guide
     integer, intent(inout) :: new_particle, electron_count
     real(kind_physics), dimension(:), intent(inout) :: CS_vector
     real(kind_physics) :: vel_mag, nu_prime
     integer :: buff_pos, ll_elem_gen, i

     ! NOTE: Currently, calculation of nu_prime involves obtaining abs_max_CS,
     !       which is the sum of all max value of cross section data of all considered
     !       reactions, disregarding the associated eV. nu_prime is then obtained
     !       by multiplying abs_max_CS with the particle's velocity magnitude and
     !       constant local density. nu_prime will thus be always larger than actual nu!
     vel_mag = sqrt(dot_product(particle%data%v,particle%data%v))
     nu_prime = abs_max_CS * vel_mag * 6545520.13889

     ! Seeding procedure for RNG (any expression that generates integer unique to the process works)
     ctr_s(1) = (my_rank + 1)*(nt - step)
     ctr_s(2:4) = CEILING(particle%x*1e5)
     key_s(1) = (my_rank + 1)*(step + 1)
     key_s(2:4) = CEILING(particle%data%v*1e17)

     ! Generating Random Number between [0,1]
     dummy = gen_norm_double_rng(ctr_s, key_s, rand_num)

    !  print *, "rand: ", rand_num(1), " expression: ", (1 - exp(-1*CS_vector(size(CS_vector))*dt))
     if (rand_num(1) < (1 - exp(-1*nu_prime*dt))) then ! type of collision determined if satisfied
       call determine_cross_sections(particle, CS_vector, CS_tables)

       ! Convert CS_vector (cross section of all reactions) to Collision freq., nu.
       CS_vector = CS_vector * vel_mag * 6545520.13889 ! test value of constant local_number_density (at 0.001Pa)

       i = 1
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

     select case(i)
     case(0) ! null collision, no update performed

     case(1) ! elastic scattering (no additional electron, no byproducts)
       ! update velocity to indicate scattering into random angle
       particle%data%v(1) = vel_mag * sin(rand_num(3)*pi) * cos(rand_num(4)*pi*2.0)
       particle%data%v(2) = vel_mag * sin(rand_num(3)*pi) * sin(rand_num(4)*pi*2.0)
       particle%data%v(3) = vel_mag * cos(rand_num(3)*pi)
       particle%data%age = 0.0_8

     case(2) ! nondissociative ionization (1 additional electron, 1 byproduct)
       call add_particle(guide, particle, vel_mag, new_particle, rand_num(3), rand_num(4), buff_pos, 0)
       guide%tmp_particles(buff_pos)%data%v = guide%tmp_particles(buff_pos)%data%v/3._8
       call add_particle(guide, particle, vel_mag, new_particle, rand_num(5), rand_num(6), buff_pos, 2)
       guide%tmp_particles(buff_pos)%data%v = guide%tmp_particles(buff_pos)%data%v/3._8

       particle%data%v(1) = vel_mag * sin(rand_num(7)*pi) * cos(rand_num(8)*pi*2.0)/3._8
       particle%data%v(2) = vel_mag * sin(rand_num(7)*pi) * sin(rand_num(8)*pi*2.0)/3._8
       particle%data%v(3) = vel_mag * cos(rand_num(7)*pi)/3._8
       particle%data%age = 0.0_8

     case(3) ! dissociative ionization (1 additional electron, 2 byproducts), Hydrogen atom is ignored!
       call add_particle(guide, particle, vel_mag, new_particle, rand_num(3), rand_num(4), buff_pos, 0)
       guide%tmp_particles(buff_pos)%data%v = guide%tmp_particles(buff_pos)%data%v*0.25_8
       call add_particle(guide, particle, vel_mag, new_particle, rand_num(5), rand_num(6), buff_pos, 1)
       guide%tmp_particles(buff_pos)%data%v = guide%tmp_particles(buff_pos)%data%v*0.25_8

       particle%data%v(1) = vel_mag * sin(rand_num(7)*pi) * cos(rand_num(8)*pi*2.0) *0.25_8
       particle%data%v(2) = vel_mag * sin(rand_num(7)*pi) * sin(rand_num(8)*pi*2.0) *0.25_8
       particle%data%v(3) = vel_mag * cos(rand_num(7)*pi) *0.25_8
       particle%data%age = 0.0_8

     end select
   end subroutine collision_update
end module
