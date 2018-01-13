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

   subroutine rotation(vector, theta_x, theta_y, theta_z)
      implicit none
      real*8, intent(inout) :: vector(:)
      real*8, intent(in) :: theta_x, theta_y, theta_z
      real*8, allocatable :: rotation_x(:, :), rotation_y(:, :), rotation_z(:, :), vec_prod(:)

      allocate (rotation_x(3, 3))
      allocate (rotation_y(3, 3))
      allocate (rotation_z(3, 3))
      allocate (vec_prod(3))

      rotation_x = 0.0_8
      rotation_y = 0.0_8
      rotation_z = 0.0_8

      rotation_x(1, 1) = 1.0_8
      rotation_x(2, 2) = cos(theta_x)
      rotation_x(2, 3) = -sin(theta_x)
      rotation_x(3, 2) = sin(theta_x)
      rotation_x(3, 3) = cos(theta_x)

      rotation_y(1, 1) = cos(theta_y)
      rotation_y(1, 3) = sin(theta_y)
      rotation_y(2, 2) = 1.0_8
      rotation_y(3, 1) = -sin(theta_y)
      rotation_y(3, 3) = cos(theta_y)

      rotation_z(1, 1) = cos(theta_z)
      rotation_z(1, 2) = -sin(theta_z)
      rotation_z(2, 1) = sin(theta_z)
      rotation_z(2, 2) = cos(theta_z)
      rotation_z(3, 3) = 1.0_8

      vec_prod = vector
      vector = matrix_vector_multiplication(rotation_x, vec_prod)
      vec_prod = vector
      vector = matrix_vector_multiplication(rotation_y, vec_prod)
      vec_prod = vector
      vector = matrix_vector_multiplication(rotation_z, vec_prod)
   end subroutine rotation

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

      ll_elem_gen = MOD(new_particle, size(guide%tmp_particles))
      buffer_pos = ll_elem_gen + 1
      ! print *, "Buffer_pos: ", buffer_pos, size(guide%tmp_particles)

      if (particle%data%age > 0.3_8) then
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

      Vm = particle%data%v + particle%results%e*q_m_ratio
      tan_w = particle%data%b*q_m_ratio
      sin_w = 2.0_8*tan_w/(1 + dot_product(tan_w, tan_w))

      V_cross = cross_product(Vm, tan_w)
      Vd = Vm + V_cross
      V_cross = cross_product(Vd, sin_w)
      Vp = Vm + V_cross

      particle%data%v = Vp + particle%results%e*q_m_ratio
   end subroutine boris_velocity_update

   subroutine particle_pusher(particle, dt)
      implicit none
      type(t_particle), intent(inout) :: particle
      real*8, intent(in) :: dt
      real*8 :: dist_vec(3)

      dist_vec = particle%data%v*dt
      particle%x = particle%x + dist_vec
      particle%data%age = particle%data%age + dot_product(dist_vec, dist_vec)

   end subroutine particle_pusher
end module
