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
!> integrator module
!>

module particles_resize
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   implicit none

   type, public :: linked_list_elem
      type(linked_list_elem), pointer :: next
      type(t_particle), allocatable :: tmp_particles(:)
   end type linked_list_elem

   type, public :: linked_list_CS
      type(linked_list_CS), pointer :: next_CS
      real(kind_physics), allocatable :: CS(:,:)
   end type linked_list_CS

contains
   subroutine allocate_ll_buffer(electron_num, temp_array)
      implicit none
      integer, intent(in) :: electron_num
      type(linked_list_elem), pointer, intent(inout) :: temp_array
      integer :: buffer_size, i

      buffer_size = 50 !electron_num

      allocate (temp_array)
      allocate (temp_array%tmp_particles(buffer_size))
      ! print *, "buffer allocated!!! ", size(temp_array%tmp_particles)
      nullify (temp_array%next)
   end subroutine allocate_ll_buffer

   subroutine deallocate_ll_buffer(temp_array)
      implicit none
      type(linked_list_elem), pointer, intent(inout) :: temp_array
      type(linked_list_elem), pointer :: temp_guide, remover

      temp_guide => temp_array
      do while (associated(temp_guide))
         remover => temp_guide
         deallocate (temp_guide%tmp_particles)
         temp_guide => temp_guide%next
         deallocate (remover)
      end do
      nullify (temp_guide)
      nullify (remover)
   end subroutine deallocate_ll_buffer

   subroutine deallocate_CS_buffer(CS_array)
      implicit none
      type(linked_list_CS), pointer, intent(inout) :: CS_array
      type(linked_list_CS), pointer :: temp_guide, remover

      temp_guide => CS_array
      do while (associated(temp_guide))
         remover => temp_guide
         deallocate (temp_guide%CS)
         temp_guide => temp_guide%next_CS
         deallocate (remover)
      end do
      nullify (temp_guide)
      nullify (remover)
   end subroutine deallocate_CS_buffer

   subroutine resize_array(array, new_size)
      implicit none
      type(t_particle), allocatable, intent(inout) :: array(:)
      type(t_particle), allocatable :: temp_array(:)
      integer, intent(in) :: new_size
      integer :: original_size

      original_size = size(array)

      allocate (temp_array(original_size))
      temp_array = array
      deallocate (array)
      allocate (array(new_size))
      array(1:original_size) = temp_array

      deallocate (temp_array)
   end subroutine resize_array

   ! NOTE: OBSOLETE SUBROUTINE
   subroutine extend_particles_list(particles, buffer, new_particles_size)
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:)
      type(linked_list_elem), pointer, intent(inout) :: buffer
      type(linked_list_elem), pointer :: temp_guide
      integer, intent(in) :: new_particles_size
      integer :: new_size, head, tail, ll_buffer_size, i, remainder, ll_elem_cnt

      ll_buffer_size = size(buffer%tmp_particles)
      ll_elem_cnt = CEILING(real(new_particles_size)/real(ll_buffer_size))
      remainder = MOD(new_particles_size, ll_buffer_size)

      new_size = size(particles) + new_particles_size
      i = 1
      head = size(particles) + 1
      tail = size(particles) + ll_buffer_size

      print *, "Linked List elements: ", ll_elem_cnt, size(particles), new_particles_size, remainder

      call resize_array(particles, new_size)

      temp_guide => buffer
      do while (associated(temp_guide))
         if (i /= ll_elem_cnt) then
            particles(head:tail) = temp_guide%tmp_particles
         else if ((i == ll_elem_cnt) .and. (remainder == 0)) then
            particles(head:new_size) = temp_guide%tmp_particles
         else
            particles(head:new_size) = temp_guide%tmp_particles(1:(remainder))
         end if

         head = tail + 1
         tail = tail + ll_buffer_size
         i = i + 1

         temp_guide => temp_guide%next
      end do
      nullify (temp_guide)

   end subroutine extend_particles_list
end module
