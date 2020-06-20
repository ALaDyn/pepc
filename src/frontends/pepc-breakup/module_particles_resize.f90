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
!> particles module
!>

module particles_resize
   use helper
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   implicit none

contains
   subroutine allocate_ll_buffer(electron_num, temp_array)
      implicit none
      integer, intent(in) :: electron_num
      type(linked_list_elem), pointer, intent(inout) :: temp_array
      integer :: buffer_size, i

      buffer_size = electron_num

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

      if (new_size < original_size) then
        array(1:new_size) = temp_array(1:new_size)
      else
        array(1:original_size) = temp_array
      end if

      deallocate (temp_array)
   end subroutine resize_array

   subroutine sort_sibling_species(particles, sibling_cnt, parent_no, unique_parents)
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, allocatable, intent(in) :: sibling_cnt(:)
     integer, intent(in) :: parent_no, unique_parents
     type(t_particle), allocatable :: sorting_buffer(:),  copied_local_buffer(:)
     integer :: i, j, k, istart, istop
     integer :: sort_width, iLeft, iRight, iEnd, sorted_i, buffer_size

     istart = 1
     do i = 2, parent_no
       istart = istart + sibling_cnt(i-1)
     end do
     if (parent_no .eq. 1) istart = 1
     istop = istart + sibling_cnt(parent_no) - 1

     ! print *, "particle parent: ", parent_no, particles(istart)%data%mp_int1, particles(istop)%data%mp_int1

     buffer_size = istop - istart + 1
     allocate(sorting_buffer(buffer_size))
     allocate(copied_local_buffer(buffer_size))

     copied_local_buffer = particles(istart:istop)

     !begin sorting
     sort_width = 1
     do while (sort_width < buffer_size)

       iLeft = 1
       do while(iLeft <= buffer_size)

         iRight = MIN(iLeft + sort_width, buffer_size)
         iEnd = MIN(iLeft + 2*sort_width-1, buffer_size)
         j = iLeft
         k = iRight
         do i = iLeft, iEnd !!!CHECK THE EXTENT, Also consider if n < largest sort_width
           if (j < iRight) then
             if (k > iEnd) then
               sorting_buffer(i) = copied_local_buffer(j)
               j = j + 1
             else if (copied_local_buffer(j)%data%species <= copied_local_buffer(k)%data%species) then
               sorting_buffer(i) = copied_local_buffer(j)
               j = j + 1
             else
               sorting_buffer(i) = copied_local_buffer(k)
               k = k + 1
             end if
           else
             sorting_buffer(i) = copied_local_buffer(k)
             k = k + 1
           end if
         end do
         iLeft = iLeft + 2*sort_width
       end do

       copied_local_buffer = sorting_buffer
       sort_width = 2*sort_width
     end do
     particles(istart:istop) = copied_local_buffer

     do i = istart, istop
       particles(i)%label = i
     end do

     deallocate(sorting_buffer)
     deallocate(copied_local_buffer)
   end subroutine

   subroutine sort_particles_by_key(particles)
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     type(t_particle), allocatable :: sorting_buffer(:)
     integer :: i, j, k
     integer :: sort_width, iLeft, iRight, iEnd, sorted_i, buffer_size

     buffer_size = size(particles)
     allocate(sorting_buffer(buffer_size))

     !begin sorting
     sort_width = 1
     do while (sort_width < buffer_size)

       iLeft = 1
       do while(iLeft <= buffer_size)

         iRight = MIN(iLeft + sort_width, buffer_size)
         iEnd = MIN(iLeft + 2*sort_width-1, buffer_size)
         j = iLeft
         k = iRight
         do i = iLeft, iEnd !!!CHECK THE EXTENT, Also consider if n < largest sort_width
           if (j < iRight) then
             if (k > iEnd) then
               sorting_buffer(i) = particles(j)
               j = j + 1
             else if (particles(j)%key <= particles(k)%key) then
               sorting_buffer(i) = particles(j)
               j = j + 1 
             else
               sorting_buffer(i) = particles(k)
               k = k + 1
             end if
           else
             sorting_buffer(i) = particles(k)
             k = k + 1
           end if
         end do
         iLeft = iLeft + 2*sort_width
       end do

       particles = sorting_buffer
       sort_width = 2*sort_width
     end do
     ! particles = copied_local_buffer
     ! do i = 1, buffer_size
     !   particles(i)%label = 0
     ! end do

     deallocate(sorting_buffer)
   end subroutine

   subroutine sort_particles_by_weight(particles)
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     type(t_particle), allocatable :: sorting_buffer(:)
     integer :: i, j, k
     integer :: sort_width, iLeft, iRight, iEnd, sorted_i, buffer_size

     buffer_size = size(particles)
     allocate(sorting_buffer(buffer_size))

     !begin sorting
     sort_width = 1
     do while (sort_width < buffer_size)

       iLeft = 1
       do while(iLeft <= buffer_size)

         iRight = MIN(iLeft + sort_width, buffer_size)
         iEnd = MIN(iLeft + 2*sort_width-1, buffer_size)
         j = iLeft
         k = iRight
         do i = iLeft, iEnd !!!CHECK THE EXTENT, Also consider if n < largest sort_width
           if (j < iRight) then
             if (k > iEnd) then
               sorting_buffer(i) = particles(j)
               j = j + 1
             else if (abs(particles(j)%data%q) <= abs(particles(k)%data%q)) then
               sorting_buffer(i) = particles(j)
               j = j + 1 
             else
               sorting_buffer(i) = particles(k)
               k = k + 1
             end if
           else
             sorting_buffer(i) = particles(k)
             k = k + 1
           end if
         end do
         iLeft = iLeft + 2*sort_width
       end do

       particles = sorting_buffer
       sort_width = 2*sort_width
     end do
     ! particles = copied_local_buffer
     ! do i = 1, buffer_size
     !   particles(i)%label = 0
     ! end do

     deallocate(sorting_buffer)
   end subroutine

   subroutine determine_siblings_at_level(particles, sibling_cnt, unique_parents, level)
     use treevars, only: maxlevel
     use module_spacefilling
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, allocatable, intent(inout) :: sibling_cnt(:)
     integer, intent(inout) :: unique_parents
     integer(kind_level), intent(in) :: level
     integer :: i, lower_bound, upper_bound
     integer(kind_level) :: level_march
     integer(kind_key) :: match_key, parent_key
     real(kind_physics) :: average_siblings

     allocate(sibling_cnt(size(particles)))
     sibling_cnt = 0
     level_march = level - maxlevel
     unique_parents = 1

     ! Compute the parent key for each particle at designated level.
     ! Since particles are sorted after every pepc_grow_tree,
     ! siblings are always contiguous!
     match_key = shift_key_by_level(particles(1)%key, level_march)
     do i = 1, size(particles)
       parent_key = shift_key_by_level(particles(i)%key, level_march)
       if (parent_key .ne. match_key) then
         unique_parents = unique_parents + 1
         match_key = parent_key
       end if
       sibling_cnt(unique_parents) = sibling_cnt(unique_parents) + 1
     end do

     average_siblings = 0.0_kind_physics
     do i = 1, unique_parents
       average_siblings = average_siblings + sibling_cnt(i)
     end do
     average_siblings = average_siblings/unique_parents

     print *, "No. of parents: ", unique_parents, " Average no. of siblings: ", average_siblings
   end subroutine

   subroutine defined_siblings_number_grouping(particles, sibling_upper_limit, sibling_cnt, unique_parents)
     use treevars, only: maxlevel
     use module_spacefilling
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     integer, intent(inout) :: sibling_upper_limit
     integer, allocatable, intent(inout) :: sibling_cnt(:)
     integer, intent(inout) :: unique_parents
     integer :: i, max_siblings
     integer(kind_level) :: level_march, level
     integer(kind_key) :: match_key, parent_key
     real(kind_physics) :: average_siblings

     allocate(sibling_cnt(size(particles)))
     level_march = -18_kind_level
     max_siblings = size(particles)

     if (sibling_upper_limit < max_siblings) then
       do while(sibling_upper_limit < max_siblings)
         sibling_cnt = 0
         level = maxlevel + level_march
         unique_parents = 1

         match_key = shift_key_by_level(particles(1)%key, level_march)
         do i = 1, size(particles)
           parent_key = shift_key_by_level(particles(i)%key, level_march)
           particles(i)%data%mp_int1 = parent_key
           if (parent_key .ne. match_key) then
             unique_parents = unique_parents + 1
             match_key = parent_key
           end if
           sibling_cnt(unique_parents) = sibling_cnt(unique_parents) + 1
         end do
         level_march = level_march + 1

         max_siblings = 0
         do i = 1, unique_parents
           if (sibling_cnt(i) > max_siblings) then
             max_siblings = sibling_cnt(i)
           end if
         end do
       end do
     else if (sibling_upper_limit > max_siblings) then
       unique_parents = 1
       level = -1
       do i = 0, size(particles)
         particles(i)%data%mp_int1 = level
       end do
       sibling_cnt(1) = max_siblings
     end if
     sibling_upper_limit = max_siblings
     print *, "Refinement level: ", level, "Max siblings: ", sibling_upper_limit, "Unique parents: ", unique_parents
   end subroutine

   subroutine gather_ll_buffers_omp(buffer, new_offset, new_particles_buffer, thread_id, thread_num)
     ! NOTE: 'buffer' refers to ll_elem that stores new particles
     !       'new_particle_buffer' that is shared across all threads must be created first!
     implicit none
     type(t_particle), allocatable, intent(inout) :: new_particles_buffer(:)
     type(linked_list_elem), pointer, intent(inout) :: buffer
     integer, allocatable, intent(in) :: new_offset(:)
     integer, intent(in) :: thread_num, thread_id
     type(linked_list_elem), pointer :: temp_guide
     integer :: head, tail, ll_buffer_size, i, remainder, ll_elem_cnt, start_i, stop_i, thread_new_p

     if (thread_id .ne. 0) then
       thread_new_p = new_offset(thread_id + 1) - new_offset(thread_id)
     else
       thread_new_p = new_offset(thread_id + 1)
     end if

     if (thread_new_p .ne. 0) then
       ll_buffer_size = size(buffer%tmp_particles)
       ll_elem_cnt = CEILING(real(thread_new_p)/real(ll_buffer_size))
       remainder = MOD(thread_new_p, ll_buffer_size)

       if (thread_id .ne. 0) then
         start_i = new_offset(thread_id) + 1
       else
         start_i = 1
       end if
       stop_i = new_offset(thread_id + 1)

       i = 1
       head = start_i
       tail = start_i + ll_buffer_size - 1

       temp_guide => buffer
       do while (associated(temp_guide))
          if (i /= ll_elem_cnt) then
             new_particles_buffer(head:tail) = temp_guide%tmp_particles
          else if ((i == ll_elem_cnt) .and. (remainder == 0)) then
             new_particles_buffer(head:stop_i) = temp_guide%tmp_particles
          else
             new_particles_buffer(head:stop_i) = temp_guide%tmp_particles(1:(remainder))
          end if

          head = tail + 1
          tail = tail + ll_buffer_size
          i = i + 1

          temp_guide => temp_guide%next
       end do
       nullify (temp_guide)
     end if
   end subroutine gather_ll_buffers_omp

   subroutine extend_particles_swap_omp(particles, new_particles_buffer, new_particles_size, swapped_cnt)
     ! NOTE: Remember to deallocate 'new_particle_buffer' outside this subroutine
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     type(t_particle), allocatable, intent(in) :: new_particles_buffer(:)
     integer, intent(in) :: new_particles_size, swapped_cnt
     integer :: new_size, head

     head = size(particles) - swapped_cnt + 1

     new_size = size(particles) + new_particles_size - swapped_cnt
     call resize_array(particles, new_size)

     if (new_particles_size .ne. 0) then
       particles(head:new_size) = new_particles_buffer
     end if
   end subroutine extend_particles_swap_omp

   subroutine extend_particles_swap_inject_omp(particles, new_particles_buffer, new_particles_size, swapped_cnt, electron_num)
     ! NOTE: Remember to deallocate 'new_particle_buffer' outside this subroutine
     implicit none
     type(t_particle), allocatable, intent(inout) :: particles(:)
     type(t_particle), allocatable, intent(in) :: new_particles_buffer(:)
     integer, intent(in) :: new_particles_size, swapped_cnt, electron_num
     integer :: new_size, head, temp_size
     real(kind_physics) :: center_pos(3), plane_orient(3)

     head = size(particles) - swapped_cnt + 1

     temp_size = size(particles) + new_particles_size - swapped_cnt
     new_size = temp_size + electron_num
     call resize_array(particles, new_size)

     if (new_particles_size .ne. 0) then
       particles(head:(temp_size)) = new_particles_buffer
     end if

     if (electron_num /= 0) then
       center_pos = 0.0
       center_pos(3) = -0.01
       plane_orient = 0.0
       plane_orient(3) = -1.0
       call injected_electrons(center_pos, plane_orient, 2, 0.025/(c*1e-12), &
                                 particles, temp_size + 1)
     end if
   end subroutine extend_particles_swap_inject_omp

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

   subroutine extend_particles_list_swap(particles, buffer, new_particles_size, swapped_cnt)
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:)
      type(linked_list_elem), pointer, intent(inout) :: buffer
      type(linked_list_elem), pointer :: temp_guide
      integer, intent(in) :: new_particles_size, swapped_cnt
      integer :: new_size, head, tail, ll_buffer_size, i, remainder, ll_elem_cnt

      ll_buffer_size = size(buffer%tmp_particles)
      ll_elem_cnt = CEILING(real(new_particles_size)/real(ll_buffer_size))
      remainder = MOD(new_particles_size, ll_buffer_size)

      new_size = size(particles) + new_particles_size - swapped_cnt
      i = 1
      head = size(particles) - swapped_cnt + 1
      tail = size(particles) - swapped_cnt + ll_buffer_size

      ! print *, "Linked List elements: ", ll_elem_cnt, size(particles), new_particles_size, remainder

      call resize_array(particles, new_size)

      if (new_particles_size /= 0) then
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
      end if
   end subroutine extend_particles_list_swap

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

     select case(geometry)
     case(1) ! square plane
       do index = starting_index, size(particles_list)
         call random_number(ran)
        !  call frand123NormDouble( state, mu, sigma2, res_mag )
         particles_list(index)%data%q = -1.0_8
         particles_list(index)%data%m = 1.0_8
         particles_list(index)%data%b = 0.0_8
         particles_list(index)%data%species = 0
         particles_list(index)%data%age = 0.0_8
         particles_list(index)%work = 1.0_8
         particles_list(index)%data%mp_int1 = 0

         particles_list(index)%x = center_pos
         particles_list(index)%data%v = 0.0
         particles_list(index)%data%f_b = 0.0_kind_physics
         particles_list(index)%data%f_e = 0.0_kind_physics
         particles_list(index)%label = 0

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
         particles_list(index)%data%b = 0.0_8
         particles_list(index)%data%species = 0
         particles_list(index)%data%age = 0.0_8
         particles_list(index)%work = 1.0_8
         particles_list(index)%data%mp_int1 = 0

         particles_list(index)%x = center_pos
         particles_list(index)%data%v = plane * magnitude
         particles_list(index)%data%f_b = 0.0_kind_physics
         particles_list(index)%data%f_e = 0.0_kind_physics
         particles_list(index)%label = 0

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

   subroutine merge_replace_particles_list(particles, buffer, merged_cnt)
      implicit none
      type(t_particle), allocatable, intent(inout) :: particles(:)
      type(linked_list_elem), pointer, intent(inout) :: buffer
      type(linked_list_elem), pointer :: temp_guide
      integer, intent(in) :: merged_cnt
      integer :: new_size, head, tail, ll_buffer_size, i, remainder, ll_elem_cnt

      ll_buffer_size = size(buffer%tmp_particles)
      ll_elem_cnt = CEILING(real(merged_cnt)/real(ll_buffer_size))
      remainder = MOD(merged_cnt, ll_buffer_size)

      new_size = merged_cnt
      i = 1
      head = 1

      if (new_size < ll_buffer_size) then
        tail = new_size
      else
        tail = ll_buffer_size
      end if

      deallocate(particles)
      allocate(particles(new_size))
      ! print *, "Linked List elements: ", ll_elem_cnt, size(particles), new_particles_size, remainder

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
   end subroutine merge_replace_particles_list

   subroutine register_density_diag_type()
     use mpi
     implicit none

     integer, parameter :: fields = 4
     integer :: mpi_err
     integer, dimension(1:fields) :: blocklengths, types
     integer(KIND=MPI_ADDRESS_KIND), dimension(1:fields) :: disp
     integer(KIND=MPI_ADDRESS_KIND), dimension(0:fields) :: address

     type(diag_vertex) :: dummy_vertex

     blocklengths(1:3)  = [3, 2, 3]
     types(1:3) = [MPI_REAL8, MPI_REAL8, MPI_REAL8]
     call MPI_GET_ADDRESS( dummy_vertex,           address(0), mpi_err)
     call MPI_GET_ADDRESS( dummy_vertex%x,         address(1), mpi_err)
     call MPI_GET_ADDRESS( dummy_vertex%q_density, address(2), mpi_err)
     call MPI_GET_ADDRESS( dummy_vertex%J_density, address(3), mpi_err)
     disp(1:3) = address(1:3) - address(0)
     call MPI_TYPE_CREATE_STRUCT( 3, blocklengths, disp, types, MPI_TYPE_density, mpi_err)
     call MPI_TYPE_COMMIT( MPI_TYPE_density, mpi_err)
   end subroutine register_density_diag_type

   subroutine free_density_diag_type()
     use mpi
     implicit none

     integer :: mpi_err
     call MPI_TYPE_FREE(MPI_TYPE_density, mpi_err)
   end subroutine free_density_diag_type
end module
