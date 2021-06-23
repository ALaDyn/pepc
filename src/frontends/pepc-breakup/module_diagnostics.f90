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
!> diagnostics module
!>

module diagnostics
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use helper
   implicit none

contains
  subroutine torus_diagnostic_xz_grid(major_radius, minor_radius, subdivisions, points)
    ! Diagnostic_grid populates domain with points at regular distance, in xz plane
    implicit none
    real(kind_physics), intent(in) :: major_radius, minor_radius
    integer, intent(in) :: subdivisions ! denotes number of subdivisions in 1 dimension
    type(t_particle), allocatable, intent(inout) :: points(:)
    real(kind_physics) :: diameter, increments, temp_dist(2), &
                          dist_sum_x, dist_sum_z, length_ratio, x_length
    integer :: total_points_x, total_points_z, total_points

    diameter = 2.0*minor_radius
    x_length = 2.0_8*(major_radius + minor_radius)
    length_ratio = x_length/diameter
    total_points_z = (2**subdivisions + 1)
    total_points_x = int(ceiling(total_points_z*length_ratio))
    total_points = total_points_x*total_points_z
    increments = diameter/(total_points_z - 1)

    allocate(points(total_points))

    temp_dist(1) = -major_radius - minor_radius - 0.25*increments
    temp_dist(2) = -minor_radius - 0.5*increments

    dist_sum_x = 0.0_kind_physics
    dist_sum_z = 0.0_kind_physics
    do i = 1, size(points)
      points(i)%x = 0.0
      points(i)%x(1) = temp_dist(1) + dist_sum_x
      points(i)%x(3) = temp_dist(2) + dist_sum_z

      dist_sum_x = dist_sum_x + increments

      if (dist_sum_x > x_length) then
        dist_sum_x = 0.0
        dist_sum_z = dist_sum_z + increments
      end if

      points(i)%label = 0
      points(i)%data%q = 0.0_8
      points(i)%data%m = 1.0_8
      points(i)%data%species = 0
      points(i)%data%age = 0.0_8
      points(i)%work = 1.0_8
      points(i)%data%v = 0.0
    end do
  end subroutine torus_diagnostic_xz_grid

  subroutine torus_diagnostic_xz_breakdown(major_radius, minor_radius, subdivisions, points)
    ! Diagnostic_grid populates domain with points at regular distance, in xz plane
    implicit none
    real(kind_physics), intent(in) :: major_radius, minor_radius
    integer, intent(in) :: subdivisions ! denotes number of subdivisions in 1 dimension
    type(t_particle), allocatable, intent(inout) :: points(:)
    real(kind_physics) :: increments, temp_dist(2), dist_sum_x, &
                          dist_sum_z, x_length
    integer :: total_points_1D, total_points

    x_length = 2.0_8*minor_radius
    total_points_1D = (2**subdivisions + 1)
    total_points = total_points_1D**2
    increments = x_length/(total_points_1D - 1)

    allocate(points(total_points))

    temp_dist(1) = major_radius - minor_radius
    temp_dist(2) = -minor_radius

    dist_sum_x = 0.0_kind_physics
    dist_sum_z = 0.0_kind_physics
    do i = 1, size(points)
      points(i)%x = 0.0
      points(i)%x(1) = temp_dist(1) + dist_sum_x
      points(i)%x(3) = temp_dist(2) + dist_sum_z

      dist_sum_x = dist_sum_x + increments

      if (MOD(i, total_points_1D) == 0) then
        dist_sum_x = 0.0
        dist_sum_z = dist_sum_z + increments
      end if

      points(i)%label = 0
      points(i)%data%q = 0.0_8
      points(i)%data%m = 1.0_8
      points(i)%data%species = 0
      points(i)%data%age = 0.0_8
      points(i)%work = 1.0_8
      points(i)%data%v = 0.0
    end do
  end subroutine torus_diagnostic_xz_breakdown

  subroutine torus_diagnostic_xy_grid(major_radius, minor_radius, subdivisions, points, z_offset)
    ! Diagnostic_grid populates domain with points at regular distance, in xy plane
    implicit none
    real(kind_physics), intent(in) :: major_radius, minor_radius
    integer, intent(in) :: subdivisions ! denotes number of subdivisions in 1 dimension
    type(t_particle), allocatable, intent(inout) :: points(:)
    real(kind_physics) :: increments, temp_dist(2), dist_sum_x, &
                          dist_sum_y, x_length, z_offset, z
    integer :: total_points_1D, total_points

    z = z_offset/(c*1e-12)
    x_length = 2.0_8*(major_radius + minor_radius)
    total_points_1D = (2**subdivisions + 1)
    total_points = total_points_1D**2
    increments = x_length/(total_points_1D - 1)

    allocate(points(total_points))

    temp_dist(1) = -major_radius - minor_radius
    temp_dist(2) = -major_radius - minor_radius

    dist_sum_x = 0.0_kind_physics
    dist_sum_y = 0.0_kind_physics
    do i = 1, size(points)
      points(i)%x = 0.0
      points(i)%x(1) = temp_dist(1) + dist_sum_x
      points(i)%x(2) = temp_dist(2) + dist_sum_y
      points(i)%x(3) = z

      dist_sum_x = dist_sum_x + increments

      if (MOD(i, total_points_1D) == 0) then
        dist_sum_x = 0.0
        dist_sum_y = dist_sum_y + increments
      end if

      points(i)%label = 0
      points(i)%data%q = 0.0_8
      points(i)%data%m = 1.0_8
      points(i)%data%species = 0
      points(i)%data%age = 0.0_8
      points(i)%work = 1.0_8
      points(i)%data%v = 0.0
    end do
  end subroutine torus_diagnostic_xy_grid

  subroutine cube_cell_vertices(cell_n, vertices)!loll, lorl, upll, uprl, lolu, loru, uplu, upru)
    implicit none
    integer(kind_particle), intent(in) :: cell_n
    integer(kind_particle), intent(out) :: vertices(8) !loll, lorl, upll, uprl, lolu, loru, uplu, upru
    integer(kind_particle) :: ny, nz, offset
    real(kind_physics) :: temp_val

    ! NOTE: cell number is taken as starting from '0', to compensate the nature of
    !       division coupled with FLOOR function. Cell numbers technically starts at '1'
    temp_val = (cell_n - 1)/(x_cell*y_cell*1.0)
    nz = FLOOR(temp_val)

    temp_val = MOD((cell_n - 1),(x_cell*y_cell))/(1.0*x_cell)
    ny = FLOOR(temp_val)
    offset = (x_cell + 1)*(y_cell + 1)

    vertices(1) = cell_n + ny + (y_cell + 1 + x_cell)*nz - 1          ! loll
    vertices(2) = vertices(1) + 1                                    ! lorl
    vertices(3) = vertices(1) + x_cell + 1                           ! upll
    vertices(4) = vertices(3) + 1                                    ! uprl
    vertices(5) = vertices(1) + offset                               ! lolu
    vertices(6) = vertices(2) + offset                               ! loru
    vertices(7) = vertices(3) + offset                               ! uplu
    vertices(8) = vertices(4) + offset                               ! upru
  end subroutine cube_cell_vertices

  subroutine construct_connectivity_cube(array)
    implicit none
    integer(kind_particle), allocatable, intent(inout) :: array(:)
    integer(kind_particle) :: ic, start_index, end_index

   !  print *, size(array)/8

    do ic = 1, size(array)/8
      start_index = (ic - 1)*8 + 1
      end_index = start_index + 7
      call cube_cell_vertices(ic, array(start_index:end_index))
     !  print *, ic, array(start_index:end_index)
    end do
  end subroutine construct_connectivity_cube

  subroutine construct_connectivity_tets(array, tets_connectivity_matrix)
    implicit none
    integer(kind_particle), allocatable, intent(inout) :: array(:)
    integer(kind_particle), allocatable, intent(in) :: tets_connectivity_matrix(:,:)
    integer(kind_particle) :: i_elem, i_cnt

    ! NOTE: VTK connectivity starts from 0 index
    do i_elem = 1, size(tets_connectivity_matrix(:,1))
      i_cnt = (i_elem - 1)*4
      array(i_cnt + 1) = tets_connectivity_matrix(i_elem,1) - 1
      array(i_cnt + 2) = tets_connectivity_matrix(i_elem,2) - 1
      array(i_cnt + 3) = tets_connectivity_matrix(i_elem,3) - 1
      array(i_cnt + 4) = tets_connectivity_matrix(i_elem,4) - 1
    end do
  end subroutine construct_connectivity_tets

  subroutine init_diagnostic_verts(vertices, min_x, min_y, min_z, &
                                   x_width, y_width, z_width)
     implicit none
     type(diag_vertex), allocatable, intent(inout) :: vertices(:)
     real(kind_physics), intent(in) :: min_x, min_y, min_z, x_width, y_width, z_width
     integer :: total_vertices, ix, iy, iz, ip
     real(kind_physics) :: s_x_width, s_y_width, s_z_width

     s_min_x = min_x/(c*1e-12)
     s_min_y = min_y/(c*1e-12)
     s_min_z = min_z/(c*1e-12)
     s_x_width = x_width/(c*1e-12)
     s_y_width = y_width/(c*1e-12)
     s_z_width = z_width/(c*1e-12)

     total_vertices = (x_cell+1) * (y_cell+1) * (z_cell+1)

     if (my_rank == 0) then
       allocate(connectivity_array(x_cell*y_cell*z_cell*8))
       call construct_connectivity_cube(connectivity_array)
       ! print *, connectivity_array
     end if

     allocate (vertices(total_vertices))
     dx = s_x_width/(x_cell - 1)
     dy = s_y_width/(y_cell - 1)
     dz = s_z_width/(z_cell - 1)

     ix = 1
     iy = 1
     iz = 1
     do ip = 1, total_vertices
       vertices(ip)%x(1) = s_min_x + (ix - 1.5)*dx
       vertices(ip)%x(2) = s_min_y + (iy - 1.5)*dy
       vertices(ip)%x(3) = s_min_z + (iz - 1.5)*dz
       vertices(ip)%q_density = 0.0
       vertices(ip)%J_density = 0.0

       ix = ix + 1
       if (ix == x_cell + 2) then
         ix = 1
         iy = iy + 1
         if (iy == y_cell + 2) then
           ix = 1
           iy = 1
           iz = iz + 1
         end if
       end if
     end do

  end subroutine init_diagnostic_verts

  subroutine density_interpolation(particle, vertices)
    implicit none
    type(t_particle), intent(in) :: particle
    type(diag_vertex), allocatable, intent(inout) :: vertices(:)
    real(kind_physics) :: min_x, min_y, min_z, inv_vol, xp, yp, zp, coeff
    real(kind_physics) :: temp_q_dens
    integer :: total_vertices, nx, ny, nz, N_cell
    integer :: loll, lorl, upll, uprl, lolu, loru, uplu, upru, offset
    ! NOTE: loll = lower left corner, lower plane
    !       lorl = lower right corner, lower plane
    !       uplu = upper left corner, upper plane
    !       upru = upper right corner, upper plane
    !       The rest can be inferred.

    total_vertices = size(vertices)

    !  NOTE: dx, dy & dz are global variables, 'min' denotes the minimum position
    !        where particles are allowed to populate.
    min_x = vertices(1)%x(1) + 0.5*dx
    min_y = vertices(1)%x(2) + 0.5*dy
    min_z = vertices(1)%x(3) + 0.5*dz
    inv_vol = 1./(dx*dy*dz)
    ! print *, dx*dy*dz, dx, dy ,dz

    offset = (x_cell + 1)*(y_cell + 1)

    xp = particle%x(1)
    yp = particle%x(2)
    zp = particle%x(3)
    coeff = particle%data%q * inv_vol**2
    ! - min_x to centre the coordinates to 0.0, easier to compute cell num.
    nx = FLOOR((xp - min_x)/dx + 0.5)
    ny = FLOOR((yp - min_y)/dy + 0.5)
    nz = FLOOR((zp - min_z)/dz + 0.5)
    N_cell = (x_cell*y_cell)*nz + (x_cell*ny + nx) + 1

    ! with cell number known, calculate the 8 vertices that bounds the cell.
    loll = N_cell + ny + (y_cell + 1 + x_cell)*nz
    lorl = loll + 1
    upll = loll + x_cell + 1
    uprl = upll + 1
    lolu = loll + offset
    loru = lorl + offset
    uplu = upll + offset
    upru = uprl + offset

    ! interpolate the particle charge to the vertices. Also computes J.
    if (particle%data%q < 0.0) then
      temp_q_dens = (vertices(lorl)%x(1) - xp) &
                    *(vertices(upll)%x(2) - yp) &
                    *(vertices(lolu)%x(3) - zp)*coeff
      vertices(loll)%q_density(1) = temp_q_dens + vertices(loll)%q_density(1)
      vertices(loll)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(loll)%J_density

      temp_q_dens = (xp - vertices(loll)%x(1)) &
                    *(vertices(uprl)%x(2) - yp) &
                    *(vertices(loru)%x(3) - zp)*coeff
      vertices(lorl)%q_density(1) = temp_q_dens + vertices(lorl)%q_density(1)
      vertices(lorl)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(lorl)%J_density

      temp_q_dens = (vertices(uprl)%x(1) - xp) &
                    *(yp - vertices(loll)%x(2)) &
                    *(vertices(uplu)%x(3) - zp)*coeff
      vertices(upll)%q_density(1) = temp_q_dens + vertices(upll)%q_density(1)
      vertices(upll)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(upll)%J_density

      temp_q_dens = (xp - vertices(upll)%x(1)) &
                    *(yp - vertices(lorl)%x(2)) &
                    *(vertices(upru)%x(3) - zp)*coeff
      vertices(uprl)%q_density(1) = temp_q_dens + vertices(uprl)%q_density(1)
      vertices(uprl)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(uprl)%J_density

      temp_q_dens = (vertices(loru)%x(1) - xp) &
                    *(vertices(uplu)%x(2) - yp) &
                    *(zp - vertices(loll)%x(3))*coeff
      vertices(lolu)%q_density(1) = temp_q_dens + vertices(lolu)%q_density(1)
      vertices(lolu)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(lolu)%J_density

      temp_q_dens = (xp - vertices(lolu)%x(1)) &
                    *(vertices(upru)%x(2) - yp) &
                    *(zp - vertices(lorl)%x(3))*coeff
      vertices(loru)%q_density(1) = temp_q_dens + vertices(loru)%q_density(1)
      vertices(loru)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(loru)%J_density

      temp_q_dens = (vertices(upru)%x(1) - xp) &
                    *(yp - vertices(lolu)%x(2)) &
                    *(zp - vertices(upll)%x(3))*coeff
      vertices(uplu)%q_density(1) = temp_q_dens + vertices(uplu)%q_density(1)
      vertices(uplu)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(uplu)%J_density

      temp_q_dens = (xp - vertices(uplu)%x(1)) &
                    *(yp - vertices(loru)%x(2)) &
                    *(zp - vertices(uprl)%x(3))*coeff
      vertices(upru)%q_density(1) = temp_q_dens + vertices(upru)%q_density(1)
      vertices(upru)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(upru)%J_density

     !  test = vertices(loll)%q_density(1) + vertices(lorl)%q_density(1) &
     !         + vertices(upll)%q_density(1) + vertices(uprl)%q_density(1) &
     !         + vertices(lolu)%q_density(1) + vertices(loru)%q_density(1) &
     !         + vertices(uplu)%q_density(1) + vertices(upru)%q_density(1)
     !  print *, test
   else if (particle%data%q > 0.0) then
     temp_q_dens = (vertices(lorl)%x(1) - xp) &
                   *(vertices(upll)%x(2) - yp) &
                   *(vertices(lolu)%x(3) - zp)*coeff
     vertices(loll)%q_density(2) = temp_q_dens + vertices(loll)%q_density(2)
     vertices(loll)%J_density = temp_q_dens*particle%data%v &
                                + vertices(loll)%J_density

     temp_q_dens = (xp - vertices(loll)%x(1)) &
                   *(vertices(uprl)%x(2) - yp) &
                   *(vertices(loru)%x(3) - zp)*coeff
     vertices(lorl)%q_density(2) = temp_q_dens + vertices(lorl)%q_density(2)
     vertices(lorl)%J_density = temp_q_dens*particle%data%v &
                                + vertices(lorl)%J_density

     temp_q_dens = (vertices(uprl)%x(1) - xp) &
                   *(yp - vertices(loll)%x(2)) &
                   *(vertices(uplu)%x(3) - zp)*coeff
     vertices(upll)%q_density(2) = temp_q_dens + vertices(upll)%q_density(2)
     vertices(upll)%J_density = temp_q_dens*particle%data%v &
                                + vertices(upll)%J_density

     temp_q_dens = (xp - vertices(upll)%x(1)) &
                   *(yp - vertices(lorl)%x(2)) &
                   *(vertices(upru)%x(3) - zp)*coeff
     vertices(uprl)%q_density(2) = temp_q_dens + vertices(uprl)%q_density(2)
     vertices(uprl)%J_density = temp_q_dens*particle%data%v &
                                + vertices(uprl)%J_density

     temp_q_dens = (vertices(loru)%x(1) - xp) &
                   *(vertices(uplu)%x(2) - yp) &
                   *(zp - vertices(loll)%x(3))*coeff
     vertices(lolu)%q_density(2) = temp_q_dens + vertices(lolu)%q_density(2)
     vertices(lolu)%J_density = temp_q_dens*particle%data%v &
                                + vertices(lolu)%J_density

     temp_q_dens = (xp - vertices(lolu)%x(1)) &
                   *(vertices(upru)%x(2) - yp) &
                   *(zp - vertices(lorl)%x(3))*coeff
     vertices(loru)%q_density(2) = temp_q_dens + vertices(loru)%q_density(2)
     vertices(loru)%J_density = temp_q_dens*particle%data%v &
                                + vertices(loru)%J_density

     temp_q_dens = (vertices(upru)%x(1) - xp) &
                   *(yp - vertices(lolu)%x(2)) &
                   *(zp - vertices(upll)%x(3))*coeff
     vertices(uplu)%q_density(2) = temp_q_dens + vertices(uplu)%q_density(2)
     vertices(uplu)%J_density = temp_q_dens*particle%data%v &
                                + vertices(uplu)%J_density

     temp_q_dens = (xp - vertices(uplu)%x(1)) &
                   *(yp - vertices(loru)%x(2)) &
                   *(zp - vertices(uprl)%x(3))*coeff
     vertices(upru)%q_density(2) = temp_q_dens + vertices(upru)%q_density(2)
     vertices(upru)%J_density = temp_q_dens*particle%data%v &
                                + vertices(upru)%J_density
    end if
  end subroutine density_interpolation

  subroutine clear_density_results(vertices)
    implicit none
    type(diag_vertex), allocatable, intent(inout) :: vertices(:)
    integer :: total_vertices, iv

    total_vertices = size(vertices)

    do iv = 1, total_vertices
      vertices(iv)%q_density = 0.0
      vertices(iv)%J_density = 0.0
    end do
  end subroutine clear_density_results

  subroutine add_neutral_points(bounding_box, nx, ny, nz, particles)
    use module_box
    implicit none
    integer, intent(in) :: nx, ny, nz
    type(t_box), intent(in) :: bounding_box
    type(t_particle), allocatable, intent(inout) :: particles(:)
    type(t_particle), allocatable :: neutral_points(:), temp_array(:)
    real(kind_physics) :: extents(3), box_min(3), delta(3)
    integer :: total_n, n, i, old_size, new_size, plane_n

    total_n = nx*ny*nz
    old_size = size(particles)
    new_size = old_size + total_n

    allocate(temp_array(old_size))
    temp_array = particles

    deallocate(particles)
    allocate(particles(new_size))
    particles(1:old_size) = temp_array
    deallocate(temp_array)

    extents = bounding_box%boxsize
    box_min = bounding_box%boxmin
    delta(1) = extents(1)/nx
    delta(2) = extents(2)/ny
    delta(3) = extents(3)/nz
    plane_n = nx*ny

    allocate(neutral_points(total_n))

    do n = 0, total_n - 1
      i = n + 1
      neutral_points(i)%x(1) = MOD(n, nx)*delta(1) + box_min(1)
      neutral_points(i)%x(2) = FLOOR(real(MOD(n, plane_n)/nx))*delta(2) + box_min(2)
      neutral_points(i)%x(3) = FLOOR(real(n/plane_n))*delta(3) + box_min(3)
      ! coordinates calculations

      ! fill out the rest of the important variables
      neutral_points(i)%data%v = 0.0_kind_physics
      neutral_points(i)%data%m = 0.0_kind_physics
      neutral_points(i)%data%q = 0.0_kind_physics
      neutral_points(i)%data%b = 0.0_kind_physics
      neutral_points(i)%data%f_e = 0.0_kind_physics
      neutral_points(i)%data%f_b = 0.0_kind_physics
      neutral_points(i)%data%species = 50.0_kind_physics
      neutral_points(i)%label = i
    end do

    particles((old_size + 1):new_size) = neutral_points
    deallocate(neutral_points)
  end subroutine

  subroutine read_msh_file(fname, vertices, connectivity, rank)
    implicit none
    type(diag_vertex), allocatable, intent(inout) :: vertices(:)
    integer(kind_particle), allocatable, intent(inout) :: connectivity(:,:)
    character(255), intent(in) :: fname
    integer, intent(in) :: rank
    character(255) :: read_string
    integer :: entries, i, j, stat, state
    integer :: int1, int2, int3, int4, int5, node_number, entity_blocks, &
               node_lines_skip, element_number, elem_dim, elem_type, &
               elem_lines_skip, resume_index
    real(kind_physics) :: flt1, flt2, flt3
    ! NOTE: state description:
    !       state -1: no action Done
    !       state 0: used in line skipping.
    !       state 1: read general statistics of $Nodes
    !       state 2: read general statistics of $Elements
    !       state 3: read statistics of $Node entity_blocks
    !       state 4: read statistics of $Elements entity_blocks
    !       state 5: begin reading and storing node coordinates
    !       state 6: begin reading and storing element connectivities

    open(11,file=fname,action='READ')
    state = -1
    node_lines_skip = 0
    elem_lines_skip = 0
    i = -1
    j = -1
    do
      read(11, '(a)', iostat=stat) read_string
      if (stat .ne. 0) EXIT
      ! print *, "State: ", state

      if (read_string(1:6) == '$Nodes') then
        state = 1
        resume_index = 0
        CYCLE
      elseif(read_string(1:9) == '$EndNodes') then
        print *, "Done reading nodes."
        node_lines_skip = 0
        state = -1
        CYCLE
      elseif (read_string(1:9) == '$Elements') then
        state = 2
        CYCLE
      elseif (read_string(1:12) == '$EndElements') then
        print *, "Done reading elements."
        state = -1
        CYCLE
      end if

      if (state == 1) then
        ! print *, read_string
        read(read_string, *) int1, int2, int3, int4
        entity_blocks = int1
        node_number = int2
        allocate(vertices(node_number))
        print *, "Rank ", rank, "Vertices array allocated!", node_number
        do int5 = 1, node_number
          vertices(int5)%x = 0.0
          vertices(int5)%q_density = 0.0
          vertices(int5)%J_density = 0.0
        end do
        state = 3
        CYCLE
      elseif (state == 2) then
        read(read_string, *) int1, int2, int3, int4
        entity_blocks = int1
        element_number = int2
        state = 4
        CYCLE
      end if

      if (state == 3) then
        ! print *, read_string
        read(read_string, *) int1, int2, int3, int4
        node_lines_skip = int4
        i = 0
        state = 0
        ! print *, node_lines_skip
        CYCLE
      elseif (state == 4) then
        ! print *, read_string
        read(read_string, *) int1, int2, int3, int4
        elem_dim = int1
        elem_type = int3
        elem_lines_skip = int4
        i = 0
        state = 0
        if ((elem_dim == 3) .and. (elem_type == 4)) then
          allocate(connectivity(elem_lines_skip, 4))
          print *, "Rank ", rank, "Connectivity array allocated!", elem_lines_skip
          state = 6
        end if
        CYCLE
      end if
      ! NOTE: Skipping lines
      if (state == 0) then
        i = i + 1
        if (i == node_lines_skip) then
          i = 0
          j = 0
          state = 5
          CYCLE
        elseif (i == elem_lines_skip) then
          i = 0
          j = 0
          state = 4
          CYCLE
        end if
      end if
      ! NOTE: Storing data to respective arrays
      if (state == 5) then
        j = j + 1
        read(read_string, *) flt1, flt2, flt3
        vertices(j+resume_index)%x(1) = flt1
        vertices(j+resume_index)%x(2) = flt2
        vertices(j+resume_index)%x(3) = flt3

        if (j == node_lines_skip) then
          resume_index = resume_index + j
          state = 3
        end if
      elseif(state == 6) then
        j = j + 1
        read(read_string, *) int1, int2, int3, int4, int5
        connectivity(j,1) = int2
        connectivity(j,2) = int3
        connectivity(j,3) = int4
        connectivity(j,4) = int5
      end if
    end do

    close(11)
  end subroutine read_msh_file

  function tet_volume_V1(point_a, point_b, point_c, point_d) result(vol)
    implicit none
    real(kind_physics), intent(in) :: point_a(3), point_b(3), point_c(3), point_d(3)
    real(kind_physics) :: vol, vec_a(3), vec_b(3), vec_c(3), temp_vec(3)

    temp_vec = 0.0
    vec_a = point_a - point_d
    vec_b = point_b - point_d
    vec_c = point_c - point_d

    temp_vec = cross_product(vec_b, vec_c)
    vol = abs(dot_product(vec_a,temp_vec))/6.0
  end function tet_volume_V1

  subroutine tet_mesh_interpolation(particle, vertices, connectivity)
    implicit none
    type(diag_vertex), allocatable, intent(inout) :: vertices(:)
    integer(kind_particle), allocatable, intent(inout) :: connectivity(:,:)
    type(t_particle), intent(in) :: particle
    real(kind_physics) :: bc1, bc2, bc3, bc4, ref_vol, temp_vol, volume_scale, eps
    real(kind_physics) :: p1(3), p2(3), p3(3), p4(3), t_p(3), check_sum, min_sum
    real(kind_physics) :: temp_q_dens
    integer :: i, j, bug_check, el1, el2, el3, el4

    bug_check = 0
    min_sum = 1E8
    eps = 1E-15
    do i = 1, size(connectivity(:,1))
      el1 = connectivity(i,1)
      el2 = connectivity(i,2)
      el3 = connectivity(i,3)
      el4 = connectivity(i,4)

      p1 = vertices(el1)%x
      p2 = vertices(el2)%x
      p3 = vertices(el3)%x
      p4 = vertices(el4)%x
      t_p = particle%x

      ref_vol = tet_volume_V1(p1, p2, p3, p4)
      volume_scale = 1.0/ref_vol

      temp_vol = tet_volume_V1(t_p, p2, p3, p4)
      bc1 = temp_vol*volume_scale

      temp_vol = tet_volume_V1(p1, t_p, p3, p4)
      bc2 = temp_vol*volume_scale

      temp_vol = tet_volume_V1(p1, p2, t_p, p4)
      bc3 = temp_vol*volume_scale

      temp_vol = tet_volume_V1(p1, p2, p3, t_p)
      bc4 = temp_vol*volume_scale

      check_sum = bc1 + bc2 + bc3 + bc4
      if (check_sum < min_sum) then
        min_sum = check_sum
      end if

      if (check_sum > 1.0 + eps) then
        CYCLE
      else
        ! print *, ref_vol!, bc1, bc2, bc3, bc4,
        j = 1
        if (particle%data%q > 0.0) then
          j = 2
        end if
        temp_q_dens = bc1*particle%data%q*volume_scale
        vertices(el1)%q_density(j) = temp_q_dens + vertices(el1)%q_density(j)
        vertices(el1)%J_density = temp_q_dens*particle%data%v + vertices(el1)%J_density

        temp_q_dens = bc2*particle%data%q*volume_scale
        vertices(el2)%q_density(j) = temp_q_dens + vertices(el2)%q_density(j)
        vertices(el2)%J_density = temp_q_dens*particle%data%v + vertices(el2)%J_density

        temp_q_dens = bc3*particle%data%q*volume_scale
        vertices(el3)%q_density(j) = temp_q_dens + vertices(el3)%q_density(j)
        vertices(el3)%J_density = temp_q_dens*particle%data%v + vertices(el3)%J_density

        temp_q_dens = bc4*particle%data%q*volume_scale
        vertices(el4)%q_density(j) = temp_q_dens + vertices(el4)%q_density(j)
        vertices(el4)%J_density = temp_q_dens*particle%data%v + vertices(el4)%J_density

        !NOTE: avoid double_counting particle on edges
        bug_check = 1
        EXIT
      end if
    end do

    if (bug_check == 0) then
      print *, "Particle is not interpolated!", t_p, min_sum
    end if
  end subroutine tet_mesh_interpolation

end module
