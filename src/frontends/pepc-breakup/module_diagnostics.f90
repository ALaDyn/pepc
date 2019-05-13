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

  subroutine cell_vertices(cell_n, vertices)!loll, lorl, upll, uprl, lolu, loru, uplu, upru)
    implicit none
    integer, intent(in) :: cell_n
    integer, intent(out) :: vertices(8) !loll, lorl, upll, uprl, lolu, loru, uplu, upru
    integer :: ny, nz, offset
    real*8 :: temp_val

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
  end subroutine cell_vertices

  subroutine construct_connectivity(array)
    implicit none
    integer, allocatable, intent(inout) :: array(:)
    integer :: ic, start_index, end_index

   !  print *, size(array)/8

    do ic = 1, size(array)/8
      start_index = (ic - 1)*8 + 1
      end_index = start_index + 7
      call cell_vertices(ic, array(start_index:end_index))
     !  print *, ic, array(start_index:end_index)
    end do
  end subroutine construct_connectivity

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
       call construct_connectivity(connectivity_array)
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

    ! interpolate the particle's charge to the vertices. Also computes J.
    if (particle%data%q < 0.0) then
      vertices(loll)%q_density(1) = (vertices(lorl)%x(1) - xp) &
                                    *(vertices(upll)%x(2) - yp) &
                                    *(vertices(lolu)%x(3) - zp) &
                                    *coeff + vertices(loll)%q_density(1)
      vertices(loll)%J_density = vertices(loll)%q_density(1) * particle%data%v &
                                 + vertices(loll)%J_density

      vertices(lorl)%q_density(1) = (xp - vertices(loll)%x(1)) &
                                    *(vertices(uprl)%x(2) - yp) &
                                    *(vertices(loru)%x(3) - zp) &
                                    *coeff + vertices(lorl)%q_density(1)
      vertices(lorl)%J_density = vertices(lorl)%q_density(1) * particle%data%v &
                                 + vertices(lorl)%J_density

      vertices(upll)%q_density(1) = (vertices(uprl)%x(1) - xp) &
                                    *(yp - vertices(loll)%x(2)) &
                                    *(vertices(uplu)%x(3) - zp) &
                                    *coeff + vertices(upll)%q_density(1)
      vertices(upll)%J_density = vertices(upll)%q_density(1) * particle%data%v &
                                 + vertices(upll)%J_density

      vertices(uprl)%q_density(1) = (xp - vertices(upll)%x(1)) &
                                    *(yp - vertices(lorl)%x(2)) &
                                    *(vertices(upru)%x(3) - zp) &
                                    *coeff + vertices(uprl)%q_density(1)
      vertices(uprl)%J_density = vertices(uprl)%q_density(1) * particle%data%v &
                                 + vertices(uprl)%J_density

      vertices(lolu)%q_density(1) = (vertices(loru)%x(1) - xp) &
                                    *(vertices(uplu)%x(2) - yp) &
                                    *(zp - vertices(loll)%x(3)) &
                                    *coeff + vertices(lolu)%q_density(1)
      vertices(lolu)%J_density = vertices(lolu)%q_density(1) * particle%data%v &
                                 + vertices(lolu)%J_density

      vertices(loru)%q_density(1) = (xp - vertices(lolu)%x(1)) &
                                    *(vertices(upru)%x(2) - yp) &
                                    *(zp - vertices(lorl)%x(3)) &
                                    *coeff + vertices(loru)%q_density(1)
      vertices(loru)%J_density = vertices(loru)%q_density(1) * particle%data%v &
                                 + vertices(loru)%J_density

      vertices(uplu)%q_density(1) = (vertices(upru)%x(1) - xp) &
                                    *(yp - vertices(lolu)%x(2)) &
                                    *(zp - vertices(upll)%x(3)) &
                                    *coeff + vertices(uplu)%q_density(1)
      vertices(uplu)%J_density = vertices(uplu)%q_density(1) * particle%data%v &
                                 + vertices(uplu)%J_density

      vertices(upru)%q_density(1) = (xp - vertices(uplu)%x(1)) &
                                    *(yp - vertices(loru)%x(2)) &
                                    *(zp - vertices(uprl)%x(3)) &
                                    *coeff + vertices(upru)%q_density(1)
      vertices(upru)%J_density = vertices(upru)%q_density(1) * particle%data%v &
                                 + vertices(upru)%J_density

     !  test = vertices(loll)%q_density(1) + vertices(lorl)%q_density(1) &
     !         + vertices(upll)%q_density(1) + vertices(uprl)%q_density(1) &
     !         + vertices(lolu)%q_density(1) + vertices(loru)%q_density(1) &
     !         + vertices(uplu)%q_density(1) + vertices(upru)%q_density(1)
     !  print *, test
   else if (particle%data%q > 0.0) then
     vertices(loll)%q_density(2) = (vertices(lorl)%x(1) - xp) &
                                   *(vertices(upll)%x(2) - yp) &
                                   *(vertices(lolu)%x(3) - zp) &
                                   *coeff + vertices(loll)%q_density(2)
     vertices(lorl)%q_density(2) = (xp - vertices(loll)%x(1)) &
                                   *(vertices(uprl)%x(2) - yp) &
                                   *(vertices(loru)%x(3) - zp) &
                                   *coeff + vertices(lorl)%q_density(2)
     vertices(upll)%q_density(2) = (vertices(uprl)%x(1) - xp) &
                                   *(yp - vertices(loll)%x(2)) &
                                   *(vertices(uplu)%x(3) - zp) &
                                   *coeff + vertices(upll)%q_density(2)
     vertices(uprl)%q_density(2) = (xp - vertices(upll)%x(1)) &
                                   *(yp - vertices(lorl)%x(2)) &
                                   *(vertices(upru)%x(3) - zp) &
                                   *coeff + vertices(uprl)%q_density(2)
     vertices(lolu)%q_density(2) = (vertices(loru)%x(1) - xp) &
                                   *(vertices(uplu)%x(2) - yp) &
                                   *(zp - vertices(loll)%x(3)) &
                                   *coeff + vertices(lolu)%q_density(2)
     vertices(loru)%q_density(2) = (xp - vertices(lolu)%x(1)) &
                                   *(vertices(upru)%x(2) - yp) &
                                   *(zp - vertices(lorl)%x(3)) &
                                   *coeff + vertices(loru)%q_density(2)
     vertices(uplu)%q_density(2) = (vertices(upru)%x(1) - xp) &
                                   *(yp - vertices(lolu)%x(2)) &
                                   *(zp - vertices(upll)%x(3)) &
                                   *coeff + vertices(uplu)%q_density(2)
     vertices(upru)%q_density(2) = (xp - vertices(uplu)%x(1)) &
                                   *(yp - vertices(loru)%x(2)) &
                                   *(zp - vertices(uprl)%x(3)) &
                                   *coeff + vertices(upru)%q_density(2)
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
end module
