! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

module module_diagnostics
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use module_helper
   use module_integrator
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
      points(i)%work = 1.0_kind_physics
      points(i)%data%species = 0
      points(i)%data%mp_int1 = 0
      points(i)%data%age = 0.0_kind_physics
      points(i)%data%v = 0.0_kind_physics
      points(i)%data%q = 0.0_kind_physics
      points(i)%data%m = 1.0_kind_physics
      points(i)%data%b = 0.0_kind_physics
      points(i)%data%f_e = 0.0_kind_physics
      points(i)%data%f_b = 0.0_kind_physics
      points(i)%results%e = 0.0_kind_physics
      points(i)%results%pot = 0.0_kind_physics
    end do
  end subroutine torus_diagnostic_xz_grid

  subroutine torus_diagnostic_xz_breakdown(major_radius, minor_radius, N_point_1D, points)
    ! Diagnostic_grid populates domain with points at regular distance, in xz plane
    implicit none
    real(kind_physics), intent(in) :: major_radius, minor_radius
    integer, intent(in) :: N_point_1D ! denotes number of subdivisions in 1 dimension
    type(t_particle), allocatable, intent(inout) :: points(:)
    real(kind_physics) :: increments, temp_dist(2), dist_sum_x, &
                          dist_sum_z, x_length
    integer :: total_points_1D, total_points

    x_length = 2.0_8*minor_radius
    total_points_1D = N_point_1D !(2**subdivisions + 1)
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
      points(i)%work = 0.0_kind_physics
      points(i)%data%species = 5
      points(i)%data%mp_int1 = 0
      points(i)%data%age = 0.0_kind_physics
      points(i)%data%v = 0.0_kind_physics
      points(i)%data%q = 0.0_kind_physics
      points(i)%data%m = 1.0_kind_physics
      points(i)%data%b = 0.0_kind_physics
      points(i)%data%f_e = 0.0_kind_physics
      points(i)%data%f_b = 0.0_kind_physics
      points(i)%results%e = 0.0_kind_physics
      points(i)%results%pot = 0.0_kind_physics
    end do
    !print *, points(1)%x(1), points(1)%x(3), points(size(points))%x(1), points(size(points))%x(3)
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

  subroutine charge_poloidal_distribution(particle_list, table, global_table, n_total, t_step)
    use mpi
    implicit none
    type(t_particle), allocatable, intent(in) :: particle_list(:)
    character(len=255) :: filename
    real(kind_physics), allocatable, intent(inout) :: table(:,:), global_table(:,:)
    integer(kind_particle), intent(in) :: n_total
    real(kind_physics) :: coords(3), z, R
    integer :: n_parts, l, local_array_size
    integer, allocatable :: offsets(:), receive_cnt(:)
    integer :: rank_size, t_step

    ! Table structure: charge 
    !                  radial distance
    !                  vertical distance (z height)

    if (root) then
      call MPI_COMM_SIZE(MPI_COMM_WORLD, rank_size, ierr)
      allocate(offsets(rank_size))
      allocate(receive_cnt(rank_size))
      allocate(global_table(3,n_total))
      offsets = 0
      receive_cnt = 0
      global_table = 0.0_kind_physics
      print *, "Starting poloidal mapping", size(global_table)
    else
      allocate(global_table(0,0))
      allocate(offsets(0))
      allocate(receive_cnt(0))
    end if

    n_parts = size(particle_list)
    allocate(table(3,n_parts))
    table = 0.0_kind_physics 

    do l = 1, n_parts
      coords = particle_list(l)%x

      ! Scaling to meters
      R = sqrt(coords(1)**2 + coords(2)**2)*(c*1e-12)
      z = coords(3)*(c*1e-12)

      table(1,l) = particle_list(l)%data%q
      table(2,l) = R
      table(3,l) = z
    end do

    local_array_size = size(table)
    call MPI_GATHER(local_array_size, 1, MPI_INT, receive_cnt, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      do l = 2, rank_size
        offsets(l) = offsets(l-1) + receive_cnt(l-1)
      end do
      do l = 1, rank_size
        print *, receive_cnt(l), offsets(l)
      end do
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(table, local_array_size, MPI_KIND_PHYSICS,&
                     global_table, receive_cnt, offsets, MPI_KIND_PHYSICS, &
                     0, MPI_COMM_WORLD, ierr) 

    if (root) then
      print *, "Start writing file"
      write(filename, '(A16,I10.10,A4)') 'charge_pol_dist_', t_step, '.txt'
      filename = trim(filename)
      open(53, file=filename, action='WRITE', position='append')
      write(53, *) 'charge    ', 'R(m)    ', 'z(m)    '
    
      do l = 1, n_total
        write(53, *) global_table(1, l), global_table(2, l), global_table(3, l)
      end do
      close(53)
    end if

    deallocate(global_table)
    deallocate(table)
    deallocate(receive_cnt)
    deallocate(offsets)
  end subroutine charge_poloidal_distribution

  subroutine define_ionisation_array(ionisation_array, ncell_r, ncell_z, resume_logic, t_step_in)
    implicit none
    real(kind_physics), allocatable, intent(inout) :: ionisation_array(:)
    integer, intent(in) :: ncell_r, ncell_z
    integer, optional, intent(in) :: resume_logic, t_step_in
    character(len=255) :: filename
    real(kind_physics) :: dummy1, dummy2
    integer :: N_cell, i_cell
    
    N_cell = ncell_r * ncell_z
    allocate(ionisation_array(N_cell))

    if (resume_logic == 1) then
      write(filename, '(A10,I10.10,A4)') 'ion_count_', t_step_in, '.txt'
      open(53, file=trim("./")//filename, action='READ')
      read(53, *) ! Skipping first line

      do i_cell = 1, N_cell
        read(53, *) dummy1, dummy2, ionisation_array(i_cell)
      end do
      close(53)
    else
      do i_cell = 1, N_cell
        ionisation_array(i_cell) = 0.0_kind_physics
      end do
    end if
  end subroutine define_ionisation_array

  subroutine record_ionisation_count(ionisation_array, ncell_r, ncell_z, start_r, start_z, r_length, z_length, particle)
    implicit none
    real(kind_physics), allocatable, intent(inout) :: ionisation_array(:)
    real(kind_physics), intent(in) :: r_length, z_length, start_r, start_z
    integer, intent(in) :: ncell_r, ncell_z
    type(t_particle), intent(in) :: particle
    real(kind_physics) :: dr, dz, r_pos, z_pos, xp, yp, top_left(2)
    integer :: r_i, z_i, cell_num

    ! start_r & start_z refers to the upper left corner of the poloidal plane grid (not the coordinate of the cell center).
    ! start_r, start_z, r_length, z_length in meters
    top_left(1) = start_r/(c*1e-12)
    top_left(2) = start_z/(c*1e-12)
    dr = r_length/(ncell_r*c*1e-12)
    dz = z_length/(ncell_z*c*1e-12)

    xp = particle%x(1)
    yp = particle%x(2)
    z_pos = particle%x(3)
    r_pos = sqrt(xp**2 + yp**2)

    ! calculate which cell in x and z direction the current particle belongs to.
    r_i = floor(abs(r_pos - top_left(1))/dr) + 1
    z_i = floor(abs(z_pos - top_left(2))/dz) + 1

    ! calculate ionisation_array cell index that particle belongs to
    cell_num = (z_i - 1)*ncell_r + r_i
    ionisation_array(cell_num) = ionisation_array(cell_num) + 1
  end subroutine record_ionisation_count

  subroutine write_ionisation_count(ionisation_array, t_step, ncell_r, ncell_z, start_r, start_z, r_length, z_length)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(in) :: ionisation_array(:)
    integer, intent(in) :: ncell_r, ncell_z, t_step
    real(kind_physics), intent(in) :: start_r, start_z, r_length, z_length
    real(kind_physics), allocatable :: output_array(:,:), temp_global_array(:)
    real(kind_physics) :: dr, dz
    character(len=255) :: filename
    integer :: N_cell, r_i, z_i, i_cell

    N_cell = ncell_r * ncell_z

    if (root) then      
      allocate(output_array(2, N_cell))
      allocate(temp_global_array(N_cell))
      dr = r_length/ncell_r
      dz = z_length/ncell_z

      do i_cell = 1, N_cell
        r_i = mod(i_cell,ncell_r) - 1
        z_i = i_cell/ncell_r
      
        output_array(1, i_cell) = start_r + r_i*dr + 0.5*dr
        output_array(2, i_cell) = start_z - z_i*dr - 0.5*dz
      end do

    else 
      allocate(output_array(1, 1))
      allocate(temp_global_array(1))
    end if

    call MPI_REDUCE(ionisation_array, temp_global_array, N_cell, MPI_KIND_PHYSICS, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      write(filename, '(A10,I10.10,A4)') 'ion_count_', t_step, '.txt'
      open(107, file=filename, action='WRITE', position='append')
      write(107, *) 'R(m)   ', 'z(m)   ', 'count'
    
      do i_cell = 1, N_cell
        write(107, *) output_array(1, i_cell), output_array(2, i_cell), temp_global_array(i_cell)
      end do
      close(107)
    end if

    deallocate(output_array)
    deallocate(temp_global_array)
  end subroutine 

  subroutine unit_vector_distribution(particle_list, table, N_theta, N_phi)
    implicit none
    type(t_particle), allocatable, intent(in) :: particle_list(:)
    integer, intent(in) :: N_theta, N_phi
    real(kind_physics), allocatable, intent(inout) :: table(:,:)
    integer :: N_parts, entry_row, entry_col, l
    real(kind_physics) :: x(3), v(3), x_mag, vel_mag, theta, phi, d_phi, d_theta, weight
    real(kind_physics) :: rot_axis(3)
    ! Table structure:
    ! columns are the range of theta (poloidal angle)
    ! rows are range of phi (toroidal angle)
    ! entry of the table denotes the counts of particles within that phi & theta

    allocate(table(N_phi, N_theta))
    d_phi = 2.*pi/N_phi
    d_theta = pi/N_theta

    table = 0.0_kind_physics
    ! theta = 0.0 aligns with +Z axis
    ! phi = 0.0 aligns with +X axis

    rot_axis = 0.0_kind_physics
    rot_axis(3) = 1.0_kind_physics
    N_parts = size(particle_list)
    do l = 1, N_parts
      x = particle_list(l)%x
      v = particle_list(l)%data%v
      weight = abs(particle_list(l)%data%q)
      x_mag = sqrt(dot_product(x,x))
      vel_mag = sqrt(dot_product(v,v))
      
      if (particle_list(l)%data%species == 0) then
      ! particle velocity NOT rotated to X-Z plane
      ! count directly added to the table entries.
        !call angles_calc(v, vel_mag, theta, phi)
        !entry_row = ceiling(phi/d_phi)
        !entry_col = ceiling(theta/d_theta)
        !table(entry_row, entry_col) = table(entry_row, entry_col) + weight

      ! particle velocity rotated to X-Z plane, particles are centered around (0, 0, 0) coord
      ! once particle velocity is aligned with respective to E field,
      ! count is then added to the table entries.
        call angles_calc(x, x_mag, theta, phi)
        v = Rodriguez_rotation(-phi, rot_axis, v)
        call angles_calc(v, vel_mag, theta, phi)
        entry_row = ceiling(phi/d_phi)
        entry_col = ceiling(theta/d_theta)
        table(entry_row, entry_col) = table(entry_row, entry_col) + weight
      end if
    end do
  end subroutine unit_vector_distribution

  subroutine toroidal_weight_distribution(particle_list, table, N_phi)
    implicit none
    type(t_particle), allocatable, intent(in) :: particle_list(:)
    integer, intent(in) :: N_phi
    real(kind_physics), allocatable, intent(inout) :: table(:,:)
    integer :: N_parts, l, idx
    real(kind_physics) :: d_phi, x(3), x_mag, theta, phi, weight
    
    ! Structure of table, column 1 is the total simulated particle,
    ! column 2 is total weight.
    allocate(table(N_phi,2))
    table = 0.0_kind_physics
    N_parts = size(particle_list)
    d_phi = 2.*pi/N_phi

    do l = 1, N_parts
      x = particle_list(l)%x
      x_mag = sqrt(dot_product(x,x))
      weight = abs(particle_list(l)%data%q)

      if (particle_list(l)%data%species == 0) then
        call angles_calc(x, x_mag, theta, phi)
        idx = ceiling(phi/d_phi)
        table(idx,1) = table(idx,1) + 1._kind_physics
        table(idx,2) = table(idx,2) + weight
      end if
    end do
  end subroutine toroidal_weight_distribution

  subroutine toroidal_max_weights(particle_list, table, table1, N_phi)
    implicit none
    type(t_particle), allocatable, intent(in) :: particle_list(:)
    integer, intent(in) :: N_phi
    real(kind_physics), allocatable, intent(inout) :: table(:), table1(:)
    integer :: N_parts, l, idx
    real(kind_physics) :: d_phi, x(3), x_mag, theta, phi, weight
    
    ! Structure of table, column 1 is the min weight,
    ! column 2 is max weight.
    allocate(table(N_phi))
    allocate(table1(N_phi))
    table = 0.0_kind_physics
    table1 = 0.0_kind_physics
    N_parts = size(particle_list)
    d_phi = 2.*pi/N_phi

    do l = 1, N_phi
      table(l) = 1e16
      table1(l) = -1e16
    end do

    do l = 1, N_parts
      x = particle_list(l)%x
      x_mag = sqrt(dot_product(x,x))
      weight = abs(particle_list(l)%data%q)

      if (particle_list(l)%data%species == 0) then
        call angles_calc(x, x_mag, theta, phi)
        idx = ceiling(phi/d_phi)
        if (table(idx) > weight) then
          table(idx) = weight
        end if
        
        if (table1(idx) < weight) then
          table1(idx) = weight
        end if 
      end if
    end do
  end subroutine toroidal_max_weights

  subroutine V_par_perp_calculation(particle_list, table)
    implicit none
    type(t_particle), allocatable, intent(in) :: particle_list(:)
    real(kind_physics), allocatable, intent(inout) :: table(:,:)
    integer :: N_parts, l
    real(kind_physics) :: ref_axis(3), sign_val
    real(kind_physics) :: x(3), x_mag, weight, Vpar(3), Vperp(3), temp_vel(3)
    real(kind_physics) :: v(3), Vsquared, direction(3)

    ref_axis = 0.0_kind_physics
    ref_axis(3) = 1.0_kind_physics

    N_parts = size(particle_list)
    allocate(table(4, N_parts))
    table = 0.0_kind_physics

    do l = 1, N_parts
      x = particle_list(l)%x
      x_mag = sqrt(dot_product(x,x))
      weight = abs(particle_list(l)%data%q)
      v = particle_list(l)%data%v
      Vsquared = dot_product(v,v)

      if (particle_list(l)%data%species == 0) then
        direction = 0.0_kind_physics
        temp_vel = cross_product(x, ref_axis)
        temp_vel = temp_vel/(sqrt(dot_product(temp_vel, temp_vel)))

        Vpar = dot_product(v,temp_vel)*temp_vel
        
        sign_val = 1.0_kind_physics
        direction = cross_product(x,Vpar)
        if (direction(3) < 0.0) then
          sign_val = -1.0_kind_physics
        end if

        Vpar = sign_val*Vpar
        Vperp = v - Vpar
        table(1, l) = sign_val*sqrt(dot_product(Vpar,Vpar))
        table(2, l) = sqrt(dot_product(Vperp, Vperp))
        table(3, l) = weight
        table(4, l) = Vsquared
      end if
    end do
  end subroutine V_par_perp_calculation

  subroutine V_par_perp_histogram(table, N, total_parts, file_ID, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: table(:,:)
    integer, intent(in) :: N, file_ID
    integer(kind_particle), intent(in) :: total_parts
    character(len=255), intent(in) :: filename
    integer :: rank_size, local_array_size, l, idx
    integer, allocatable :: receive_cnt(:), offsets(:)
    real(kind_physics), allocatable :: global_table(:,:), bin_count(:,:)
    real(kind_physics) :: total_weight, min_limits(4), max_limits(4), Vpar_width, Vperp_width, Vsquared_width
    real(kind_physics) :: Vpar_diff, Vperp_diff, Vpar_mean, Vperp_mean, total_KE, Vsquared, Vsquared_diff

    ! Structure of bin_counts: Vpar_bin_mid, Vpar_cnt, Vperp_bin_mid, Vperp_cnt, V^2, V^2 cnt
    ! min_limits and max_limits: Vpar, Vperp, weight, V^2
    if (root) then
      call MPI_COMM_SIZE(MPI_COMM_WORLD, rank_size, ierr)
      allocate(receive_cnt(rank_size))
      allocate(offsets(rank_size))
      allocate(global_table(4, total_parts))
      allocate(bin_count(6, N))
      receive_cnt = 0
      offsets = 0
      global_table = 0.0_kind_physics
      bin_count = 0.0_kind_physics
    end if

    if (.not. root) then
      allocate(receive_cnt(0))
      allocate(global_table(0,0))
      allocate(offsets(0))
    end if
    call MPI_COMM_SIZE(MPI_COMM_WORLD, rank_size, ierr)
    local_array_size = size(table)
    call MPI_GATHER(local_array_size, 1, MPI_INT, receive_cnt, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      do l = 2, rank_size
        offsets(l) = offsets(l-1) + receive_cnt(l-1)
      end do
    end if

    call MPI_GATHERV(table, local_array_size, MPI_KIND_PHYSICS,&
                     global_table, receive_cnt, offsets, MPI_KIND_PHYSICS, &
                     0, MPI_COMM_WORLD, ierr)

    if (root) then
      min_limits = minval(global_table, 2)
      max_limits = maxval(global_table, 2)

      Vpar_width = (max_limits(1) - min_limits(1))/N
      Vperp_width = (max_limits(2) - min_limits(2))/N
      Vsquared_width = (max_limits(4) - min_limits(4))/N 

      bin_count(1, 1) = min_limits(1) + Vpar_width*0.5
      bin_count(3, 1) = min_limits(2) + Vperp_width*0.5
      bin_count(5, 1) = min_limits(4) + Vsquared_width*0.5
      do l = 2, N
        bin_count(1, l) = bin_count(1, l-1) + Vpar_width
        bin_count(3, l) = bin_count(3, l-1) + Vperp_width
        bin_count(5, l) = bin_count(5, l-1) + Vsquared_width
      end do

      Vpar_mean = 0.0_kind_physics
      Vperp_mean = 0.0_kind_physics
      total_weight = 0.0_kind_physics
      total_KE = 0.0_kind_physics
      do l = 1, total_parts
        Vpar_diff = global_table(1,l) - min_limits(1)
        idx = ceiling(Vpar_diff/Vpar_width)
        if (idx < 1) then
          idx = 1
        else if (idx > N) then
          idx = N - 1
        end if

        bin_count(2, idx) = bin_count(2, idx) + global_table(3,l)

        Vperp_diff = global_table(2,l) - min_limits(2)
        idx = ceiling(Vperp_diff/Vperp_width)
        if (idx < 1) then
          idx = 1
        else if (idx > N) then
          idx = N - 1
        end if

        Vpar_mean = Vpar_mean + global_table(1,l)*global_table(3,l)
        Vperp_mean = Vperp_mean + global_table(2,l)*global_table(3,l)
        bin_count(4, idx) = bin_count(4, idx) + global_table(3,l)
        total_weight = total_weight + global_table(3,l)

        Vsquared_diff = global_table(4,l) - min_limits(4)
        idx = ceiling(Vsquared_diff/Vsquared_width)
        if (idx < 1) then
          idx = 1
        else if (idx > N) then
          idx = N - 1
        end if

        bin_count(6, idx) = bin_count(6, idx) + global_table(3,l)

        total_KE = total_KE + global_table(4,l)*global_table(3,l)
      end do
      Vpar_mean = Vpar_mean/total_weight
      Vperp_mean = Vperp_mean/total_weight
      total_KE = total_KE/total_weight

      open(file_ID, file=filename, action='WRITE', position='append')
      write(file_ID, *) 'Vpar_mid    ', 'Vpar_cnt   ', 'Vperp_mid    ', 'Vperp_cnt    ', 'Vsquared_mid    ', 'Vsquared_cnt'
    
      write(file_ID, *) 'Mean: ', Vpar_mean, Vperp_mean, total_KE
      do l = 1, N
        total_weight = total_weight + bin_count(2, l)
        write(file_ID, *) bin_count(1, l), bin_count(2, l), bin_count(3, l), bin_count(4, l), bin_count(5, l), bin_count(6, l)
      end do
      close(file_ID)

      deallocate(bin_count)
    end if
    deallocate(global_table)
    deallocate(offsets)
    deallocate(receive_cnt)

    deallocate(table)
  end subroutine V_par_perp_histogram

  subroutine V_par_perp_histogram_fix(table, N, total_parts, file_ID, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: table(:,:)
    integer, intent(in) :: N, file_ID
    integer(kind_particle), intent(in) :: total_parts
    character(len=255), intent(in) :: filename
    integer :: rank_size, local_array_size, l, idx
    integer, allocatable :: receive_cnt(:), offsets(:)
    real(kind_physics), allocatable :: global_table(:,:), bin_count(:,:)
    real(kind_physics) :: total_weight, min_limits(4), max_limits(4), Vpar_width, Vperp_width, Vsquared_width
    real(kind_physics) :: Vpar_diff, Vperp_diff, Vpar_mean, Vperp_mean, total_KE, Vsquared, Vsquared_diff

    ! Structure of bin_counts: Vpar_bin_mid, Vpar_cnt, Vperp_bin_mid, Vperp_cnt, V^2, V^2 cnt
    ! min_limits and max_limits: Vpar, Vperp, weight, V^2
    if (root) then
      call MPI_COMM_SIZE(MPI_COMM_WORLD, rank_size, ierr)
      allocate(receive_cnt(rank_size))
      allocate(offsets(rank_size))
      allocate(global_table(4, total_parts))
      allocate(bin_count(6, N))
      receive_cnt = 0
      offsets = 0
      global_table = 0.0_kind_physics
      bin_count = 0.0_kind_physics
    end if

    if (.not. root) then
      allocate(receive_cnt(0))
      allocate(global_table(0,0))
      allocate(offsets(0))
    end if
    call MPI_COMM_SIZE(MPI_COMM_WORLD, rank_size, ierr)
    local_array_size = size(table)
    call MPI_GATHER(local_array_size, 1, MPI_INT, receive_cnt, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      do l = 2, rank_size
        offsets(l) = offsets(l-1) + receive_cnt(l-1)
      end do
    end if

    call MPI_GATHERV(table, local_array_size, MPI_KIND_PHYSICS,&
                     global_table, receive_cnt, offsets, MPI_KIND_PHYSICS, &
                     0, MPI_COMM_WORLD, ierr)

    if (root) then
      min_limits = minval(global_table, 2)
      max_limits = maxval(global_table, 2)

      Vpar_width = (max_limits(1) - min_limits(1))/N
      Vperp_width = (max_limits(2) - min_limits(2))/N
      Vsquared_width = 1.9208321467648882E-002_kind_physics/N 
      print *, max_limits(4), min_limits(4) 

      bin_count(1, 1) = min_limits(1) + Vpar_width*0.5
      bin_count(3, 1) = min_limits(2) + Vperp_width*0.5
      bin_count(5, 1) = min_limits(4) + Vsquared_width*0.5
      do l = 2, N
        bin_count(1, l) = bin_count(1, l-1) + Vpar_width
        bin_count(3, l) = bin_count(3, l-1) + Vperp_width
        bin_count(5, l) = bin_count(5, l-1) + Vsquared_width
      end do

      Vpar_mean = 0.0_kind_physics
      Vperp_mean = 0.0_kind_physics
      total_weight = 0.0_kind_physics
      total_KE = 0.0_kind_physics
      do l = 1, total_parts
        Vpar_diff = global_table(1,l) - min_limits(1)
        idx = ceiling(Vpar_diff/Vpar_width)
        if (idx < 1) then
          idx = 1
        else if (idx > N) then
          idx = N - 1
        end if

        bin_count(2, idx) = bin_count(2, idx) + global_table(3,l)

        Vperp_diff = global_table(2,l) - min_limits(2)
        idx = ceiling(Vperp_diff/Vperp_width)
        if (idx < 1) then
          idx = 1
        else if (idx > N) then
          idx = N - 1
        end if

        Vpar_mean = Vpar_mean + global_table(1,l)*global_table(3,l)
        Vperp_mean = Vperp_mean + global_table(2,l)*global_table(3,l)
        bin_count(4, idx) = bin_count(4, idx) + global_table(3,l)
        total_weight = total_weight + global_table(3,l)

        Vsquared_diff = global_table(4,l)! - min_limits(4)
        idx = ceiling(Vsquared_diff/Vsquared_width)
        if (idx < 1) then
          idx = 1
        else if (idx > N) then
          idx = N - 1
        end if

        bin_count(6, idx) = bin_count(6, idx) + global_table(3,l)

        total_KE = total_KE + global_table(4,l)*global_table(3,l)
      end do
      Vpar_mean = Vpar_mean/total_weight
      Vperp_mean = Vperp_mean/total_weight
      total_KE = total_KE/total_weight

      open(file_ID, file=filename, action='WRITE', position='append')
      write(file_ID, *) 'Vpar_mid    ', 'Vpar_cnt   ', 'Vperp_mid    ', 'Vperp_cnt    ', 'Vsquared_mid    ', 'Vsquared_cnt'
    
      write(file_ID, *) 'Mean: ', Vpar_mean, Vperp_mean, total_KE
      do l = 1, N
        total_weight = total_weight + bin_count(2, l)
        write(file_ID, *) bin_count(1, l), bin_count(2, l), bin_count(3, l), bin_count(4, l), bin_count(5, l), bin_count(6, l)
      end do
      close(file_ID)

      deallocate(bin_count)
    end if
    deallocate(global_table)
    deallocate(offsets)
    deallocate(receive_cnt)

    deallocate(table)
  end subroutine V_par_perp_histogram_fix

  subroutine V_mean_phi_distribution(particle_list, table, N_phi)
    implicit none
    type(t_particle), allocatable, intent(in) :: particle_list(:)
    integer, intent(in) :: N_phi
    real(kind_physics), allocatable, intent(inout) :: table(:,:)
    integer :: N_parts, l, idx
    real(kind_physics) :: d_phi, x(3), x_mag, theta, phi, Vpar(3), Vperp(3), v(3), weight
    real(kind_physics) :: Vsquared, ref_axis(3), direction(3), temp_vel(3), sign_val

    ! Structure of table, column 1 is the total simulated particle,
    ! column 2 is total weight.
    allocate(table(4, N_phi))
    table = 0.0_kind_physics
    N_parts = size(particle_list)
    d_phi = 2.*pi/N_phi

    ref_axis = 0.0_kind_physics
    ref_axis(3) = 1.0_kind_physics

    do l = 1, N_parts
      x = particle_list(l)%x
      x_mag = sqrt(dot_product(x,x))
      weight = abs(particle_list(l)%data%q)
      v = particle_list(l)%data%v
      Vsquared = dot_product(v,v)

      if (particle_list(l)%data%species == 0) then
        direction = 0.0_kind_physics
        temp_vel = cross_product(x, ref_axis)
        temp_vel = temp_vel/(sqrt(dot_product(temp_vel, temp_vel)))

        Vpar = dot_product(v,temp_vel)*temp_vel
        
        sign_val = 1.0_kind_physics
        direction = cross_product(x,Vpar)
        if (direction(3) < 0.0) then
          sign_val = -1.0_kind_physics
        end if

        Vpar = sign_val*Vpar
        Vperp = v - Vpar

        call angles_calc(x, x_mag, theta, phi)
        idx = ceiling(phi/d_phi)
        table(1, idx) = table(1, idx) + 1._kind_physics
        table(2, idx) = table(2, idx) + weight
        table(3, idx) = table(3, idx) + sign_val*sqrt(dot_product(Vpar,Vpar))
        table(4, idx) = table(4, idx) + sqrt(dot_product(Vperp, Vperp))
      end if
    end do
  end subroutine V_mean_phi_distribution

  subroutine V_mean_phi_distribution_gather(table, global_table, file_ID, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: table(:,:), global_table(:,:)
    integer, intent(in) :: file_ID
    character(len=255), intent(in) :: filename
    integer :: N_entries, m
    integer, allocatable :: N_counts(:)
    real(kind_physics) :: d_phi, phi, total_electron

    allocate(N_counts(size(shape(table))))
    N_counts = shape(table)
    N_entries = size(table)

    if (root) then
      allocate(global_table(N_counts(1), N_counts(2)))
      global_table = 0.0_kind_physics
    end if
    call MPI_REDUCE(table, global_table, N_entries, MPI_KIND_PHYSICS, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      d_phi = 2.0*pi/N_counts(2)
      open(file_ID, file=filename, action='WRITE', position='append')
      write(file_ID, *) 'phi    ', 'Vpar_mean  ', 'Vperp_mean    '
    
      total_electron = 0.0_kind_physics
      do m = 1, N_counts(2)
        phi = m*d_phi - d_phi*0.5
        write(file_ID, *) phi, global_table(3, m)/global_table(1, m), global_table(4, m)/global_table(1, m)
        total_electron = total_electron + global_table(2,m)
      end do
      print *, "total counted electrons: ", total_electron
      close(file_ID)
      deallocate(global_table)
    end if

    deallocate(table)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine V_mean_phi_distribution_gather

  subroutine gather_weights_tables(table, global_table, file_ID, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: table(:,:), global_table(:,:)
    integer, intent(in) :: file_ID
    character(len=255), intent(in) :: filename
    integer :: N_entries, m
    integer, allocatable :: N_counts(:)

    allocate(N_counts(size(shape(table))))
    N_counts = shape(table)
    N_entries = size(table)

    if (root) then
      allocate(global_table(N_counts(1), N_counts(2)))
      global_table = 0.0_kind_physics
    end if
    call MPI_REDUCE(table, global_table, N_entries, MPI_KIND_PHYSICS, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      open(file_ID, file=filename, action='WRITE', position='append')
      write(file_ID, *) 'Counts    ', 'total weight   '
    
      do m = 1, N_counts(1)
        write(file_ID, *) global_table(m,1), global_table(m,2)
      end do
      close(file_ID)
      deallocate(global_table)
    end if

    deallocate(table)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine gather_weights_tables

  subroutine gather_minmaxWeights_tables(table, table1, global_table, global_table1, file_ID, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: table(:), global_table(:)
    real(kind_physics), allocatable, intent(inout) :: table1(:), global_table1(:)
    integer, intent(in) :: file_ID
    character(len=255), intent(in) :: filename
    integer :: N_entries, m

    N_entries = size(table)

    if (root) then
      allocate(global_table(N_entries))
      allocate(global_table1(N_entries))
      global_table = 1e16_kind_physics
      global_table1 = -1e16_kind_physics
    end if
    call MPI_REDUCE(table, global_table, N_entries, MPI_KIND_PHYSICS, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(table1, global_table1, N_entries, MPI_KIND_PHYSICS, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      open(file_ID, file=filename, action='WRITE', position='append')
      write(file_ID, *) 'min_weight    ', 'max_weight   '
    
      do m = 1, N_entries
        write(file_ID, *) global_table(m), global_table1(m)
      end do
      close(file_ID)
      deallocate(global_table)
      deallocate(global_table1)
    end if

    deallocate(table)
    deallocate(table1)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine gather_minmaxWeights_tables

  subroutine gather_spherical_angle_tables(table, global_table, file_ID, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: table(:,:), global_table(:,:)
    integer, intent(in) :: file_ID
    character(len=255), intent(in) :: filename
    integer, allocatable :: N_counts(:)
    integer :: N_entries, m, n
    real(kind_physics) :: phi, theta, d_phi, d_theta

    allocate(N_counts(size(shape(table))))
    N_counts = shape(table)
    N_entries = size(table)

    if (root) then
      allocate(global_table(N_counts(1), N_counts(2)))
      global_table = 0.0_kind_physics
    end if
    call MPI_REDUCE(table, global_table, N_entries, MPI_KIND_PHYSICS, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (root) then
      d_phi = 2.0*pi/N_counts(1)
      d_theta = pi/N_counts(2)

      open(file_ID, file=filename, action='WRITE', position='append')
      write(file_ID, *) 'PHI    ', 'THETA    ', 'VALUE   '
    
      do m = 1, N_counts(1)
        phi = m*d_phi - d_phi*0.5
        do n = 1, N_counts(2)
          theta = n*d_theta - d_theta*0.5
          write(file_ID, *) phi, theta, global_table(m,n)
        end do
      end do
      close(file_ID)
      deallocate(global_table)
    end if

    deallocate(N_counts)
    deallocate(table)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine gather_spherical_angle_tables

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

  subroutine Sum_Vde(particle, vertices)
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
      temp_q_dens = 1.0
      vertices(loll)%q_density(1) = temp_q_dens + vertices(loll)%q_density(1)
      vertices(loll)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(loll)%J_density

      vertices(lorl)%q_density(1) = temp_q_dens + vertices(lorl)%q_density(1)
      vertices(lorl)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(lorl)%J_density

      vertices(upll)%q_density(1) = temp_q_dens + vertices(upll)%q_density(1)
      vertices(upll)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(upll)%J_density

      vertices(uprl)%q_density(1) = temp_q_dens + vertices(uprl)%q_density(1)
      vertices(uprl)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(uprl)%J_density

      vertices(lolu)%q_density(1) = temp_q_dens + vertices(lolu)%q_density(1)
      vertices(lolu)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(lolu)%J_density

      vertices(loru)%q_density(1) = temp_q_dens + vertices(loru)%q_density(1)
      vertices(loru)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(loru)%J_density

      vertices(uplu)%q_density(1) = temp_q_dens + vertices(uplu)%q_density(1)
      vertices(uplu)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(uplu)%J_density

      vertices(upru)%q_density(1) = temp_q_dens + vertices(upru)%q_density(1)
      vertices(upru)%J_density = temp_q_dens*particle%data%v &
                                 + vertices(upru)%J_density

     !  test = vertices(loll)%q_density(1) + vertices(lorl)%q_density(1) &
     !         + vertices(upll)%q_density(1) + vertices(uprl)%q_density(1) &
     !         + vertices(lolu)%q_density(1) + vertices(loru)%q_density(1) &
     !         + vertices(uplu)%q_density(1) + vertices(upru)%q_density(1)
     !  print *, test
   end if
  end subroutine Sum_Vde

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

  subroutine slice_particles(particles_list, xmin, xmax, ymin, ymax, slab_particles)
    implicit none
    type(t_particle), allocatable, intent(inout) :: particles_list(:)
    type(linked_list_elem), pointer, intent(inout) :: slab_particles
    real(kind_physics), intent(in) :: xmin, xmax, ymin, ymax
    type(linked_list_elem), pointer :: slab_guide
    integer :: n_local_slab, k, l, ll_buffer_size, ll_elem_cnt, remainder, head, tail
    real(kind_physics) :: x_sign_check, y_sign_check, s_xmin, s_xmax, s_ymin, s_ymax

    ll_buffer_size = 50
    n_local_slab = 0
    l = 0

    s_xmin = xmin/(c*1e-12)
    s_xmax = xmax/(c*1e-12)
    s_ymin = ymin/(c*1e-12)
    s_ymax = ymax/(c*1e-12)

    ! Record all particles that lies within the x and y range into a linked list element
    do k = 1, size(particles_list)
      x_sign_check = (particles_list(k)%x(1) - s_xmin)*(particles_list(k)%x(1) - s_xmax)
      y_sign_check = (particles_list(k)%x(2) - s_ymin)*(particles_list(k)%x(2) - s_ymax)
      if (x_sign_check < 0.0_kind_physics .and. y_sign_check < 0.0_kind_physics) then
        if (n_local_slab == 0) then
          call allocate_ll_buffer(ll_buffer_size, slab_particles)
          slab_guide => slab_particles
          nullify(slab_guide%next)
        end if

        if ((n_local_slab .ne. 0) .and. (mod(n_local_slab,ll_buffer_size) == 0)) then
          call allocate_ll_buffer(ll_buffer_size, slab_guide%next)
          slab_guide => slab_guide%next
          nullify(slab_guide%next)
          l = 0
        end if

        n_local_slab = n_local_slab + 1
        l = l + 1
        slab_guide%tmp_particles(l) = particles_list(k)   
      end if
    end do
    
    deallocate(particles_list)
    allocate(particles_list(n_local_slab))

    if (n_local_slab .ne. 0) then
      head = 1
      
      ll_elem_cnt = CEILING(real(n_local_slab)/real(ll_buffer_size))
      remainder = MOD(n_local_slab, ll_buffer_size)
      tail = ll_buffer_size 

      i = 1
      slab_guide => slab_particles
      do while (associated(slab_guide))
        if (i /= ll_elem_cnt) then
          particles_list(head:tail) = slab_guide%tmp_particles
        else if ((i == ll_elem_cnt) .and. (remainder == 0)) then
          particles_list(head:n_local_slab) = slab_guide%tmp_particles
        else
          particles_list(head:n_local_slab) = slab_guide%tmp_particles(1:(remainder))
        end if

        head = tail + 1
        tail = tail + ll_buffer_size
        i = i + 1

        slab_guide => slab_guide%next
      end do
      nullify (slab_guide)
    end if
  end subroutine slice_particles

  subroutine torus_cut(particles_list, torus_r, torus_z, minor_r, remain_particles)
    implicit none
    type(t_particle), allocatable, intent(inout) :: particles_list(:)
    type(linked_list_elem), pointer, intent(inout) :: remain_particles
    real(kind_physics), intent(in) :: torus_r, torus_z, minor_r
    type(linked_list_elem), pointer :: slab_guide
    integer :: n_local_slab, k, l, ll_buffer_size, ll_elem_cnt, remainder, head, tail
    real(kind_physics) :: s_torus_r, s_torus_z, s_minor_r, particle_r, delta_r, delta_z
    real(kind_physics) :: delta_radius

    ll_buffer_size = 50
    n_local_slab = 0
    l = 0

    s_torus_r = torus_r/(c*1e-12)
    s_torus_z = torus_z/(c*1e-12)
    s_minor_r = minor_r/(c*1e-12)

    ! Record all particles whose delta_radius is less than the given minor_r into a linked list element
    do k = 1, size(particles_list)
      particle_r = sqrt(particles_list(k)%x(1)**2 + particles_list(k)%x(2)**2)
      delta_r = particle_r - s_torus_r
      delta_z = particles_list(k)%x(3) - s_torus_z
      delta_radius = sqrt(delta_r**2 + delta_z**2)

      if (delta_radius < s_minor_r) then
        if (n_local_slab == 0) then
          call allocate_ll_buffer(ll_buffer_size, slab_particles)
          slab_guide => slab_particles
          nullify(slab_guide%next)
        end if

        if ((n_local_slab .ne. 0) .and. (mod(n_local_slab,ll_buffer_size) == 0)) then
          call allocate_ll_buffer(ll_buffer_size, slab_guide%next)
          slab_guide => slab_guide%next
          nullify(slab_guide%next)
          l = 0
        end if

        n_local_slab = n_local_slab + 1
        l = l + 1
        slab_guide%tmp_particles(l) = particles_list(k)   
      end if
    end do
    
    deallocate(particles_list)
    allocate(particles_list(n_local_slab))

    if (n_local_slab .ne. 0) then
      head = 1
      
      ll_elem_cnt = CEILING(real(n_local_slab)/real(ll_buffer_size))
      remainder = MOD(n_local_slab, ll_buffer_size)
      tail = ll_buffer_size 

      i = 1
      slab_guide => slab_particles
      do while (associated(slab_guide))
        if (i /= ll_elem_cnt) then
          particles_list(head:tail) = slab_guide%tmp_particles
        else if ((i == ll_elem_cnt) .and. (remainder == 0)) then
          particles_list(head:n_local_slab) = slab_guide%tmp_particles
        else
          particles_list(head:n_local_slab) = slab_guide%tmp_particles(1:(remainder))
        end if

        head = tail + 1
        tail = tail + ll_buffer_size
        i = i + 1

        slab_guide => slab_guide%next
      end do
      nullify (slab_guide)
    end if
  end subroutine torus_cut

  subroutine streakline_integral(particle, init_dx, gradB_lim, iters, major_radius, minor_radius, tok_origin, conn_length)
    implicit none
    type(t_particle), intent(inout) :: particle
    real(kind_physics), intent(in) :: init_dx, gradB_lim, major_radius, minor_radius
    real(kind_physics), intent(in) :: tok_origin(3)
    real(kind_physics), intent(out) :: conn_length
    real(kind_physics) :: gradB, B_mag, new_B_mag, new_pos(3), B(3), new_B(3), pos_a
    real(kind_physics) :: r_xy, r_x, r_y, r_z, init_pos(3), dx, length_lim, total_length
    type(t_particle) :: particle_tmp
    integer(kind_particle) :: iters, i, j

    particle_tmp%x = 0.0_kind_physics
    particle_tmp%data%q = 0.0_kind_physics
    particle_tmp%data%v = 0.0_kind_physics
    particle_tmp%data%m = 0.0_kind_physics
    particle_tmp%data%b = 0.0_kind_physics
    particle_tmp%data%f_e = 0.0_kind_physics
    particle_tmp%data%f_b = 0.0_kind_physics
    particle_tmp%results%e = 0.0_kind_physics
    particle_tmp%results%pot = 0.0_kind_physics

    dx = init_dx
    pos_a = 0.0
    total_length = 0.0
    length_lim = 80000./(c*1e-12) 
    i = 0
    dx = 0.003335_kind_physics! 0.0001335_kind_physics
    do while(pos_a < minor_radius .and. total_length < length_lim)
      call particle_EB_field(particle, external_e, B_pol_grid)
      B_mag = sqrt(dot_product(particle%data%b, particle%data%b))
      B = 1.0_kind_physics*particle%data%b/B_mag

      new_pos = particle%x + B*dx
      particle_tmp%x = new_pos
      call particle_EB_field(particle_tmp, external_e, B_pol_grid)
      new_B_mag = sqrt(dot_product(particle_tmp%data%b, particle_tmp%data%b))
      gradB = (new_B_mag - B_mag)/dx

      ! j = 0
      ! do while(abs(gradB) > gradB_lim .and. dx > 1e-4_kind_physics)
      !   dx = dx*0.5_kind_physics
      !   new_pos = particle%x + B*dx
      !   particle_tmp%x = new_pos
      !   call particle_EB_field(particle_tmp, external_e, B_pol_grid)
      !   new_B_mag = sqrt(dot_product(particle_tmp%data%b, particle_tmp%data%b))
      !   gradB = (new_B_mag - B_mag)/dx

      !   j = j + 1
      ! end do
      total_length = total_length + dx

      particle%x = new_pos
      r_x = particle%x(1) - tok_origin(1)
      r_y = particle%x(2) - tok_origin(2)
      r_z = particle%x(3) - tok_origin(3)
      r_xy = sqrt(r_x**2 + r_y**2)
      r_xy = r_xy - major_radius

      pos_a = sqrt(r_xy**2 + r_z**2)
      i = i + 1
      if (mod(i, 200000000) == 0) i = 0 !print *, "pos: ", i, B, "dx: ", dx, "B_mag: ", B_mag
    end do
    conn_length = total_length
    particle%data%f_b(1) = conn_length
  end subroutine streakline_integral

  subroutine streakline_integral_timing(particle, major_radius, minor_radius, tok_origin, conn_length, time_start, rate, end_time)
    implicit none
    type(t_particle), intent(inout) :: particle
    integer(kind_particle), intent(in) :: time_start
    real(kind_physics), intent(in) :: rate, end_time
    real(kind_physics), intent(in) :: major_radius, minor_radius
    real(kind_physics), intent(in) :: tok_origin(3)
    real(kind_physics), intent(out) :: conn_length
    real(kind_physics) :: B_mag, new_pos(3), B(3), pos_a, elapsed_time
    real(kind_physics) :: r_xy, r_pos(3), r_pos2(3), dx, total_length
    integer(kind_particle) :: time_finish

    total_length = particle%data%f_b(1)
    dx = 0.003335_kind_physics! 0.0001335_kind_physics

    r_pos = particle%x - tok_origin
    r_pos2 = r_pos**2

    r_xy = sqrt(r_pos2(1) + r_pos2(2)) - major_radius
    pos_a = sqrt(r_xy**2 + r_pos2(3))

    do while(pos_a < minor_radius .and. elapsed_time < end_time)
      call particle_EB_field(particle, external_e, B_pol_grid)
      B_mag = sqrt(dot_product(particle%data%b, particle%data%b))
      B = 1.0_kind_physics*particle%data%b/B_mag

      new_pos = particle%x + B*dx
      total_length = total_length + dx

      particle%x = new_pos
      r_pos = particle%x - tok_origin
      r_pos2 = r_pos**2
      r_xy = sqrt(r_pos2(1) + r_pos2(2)) - major_radius

      pos_a = sqrt(r_xy**2 + r_pos2(3))
      
      call system_clock(time_finish)
      elapsed_time = (time_finish - time_start)/rate
    end do
    conn_length = total_length
    particle%data%f_b(1) = conn_length
  end subroutine streakline_integral_timing

  subroutine connection_length_output(global_table, file_ID)!, filename)
    use mpi
    implicit none
    real(kind_physics), allocatable, intent(inout) :: global_table(:,:)
    integer, intent(in) :: file_ID
    character(len=255) :: filename
    integer :: N_entries, m
    integer, allocatable :: N_counts(:)    

    ! local_table2(1:3,:) = particle's initial coordinate.
    ! local_table2(4,:)   = connection length
    ! local_table2(5,:)   = particle's init B_mag
    ! local_table2(6,:)   = particle's final B_mag

    allocate(N_counts(size(shape(global_table))))
    N_counts = shape(global_table)
    N_entries = size(global_table)

    if (root) then
      print *, "N_counts", N_counts
      write (filename, '(a,"_",i2.2)') "test_sce3C", itime_in
      open(file_ID, file=filename, action='WRITE', position='append')
      print *, "file opened"
      write(file_ID, *) 'init_x    ', 'init_y   ', 'init_z   ', 'length   ', 'init_Bmag   ', 'final_Bmag   '
      print *, "Start looping"
      do m = 1, N_counts(2)
        write(file_ID, *) global_table(1,m), global_table(2,m), global_table(3,m), global_table(4,m), global_table(5,m), global_table(6,m)
      end do
      close(file_ID)
    end if
    deallocate(N_counts)
  end subroutine connection_length_output

  subroutine circle_points(particle, x_offset, z_offset, radius, n_points)
    implicit none
    type(t_particle), allocatable, intent(inout) :: particle(:)
    real(kind_physics), intent(in) :: x_offset, z_offset, radius
    integer, intent(in) :: n_points
    real(kind_physics) :: theta, x_scaled, z_scaled, radius_scaled, delta_theta
    integer :: l

    allocate(particle(n_points))
    x_scaled = x_offset/(c*1e-12) 
    z_scaled = z_offset/(c*1e-12) 
    radius_scaled = radius/(c*1e-12) 

    delta_theta = 2*pi/n_points
    do l = 1, n_points
      particle(l)%x(1) = x_scaled - radius_scaled*cos(delta_theta*(l-1))
      particle(l)%x(2) = 0.0_kind_physics
      particle(l)%x(3) = z_scaled - radius_scaled*sin(delta_theta*(l-1))

      particle(l)%data%species = 0
      particle(l)%data%mp_int1 = 0
      particle(l)%data%age = 0.0_kind_physics
      particle(l)%work = 1.0_Kind_physics
      particle(l)%data%q = 0.0_kind_physics
      particle(l)%data%v = 0.0_kind_physics
      particle(l)%data%m = 0.0_kind_physics
      particle(l)%data%b = 0.0_kind_physics
      particle(l)%data%f_e = 0.0_kind_physics
      particle(l)%data%f_b = 0.0_kind_physics
      particle(l)%results%e = 0.0_kind_physics
      particle(l)%results%pot = 0.0_kind_physics
    end do
!    particle(n_points+1)%x(1) = x_scaled
!    particle(n_points+1)%x(2) = 0.0_kind_physics
!    particle(n_points+1)%x(3) = z_scaled
!    
!    particle(n_points+1)%data%species = 0
!    particle(n_points+1)%data%age = 0.0_kind_physics
!    particle(n_points+1)%work = 1.0_Kind_physics
!    particle(n_points+1)%data%q = 0.0_kind_physics
!    particle(n_points+1)%data%v = 0.0_kind_physics
!    particle(n_points+1)%data%m = 0.0_kind_physics
!    particle(n_points+1)%data%b = 0.0_kind_physics
!    particle(n_points+1)%data%f_e = 0.0_kind_physics
!    particle(n_points+1)%data%f_b = 0.0_kind_physics
!    particle(n_points+1)%results%e = 0.0_kind_physics
!    particle(n_points+1)%results%pot = 0.0_kind_physics

  end subroutine circle_points
end module
