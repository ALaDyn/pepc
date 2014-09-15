! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

program pepc
  ! pepc modules
  use module_pepc
  use module_pepc_kinds
  use module_pepc_types
  use module_directsum
  use module_timings
  use module_debug
  use module_interaction_specific, only : theta2

  ! frontend helper routines
  use helper
  implicit none
  include 'mpif.h'

  ! particle and error related variables
  type(t_particle), allocatable   :: particles(:)
  type(t_particle_results), allocatable :: direct_results(:)
  integer(kind_particle) :: ip
  integer(kind_particle), allocatable :: tindx(:)
  real(kind_physics) :: mean_relerrs(2)
  real(kind_physics), allocatable :: relerrs(:,:)
  real(kind_physics), allocatable :: gathered_relerrs(:,:)

  ! initialize pepc library and MPI
  call pepc_initialize("pepc-cellcell", my_rank, n_ranks, .true.)

  root = my_rank.eq.0

  call timer_start(t_user_total)
  call timer_start(t_user_init)
  call set_parameter()
  call init_particles(particles)
  call timer_stop(t_user_init)

  if(root) then
    write(*,'(a,es12.4)') " === init time [s]: ", timer_read(t_user_init)
    write(*,*) " "
  end if

  ! calculate the exact solution
  allocate(tindx(size(particles)), relerrs(size(particles),2), gathered_relerrs(tnp, 2))
  mean_relerrs = 0
  relerrs = 0
  tindx = (/ (ip, ip=1,size(particles)) /)
  call timer_start(t_user_directsum)
  call directforce(particles, tindx, size(tindx, 1, kind_particle), direct_results, MPI_COMM_WORLD)
  call timer_stop(t_user_directsum)

  if(root) then
    write(*,'(a,es12.4)') " == direct summation time [s]: ", timer_read(t_user_directsum)
    write(*,*) " "
  end if

  ! calculate approximate solution
  do while (theta2 < 1)
    call pepc_particleresults_clear(particles)

    select case (method)
      case (0)
        call grow_and_traverse(particles)
        call calculate_errors(particles, direct_results, tindx, 2, mean_relerrs, relerrs)
      case (1)
        call calculate_internal(particles)
        call calculate_errors(particles, direct_results, tindx, 1, mean_relerrs, relerrs)
      case (2)
        call calculate_internal(particles)
        call calculate_errors(particles, direct_results, tindx, 1, mean_relerrs, relerrs)

        call pepc_particleresults_clear(particles)

        call grow_and_traverse(particles)
        call calculate_errors(particles, direct_results, tindx, 2, mean_relerrs, relerrs)
      case default
        if(root) write(*,*) "=== Invalid  method specification: ", method
    end select

    call gather_results(relerrs, gathered_relerrs)
    if (root) call write_results(mean_relerrs, gathered_relerrs)

    theta2 = (sqrt(theta2) + 0.2)**2
  end do

  deallocate(direct_results, tindx, relerrs, gathered_relerrs)
  call finalize(particles)

end program pepc

