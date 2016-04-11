! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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
  use module_pepc_types
  use module_pepc_kinds
  use module_timings
  use module_debug

  ! frontend helper routines
  use module_globals
  use module_math
  use helper
  use encap
  use field_helper, only: compute_field,setup_field_grid,pepc_setup
  !use module_diagnostic

  implicit none

  ! timing variables
  real*8 :: timer(5)
  real*8 :: t1, t2!,t
  type(t_particle), allocatable   :: q(:),r(:)
  type(physics_pars_t)            :: physics_pars
  type(field_grid_t)              :: field_grid
  type(pepc_pars_t)               :: pepc_pars
  integer(kind_particle)          :: ii,flag
  real(kind_particle),allocatable:: divergence(:,:),prova
  ! control variable
  logical :: doDiag!

  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-2dd", my_rank, n_ranks, .true., comm=communicator)
  !call pepc_initialize('pepc-2dd', init_mpi=.false., db_level_in=DBG_STATUS, comm=communicator)

  root = my_rank.eq.0

  timer(1) = get_time()

  call set_parameter()

  select case (initial_setup)
          case (2)  !  Default initial setup - 2D random particles with thermal velocity
              call init_particles(particles)
          case (3)  !  2D random particles in a square
              call init_particles_square(particles)
          case (4)  ! 2D random particles in a ring
              call init_particles_ring_2D(particles)
          case (5)  ! wire random distribution
              call init_particles_wire(particles)
          case (6)  ! stripe random distribution
              call init_particles_stripe(particles)
          case (7)  ! disc random distribution
              call init_particles_disc(particles)
          case (8)  ! 2D Two wires distribution
              call init_particles_two_wires(particles)
          case (9)  ! Default initial setup - 3D random particles with thermal velocity
              call init_particles3D(particles)
          case (10)  ! 3D Solenoid
              call init_particles_solenoid(particles)
          case (11)  ! two bodies problem - grav
              call two_bodies_grav(particles

          case default

              call init_particles(particles)

  end select

  call pepc_setup(pepc_pars)
  call setup_field_grid(field_grid, pepc_pars%pepc_comm)

  timer(2) = get_time()

  if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)

  do step=0, nt
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step  :", step
      write(*,'(a,es12.4)') " ====== simulation time :", step * dt
    end if

    timer(3) = get_time()

!    doDiag = MOD(step, diag_interval) .eq. 0

    call pepc_particleresults_clear(particles)
    timer(1) = get_time()
    t1 = get_time()

    call pepc_grow_tree(particles)

    !print *, global_tree%nodes(global_tree%node_root)%interaction_data

    np = size(particles, kind=kind(np))
    t2 = get_time()
    if(root) write(*,'(a,es12.4)') " ====== tree grow time  :", t2-t1
    t1 = get_time()

    call pepc_traverse_tree(particles)
    t2 = get_time()
    if(root) write(*,'(a,es12.4)') " ====== tree walk time  :", t2-t1

    call compute_field(pepc_pars, field_grid, particles)
!    write(*,*) "====== Debug on B   "
!    call test_fields(particles,field_grid,pepc_pars,initial_setup)


!    call write_field_rectilinear_grid(particles)

    !if(doDiag .and. domain_output) call write_domain(particles)

!    if (dbg(DBG_STATS)) call pepc_statistics(step)
    call pepc_timber_tree()

    !call pepc_restore_particles(np, particles)

    !if(doDiag .and. particle_test) call test_particles()

!    if(root .and. (step == 0 ) ) call reordering(particles,q)
!    if(root) call write_particles(q)
!    if(root .and. (step == 0 ) ) then
!        call write_field(particles)
!    endif
    write(*,*) "====== writing field on grid   "
    call write_field_on_grid_ascii(field_grid)
    !call write_field_on_grid(field_grid)
    write(*,*) "====== particles   "
!    call write_particles(particles)
    call write_particle_ascii(particles)
    flag = 2
!    call div( field_grid, divergence, flag )
!    call curl( field_grid, divergence,divergence,divergence, flag )

    !prova  = real(doDiag,8)!*10.0 - (1- doDiag )*20
!    write(*,*) "prova value=======", prova
!    call test_particles()
!    call push_particles(particles)
!
!    !if(particle_filter) call filter_particles(particles)
    call energy(particles)
    timer(4) = get_time()
    if(root) write(*,'(a,es12.4)') " == time in step [s]                              : ", timer(4) - timer(3)

    call timings_GatherAndOutput(step, 0)

  end do

  deallocate(particles)

  timer(5) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time [s]: ",timer(5) - timer(1)
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

