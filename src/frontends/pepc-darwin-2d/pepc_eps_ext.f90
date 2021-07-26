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
!  real*8 :: timer(5)
!  real*8 :: t1, t2!,t
!  type(t_particle), allocatable   :: q(:),r(:)
!  type(physics_pars_t)            :: physics_pars
!  type(field_grid_t)              :: field_grid
!  type(pepc_pars_t)               :: pepc_pars
!  integer(kind_particle)          :: ii,flag
!  real(kind_particle),allocatable:: divergence(:,:),prova
!  ! control variable
!  logical :: doDiag!

  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-darwin-2d", my_rank, n_ranks, .true., comm=communicator)
  !call pepc_initialize('pepc-darwin-2d', init_mpi=.false., db_level_in=DBG_STATUS, comm=communicator)

  root = my_rank.eq.0


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
          case (8)  ! Solenoidal distribution
              call init_particles_solenoidal(particles)

          case default

              call init_particles(particles)

  end select

  call pepc_particleresults_clear(particles)
  call pepc_grow_tree(particles)
  call pepc_traverse_tree(particles)
  call pepc_timber_tree()
  call estimation_eps2(particles)
  deallocate(particles)


  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

