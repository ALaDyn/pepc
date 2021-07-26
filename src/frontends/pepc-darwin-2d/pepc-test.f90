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
  use module_pepc_types
  use module_pepc_kinds
  use module_timings
  use module_debug

  ! frontend helper routines
  !use module_globals
  use helper
  use encap
  use module_interaction_specific, only : theta2
  use field_helper               , only : compute_field,setup_field_grid,pepc_setup
  !use module_diagnostic

  implicit none

  ! timing variables
  real*8 :: timer(5)
  real*8 :: t1, t2,delta,t
  type(t_particle), allocatable   :: q(:),r(:)
  type(field_grid_t)              :: field_grid
  type(pepc_pars_t)               :: pepc_pars
  integer(kind_particle)          :: ii
  ! control variable
  logical :: doDiag!

  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-darwin-2d", my_rank, n_ranks, .true., comm=communicator)
  !call pepc_initialize('pepc-darwin-2d', init_mpi=.false., db_level_in=DBG_STATUS, comm=communicator)

  root = my_rank.eq.0

  timer(1) = get_time()

  call set_parameter()
  np    = size(particles, kind=kind_particle)

  call pepc_setup(pepc_pars)
  call setup_field_grid(field_grid, pepc_pars%pepc_comm)
  call init_particles(particles,field_grid)

  timer(2) = get_time()
  t        = 0
  

  if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)

  delta = one/real( nt - 1 , kind=kind_particle)     
  
  
  do step=1, nt
      
    theta2 =  real( (step-1) , kind=kind_particle )*delta  
    

    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step   :", step
      write(*,'(a,es12.4)') " ====== simulation time  :", step * dt
      write(*,'(a,es12.4)') " ====== Theta squared    :", theta2
    end if

    timer(3) = get_time()
    
    
    call pepc_particleresults_clear(particles)
    timer(1) = get_time()
    t1 = get_time()

    call pepc_grow_tree(particles)

    t2 = get_time()
    if(root) write(*,'(a,es12.4)') " ====== tree grow time  :", t2-t1
    t1 = get_time()
    call pepc_traverse_tree(particles)
    t2 = get_time()
    if(root) write(*,'(a,es12.4)') " ====== tree walk time  :", t2-t1
    call pepc_restore_particles(particles)
    call pepc_timber_tree()

    
    timer(5) = get_time()
    t        = t + timer(5) - timer(1)
    
!    call normalize(np, particles)
    call test_particles()
    
  end do

  deallocate(particles)

  timer(5) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time [s]: ", timer(5) - timer(2) 
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

