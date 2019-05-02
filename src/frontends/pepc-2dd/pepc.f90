! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2019 Juelich Supercomputing Centre,
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
  use module_timings
  use module_debug
  use module_globals             , only : adv
  use module_interaction_specific, only : force_law
  use module_tool                , only : copy_particle
  use field_helper               , only : compute_field, setup_field_grid, pepc_setup,write_field_on_grid,write_field_on_grid3D

  use module_integration
  use zufall
  use encap
  use module_diagnostic

  ! frontend helper routines
  use helper
!  use module_mpi_io
  use module_shortcut
  implicit none

  ! timing variables
  real(kind_particle)                        :: timer(5)
  real(kind_particle)                        :: t1, t2,t!,tau,I
!  type(t_particle), allocatable              :: q(:),r(:)
!  type(physics_pars_t)                       :: physics_pars
  type(field_grid_t)                         :: field_grid
  type(pepc_pars_t)                          :: pepc_pars
!  type(pepc_comm_t)                          :: pepc_comm
  integer(kind_particle)                     :: rc = 101
!  integer(kind_particle)                     :: step

  ! control variable
  logical :: dorestart,dointerp, explicit =.false. !dogrid,

  !!! initialize pepc library and MPI

  call pepc_initialize("pepc-2dd", my_rank, n_ranks, .true., comm=communicator)
!  call pepc_initialize("pepc-2dd", my_rank, n_ranks, .true.)

  root = my_rank.eq.0
  timer(1) = get_time()
  call set_parameter()
  np    = size(particles, kind=kind_particle)

  call pepc_setup(pepc_pars)
  call setup_field_grid(field_grid, pepc_pars%pepc_comm)
  call init_particles(particles,field_grid)

  timer(2) = get_time()
  t        = 0
  

  if(root) then
      write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)
      write(*,'(a,es12.4)') " ====== Plasma frequency :", wpe
    end if
 
    
  do step=0, nt
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step   :", step
      write(*,'(a,es12.4)') " ====== simulation time  :", step * dt
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
    
!    np = size(particles, kind=kind_particle)
!    if (root) write(*,*) "====== time  integration"

    np = size(particles, kind=kind_particle)
!    if (periodicity_particles) call periodic_particles(np,particles)
    dointerp = mod( step , diag_interval) .eq. 0 .and. step .ne. 0
    
    call normalize(np, particles)
    
    call hamiltonian(np,particles,pold,real(step, kind = kind_particle)*dt)

   
!    call write_particles_ascii(step, particles)
    
    call beam_rnv(tnp,particles,real(step, kind = kind_particle)*dt)
    call densities(np,particles,real(step, kind = kind_particle)*dt)

    dorestart = (mod( step , restart_step) .eq. 0).and.(step.gt.0)

    if (dorestart)   call write_restart(particles)
    
    if ( (mod( step , diag_interval) .eq. 0) )   then 
!!    if ( ( step .eq.0 ) .or.( ( step .gt. 250).and.(mod( step , 10) .eq. 0) ) ) then
        call compute_field(pepc_pars, field_grid, particles)
        
        if (force_law .eq. 3) then
            
            if (root)  call write_field_on_grid3D(field_grid, step, step*dt)
            
        else if (force_law .eq. 2) then
            
            call write_field_on_grid(pepc_pars%pepc_comm, step, field_grid)
            
        end if
        
        call write_particles_ascii(step, particles)
!        
        
        call write_particles_vtk(particles, step, dt*step)
    end if

    !    call copy_particle(pold,poldold,np)
    call copy_particle(particles,pold,np)
    call march(np,dt,particles,ischeme,adv)

    

    timer(4) = get_time()
    if(root) write(*,'(a,es12.4)') " == time in step [s]                              : ", timer(4) - timer(1)

  end do


  deallocate(particles)
  if (allocated(pold)) deallocate(pold)
  if (allocated(poldold)) deallocate(poldold)

  timer(5) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time [s]: ",timer(5) - timer(2)
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

