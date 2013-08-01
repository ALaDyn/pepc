! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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
  ! frontend helper routines
  use pepca_helper
  use pepca_integrator
  use pepca_diagnostics
  use pepca_globals

  use pfasst
  use pf_helper
  
  implicit none
  
  integer :: step
  
  ! variables for pfasst
  integer(kind_default):: MPI_COMM_SPACE, MPI_COMM_TIME
  type(pf_nml_t) :: pf_nml
  

  !integer(kind_particle), parameter :: numparts = 296568
  integer(kind_particle), parameter :: numparts = 2500
  
  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)

  ! Take care of communication stuff, set up PFASST and PMG
  call init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)
  ! initialize pepc library and MPI
  call pepc_initialize("pepc-a-pfasst", my_rank, n_ranks, .false., comm=MPI_COMM_SPACE)

  root = my_rank.eq.0

  call timer_start(t_user_total)
  call timer_start(t_user_init)
  call set_parameter()
  call read_particles(particles, 'E_phase_space.dat', numparts, 'I_phase_space.dat', numparts)
  call timer_stop(t_user_init)

  if(root) write(*,'(a,f12.4," s")') ' === init time [s]: ', timer_read(t_user_init)
 
  do step=0, nt - 1
    if(root) then
      write(*,*) " "
      write(*,'(a,i12,"/",i0)')   ' ====== computing step       :', step, nt-1
      write(*,'(a,f12.4)') ' ====== simulation time (fs) :', step*dt*unit_time_fs_per_simunit
    end if
    
    call timer_start(t_user_step)
    
    call pepc_particleresults_clear(particles)
    call pepc_grow_tree(particles)
    if(root) write(*,'(a,f12.2," s")') ' ====== tree grow time   :', timer_read(t_fields_tree)
    call pepc_traverse_tree(particles)
    if(root) write(*,'(a,f12.2," s")') ' ====== tree walk time   :', timer_read(t_fields_passes)

    if (dbg(DBG_STATS)) call pepc_statistics(step)
    call pepc_timber_tree()
    
    if (step > 0) then
      ! first half step to synchronize velocities with positions
      call update_velocities(particles, dt/2.)
    end if
    
    ! do diagnostics etc here
    call timer_start(t_user_diag)
    if ((particle_output_interval>0) .and. ((mod(step, particle_output_interval)==0) .or. (step==nt-1))) then
      call write_particles_vtk(particles, step, dt*step*unit_time_fs_per_simunit)
      call write_particles_ascii(my_rank, step, particles)
      call gather_and_write_densities(particles, step, dt*step*unit_time_fs_per_simunit)
    endif
    if ((domain_output_interval  >0) .and. ((mod(step, domain_output_interval)  ==0) .or. (step==nt-1))) call write_domain(   particles, step, dt*step*unit_time_fs_per_simunit)
    call diagnose_energy(particles, step, dt*step*unit_time_fs_per_simunit)
    call timer_stop(t_user_diag)
    if(root) write(*,'(a,f12.2," s")') ' ====== diagnostics time :', timer_read(t_user_diag)
    
    ! second half step: velocities are one half step ahead again
    call update_velocities(particles, dt/2.)
    ! full step for positions: now positions are one half step ahead
    call push_particles(particles, dt)    
    
    call timer_stop(t_user_step)
    if(root) write(*,'(a,f12.2," s")') ' == time in step: ', timer_read(t_user_step)

    call timings_GatherAndOutput(step, 0, 0 == step)
    
  end do 
 
  deallocate(particles)

  call timer_stop(t_user_total)

  if(root) then
    write(*,*)            ' '
    write(*,'(a)')        ' ===== finished pepc simulation'
    write(*,'(a,es12.4, " s")') ' ===== total run time    : ', timer_read(t_user_total)
  end if

  ! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

