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


! ==============================================================
!
!
!                  PEPC-V
!
!    Parallel Efficient Parallel Coulomb-solver: Vortex particles
!
!  ==============================================================

program pepcv

  use physvars
  use manipulate_particles
  use module_pepc
  use module_timings
  use files
  use diagnostics
  use module_interaction_specific, only: theta2
  use pfasst_helper_module
  use pfasst_calc_module
  use pfasst_run_module
  use module_pepc_types
  implicit none

  integer :: i, ierr
  real :: trun                     ! total run time including restarts and offset
  integer :: itime, stage, t_flag
  integer, parameter :: t_remesh = t_userdefined_first + 1

  ! Allocate array space for tree

  call init_communication()

  call pepc_initialize("pepc-v" ,my_rank_space, n_cpu_space, .false., 0, comm=MPI_COMM_SPACE)
  call pepc_read_parameters_from_first_argument()

  ! Set up O/P files
  call openfiles

  ! Each CPU gets copy of initial data
  call pepc_setup(itime, trun)

  ! Set up particles
  call special_start()

  call pepc_prepare(3_kind_dim)

  call init_pfasst(vortex_particles(1:np),np)

  call dump(0,ts)

  if (parallel == 1) then
    call run_parallel(y0, 1.0D00*dt, 1.0D00*te)
  else
    !call run_rk3(y0, 1.0D00*dt, nt)
    !call run_rk4(y0, 1.0D00*dt, nt)
    call run_serial(y0, 1.0D00*dt, nt)
  end if

  !call remeshing()
  !if (my_rank_time == 0) write(*,*) n,'particles now in use'

  if (my_rank_time == n_cpu_time-1) call dump_results()

  call MPI_BARRIER(MPI_COMM_TIME,ierr)
  call MPI_BARRIER(MPI_COMM_SPACE,ierr)

  ! Loop over all timesteps
!  do while (itime < nt)
!
!     itime = itime + 1
!     !if (my_rank==0 ) write(*,'(a5,i8,a3,i8,a7,i8,a)') 'Step',itime,' of',nt,', using',n,' particles -------------'
!
!     ! Runge-Kutta time integration
!     do stage = 1,rk_stages
!
!        !if (my_rank == 0) write(*,*) 'Time:',trun,'/',te
!
!        call timer_start(t_tot)
!
!        call pepc_particleresults_clear(vortex_particles, np)
!
!        call pepc_grow_and_traverse(np, n, vortex_particles, itime, .false., .false.)
!
!        do i=1,np
!          vortex_particles(i)%results%u( 1:3) = vortex_particles(i)%results%u( 1:3) * force_const
!          vortex_particles(i)%results%af(1:3) = vortex_particles(i)%results%af(1:3) * force_const
!          vortex_particles(i)%results%div     = vortex_particles(i)%results%div * force_const
!        end do
!
!        call push_rk2(stage)
!
!        if (stage .lt. rk_stages) then
!            call timer_stop(t_tot)   ! total loop time without diags
!            call timings_LocalOutput(itime,stage)
!            call timings_GatherAndOutput(itime,stage)
!        end if
!
!        flush(6)
!
!        trun  = trun  + dt/rk_stages
!
!     end do
!
!     ! dump, if needed (need to do this before remeshing or we will loose velocity information)
!     call dump(itime, trun)
!
!     ! if remeshing is requested at all and if it is time right now, do it!
!     if ((rem_freq .gt. 0) .and. (mod(itime,rem_freq)==0)) then
!
!        !if (my_rank==0) write(*,'("PEPC-V | ", a)') 'Starting remeshing...'
!        call timer_start(t_remesh)
!
!        call remeshing()
!
!        call timer_stop(t_remesh)
!        !if (my_rank==0) write(*,'("PEPC-V | ", a,f12.8,a)') 'Finished remeshing after ',timer_read(t_remesh),' seconds'
!        t_flag = -rk_stages
!
!     else
!
!        t_flag = rk_stages
!
!     end if
!
!     call timer_stop(t_tot)   ! total loop time incl. remeshing if requested
!     call timings_LocalOutput(itime,t_flag)
!     call timings_GatherAndOutput(itime,t_flag)
!
!     ! Some linear diagnostics
!     call linear_diagnostics(itime,trun)
!     call divergence_diag(itime,trun)
!
!  end do

  ! deallocate array space for particles
  call cleanup()
  
  call finish_pfasst()

  ! Tidy up O/P files
  call closefiles

  ! cleanup of lpepc static data
  call pepc_finalize()

end program pepcv
