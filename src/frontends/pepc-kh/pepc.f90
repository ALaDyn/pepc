! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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
  
   use encap
   use pepc_helper
   use time_helper
   use field_helper

   use mpi

   implicit none

   type(pepc_comm_t) :: pepc_comm
   type(pepc_nml_t)  :: pepc_nml
   type(pepc_pars_t) :: pepc_pars
   type(time_pars_t) :: time_pars
   type(t_particle), dimension(:), allocatable :: p

   integer :: mpi_err, MPI_COMM_SPACE
   integer :: step
   real(kind=8) :: timer(5)

   logical :: root, do_pdump, do_fdump

   ! global MPI initialization
   call MPI_Init( mpi_err )
   call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_SPACE)

   timer(1) = get_time()

   ! initialize pepc library and MPI, setup particle data
   call init_pepc(pepc_comm, pepc_nml, MPI_COMM_SPACE)
   call pepc_setup(p, pepc_pars, pepc_comm, pepc_nml)
   call setup_time(time_parspepc_pars%pepc_comm)

   root = pepc_pars%pepc_comm%mpi_rank.eq.0

   timer(2) = get_time()

   do step = 0, time_pars%nsteps
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step  :", step
      write(*,'(a,es12.4)') " ====== simulation time :", step * time_pars%dt
    end if
    
    timer(3) = get_time()
    
    do_pdump = MOD(step, pepc_pars%pdump) .eq. 0
    do_fdump = MOD(step, pepc_pars%fdump) .eq. 0
    
    call pepc_particleresults_clear(p, pepc_pars%npp)

    call pepc_grow_tree(pepc_pars%npp, pepc_pars%np, p)
    call pepc_traverse_tree(pepc_pars%npp, p)

    call apply_external_field()

    if(do_fdump) call compute_field()
    
    call pepc_timber_tree()
    call pepc_restore_particles(pepc_pars%npp, p)

    if(do_fdump) call write_particles(particles)
        
    call push_particles(p, pepc_pars)

    timer(4) = get_time()
    if(root) write(*,'(a,es12.4)') " == time in step [s]                              : ", timer(4)-timer(3)

    call timings_GatherAndOutput(step, 0)
    
  end do 
 
  deallocate(particles)

  timer(5) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time without setup [s]: ", timer(5) - timer(2)
    write(*,'(a,es12.4)') " ===== total run time wit setup [s]:     ", timer(5) - timer(1)
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

