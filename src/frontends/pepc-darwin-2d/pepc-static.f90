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

program pepc

   ! pepc modules
   use module_pepc
   use module_pepc_types
   use module_timings
   use module_debug
   use module_globals, only: adv
   use newton_krylov
   use module_integration
   use module_tool, only: copy_particle, gyrofrequency, integrate
   use zufall
   use encap
   use module_diagnostic
   use field_helper, only: compute_field, setup_field_grid, pepc_setup, prepare_grid, field_grid_update_grid, write_field_on_grid

   ! frontend helper routines
   use helper
!  use module_mpi_io
   use module_shortcut
   use mpi
   implicit none
   ! timing variables
   real(kind_particle)                        :: timer(5)
   real(kind_particle)                        :: t1, t2, t, tau, I
!  type(t_particle), allocatable              :: q(:),r(:)
!  type(physics_pars_t)                       :: physics_pars
   type(field_grid_t)                         :: field_grid
   type(time_pars_t)                          :: time_pars
   type(pepc_pars_t)                          :: pepc_pars
   type(pepc_comm_t)                          :: pepc_comm
   real(kind_particle)                        :: new_extent(3), new_offset(3), vxmin, vymin, vzmin, vxmax, vymax, vzmax
   integer(kind_particle)                     :: ip, jp, rc = 1

   ! control variable
   logical :: dogrid, dorestart, dointerp, explicit = .false.

  !!! initialize pepc library and MPI

   call pepc_initialize("pepc-darwin-2d", my_rank, n_ranks, .true., comm=communicator)
!  call pepc_initialize("pepc-darwin-2d", my_rank, n_ranks, .true.)

   root = my_rank .eq. 0
   timer(1) = get_time()
   call set_parameter()
   np = size(particles, kind=kind_particle)

   call pepc_setup(pepc_pars)
   call setup_field_grid(field_grid, pepc_pars%pepc_comm)
   call init_particles(particles, field_grid)

   timer(2) = get_time()
   t = 0

   if (root) then
      write (*, '(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)
      write (*, '(a,es12.4)') " ====== Plasma frequency :", wpe
   end if

   timer(3) = get_time()

   call pepc_particleresults_clear(particles)
   timer(1) = get_time()
   t1 = get_time()

   call pepc_grow_tree(particles)

   t2 = get_time()
   if (root) write (*, '(a,es12.4)') " ====== tree grow time  :", t2 - t1
   t1 = get_time()
   call pepc_traverse_tree(particles)
   t2 = get_time()
   if (root) write (*, '(a,es12.4)') " ====== tree walk time  :", t2 - t1
   call pepc_restore_particles(particles)
   call pepc_timber_tree()

   timer(5) = get_time()
   t = t + timer(5) - timer(1)

   np = size(particles, kind=kind_particle)

   call normalize(np, particles)

   call time_step(particles)

   call compute_field(pepc_pars, field_grid, particles)
!    call write_field_on_grid(pepc_comm, 0, field_grid)
!    call write_particles(particles)
!    call write_particles_vtk(particles, 0, 0.0_8)
   call write_particles_ascii(0, particles)
   call write_field_on_grid_ascii(field_grid, 0)
!    call hamiltonian_weibel(np,particles,particles,real(0, kind = kind_particle)*dt)
!    call write_particles_ascii( int(step, kind = kind_particle) , particles )

!    call copy_particle(particles,pold,np)
!    call march(np,dt,particles,ischeme,adv)
!    call hamiltonian_weibel(np,particles,particles,real(1, kind = kind_particle)*dt)
!--------------------------------------------------------------------------------------------
!    Is F a contraction?
!    write(*,*) "Distances between particles electrons"
!    jp   = 1
!    ip   = 33333
!    write(*,*)  "labels =,",pold(ip)%label,pold(jp)%label
!    write(*,*)  "|| x - y || =,", sqrt( sum( ( pold(ip)%x(1:2) - ( pold(jp)%x(1:2) ) )**2 ) &
!                + sum( ( pold(ip)%data%v(1:3) - ( pold(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*)  "T(|| x - y ||) =,", sqrt( sum( ( particles(ip)%x(1:2) - ( particles(jp)%x(1:2) ) )**2 ) &
!                + sum( ( particles(ip)%data%v(1:3) - ( particles(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*) "Distances between particles protons"
!    jp   = 22002
!    ip   = 404
!    write(*,*)  "labels =,",pold(ip)%label,pold(jp)%label
!    write(*,*)  "|| x - y || =", sqrt( sum( ( pold(ip)%x(1:2) - ( pold(jp)%x(1:2) ) )**2 ) &
!                + sum( ( pold(ip)%data%v(1:3) - ( pold(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*)  "T(|| x - y ||) =", sqrt( sum( ( particles(ip)%x(1:2) - ( particles(jp)%x(1:2) ) )**2 ) &
!                + sum( ( particles(ip)%data%v(1:3) - ( particles(jp)%data%v(1:3) ) )**2 ) )
!
!    write(*,*) "Distances between particles electrons beam"
!    jp   = 10
!    ip   = 8000
!    write(*,*)  "labels =,",pold(ip)%label,pold(jp)%label
!    write(*,*)  "|| x - y || =", sqrt( sum( ( pold(ip)%x(1:2) - ( pold(jp)%x(1:2) ) )**2 ) &
!                + sum( ( pold(ip)%data%v(1:3) - ( pold(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*)  "T(|| x - y ||) =", sqrt( sum( ( particles(ip)%x(1:2) - ( particles(jp)%x(1:2) ) )**2 ) &
!                + sum( ( particles(ip)%data%v(1:3) - ( particles(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*) "Distances between different type of electron-beam"
!    jp   = 10000
!    ip   = 78004
!    write(*,*)  "labels =",pold(ip)%label,pold(jp)%label
!    write(*,*)  "|| x - y || =", sqrt( sum( ( pold(ip)%x(1:2) - ( pold(jp)%x(1:2) ) )**2 ) &
!                + sum( ( pold(ip)%data%v(1:3) - ( pold(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*)  "T(|| x - y ||) =", sqrt( sum( ( particles(ip)%x(1:2) - ( particles(jp)%x(1:2) ) )**2 ) &
!                + sum( ( particles(ip)%data%v(1:3) - ( particles(jp)%data%v(1:3) ) )**2 ) )
!
!
!    write(*,*) "Distances between different type of proton-beam"
!    jp   = 10000
!    ip   = 41
!    write(*,*)  "labels =,",pold(ip)%label,pold(jp)%label
!    write(*,*)  "|| x - y || =,", sqrt( sum( ( pold(ip)%x(1:2) - ( pold(jp)%x(1:2) ) )**2 ) &
!                + sum( ( pold(ip)%data%v(1:3) - ( pold(jp)%data%v(1:3) ) ) ) )
!
!
!    write(*,*)  "T(|| x - y ||) =,", sqrt( sum( ( particles(ip)%x(1:2) - ( particles(jp)%x(1:2) ) )**2 ) &
!                + sum( ( particles(ip)%data%v(1:3) - ( particles(jp)%data%v(1:3) ) )**2 ) )

!--------------------------------------------------------------------------------------------
   timer(4) = get_time()
   if (root) write (*, '(a,es12.4)') " == time in step [s]                              : ", timer(4) - timer(1)

   deallocate (particles)
   if (allocated(pold)) deallocate (pold)
   if (allocated(poldold)) deallocate (poldold)

   timer(5) = get_time()

   if (root) then
      write (*, *) " "
      write (*, '(a)') " ===== finished pepc simulation"
      write (*, '(a,es12.4)') " ===== total run time [s]: ", timer(5) - timer(2)
   end if

  !!! cleanup pepc and MPI
   call pepc_finalize()

end program pepc

