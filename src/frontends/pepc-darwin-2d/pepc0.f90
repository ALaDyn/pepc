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
  use module_timings
  use module_debug
  use newton_krylov

  ! frontend helper routines
  use helper
  implicit none


  ! timing variables

  real*8 :: timer(5)
  real*8 :: t1, t2,errtol
  integer*8 :: itsub,iter,rc,ierr,i
  type(t_particle), allocatable :: particles0(:),particles1(:)

  ! control variable
  logical :: doDiag
  itsub = 100
  errtol = 1.0d-8
  !allocate(particles0(10), stat=rc )
  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-darwin-2d", my_rank, n_ranks, .true.)

  root = my_rank.eq.0

  call set_parameter()
  call init_particles_square(particles)


  call pepc_particleresults_clear(particles)

  call pepc_grow_tree(particles)

  call pepc_traverse_tree(particles)
  call pepc_timber_tree()

!  call linear(particles(1),particles0(1))
!  call gmres( particles(1), linear, particles(2), &
!         errtol,itsub,particles0(3),iter)

  call nsolgm(np,particles,linear,errtol, errtol,itsub,itsub,errtol, particles0)
!  call dirder(np,particles,particles,linear,particles,particles1)
!  deallocate(particles)
  !deallocate(particles0)
  !deallocate(particles1)
  do i = 1,np
    write(*,*) particles0(i)%x
    write(*,*) particles0(i)%data%v
    write(*,*)
 enddo


  call pepc_finalize()

end program pepc

