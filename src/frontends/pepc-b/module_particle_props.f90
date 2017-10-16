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


module module_particle_props

  type particle
     real*8 :: x    ! coords
     real*8 :: y
     real*8 :: z
     real*8 :: ux    ! momenta
     real*8 :: uy
     real*8 :: uz 
     real*8 :: q     ! charge
     real*8 :: m     ! mass
     real*8 :: work  ! work load from force sum
     integer :: label    ! label
  end type particle
  
  integer :: mpi_type_particle


! particle arrays
  real*8, allocatable, public :: x(:),  y(:),  z(:), &     ! position
                      ux(:), uy(:), uz(:), &     ! velocity
                              q(:),  m(:), &     ! charge and mass
			Ex(:), Ey(:), Ez(:), &   ! E-field
			Bx(:), By(:), Bz(:), &   ! B-field
			pot(:), &	         ! scalar potential
			Ax(:), Ay(:), Az(:), &   ! vector potential
			work(:)	         ! work load (interaction list length)

  integer, allocatable, public ::  pelabel(:), pepid(:)  ! particle label

  public setup_particle_arrays

contains

  subroutine setup_particle_arrays(nppm)

  use module_physvars
  implicit none
  integer, intent(in) :: nppm   

  allocate( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), q(nppm), m(nppm), &
	Ex(nppm), Ey(nppm), Ez(nppm), Bx(nppm), By(nppm), Bz(nppm), &
	pot(nppm), Ax(nppm), Ay(nppm), Az(nppm), work(nppm) )

  allocate( pelabel(nppm), pepid(nppm) )

  end subroutine setup_particle_arrays

 
end module module_particle_props

