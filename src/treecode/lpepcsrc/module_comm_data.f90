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

module module_comm_data
  implicit none
  private

  type, public :: t_comm_data
    integer :: comm !< an MPI communicator
    integer :: size !< size of the communicator
    integer :: rank !< rank within the communicator
    logical :: first !< whether this process is the first process
    logical :: last !< whether this process is the last process
  end type t_comm_data

  public comm_data_create
  public comm_data_dup
  public comm_data_destroy

  contains

  !>
  !> initialize a t_comm_data
  !>
  subroutine comm_data_create(c, comm)
    use mpi
    implicit none

    type(t_comm_data), intent(out) :: c
    integer, intent(in) :: comm
    integer :: ierr

    call mpi_comm_dup(comm, c%comm, ierr)
    call mpi_comm_size(c%comm, c%size, ierr)
    call mpi_comm_rank(c%comm, c%rank, ierr)
    c%first = c%rank == 0
    c%last = c%rank == (c%size - 1)
  end subroutine comm_data_create


  !>
  !> duplicate a t_comm_data
  !>
  subroutine comm_data_dup(c1, c2)
    implicit none

    type(t_comm_data), intent(in) :: c1
    type(t_comm_data), intent(out) :: c2

    call comm_data_create(c2, c1%comm)
  end subroutine comm_data_dup


  !>
  !> destroy t_comm_data
  !>
  subroutine comm_data_destroy(c)
    use mpi
    implicit none

    type(t_comm_data), intent(inout) :: c
    integer :: ierr
    
    call mpi_comm_free(c%comm, ierr)
    c%size = 0
    c%rank = 0
    c%first = .false.
    c%last = .false.
  end subroutine comm_data_destroy

end module module_comm_data
