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

!>
!> Defines a derived type that describes a communication environment and
!> associated procedures.
!>
module module_comm_env
  use module_pepc_types
  implicit none
  private

  !>
  !> Derived type representing a communication environment.
  !>
  type, public :: t_comm_env
    integer(kind_default) :: comm !< an MPI communicator
    integer(kind_pe) :: size !< size of the communicator
    integer(kind_pe) :: rank !< rank within the communicator
    logical :: first !< whether this process is the first process
    logical :: last !< whether this process is the last process
    logical :: mirror !< whether this `t_comm_env`'s `%comm` is a mirror of another
  end type t_comm_env

  interface comm_env_mirror
    module procedure comm_env_mirror
    module procedure comm_env_mirror_from_mpi
  end interface comm_env_mirror

  interface comm_env_dup
    module procedure comm_env_dup
    module procedure comm_env_dup_from_mpi
  end interface comm_env_dup

  public comm_env_mirror
  public comm_env_dup
  public comm_env_destroy
  public comm_env_world
  public comm_env_self

  contains

  !>
  !> Initialize a communication environment `c` from an MPI communicator `comm`.
  !>
  subroutine comm_env_init_from_mpi(comm, c)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: comm !< the MPI communicator used for communication
    type(t_comm_env), intent(out) :: c !< the communication environment to initialize
    integer(kind_default) :: ierr

    c%comm = comm
    call MPI_COMM_SIZE(c%comm, c%size, ierr)
    call MPI_COMM_RANK(c%comm, c%rank, ierr)
    c%first = c%rank == 0
    c%last = c%rank == (c%size - 1)
  end subroutine comm_env_init_from_mpi


  !>
  !> Mirror an MPI communicator `comm`.
  !>
  !> This does not `mpi_comm_dup` `comm`, instead `comm` will be encapsulated.
  !>
  function comm_env_mirror_from_mpi(comm) result(res)
    implicit none
    include 'mpif.h'

    type(t_comm_env) :: res

    integer, intent(in) :: comm !< the MPI communicator used for communication

    call comm_env_init_from_mpi(comm, res)
    res%mirror = .true.
  end function comm_env_mirror_from_mpi


  !>
  !> Mirrors a communication environment `c` to another.
  !>
  !> Both environments share the same communicator.
  !>
  function comm_env_mirror(c) result(res)
    implicit none

    type(t_comm_env) :: res

    type(t_comm_env), intent(in) :: c !< environment to mirror

    res = comm_env_mirror_from_mpi(c%comm)
  end function comm_env_mirror


  !>
  !> Duplicate MPI communicator `comm` and encapsulate.
  !>
  function comm_env_dup_from_mpi(comm) result(res)
    implicit none

    type(t_comm_env) :: res

    integer, intent(in) :: comm

    integer(kind_default) :: cdup, ierr

    call mpi_comm_dup(comm, cdup, ierr)
    call comm_env_init_from_mpi(cdup, res)
    res%mirror = .false.
  end function comm_env_dup_from_mpi


  !>
  !> Duplicate an existing communication environment.
  !>
  !> The underlying MPI communicator is duplicated via `mpi_comm_dup`.
  !>
  function comm_env_dup(c) result(res)
    implicit none

    type(t_comm_env) :: res

    type(t_comm_env), intent(in) :: c !< environment to duplicate

    res = comm_env_dup_from_mpi(c%comm)
  end function comm_env_dup


  !>
  !> Destroy a communication environment.
  !>
  !> Calls `mpi_comm_free` on the MPI communicator contained in `c`.
  !>
  subroutine comm_env_destroy(c)
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_comm_env), intent(inout) :: c !< environment to destroy
    integer(kind_default) :: ierr

    DEBUG_ASSERT(c%comm /= MPI_COMM_NULL)

    if (.not. c%mirror) then
      call mpi_comm_free(c%comm, ierr)
    else
      c%comm = MPI_COMM_NULL
    end if

    c%size = 0
    c%rank = 0
    c%first = .false.
    c%last = .false.
    c%mirror = .false.
  end subroutine comm_env_destroy


  !>
  !> Returns a `t_comm_env` that is a mirror of `MPI_COMM_WORLD`.
  !>
  function comm_env_world()
    implicit none
    include 'mpif.h'

    type(t_comm_env) :: comm_env_world

    comm_env_world = comm_env_mirror(MPI_COMM_WORLD)
  end function


  !>
  !> Returns a `t_comm_env` that is a mirror of `MPI_COMM_SELF`.
  !>
  function comm_env_self()
    implicit none
    include 'mpif.h'

    type(t_comm_env) :: comm_env_self

    comm_env_self = comm_env_mirror(MPI_COMM_SELF)
  end function

end module module_comm_env
