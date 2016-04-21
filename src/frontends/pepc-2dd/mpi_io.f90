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


!!!!!!!!!!!!!!   MPI Input/Output Routines      !!!!!!!!!!!!!!!

module module_mpi_io
  use module_pepc_kinds
  use module_pepc_types
  implicit none
  include 'mpif.h'
  save

  public read_ascii
  public write_ascii

  contains

    subroutine read_ascii(np,p,filename,pepc_comm)
    use encap
    implicit none

    type(pepc_comm_t)                        , intent(in)  :: pepc_comm
    integer(kind_particle)                   , intent(in)  :: np
    character(255)                           , intent(in)  :: filename
    real(kind_particle),allocatable          , intent(out) :: p(:)
!    type(t_particle),allocatable             , intent(out) :: p(:)
    integer(kind_particle)                                 :: rc,ierr,fh,mpi_rank
    integer(kind_particle)                                 :: status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND)                          :: offset, empty

    if (allocated(p)) deallocate(p)
    allocate(p(np), stat = rc)


    mpi_rank    = int(pepc_comm%mpi_rank, kind = kind(mpi_rank))
    empty       = 0
    fh          = 0
    offset      = np*mpi_rank
    write(*,*) "my rank is", mpi_rank

    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FILE_SET_VIEW(fh, 0, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)
    call MPI_FILE_READ_AT(fh, offset, p, np, MPI_REAL4, status, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FILE_CLOSE(fh, ierr)


    end subroutine read_ascii


    subroutine write_ascii(np,p,filename,pepc_comm)
    use encap
    implicit none

    type(pepc_comm_t)                        , intent(in)  :: pepc_comm
    integer(kind_particle)                   , intent(in)  :: np
    character(*)                             , intent(in)  :: filename
    real(kind_particle)                      , intent(in)  :: p(np)
!    real(kind_particle),allocatable          , intent(in)  :: p(:)
!    type(t_particle),allocatable             , intent(out) :: p(:)
    integer(kind_particle)                                 :: rc,ierr,fh,mpi_rank
    integer(kind_particle)                                 :: status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND)                          :: offset, empty


    mpi_rank    = int(pepc_comm%mpi_rank, kind = kind(mpi_rank))
    empty       = 0
    fh          = 0
    offset      = np!*mpi_rank
!    write(*,*) "my rank is", mpi_rank

!    if (mpi_rank.eq.0)  call MPI_File_delete(filename, MPI_INFO_NULL, ierr)


    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FILE_SET_VIEW(fh, empty, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE_AT(fh, offset, p, np, MPI_DOUBLE_PRECISION, status, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FILE_CLOSE(fh, ierr)


    end subroutine write_ascii


end module module_mpi_io
