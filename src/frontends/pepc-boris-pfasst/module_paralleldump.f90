! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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
!> module for asynchronous dumping into fort.XXX files
!>
!> upon initialization, a copy of MPI_COMM_WORLD is created
!>
!> if paralleldump_dump() is called, the given character string is sent
!> to the master node (rank==0). the message tag contains the unit number to dump to
!>
!> actual output is performed upon paralleldump_flush() only by rank 0
!> if some other rank calls that function, nothing will happen
!>
!> in this implementation, MPI-internal buffers are used for buffering the output
!> lets hope they do not overflow :-)
module pepcboris_paralleldump
  use module_pepc_kinds
  implicit none
  include 'mpif.h'
  private

  public paralleldump_init
  public paralleldump_cleanup
  public paralleldump_dump
  public paralleldump_flush

  #undef PARALLELDUMP_OWN_BSEND_BUFFER

  integer, public, parameter :: PARALLELDUMP_MAXLEN = 1024 ! maximum length of a single line sent with this method
  #ifdef PARALLELDUMP_OWN_BSEND_BUFFER
    integer, parameter :: PARALLELDUMP_BUFFLEN = 64 ! number of lines in the used send-buffer
    character, dimension(:), allocatable :: paralleldump_bsend_buffer
  #endif

  integer(kind_default) :: paralleldump_comm = MPI_COMM_NULL
  integer(kind_pe) :: paralleldump_rank = -1
  integer(kind_pe), parameter :: PARALLELDUMP_ROOT = 0

  logical, public :: paralleldump_autoflush = .true. !< if set to .true. every call to paralleldump_dump() automatically also calls paralleldump_flush()

  contains


  subroutine paralleldump_dump(istream, message)
    implicit none
    integer(kind_default), intent(in) :: istream
    character(len=PARALLELDUMP_MAXLEN), intent(in) :: message
    integer(kind_default) :: ierr, msg_len
    msg_len = len(trim(message))
    call MPI_BSEND(message, msg_len, MPI_CHARACTER, PARALLELDUMP_ROOT, istream, paralleldump_comm, ierr)
    if (paralleldump_autoflush) call paralleldump_flush()
  end subroutine


  subroutine paralleldump_flush()
    implicit none
    logical :: msg_avail
    integer(kind_default) :: stat(MPI_STATUS_SIZE), ierr
    character(len=PARALLELDUMP_MAXLEN) :: message
    integer(kind_pe) :: sender
    integer(kind_default) :: msg_tag, msg_len

    if (paralleldump_rank == PARALLELDUMP_ROOT) then
      do
        call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, paralleldump_comm, msg_avail, stat, ierr)

        if (.not. msg_avail) exit

        sender  = stat(MPI_SOURCE)
        msg_tag = stat(MPI_TAG)
        call MPI_GET_COUNT(stat, MPI_CHARACTER, msg_len, ierr)
        call MPI_RECV(message, msg_len, MPI_CHARACTER, sender, msg_tag, paralleldump_comm, MPI_STATUS_IGNORE, ierr)

        write(msg_tag, *) message(1:msg_len)
      end do
    endif

  end subroutine


  subroutine paralleldump_init()
    implicit none
    integer(kind_default) :: ierr
    call MPI_COMM_DUP(MPI_COMM_WORLD, paralleldump_comm, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, paralleldump_rank, ierr)
    #ifdef PARALLELDUMP_OWN_BSEND_BUFFER
      allocate(paralleldump_bsend_buffer(PARALLELDUMP_MAXLEN*PARALLELDUMP_BUFFLEN))
      call MPI_BUFFER_ATTACH(paralleldump_bsend_buffer, PARALLELDUMP_MAXLEN*PARALLELDUMP_BUFFLEN, ierr)
    #endif
  end subroutine


  subroutine paralleldump_cleanup()
    implicit none
    integer(kind_default) :: ierr
    #ifdef PARALLELDUMP_OWN_BSEND_BUFFER
      integer(kind_default) :: buffsize
    #endif
    ! make sure that nobody wants to send any messages
    call MPI_BARRIER(paralleldump_comm, ierr)
    ! flush anz pending messages
    call paralleldump_flush()
    #ifdef PARALLELDUMP_OWN_BSEND_BUFFER
      call MPI_BUFFER_DETACH(paralleldump_bsend_buffer, buffsize, ierr)
      deallocate(paralleldump_bsend_buffer)
    #endif
    call MPI_COMM_FREE(paralleldump_comm, ierr)
  end subroutine



end module
