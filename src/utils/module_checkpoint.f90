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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Helper functions for checkpointing and restarting purposes
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_checkpoint
      implicit none
      include 'mpif.h'

      private

      integer, parameter :: filehandle = 91
      character(12), parameter :: directory = './particles'
      integer(KIND=MPI_OFFSET_KIND) :: max_header_size = 128

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public write_particles_ascii
      public write_particles_binary
      public write_particles_mpiio
      public read_particles_mpiio

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_ascii(my_rank, itime, np_local, dp, filename)
            use module_pepc_types
            use module_utils
            use module_pepc, only : pepc_write_parameters
            use module_debug, only : pepc_status
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            type(t_particle), intent(in), dimension(np_local) :: dp
            character(*), intent(out) :: filename
            logical :: firstcall  = .true.
            character(50) :: dir
            integer :: i

            dir = trim(directory)//"/ascii/"
            write(filename,'(a,"particle_",i6.6,"_",i6.6,".dat")') trim(dir), itime, my_rank
            call pepc_status("DUMP PARTICLES ASCII: "//filename)

            if (firstcall) then
              call create_directory(trim(directory))
              call create_directory(trim(dir))
              firstcall = .false.
            endif

            open(filehandle, file=trim(filename), STATUS='REPLACE')
            do i=1, np_local
              write(filehandle,*) dp(i)
            end do
            close(filehandle)

            filename = trim(filename)//".h"
            if (my_rank == 0) then
              open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
              call pepc_write_parameters(filehandle)
              close(filehandle)
            endif
          end subroutine


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_binary(my_rank, itime, np_local, dp, filename)
            use module_pepc_types
            use module_utils
            use module_pepc, only : pepc_write_parameters
            use module_debug, only : pepc_status
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            type(t_particle), intent(in), dimension(np_local) :: dp
            character(*), intent(out) :: filename
            logical :: firstcall  = .true.

            character(50) :: dir

            dir = trim(directory)//"/binary/"
            write(filename,'(a,"particle_",i6.6,"_",i6.6,".bin")') trim(dir), itime, my_rank
            call pepc_status("DUMP PARTICLES BINARY: "//filename)

            if (firstcall) then
              call create_directory(trim(directory))
              call create_directory(trim(dir))
              firstcall = .false.
            endif

            open(filehandle, file=trim(filename), STATUS='REPLACE', ACCESS="STREAM")
            write(filehandle) dp(1:np_local)
            close(filehandle)

            filename = trim(filename)//".h"
            if (my_rank == 0) then
              open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
              call pepc_write_parameters(filehandle)
              close(filehandle)
            endif
          end subroutine


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_mpiio(comm, my_rank, itime, np_local, n_total, dp, filename)
            use module_pepc_types
            use module_debug
            use module_utils
            use module_pepc, only : pepc_write_parameters
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            integer*8, intent(in) :: n_total
            integer, intent(in) :: comm
            type(t_particle), intent(in) :: dp(np_local)
            character(*), intent(out) :: filename
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer(KIND=MPI_OFFSET_KIND) :: disp
            logical :: firstcall  = .true.
            character(50) :: dir

            dir = trim(directory)//"/mpi/"
            write(filename,'(a,"particle_",i6.6,".mpi")') trim(dir), itime
            call pepc_status("DUMP PARTICLES MPI: "//filename)

            if (firstcall) then
              call create_directory(trim(directory))
              call create_directory(trim(dir))
              firstcall = .false.
            endif

            call MPI_FILE_OPEN(comm,filename,IOR(MPI_MODE_RDWR,MPI_MODE_CREATE),MPI_INFO_NULL,fh,ierr)

            if (ierr .ne. MPI_SUCCESS) then
              DEBUG_ERROR(*, 'write_particles_mpiio(): file open failed', my_rank, ierr, filename)
            end if

            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            if (my_rank == 0) then
              call MPI_FILE_WRITE(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
              call MPI_FILE_WRITE(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)
              call MPI_FILE_GET_POSITION(fh, disp, ierr)

              if (disp > max_header_size) then
                DEBUG_ERROR(*, "header_size is too small: ", max_header_size, "<", disp)
              end if
            end if

            ! Redefine file view, now with our custom type
            call MPI_FILE_SET_VIEW(fh, max_header_size, MPI_TYPE_particle, MPI_TYPE_particle, 'native', MPI_INFO_NULL, ierr)
            ! Write particle data
            call MPI_FILE_WRITE_ORDERED(fh, dp(:), np_local, MPI_TYPE_particle, status, ierr)
            ! Take care before closing
            call MPI_FILE_SYNC(fh,ierr)
            call MPI_FILE_CLOSE(fh,ierr)

            filename = trim(filename)//".h"
            if (my_rank == 0) then
              open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
              call pepc_write_parameters(filehandle)
              close(filehandle)
            endif

          end subroutine



          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine read_particles_mpiio(itime_in, comm, my_rank, n_cpu, itime, np_local, n_total, dp, filename)
            use module_pepc_types
            use module_debug
            use module_pepc, only : pepc_read_parameters
            implicit none
            integer, intent(in) :: my_rank, n_cpu, itime_in
            integer, intent(out) :: itime, np_local
            integer*8, intent(out) :: n_total
            integer, intent(in) :: comm
            character(*), intent(out) :: filename
            type(t_particle), allocatable :: dp(:)
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer*8 :: remain

            character(50) :: dir

            dir = trim(directory)//"/mpi/"
            write(filename,'(a,"particle_",i6.6,".mpi")') trim(dir), itime_in
            call pepc_status("READ PARTICLES MPI: "//filename)

            call MPI_FILE_OPEN(comm,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

            if (ierr .ne. MPI_SUCCESS) then
              DEBUG_ERROR(*,'read_particles_mpiio(): file open failed', my_rank, ierr, filename)
            end if

            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_READ(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
            call MPI_FILE_READ(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)

            np_local = int(n_total/n_cpu, kind(np_local))
            remain = n_total-np_local*n_cpu
            if ((remain > 0) .and. (my_rank < remain)) np_local = np_local+1

            if (allocated(dp)) deallocate(dp)
            allocate(dp(1:np_local))

            ! Redefine file view, now with our custom type
            call MPI_FILE_SET_VIEW(fh, max_header_size, MPI_TYPE_particle, MPI_TYPE_particle, 'native', MPI_INFO_NULL, ierr)
            ! Read particle data
            call MPI_FILE_READ_ORDERED(fh, dp, np_local, MPI_TYPE_particle, status, ierr)
            ! Close file
            call MPI_FILE_CLOSE(fh,ierr)

            filename = trim(filename)//".h"
            open(filehandle, file=trim(filename),action='read')
            call pepc_read_parameters(filehandle)
            close(filehandle)

          end subroutine


end module module_checkpoint
