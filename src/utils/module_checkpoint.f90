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

!>
!> Helper functions for checkpointing and restarting purposes
!>
module module_checkpoint
      implicit none
      include 'mpif.h'

      private

      integer, parameter :: filehandle = 91
      character(12), parameter :: directory = './particles'
      integer(KIND=MPI_OFFSET_KIND) :: max_header_size = 128

      public write_particles_ascii
      public write_particles_binary
      public write_particles_mpiio
      public read_particles_mpiio
      public read_particles_mpiio_from_filename

      contains

          subroutine write_particles_ascii(my_rank, itime, dp, filename_out)
            use module_pepc_types
            use module_utils
            use module_pepc, only : pepc_write_parameters
            use module_debug, only : pepc_status
            implicit none
            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in) :: itime
            type(t_particle), intent(in), dimension(:) :: dp
            character(*), optional, intent(out) :: filename_out
            logical :: firstcall  = .true.
            character(50) :: dir
            integer(kind_particle) :: i
            character(100) :: filename

            dir = trim(directory)//"/ascii/"
            write(filename,'(a,"particle_",i6.6,"_",i6.6,".dat")') trim(dir), itime, my_rank
            call pepc_status("DUMP PARTICLES ASCII: "//trim(filename))

            if (firstcall) then
              call create_directory(trim(directory))
              call create_directory(trim(dir))
              firstcall = .false.
            endif

            open(filehandle, file=trim(filename), STATUS='REPLACE')
            do i=1, size(dp,kind=kind(i))
              write(filehandle,*) dp(i)
            end do
            close(filehandle)

            filename = trim(filename)//".h"
            if (my_rank == 0) then
              open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
              call pepc_write_parameters(filehandle)
              close(filehandle)
            endif
            
            if (present(filename_out)) then
              filename_out = trim(filename)
            endif
          end subroutine


          subroutine write_particles_binary(my_rank, itime, dp, filename)
            use module_pepc_types
            use module_utils
            use module_pepc, only : pepc_write_parameters
            use module_debug, only : pepc_status
            implicit none
            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in) :: itime
            type(t_particle), intent(in), dimension(:) :: dp
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
            write(filehandle) dp(1:size(dp))
            close(filehandle)

            filename = trim(filename)//".h"
            if (my_rank == 0) then
              open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
              call pepc_write_parameters(filehandle)
              close(filehandle)
            endif
          end subroutine


          subroutine write_particles_mpiio(comm, itime, n_total, dp, filename)
            use module_pepc_types
            use module_debug
            use module_utils
            use module_pepc, only : pepc_write_parameters
            implicit none
            integer(kind_pe) :: my_rank
            integer(kind_default), intent(in) :: itime
            integer*8, intent(in) :: n_total
            integer, intent(in) :: comm
            type(t_particle), intent(in) :: dp(:)
            character(*), intent(out) :: filename
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer(KIND=MPI_OFFSET_KIND) :: disp
            logical :: firstcall  = .true.
            character(50) :: dir
            integer(kind_particle) :: n_totsum

            dir = trim(directory)//"/mpi/"
            write(filename,'(a,"particle_",i6.6,".mpi")') trim(dir), itime
            call pepc_status("DUMP PARTICLES MPI: "//filename)
            
            call MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, ierr )

            if (firstcall) then
              call create_directory(trim(directory))
              call create_directory(trim(dir))
              firstcall = .false.
            endif

            n_totsum = size(dp, kind=kind(n_totsum))
            call MPI_ALLREDUCE(MPI_IN_PLACE, n_totsum,  1, MPI_KIND_PARTICLE, MPI_SUM, comm, ierr)
            if (n_totsum .ne. n_total) then
              DEBUG_ERROR('("Invalid total particle number: sum(nparticles_local) = ", I0, " but n_total = ", I0)', n_totsum, n_total)
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
            call MPI_FILE_WRITE_ORDERED(fh, dp(:), size(dp,kind=kind_default), MPI_TYPE_particle, status, ierr)
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


          subroutine read_particles_mpiio(itime_in, comm, itime, n_total, dp, filename, nparticles_local, file_exists)
            use module_pepc_types
            use module_debug
            use module_pepc, only : pepc_read_parameters
            implicit none
            integer(kind_default) :: itime_in
            integer(kind_default), intent(out) :: itime
            integer(kind_particle), intent(out) :: n_total
            integer, intent(in) :: comm
            character(*), intent(out) :: filename
            type(t_particle), allocatable :: dp(:)
            integer(kind_default), optional, intent(in) :: nparticles_local ! this has to be kind_default - not kind_particle as it is parameter to MPI functions
            logical, optional, intent(out) :: file_exists

            character(50) :: dir

            dir = trim(directory)//"/mpi/"
            write(filename,'(a,"particle_",i6.6,".mpi")') trim(dir), itime_in

            call read_particles_mpiio_from_filename(comm, itime, n_total, dp, filename, nparticles_local, file_exists)

            filename = trim(filename)//".h"
          end subroutine


          subroutine read_particles_mpiio_from_filename(comm, itime, n_total, dp, filename, nparticles_local, file_exists)
            use module_pepc_types
            use module_debug
            use module_pepc, only : pepc_read_parameters
            use module_utils, only : utils_file_exists => file_exists
            implicit none
            integer(kind_pe) :: my_rank, n_cpu
            integer(kind_default), intent(out) :: itime
            integer(kind_default) :: np_local ! this has to be default - not kind_particle as it is parameter to MPI functions
            integer(kind_particle), intent(out) :: n_total
            integer(kind_default), intent(in) :: comm
            character(*), intent(in) :: filename
            logical, optional, intent(out) :: file_exists

            character(255) :: filename2
            type(t_particle), allocatable :: dp(:)
            integer(kind_default) :: fh, ierr, status(MPI_STATUS_SIZE)
            integer(kind_particle) :: remain
            integer(kind_default), optional, intent(in) :: nparticles_local ! this has to be kind_default - not kind_particle as it is parameter to MPI functions
            integer(kind_particle) :: n_totsum


            call pepc_status("READ PARTICLES MPI: "//filename)
            
            call MPI_COMM_SIZE( MPI_COMM_WORLD, n_cpu,   ierr )
            call MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, ierr )
            
            if (utils_file_exists(trim(filename))) then
            
              if (present(file_exists)) file_exists = .true.

              call MPI_FILE_OPEN(comm,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

              if (ierr .ne. MPI_SUCCESS) then
                DEBUG_ERROR(*,'read_particles_mpiio_from_file(): file open failed', my_rank, ierr, filename)
              end if

              ! Set file view to BYTE for header, only rank 0 writes it
              call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
              call MPI_FILE_READ(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
              call MPI_FILE_READ(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)

              if (present(nparticles_local)) then
                np_local = nparticles_local
                n_totsum = int(np_local, kind=kind(n_totsum))
                call MPI_ALLREDUCE(MPI_IN_PLACE, n_totsum,  1, MPI_KIND_PARTICLE, MPI_SUM, comm, ierr)

                if (n_totsum .ne. n_total) then
                  DEBUG_ERROR('("Invalid total particle number: sum(nparticles_local) = ", I0, " but the data file says n_total = ", I0)', n_totsum, n_total)
                endif
              else
                np_local = int(n_total/n_cpu, kind(np_local))
                remain = n_total-np_local*n_cpu
                if ((remain > 0) .and. (my_rank < remain)) np_local = np_local+1
              endif

              if (allocated(dp)) deallocate(dp)
              allocate(dp(1:np_local))

              ! Redefine file view, now with our custom type
              call MPI_FILE_SET_VIEW(fh, max_header_size, MPI_TYPE_particle, MPI_TYPE_particle, 'native', MPI_INFO_NULL, ierr)
              ! Read particle data
              call MPI_FILE_READ_ORDERED(fh, dp, np_local, MPI_TYPE_particle, status, ierr)
              ! Close file
              call MPI_FILE_CLOSE(fh,ierr)

              filename2 = trim(filename)//".h"
              open(filehandle, file=trim(filename2),action='read')
              call pepc_read_parameters(filehandle)
              close(filehandle)
              
            else
              if (present(file_exists)) then
                file_exists = .false.
              else
                DEBUG_ERROR('("read_particles_mpiio_from_filename(filename=", a, ") - file does not exist")', trim(filename))
              endif
            endif
          end subroutine
end module module_checkpoint
