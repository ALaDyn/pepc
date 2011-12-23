!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Helper functions for checkpointing and restarting purposes
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "pepc_debug.h"
module module_checkpoint
      implicit none
      include 'mpif.h'

      private

      integer, parameter :: filehandle = 91
      character(12), parameter :: directory = './particles/'
      integer(KIND=MPI_OFFSET_KIND) :: max_header_size = 1024

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
          subroutine write_particles_ascii(my_rank, itime, np_local, dp)
            use module_pepc_types
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            type(t_particle), intent(in), dimension(np_local) :: dp
            character(50) :: filename

            call system("mkdir -p " // trim(directory))
            write(filename,'(a,"particle_",i6.6,"_",i6.6,".dat")') trim(directory), itime, my_rank
            if(my_rank == 0) write(*,*) "write particles in text mode to file ", filename
            open(filehandle, file=trim(filename), STATUS='REPLACE')

            write(filehandle,*) "# particle positions at timestep ", itime, ": x y z vx vy vz q m pot ex ey ez pelabel"
            write(filehandle,'(12E14.4E2," ",I16.12)') dp(1:np_local)

            close(filehandle)
          end subroutine


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_binary(my_rank, itime, np_local, dp)
            use module_pepc_types
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            type(t_particle), intent(in), dimension(np_local) :: dp
            character(50) :: filename

            call system("mkdir -p " // trim(directory))
            write(filename,'(a,"particle_",i6.6,"_",i6.6,".bin")') trim(directory), itime, my_rank
            if(my_rank == 0) write(*,*) "write particles in binary mode to file ", filename
            open(filehandle, file=trim(filename), STATUS='REPLACE', ACCESS="STREAM")

             write(filehandle) dp(1:np_local)

            close(filehandle)
          end subroutine


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_mpiio(comm, my_rank, itime, trun, np_local, n_total, dp)
            use module_pepc_types
            use module_debug
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            integer*8, intent(in) :: n_total
            integer, intent(in) :: comm
            real*8, intent(in) :: trun
            type(t_particle), intent(in) :: dp(np_local)
            character(50) :: filename
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer(KIND=MPI_OFFSET_KIND) :: disp

            call system("mkdir -p " // trim(directory))
            write(filename,'(a,"particle_",i6.6,".mpi")') trim(directory), itime
            if(my_rank == 0) write(*,*) "write particles using MPI-IO to file ", filename

            call MPI_FILE_OPEN(comm,filename,IOR(MPI_MODE_RDWR,MPI_MODE_CREATE),MPI_INFO_NULL,fh,ierr)

            if (ierr .ne. MPI_SUCCESS) then
              DEBUG_ERROR(*, 'write_particles_mpiio(): file open failed', my_rank, ierr, filename)
            end if

            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            if (my_rank == 0) then
              call MPI_FILE_WRITE(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
              call MPI_FILE_WRITE(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)
              call MPI_FILE_WRITE(fh, trun,    1, MPI_REAL,     status, ierr) ! current simulation time
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

          end subroutine



          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine read_particles_mpiio(itime_in, comm, my_rank, n_cpu, itime, trun, np_local, n_total, dp)
            use module_pepc_types
            use module_debug
            implicit none
            integer, intent(in) :: my_rank, n_cpu, itime_in
            integer, intent(out) :: itime, np_local
            integer*8, intent(out) :: n_total
            real*8, intent(out) :: trun
            integer, intent(in) :: comm
            type(t_particle), allocatable :: dp(:)
            character(50) :: filename
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer*8 :: remain

            write(filename,'(a,"particle_",i6.6,".mpi")') trim(directory), itime_in
            if(my_rank == 0) write(*,*) "loading particles using MPI-IO from file", filename

            call MPI_FILE_OPEN(comm,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

            if (ierr .ne. MPI_SUCCESS) then
              DEBUG_ERROR(*,'read_particles_mpiio(): file open failed', my_rank, ierr, filename)
            end if

            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_READ(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
            call MPI_FILE_READ(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)
            call MPI_FILE_READ(fh, trun,    1, MPI_REAL,     status, ierr) ! current simulation time

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

          end subroutine


end module module_checkpoint
