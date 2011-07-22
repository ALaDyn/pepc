!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Helper functions for checkpointing and restarting purposes
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_checkpoint
      implicit none
      include 'mpif.h'

      type, public :: t_dumpparticle
        sequence
        real*8 :: x, y, z, ux, uy, uz, q, m, pot, ex, ey, ez
        integer :: pelabel
      end type t_dumpparticle

      private

      integer, parameter :: filehandle = 91
      character(12), parameter :: directory = './particles/'
      integer(KIND=MPI_OFFSET_KIND) :: max_header_size = 1024

      integer, parameter :: nprops_dumpparticle = 13


      integer :: MPI_TYPE_DUMPPARTICLE
      logical :: registered_dumpparticle = .false.

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
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            type(t_dumpparticle), intent(in), dimension(np_local) :: dp
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
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            type(t_dumpparticle), intent(in), dimension(np_local) :: dp
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
          subroutine register_dumpparticle_type()
            implicit none
            integer :: ierr
            integer(KIND=MPI_ADDRESS_KIND), dimension(nprops_dumpparticle) :: displacements
            integer, dimension(nprops_dumpparticle) :: blocklengths, types
            integer(KIND=MPI_ADDRESS_KIND) :: base
            type(t_dumpparticle) :: dummy

            if (.not. registered_dumpparticle) then
              types = [ MPI_REAL8,  MPI_REAL8,  MPI_REAL8,  MPI_REAL8, MPI_REAL8,  MPI_REAL8,  MPI_REAL8,  MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_INTEGER ]
              blocklengths = 1

              call MPI_GET_ADDRESS(dummy, base, ierr)
              call MPI_GET_ADDRESS(dummy%x,  displacements(1),  ierr)
              call MPI_GET_ADDRESS(dummy%y,  displacements(2),  ierr)
              call MPI_GET_ADDRESS(dummy%z,  displacements(3),  ierr)
              call MPI_GET_ADDRESS(dummy%ux, displacements(4),  ierr)
              call MPI_GET_ADDRESS(dummy%uy, displacements(5),  ierr)
              call MPI_GET_ADDRESS(dummy%uz, displacements(6),  ierr)
              call MPI_GET_ADDRESS(dummy%q,  displacements(7),  ierr)
              call MPI_GET_ADDRESS(dummy%m,  displacements(8),  ierr)
              call MPI_GET_ADDRESS(dummy%pot,displacements(9),  ierr)
              call MPI_GET_ADDRESS(dummy%ex, displacements(10), ierr)
              call MPI_GET_ADDRESS(dummy%ey, displacements(11), ierr)
              call MPI_GET_ADDRESS(dummy%ez, displacements(12), ierr)
              call MPI_GET_ADDRESS(dummy%pelabel,  displacements(13), ierr)

              displacements = displacements - base

              call MPI_TYPE_CREATE_STRUCT(nprops_dumpparticle, blocklengths, displacements, types, MPI_TYPE_DUMPPARTICLE, ierr)
              registered_dumpparticle = .true.
            endif

          end subroutine


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_mpiio(comm, my_rank, itime, trun, np_local, n_total, dp)
            implicit none
            integer, intent(in) :: my_rank, itime, np_local
            integer*8, intent(in) :: n_total
            integer, intent(in) :: comm
            real*8, intent(in) :: trun
            type(t_dumpparticle), intent(in) :: dp(np_local)
            character(50) :: filename
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer(KIND=MPI_OFFSET_KIND) :: disp

            call system("mkdir -p " // trim(directory))
            write(filename,'(a,"particle_",i6.6,".mpi")') trim(directory), itime
            if(my_rank == 0) write(*,*) "write particles using MPI-IO to file ", filename

            call MPI_FILE_OPEN(comm,filename,IOR(MPI_MODE_RDWR,MPI_MODE_CREATE),MPI_INFO_NULL,fh,ierr)

            if (ierr .ne. MPI_SUCCESS) then
              write(*,*) 'write_particles_mpiio(): file open failed', my_rank, ierr, filename
              call MPI_ABORT(comm,1,ierr)
            end if

            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            if (my_rank == 0) then
              call MPI_FILE_WRITE(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
              call MPI_FILE_WRITE(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)
              call MPI_FILE_WRITE(fh, trun,    1, MPI_REAL,     status, ierr) ! current simulation time
              call MPI_FILE_GET_POSITION(fh, disp, ierr)

              if (disp > max_header_size) then
                write(*,*) "header_size is too small: ", max_header_size, "<", disp
                call MPI_ABORT(comm, 1,ierr)
              end if
            end if

            ! register dumpparticle type (if not already done)
            call register_dumpparticle_type()

            ! Redefine file view, now with our custom type
            call MPI_FILE_SET_VIEW(fh, max_header_size, MPI_TYPE_DUMPPARTICLE, MPI_TYPE_DUMPPARTICLE, 'native', MPI_INFO_NULL, ierr)
            ! Write particle data
            call MPI_FILE_WRITE_ORDERED(fh, dp(:), np_local, MPI_TYPE_DUMPPARTICLE, status, ierr)
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
            implicit none
            integer, intent(in) :: my_rank, n_cpu, itime_in
            integer, intent(out) :: itime, np_local
            integer*8, intent(out) :: n_total
            real*8, intent(out) :: trun
            integer, intent(in) :: comm
            type(t_dumpparticle), allocatable :: dp(:)
            character(50) :: filename
            integer :: fh, ierr, status(MPI_STATUS_SIZE)
            integer*8 :: remain

            write(filename,'(a,"particle_",i6.6,".mpi")') trim(directory), itime_in
            if(my_rank == 0) write(*,*) "loading particles using MPI-IO from file", filename

            call MPI_FILE_OPEN(comm,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

            if (ierr .ne. MPI_SUCCESS) then
              write(*,*) 'read_particles_mpiio(): file open failed', my_rank, ierr, filename
              call MPI_ABORT(comm,1,ierr)
            end if

            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_READ(fh, n_total, 1, MPI_INTEGER8, status, ierr) ! # particles
            call MPI_FILE_READ(fh, itime,   1, MPI_INTEGER,  status, ierr) ! Last successful timestep (new ts)
            call MPI_FILE_READ(fh, trun,    1, MPI_REAL,     status, ierr) ! current simulation time

            np_local = int(n_total/n_cpu, kind(np_local))
            remain = n_total-np_local*n_cpu
            if ((remain > 0) .and. (my_rank < remain)) np_local = np_local+1

            ! register dumpparticle type (if not already done)
            call register_dumpparticle_type()

            if (allocated(dp)) deallocate(dp)
            allocate(dp(1:np_local))

            ! Redefine file view, now with our custom type
            call MPI_FILE_SET_VIEW(fh, max_header_size, MPI_TYPE_DUMPPARTICLE, MPI_TYPE_DUMPPARTICLE, 'native', MPI_INFO_NULL, ierr)
            ! Read particle data
            call MPI_FILE_READ_ORDERED(fh, dp, np_local, MPI_TYPE_DUMPPARTICLE, status, ierr)
            ! Close file
            call MPI_FILE_CLOSE(fh,ierr)

          end subroutine


end module module_checkpoint
