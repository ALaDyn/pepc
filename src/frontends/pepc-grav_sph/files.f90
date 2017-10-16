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

module files
  implicit none
  private

  logical, public :: IO_debug = .false.
  integer, public :: SPH_ASCII02_columnwidth = 12

  integer, private :: num_known_output_formats = 3

  character(255), public :: parameter_file_name
  character(255), public :: data_file_name

  logical, public :: parameter_file_available
  logical, public :: data_file_available
  
  integer, public :: output_type
  integer, public :: run_type

  public openfiles
  public closefiles
  public write_particles
  public read_in_checkpoint
  public write_particles_to_vtk
  public process_commandline_arguments
  public read_particles

contains



  subroutine openfiles
    use physvars

    if (my_rank == 0) then
       !  master diagnostics output
       open(15,file='run.out')
       open(70,file='domains.dat')
    endif

  end subroutine openfiles




  subroutine closefiles
    use physvars

    if (my_rank == 0) then
       close(15)
       close(70)
    endif

    close(80)  ! initial particle data

  end subroutine closefiles


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>   Dump VTK or checkpoint
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_particles(i,simtime)

    use physvars

    use module_pepc_types

    implicit none
    include 'mpif.h'

    integer, intent(in) :: i
    real, intent(in) :: simtime
    integer :: p, fh, ierr, err, status(MPI_STATUS_SIZE)
    integer(KIND=MPI_OFFSET_KIND) :: disp, header_disp=1024

    character(50) :: cfile

    character(100) :: filename

    if ((dump_time.ne.0).and.(mod(i,dump_time) == 0)) then

       call write_particles_to_vtk(i,simtime)
       
    end if
    
    if ((cp_time.ne.0).and.(mod(i,cp_time) == 0)) then
       ! Open new file for i-th timestep
       write(mpifile,'(a,i6.6,a)') "part_data/particle_", i,".mpi"
       call MPI_FILE_OPEN(MPI_COMM_WORLD,mpifile,IOR(MPI_MODE_RDWR,MPI_MODE_CREATE),MPI_INFO_NULL,fh,ierr)
       if (ierr .ne. MPI_SUCCESS) then
          write(*,*) 'something is wrong here: file open failed',my_rank,ierr,cfile
          call MPI_ABORT(MPI_COMM_WORLD,err,ierr)
          stop
       end if
       ! Set file view to BYTE for header, only rank 0 writes it
       call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
       if (my_rank == 0) then
          call MPI_FILE_WRITE(fh,npart_total,  1,MPI_INTEGER,status,ierr)    ! # particles
          call MPI_FILE_WRITE(fh,dt,           1,MPI_REAL,status,ierr)          ! timestep
          call MPI_FILE_WRITE(fh,i,            1,MPI_INTEGER,status,ierr)        ! Last successful timestep (number)
          call MPI_FILE_GET_POSITION(fh, disp, ierr)
          if (disp .gt. header_disp) then
             write(*,*) "header_size is too small: ", header_disp, "<", disp
             call MPI_ABORT(MPI_COMM_WORLD,err,ierr)
             stop
          end if
       end if

       ! Redefine file view, now with our custom type
       call MPI_FILE_SET_VIEW(fh, header_disp, MPI_TYPE_PARTICLE, MPI_TYPE_PARTICLE, 'native', MPI_INFO_NULL, ierr)
       ! Write particle data
       call MPI_FILE_WRITE_ORDERED(fh, particles(1:np_local), np_local, MPI_TYPE_PARTICLE, status, ierr)
       ! Take care before closing
       call MPI_FILE_SYNC(fh,ierr)
       call MPI_FILE_CLOSE(fh,ierr)
    end if


    write(filename, '(a,I6.6)' ) 'particles_', i          ! all processes start write subroutine with the correct filename so BCAST before
    call MPI_BCAST(filename, 100, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )     ! just to make sure all write to same file
    call write_particles_sph_ascii02(filename,i)


  end subroutine write_particles



  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>   Read in data from MPI checkpoint file
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_in_checkpoint
    
    use physvars
   
    use module_pepc_types
 
    implicit none
    include 'mpif.h'

    integer :: ierr, fh
    integer(KIND=MPI_OFFSET_KIND) :: header_disp=1024
    integer :: status(MPI_STATUS_SIZE)
    
    ! Open new file for i-th timestep
    call MPI_FILE_OPEN(MPI_COMM_WORLD,mpifile,IOR(MPI_MODE_RDWR,MPI_MODE_CREATE),MPI_INFO_NULL,fh,ierr)
    ! Redefine file view, now with our custom type 
    call MPI_FILE_SET_VIEW(fh, header_disp, MPI_TYPE_PARTICLE, MPI_TYPE_PARTICLE, 'native', MPI_INFO_NULL, ierr)
    ! Read particle data
    call MPI_FILE_READ_ORDERED(fh, particles(1:np_local), np_local, MPI_TYPE_PARTICLE, status, ierr)
    call MPI_FILE_CLOSE(fh,ierr)

  end subroutine read_in_checkpoint


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>   Dump particles to binary parallel VTK
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_particles_to_vtk(step,time)

    use physvars, only: &
         np => np_local, &
         particles, &
         n_cpu, &
         my_rank, &
         nt, &
         dump_time

    use module_interaction_specific

    use module_pepc_types
    
    use module_vtk


    implicit none

    real, intent(in) :: time
    integer, intent(in) :: step
    integer :: i
    type(vtkfile_unstructured_grid) :: vtk
    integer :: vtk_step

    if (step .eq. 0) then
       vtk_step = VTK_STEP_FIRST
    else if (step .eq. (nt - mod(nt, dump_time)  ) ) then ! last call of write_particles_to_vtk according to dump_time
       vtk_step = VTK_STEP_LAST
    else
       vtk_step = VTK_STEP_NORMAL
    endif
    
    call vtk%create_parallel("particles", step, my_rank, n_cpu, 0.1D01*time, vtk_step)
    call vtk%write_headers(np,0)
    call vtk%startpoints()
    call vtk%write_data_array("xyz",              particles(1:np)%x(1), particles(1:np)%x(2), particles(1:np)%x(3) )
    call vtk%finishpoints()
    call vtk%startpointdata()
    call vtk%write_data_array("velocity",         particles(1:np)%data%v(1), particles(1:np)%data%v(2), particles(1:np)%data%v(3) )
    call vtk%write_data_array("mass",             particles(1:np)%data%q)
    call vtk%write_data_array("work",             particles(1:np)%work)
    call vtk%write_data_array("label",            particles(1:np)%label)
    call vtk%write_data_array("temperature",      particles(1:np)%data%temperature)
    call vtk%write_data_array("smoothing-length", particles(1:np)%results%h)
    call vtk%write_data_array("density",          particles(1:np)%results%rho)
    call vtk%write_data_array("pid",              np, my_rank)
    call vtk%write_data_array("sph-force",        particles(1:np)%results%sph_force(1), particles(1:np)%results%sph_force(2), particles(1:np)%results%sph_force(3) )
    call vtk%write_data_array("particle-type",    particles(1:np)%data%type)
    call vtk%finishpointdata()
    call vtk%dont_write_cells()
    call vtk%write_final()
    call vtk%close()

  end subroutine write_particles_to_vtk




  subroutine write_particles_sph_ascii02(filename, step)
    ! write particle data as ascii using MPI-IO with own SPH_ASCII02 format.
    ! SPH_ASCII02 format allows the specification of the width of the columns.
    ! This output method is probably very slow due to the conversion.
    ! It is not intended to be used for real simulations, only for tests with low particle numbers.
    
    use physvars, only: &
         my_rank, &
         np_local, &
         npart_total, &
         particles
    
!     use utils
 
     implicit none
 
     include 'mpif.h'
 
     character(*), intent(in) :: filename             ! all processes start this subroutine with the correct filename so BCAST before
     integer, intent(in) :: step

     integer :: width                     ! width of the data columns in the file. columns are seperated by a blank

     integer :: p, ierr
     real*8 :: writesize, t1, t2, t3, t4
 
 
     logical :: file_exists
 
 
     integer :: part_including_mine
     integer :: part_before_me
     integer :: fh      ! filehandle
     integer :: bufferlength
     integer*8 :: offset, current_offset
     character :: buffer*400
     integer :: headerlength
     character :: header*30
     integer :: extended_header_length = 5000
     character :: extended_header*5000
     integer, dimension(MPI_STATUS_SIZE) :: status
     integer :: i
     integer :: all_part


 
     character(100) :: formatstring

     !TODO: should be set for inputfile
     SPH_ASCII02_columnwidth = 16

     if(my_rank == 0) write(*,*) "IO: write particles in mpi-io ascii 02 mode. Step:", step
!     flush(6)

     width = SPH_ASCII02_columnwidth

     if(width<12) then
        write(*,*) 'width of SPH_ASCII02 datacolumns is below 12. Set to a value >=12.'
        ! width < 12 makes no sense and is therefore not implemented.
        call MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if


 
     call MPI_SCAN(np_local, part_including_mine, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     part_before_me = part_including_mine - np_local
 
     all_part = 0
     call MPI_REDUCE(np_local, all_part, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

     if(IO_debug) then
        write(*,*) 'opening file'
!        flush(6)
     end if

 

     call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
 
 
     if( my_rank == 0 ) then

        if(IO_debug) then
           write(*,*) 'writing header', all_part, step
!           flush(6)
        end if
 
        write(header, '(a,I10,a)') '# SPH_ASCII02 ', all_part, new_line('A')
        headerlength = len_trim(header)
 
        
        if(IO_debug) then
           write(*,*) 'writing header 1 ', headerlength, trim(header)
!           flush(6)
        end if
 
        offset = 0_8
        call MPI_FILE_WRITE_AT(fh, offset, header, headerlength, MPI_CHARACTER, status, ierr)
 
        if(IO_debug) then
           write(*,*) 'writing header 2 ', status
!           flush(6)
        end if

        write( formatstring, '(a,I2,a,I2,a)' ) '(a,I6,2a,I2,2a,a', width-2, ',x,15(a', width, ',x),a)'
        
        if(IO_debug) then
           write (*,*) 'format string for header:', formatstring
        end if


        write(extended_header, formatstring ) "# particle properties for timestep ", step, new_line('A'), &
             "# width of datacolumns ", width, new_line('A'), &
             "# ", &
             "pelabel", "type", &
             "x", "y", "z", &
             "vx", "vy", "vz", &
             "m", "temperature", &
             "h", "rho", &
             "ex", "ey", "ez", "temp_change", &
             new_line('A')
 
        extended_header_length = len_trim(extended_header)
 
        if(IO_debug) then
           write(*,*) 'writing header 3'
!           flush(6)
        end if

        offset = int(headerlength, 8)
        call MPI_FILE_WRITE_AT(fh, offset, extended_header, extended_header_length, MPI_CHARACTER, status, ierr)
 
        if(IO_debug) then
           write(*,*) 'writing header done'
!           flush(6)
        end if
        
        current_offset = int(headerlength, 8) + int(extended_header_length, 8)

     end if
 
 
     call MPI_BCAST(current_offset, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)


     bufferlength = 16*(width+1)
 

     if( width-7 > 9) then
        write( formatstring, '(a,I2,a,I2,a,I2,a,I2,a)' ) '(I', width, ',x,I', width, ',14(x,E', width, '.', width-7, 'E2),a)'
     else
        write( formatstring, '(a,I2,a,I2,a,I2,a,I1,a)' ) '(I', width, ',x,I', width, ',14(x,E', width, '.', width-7, 'E2),a)'
     end if
     
     if(IO_debug) then
        write (*,*) 'format string for ascii conversion:', formatstring
     end if

     do p =1, np_local
        
        !if the number of written values is changed, the formatstring and bufferlength has to be changed, too.
        write(buffer, formatstring ) particles(p)%label, particles(p)%data%type, particles(p)%x(1), particles(p)%x(2), particles(p)%x(3), &
             particles(p)%data%v(1), particles(p)%data%v(2), particles(p)%data%v(3), particles(p)%data%q, particles(p)%data%temperature, &
             particles(p)%results%h, particles(p)%results%rho, &
             particles(p)%results%sph_force(1), particles(p)%results%sph_force(2), particles(p)%results%sph_force(3), &
             particles(p)%results%temperature_change, new_line('A')
        offset = current_offset + ( (part_before_me + p -1 ) * bufferlength )
        call MPI_FILE_WRITE_AT(fh, offset, buffer, bufferlength, MPI_CHARACTER, status, ierr)
     end do
     
     
     call MPI_FILE_CLOSE(fh, ierr)
 
     !\bug ab: update writesize
     writesize = real(current_offset + npart_total * bufferlength )/1024.0/1024.0
 
   end subroutine write_particles_sph_ascii02
 


   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! >
   ! > Process the commandline arguments
   ! >
   ! > Gets name of parameterfile (if available) and passes it to pepc_read_parameters_from_file_name.
   ! > Gets name of datafile (if available) and returns it.
   ! > Has same functionality for pepc as pepc_read_parameters_from_first_argument.
   ! >
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine process_commandline_arguments()
     
     use physvars, only: &
          my_rank, &
          dump_time, &
          cp_time, &
          ispecial

     use module_debug, only: &
          pepc_status

     use module_pepc, only: &
          pepc_read_parameters_from_file_name
     
     implicit none
     include 'mpif.h'
     
     logical            :: available
     character(len=255) :: file_name
     logical            :: data_available
     character(len=255) :: datafile
!     integer            :: run_type         ! what is performed during the run.

     character(len=255) :: argument1
     character(len=255) :: argument2
     character(len=255) :: argument3
     integer :: ierr
     logical :: syntax_error = .false.
     integer :: output_type
     

     available      = .false.
     file_name      = ''
     data_available = .false.
     datafile       = ''
     run_type       = 1


 
     ! rank 0 reads in first command line argument
     if (my_rank .eq. 0) then
        if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
           call GET_COMMAND_ARGUMENT(1, argument1)
           write(*,*) "argument1:", argument1
           if( argument1(1:7) .eq. '-resume' ) then
              if( COMMAND_ARGUMENT_COUNT() > 2 ) then
                 call GET_COMMAND_ARGUMENT(2, argument2)
                 call GET_COMMAND_ARGUMENT(3, argument3)
                 available = .true.
                 file_name = trim(argument2)
                 datafile = trim(argument3)
                 run_type = 2
              else
                 write(*,*) "======================================================================="
                 write(*,*) "Error: To few arguments for -resume. Need parameter file and data file."
                 write(*,*) ""
                 syntax_error = .true.
              end if
              
           else if( argument1(1:8) .eq. '-convert' ) then
              write(*,*) "converting"
              if( COMMAND_ARGUMENT_COUNT() > 2 ) then
                 call GET_COMMAND_ARGUMENT(2, argument2)
                 read(argument2,*) output_type
                 if( output_type > 0 .and. output_type < num_known_output_formats) then
                    call GET_COMMAND_ARGUMENT(3, argument3)
                    data_available = .true.
                    datafile = trim(argument3)
                    run_type = 3
                    ispecial = -1
                    ! vtk dump every dump_time steps
                    dump_time = 1
                    
                    ! checkpoint files every cp_time steps
                    cp_time = 1
                    write (*,*) "set runtype to 3"
                 else
                    write(*,*) "======================================================================="
                    write(*,*) "Error: Unknown output_type."
                    write(*,*) ""
!                    call flush()
                    syntax_error = .true.
                 end if
              else
                 write(*,*) "======================================================================="
                 write(*,*) "Error: To few arguments for -convert. Need target type and data file."
                 write(*,*) ""
!                 call flush()
                 syntax_error = .true.
              end if
              
           else
              !if( argument1(1:5) .eq. 'run.h' ) then
              available = .true.
              file_name = trim(argument1)
              run_type = 1

              if( COMMAND_ARGUMENT_COUNT() > 1 ) then
                 call GET_COMMAND_ARGUMENT(2, argument2)
                 datafile = trim(argument2)
              end if

           end if

        else
           ! default run

        end if
     
     end if
  
 
     if( syntax_error ) then 
        
        if(my_rank .eq. 0) then
           write(*,*) "Usage: pepc-...                                   To run the standard setup."
           write(*,*) "       pepc-... run.h                             To run a user specified setup."
           write(*,*) "       pepc-... -resume run.h datafile            To read physical conditions from run.h and particle data from datafile."
           write(*,*) "       pepc-... -convert outputformat datafile    To convert datafile into outputformat."
           write(*,*) ""
           write(*,*) "Allowed outputformats are:"
           write(*,*) "1      for ASCII"
           write(*,*) "2      for VTK BINARY"
           write(*,*) "3      for MPIIO"
           write(*,*) ""
           write(*,*) "Exiting."
        end if
        
        call MPI_ABORT(MPI_COMM_WORLD, ierr)
        
     end if

        

     
        
        
     call MPI_BCAST( available, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
     
     if( available ) then
        
        call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
        
        if(my_rank .eq. 0) write(*,*) "found parameter file: ", file_name
        
     end if
        
        
     call MPI_BCAST( data_available, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
     
     if( data_available ) then
        
        call MPI_BCAST( datafile, 255, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
        
        write(*,*) "found data file: ", datafile
     end if

     call MPI_BCAST( run_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

     if( available ) then
        call pepc_status("INIT FROM FILE "//file_name)
        call pepc_read_parameters_from_file_name(file_name)
     else
        call pepc_status("INIT WITH DEFAULT PARAMETERS")
     end if

     parameter_file_available = available
     parameter_file_name      = file_name
     data_file_available      = data_available
     data_file_name           = datafile

   end subroutine process_commandline_arguments



   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! >
   ! >
   ! >
   ! >
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_particles()
     
     use physvars

     implicit none
     include 'mpif.h'
     
     integer :: ierr
     integer :: filetype = -1    ! 1 for ascii02, 2 for binary vtk, 3 for mpiio, 4 for ascii01, -1 for unknown
     character(50) ::  header

     
     if( my_rank .eq. 0 ) write(*,*) "Reading particles from file", trim(data_file_name)
     
     if( my_rank .eq. 0 ) then
        
        open(33, file=trim(data_file_name)) 
        read(33, '(a)' ) header
        
        if( header(1:14) .eq. '# SPH_ASCII01' ) then
           
           filetype = 4
           
!        else if( header(1:27) .eq. '# vtk DataFile Version 3.0' ) then
!           
!           filetype = -1
           
        else
           
           write(*,*) 'unknown data file format', header   
           call MPI_ABORT(MPI_COMM_WORLD, 12, ierr)
           
        end if
        
        close(33)
        
     end if
     
     
     call MPI_BCAST( filetype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
     
     if( filetype .eq. 4 ) then
        
        call read_sph_ascii01(data_file_name)
        
!     else if( filetype .eq. 1 ) then 
        
!        call read_sph_vtk(datafilename, resume_time)
        
     else
        
        write(*,*) 'not yet implemented'
        call MPI_ABORT(MPI_COMM_WORLD, 12, ierr)
        
     end if

     
   end subroutine read_particles


   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! > 
   ! > read particle data from SPH_ASCII01 file using MPI-IO
   ! >
   ! >
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_sph_ascii01(filename)

     use physvars
!, only: &
!          my_rank, &
!!          np_local, &
!          npart_total, &
!          particles

    implicit none
    include 'mpif.h'

    character(*), intent(in) :: filename               ! all processes start this subroutine with the correct filename so BCAST before

    integer :: step
    integer :: p, ierr, tmp_int

    integer :: part_including_mine
    integer :: part_before_me
    integer :: fh      ! filehandle
    integer :: bufferlength
    integer*8 :: offset
    character(5000) :: buffer
    integer :: headerlength
    character(5000) ::  header
    character(5000) ::  header2
    integer :: namelistlength = 5000
    character :: namelist*5000
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: input_particles
    integer :: current_pos

    integer :: resume_time


    integer :: ftell  ! FIXME: ftell is a GNU extension, not available in xlf, see http://gcc.gnu.org/onlinedocs/gfortran/FTELL.html


    if(my_rank .eq. 0) write(*,*) "read particles from SPH_ASCII01 file with MPIIO: ", TRIM(filename)

    call MPI_SCAN(np_local, part_including_mine, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    part_before_me = part_including_mine - np_local

    if( my_rank .eq. 0) then 

       open(33, file=trim(filename)) 
       read(33, '(a)' ) header

       if( header(1:14) .eq. '# SPH_ASCII01' ) then

          read(header(15:25), '(I10.10)')  input_particles

          if( input_particles .ne. npart_total) then
             
             if(run_type .eq. 3) then ! converting data
                npart_total = input_particles
                
                np_local = npart_total/n_cpu
                if (np_local*n_cpu .ne. npart_total .and. mod(npart_total,n_cpu) > my_rank)  np_local = npart_total/n_cpu+1
                
                deallocate(particles)
                allocate( particles(np_local) )
                
             else
                write(69,*) 'particle number in data file not equal particle number in setup file', input_particles, npart_total
                close(69)
                
                call MPI_ABORT(MPI_COMM_WORLD, 12, ierr)

             end if
          end if

       else

          write(*,*) 'unknown data file format', header   
          call MPI_ABORT(MPI_COMM_WORLD, 12, ierr)
       end if

       read(33, '(a)' ) header

       headerlength = len_trim(header)

       read(header(headerlength-5:headerlength) , '(I6)' ) resume_time

       itime = resume_time

       close(33)

       ! find first character of first not comment line
       open(33, file=trim(filename), action='read') 
       read(33, '(a)') header
       current_pos = ftell(33) ! ftell is a intrinsic function
       headerlength = 0
       do while(header(1:1) .eq. '#')
          read(33 , '(a)' ) header
          headerlength = current_pos
          current_pos = ftell(33)
       end do
       close(33)
       ! headerlength now points to first byte of first not comment line

    end if

    call MPI_BCAST( resume_time, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    call MPI_BCAST( headerlength, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

    call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)


    bufferlength = 10+1+10*(12+1)+6+1+3*(12+1)+12+1

    do p =1, np_local
       !if the format string is changed, the bufferlength has to be changed, too.
       offset = int(headerlength, 8) + ( (part_before_me + p -1 ) * bufferlength )
       call MPI_FILE_READ_AT(fh, offset, buffer, bufferlength, MPI_CHARACTER, status, ierr)
       read(buffer,'(I10.10,x,10(E12.5E2,x),I6.6,x,3(E12.5E2,x),E12.5E2,a))') particles(p)%label, particles(p)%x(1), particles(p)%x(2), particles(p)%x(3), &
            particles(p)%data%v(1), particles(p)%data%v(2), particles(p)%data%v(3), particles(p)%data%q, particles(p)%results%h, particles(p)%results%rho, &
            particles(p)%data%temperature, tmp_int,  particles(p)%results%e(1), particles(p)%results%e(2), particles(p)%results%e(3), particles(p)%results%temperature_change
       particles(p)%results%h = particles(p)%results%h/2. ! SPH_ASCII01 files contain r_nn = 2h
    end do

    call MPI_FILE_CLOSE(fh, ierr)

  end subroutine read_sph_ascii01


end module files
