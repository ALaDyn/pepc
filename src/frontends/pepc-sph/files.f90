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


  public openfiles
  public closefiles
  public write_particles
  public read_in_checkpoint
  public write_particles_to_vtk


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
     flush(6)

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
        flush(6)
     end if

 

     call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
 
 
     if( my_rank == 0 ) then

        if(IO_debug) then
           write(*,*) 'writing header', all_part, step
           flush(6)
        end if
 
        write(header, '(a,I10,a)') '# SPH_ASCII02 ', all_part, new_line('A')
        headerlength = len_trim(header)
 
        
        if(IO_debug) then
           write(*,*) 'writing header 1 ', headerlength, trim(header)
           flush(6)
        end if
 
        offset = 0_8
        call MPI_FILE_WRITE_AT(fh, offset, header, headerlength, MPI_CHARACTER, status, ierr)
 
        if(IO_debug) then
           write(*,*) 'writing header 2 ', status
           flush(6)
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
           flush(6)
        end if

        offset = int(headerlength, 8)
        call MPI_FILE_WRITE_AT(fh, offset, extended_header, extended_header_length, MPI_CHARACTER, status, ierr)
 
        if(IO_debug) then
           write(*,*) 'writing header done'
           flush(6)
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
 
 


end module files
