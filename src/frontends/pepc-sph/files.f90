module files
  implicit none
  private

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
         nt

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
    else if (step .eq. nt) then
       vtk_step = VTK_STEP_LAST
    else
       vtk_step = VTK_STEP_NORMAL
    endif
    
    call vtk%create_parallel("particles", step, my_rank, n_cpu, 0.1D01*time, vtk_step)
    call vtk%write_headers(np,0)
    call vtk%startpoints()
    call vtk%write_data_array("xyz",              np, particles(1:np)%x(1), particles(1:np)%x(2), particles(1:np)%x(3) )
    call vtk%finishpoints()
    call vtk%startpointdata()
    call vtk%write_data_array("velocity",         np, particles(1:np)%data%v(1), particles(1:np)%data%v(2), particles(1:np)%data%v(3) )
    call vtk%write_data_array("mass",             np, particles(1:np)%data%q)
    call vtk%write_data_array("work",             np, particles(1:np)%work)
    call vtk%write_data_array("label",            np, particles(1:np)%label)
    call vtk%write_data_array("temperature",      np, particles(1:np)%data%temperature)
    call vtk%write_data_array("smoothing-length", np, particles(1:np)%results%h)
    call vtk%write_data_array("density",          np, particles(1:np)%results%rho)
    call vtk%write_data_array("pid",              np, [(my_rank,i=1,np)])
    call vtk%finishpointdata()
    call vtk%dont_write_cells()
    call vtk%write_final()
    call vtk%close()

  end subroutine write_particles_to_vtk






end module files
