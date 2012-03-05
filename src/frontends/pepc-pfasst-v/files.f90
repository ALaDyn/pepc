module files
  implicit none


    public openfiles
    public closefiles

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Open plain text output files
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine openfiles
      use physvars
      use module_utils


      if (my_rank_space == 0) then
         !  master diagnostics output
         open(15,file='run.out')
         open(70,file='domains.dat')
         open(66,file='linear_diag.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
     endif

     ! for MPI I/O
     call create_directory("part_data")

    end subroutine openfiles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Close plain text output files
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine closefiles
      use physvars

      if (my_rank_space == 0) then
         close(15)
         close(70)
         close(66)
      endif

      !close(20)
      close(80)  ! initial particle data

    end subroutine closefiles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Dump VTK or checkpoint
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dump(i,simtime)

        use physvars
        implicit none
        include 'mpif.h'

        integer, intent(in) :: i
        real, intent(in) :: simtime
        integer :: fh, ierr, err, status(MPI_STATUS_SIZE)
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
                write(*,*) 'something is wrong here: file open failed',my_rank_space,ierr,cfile
                call MPI_ABORT(MPI_COMM_WORLD,err,ierr)
            end if
            ! Set file view to BYTE for header, only rank 0 writes it
            call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
            if (my_rank_space == 0) then
                call MPI_FILE_WRITE(fh,n,1,MPI_INTEGER,status,ierr)    ! # particles
                call MPI_FILE_WRITE(fh,dt,1,MPI_REAL,status,ierr)          ! timestep
                call MPI_FILE_WRITE(fh,ts,1,MPI_REAL,status,ierr)          ! Starting time
                call MPI_FILE_WRITE(fh,i,1,MPI_INTEGER,status,ierr)        ! Last successful timestep (number)
                call MPI_FILE_WRITE(fh,te,1,MPI_REAL,status,ierr)          ! Final time
                call MPI_FILE_WRITE(fh,nu,1,MPI_REAL,status,ierr)          ! Viscousity
                call MPI_FILE_WRITE(fh,h,1,MPI_REAL,status,ierr)           ! Original particle distance
                call MPI_FILE_WRITE(fh,m_h,1,MPI_REAL,status,ierr)         ! Remeshing distance
                call MPI_FILE_WRITE(fh,rem_freq,1,MPI_INTEGER,status,ierr) ! Remeshing frequence
                call MPI_FILE_WRITE(fh,thresh,1,MPI_REAL8,status,ierr)     ! threshold for pop. control
                call MPI_FILE_WRITE(fh,eps,1,MPI_REAL8,status,ierr)         ! core size
                call MPI_FILE_GET_POSITION(fh, disp, ierr)
                if (disp .gt. header_disp) then
                    write(*,*) "header_size is too small: ", header_disp, "<", disp
                    call MPI_ABORT(MPI_COMM_WORLD,err,ierr)
                end if
            end if

            ! Redefine file view, now with our custom type
            call MPI_FILE_SET_VIEW(fh, header_disp, MPI_TYPE_PARTICLE, MPI_TYPE_PARTICLE, 'native', MPI_INFO_NULL, ierr)
            ! Write particle data
            call MPI_FILE_WRITE_ORDERED(fh, vortex_particles(1:np), np, MPI_TYPE_PARTICLE, status, ierr)
            ! Take care before closing
            call MPI_FILE_SYNC(fh,ierr)
            call MPI_FILE_CLOSE(fh,ierr)
        end if

    end subroutine dump


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Read in data from MPI checkpoint file
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_in_checkpoint

        use physvars
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
        call MPI_FILE_READ_ORDERED(fh, vortex_particles(1:np), np, MPI_TYPE_PARTICLE, status, ierr)
        call MPI_FILE_CLOSE(fh,ierr)

    end subroutine read_in_checkpoint


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Dump particles to binary parallel VTK
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles_to_vtk(step,time)

        use physvars
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

        call vtk%create_parallel("particles", step, my_rank_space, n_cpu_space, 0.1D01*time, vtk_step, MPI_COMM_SPACE)
        call vtk%write_headers(np,0)
        call vtk%startpoints()
            call vtk%write_data_array("xyz", np, vortex_particles(1:np)%x(1), vortex_particles(1:np)%x(2), vortex_particles(1:np)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
            call vtk%write_data_array("velocity", np, vortex_particles(1:np)%results%u(1), vortex_particles(1:np)%results%u(2), vortex_particles(1:np)%results%u(3))
            call vtk%write_data_array("vorticity", np, vortex_particles(1:np)%data%alpha(1), vortex_particles(1:np)%data%alpha(2), vortex_particles(1:np)%data%alpha(3))
            call vtk%write_data_array("work", np, vortex_particles(1:np)%work)
            call vtk%write_data_array("label", np, vortex_particles(1:np)%label)
            call vtk%write_data_array("pid", np, [(my_rank_space,i=1,np)])
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

    end subroutine write_particles_to_vtk

    subroutine dump_results()

        use physvars
        implicit none

        integer :: i
        character(50) :: resfile

        write(resfile,'(a,i6.6,a)') "part_data/results_", my_rank_space,".dat"

        open(11,file=resfile)

        do i = 1, np
            write(11,*) my_rank_space, i, vortex_particles(i)%label, vortex_particles(i)%x(1:3), vortex_particles(i)%data%alpha(1:3)
        end do

        close(11)

    end subroutine dump_results

end module files
