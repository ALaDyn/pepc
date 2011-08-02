!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_diagnostics
   implicit none
   private



     public write_total_momentum
     public write_particles
     public read_particles



   contains

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles(allowcheckpoint)
            use physvars
            implicit none
            logical, intent(in) :: allowcheckpoint
            logical :: bin, asc, check, vtk
            integer, save :: lasttimestep = -1

            ! avoid calling this function several times per timestep
            if (lasttimestep .ne. itime) then

              bin = (idump_binary       > 0)
              if (bin) bin = ((mod(itime, idump_binary        ) == 0) .or. (itime == nt))
              asc = (idump              > 0)
              if (asc) asc = ((mod(itime, idump               ) == 0) .or. (itime == nt))
              check = (idump_checkpoint > 0) .and. (allowcheckpoint)
              if (check) check = ((mod(itime, idump_checkpoint) == 0) .or. (itime == nt))
              vtk = (idump_vtk          > 0)
              if (vtk) vtk = ((mod(itime, idump_vtk           ) == 0) .or. (itime == nt))

              call write_particles_type( bin, asc, check, vtk , itime==nt)
            endif

            lasttimestep = itime

          end subroutine write_particles


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine write_particles_type(binary, ascii, mpiio, vtk, final)
            use physvars
            use module_checkpoint
            implicit none
            include 'mpif.h'
            logical, intent(in) :: binary, ascii, mpiio, vtk, final
            integer :: p
            type(t_dumpparticle), allocatable :: dp(:)
            integer*8 :: npart

            if (binary .or. ascii .or. mpiio) then

              npart = npart_total ! TODO: conversion real*4 --> real*8 will be unneccessary soon

              allocate(dp(np_local))

              ! pack our data
              do p=1,np_local
                dp(p) = t_dumpparticle( x(p), y(p), z(p), ux(p), uy(p), uz(p), q(p), m(p), pot(p), ex(p), ey(p), ez(p), pelabel(p))
              end do

              !!! write particle date as a binary file
              if (binary) call write_particles_binary(my_rank, itime, np_local, dp)

              !!! write particle date as a text file
              if (ascii) call write_particles_ascii(my_rank, itime, np_local, dp)

              !!! write particle checkpoint data using mpi-io
              if (mpiio) call write_particles_mpiio(MPI_COMM_WORLD, my_rank, itime, trun, np_local, npart, dp)

              deallocate(dp)
            endif


            if (vtk) call write_particles_to_vtk(itime, final)

        end subroutine write_particles_type


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine read_particles(itime_in_)
            implicit none
            integer, intent(in) :: itime_in_

            call read_particles_type(itime_in_, .false., .false., .true.)

          end subroutine read_particles



          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine read_particles_type(itime_in_, binary, ascii, mpiio)
            use physvars
            use module_checkpoint
            implicit none
            include 'mpif.h'
            logical, intent(in) :: binary, ascii, mpiio
            integer, intent(in) :: itime_in_
            integer :: p
            type(t_dumpparticle), allocatable :: dp(:)
            integer*8 :: npart

            if (binary .or. ascii .or. mpiio) then

              !!! read particle date as a binary file
              if (binary) write(*,*) "read_particles(): binary mode unsupported" !call read_particles_binary(my_rank, itime, np_local, dp)

              !!! read particle date as a text file
              if (ascii)  write(*,*) "read_particles(): ascii mode unsupported" !call read_particles_ascii(my_rank, itime, np_local, dp)

              !!! read particle checkpoint data using mpi-io
              if (mpiio) call read_particles_mpiio(itime_in_, MPI_COMM_WORLD, my_rank, n_cpu, itime, trun, np_local, npart, dp)

              npart_total = npart

              ! unpack our data
              do p=1,np_local
                 x(p) = dp(p)%x
                 y(p) = dp(p)%y
                 z(p) = dp(p)%z
                ux(p) = dp(p)%ux
                uy(p) = dp(p)%uy
                uz(p) = dp(p)%uz
                 q(p) = dp(p)%q
                 m(p) = dp(p)%m
               pot(p) = dp(p)%pot
                ex(p) = dp(p)%ex
                ey(p) = dp(p)%ey
                ez(p) = dp(p)%ez
                pelabel(p) = dp(p)%pelabel
              end do

              if (allocated(dp)) deallocate(dp)
            endif

        end subroutine read_particles_type



          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_particles_to_vtk(step, final)
        use physvars
        use module_vtk
        use module_units
        implicit none
        integer, intent(in) :: step
        logical, intent(in) :: final
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

        call vtk%create_parallel("particles", step, my_rank, n_cpu, trun*unit_t0_in_fs, vtk_step)
          call vtk%write_headers(np_local, 0)
            call vtk%startpoints()
              call vtk%write_data_array("xyz", np_local, x, y, z)
            call vtk%finishpoints()
            call vtk%startpointdata()
              call vtk%write_data_array("velocity", np_local, ux, uy, uz)
              call vtk%write_data_array("el_field", np_local, ex, ey, ez)
              call vtk%write_data_array("el_pot", np_local, pot)
              call vtk%write_data_array("charge", np_local, q)
              call vtk%write_data_array("mass", np_local, m)
              call vtk%write_data_array("pelabel", np_local, pelabel)
              call vtk%write_data_array("local index", np_local, [(i,i=1,np_local)])
              call vtk%write_data_array("processor", np_local, [(my_rank,i=1,np_local)])
            call vtk%finishpointdata()
          call vtk%dont_write_cells()
          call vtk%write_final()
        call vtk%close()

      end subroutine write_particles_to_vtk


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !>
          !>
          !>
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_total_momentum(filename, itime_, trun_, selection, mom)
        use physvars
        implicit none
        include 'mpif.h'

        integer, intent(in) :: itime_
        character(*), intent(in) :: filename
        real*8, intent(in) :: trun_
        logical, intent(in) :: selection(1:np_local)
        real*8, intent(out) :: mom(4)
        real*8 :: tmp(4)
        real*8 :: r(4)
        integer :: p, ierr

        tmp = 0.

        do p = 1,np_local
          if (selection(p)) then
            r = [ux(p), uy(p), uz(p), 0._8]
            r(4) = sqrt(dot_product(r,r))
            tmp = tmp + r
          endif
        end do

        call MPI_ALLREDUCE(tmp, mom, 4, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr )

        if (my_rank == 0) then
          if (itime_ <= 1) then
             open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
          else
             open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
           endif
           write(87,'(i10,5g25.12)') itime_, trun_, mom
           close(87)
        endif
      end subroutine write_total_momentum

end module module_diagnostics
