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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_diagnostics
    use module_pepc_kinds
    implicit none
    private

    public write_particles
    public write_particles_vtk
    public read_particles
    public compute_force_direct

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compute_force_direct(particles)
        use physvars, only : MPI_COMM_PEPC
        use module_pepc_types
        use module_timings
        use module_directsum
        use module_debug, only : pepc_status
        use module_interaction_specific_types, only: t_particle_results
        implicit none
        type(t_particle), intent(inout) :: particles(:) !< input particle data, initializes %x, %data appropriately (and optionally set %label) before calling this function

        integer(kind_particle) :: i
        type(t_particle_results), allocatable :: directresults(:)

        call pepc_status('PEPC-MW: DIRECTSUM')

        call timer_start(t_all)

        allocate(directresults(size(particles,kind=kind_particle)))

        call directforce(particles, [(i,i=1,size(particles,kind=kind_particle))], size(particles,kind=kind_particle), directresults, MPI_COMM_PEPC)
        particles(1:size(particles,kind=kind_particle))%results = directresults(1:size(particles,kind=kind_particle))

        deallocate(directresults)

        call timer_stop(t_all)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles(particles, allowcheckpoint)
        use module_pepc_types
        use physvars
        implicit none
        type(t_particle), intent(in) :: particles(:)
        logical, intent(in) :: allowcheckpoint
        logical :: bin, asc, check
        integer, save :: lasttimestep = -1

        ! avoid calling this function several times per timestep
        if (lasttimestep .ne. itime) then

            bin = (idump_binary       > 0)
            if (bin) bin = ((mod(itime, idump_binary        ) == 0) .or. (itime == nt))
            asc = (idump              > 0)
            if (asc) asc = ((mod(itime, idump               ) == 0) .or. (itime == nt))
            check = (idump_checkpoint > 0) .and. (allowcheckpoint)
            if (check) check = ((mod(itime, idump_checkpoint) == 0) .or. (itime == nt))

            call write_particles_type(particles, bin, asc, check , itime==nt)
        endif

        lasttimestep = itime

    end subroutine write_particles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles_type(particles, binary, ascii, mpiio, final)
        use module_pepc_types
        use physvars
        use module_checkpoint
        use module_namelist
        implicit none
        include 'mpif.h'
        type(t_particle), intent(in) :: particles(:)
        logical, intent(in) :: binary, ascii, mpiio, final
        character(255) :: filename

        if (binary .or. ascii .or. mpiio) then

            !!! write particle date as a binary file
            if (binary) then
                call write_particles_binary(my_rank, itime, particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

            !!! write particle date as a text file
            if (ascii) then
                call write_particles_ascii(my_rank, itime, particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

            !!! write particle checkpoint data using mpi-io
            if (mpiio) then
                call write_particles_mpiio(MPI_COMM_WORLD, itime, npart_total, particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

        endif

    end subroutine write_particles_type


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_particles(particles, itime_in_)
        use module_pepc_types
        implicit none
        integer, intent(in) :: itime_in_
        type(t_particle), allocatable, intent(out) :: particles(:)

        call read_particles_type(particles, itime_in_, .false., .false., .true.)

    end subroutine read_particles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function vtk_step_of_step(step) result(vtk_step)
      use module_vtk
      use physvars
      implicit none

      integer, intent(in) :: step

      if (step .eq. 0) then
        vtk_step = VTK_STEP_FIRST
      else if (trun+idump_vtk*dt>tend ) then
        vtk_step = VTK_STEP_LAST
      else
        vtk_step = VTK_STEP_NORMAL
      endif
    end function vtk_step_of_step


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles_vtk(p)
      use module_pepc_types
      use module_vtk_helpers
      use physvars
      use module_units
      implicit none

      include 'mpif.h'

      type(t_particle), intent(in) :: p(:)

      integer :: vtk_step

      vtk_step = vtk_step_of_step(itime)
      call vtk_write_particles("particles", MPI_COMM_WORLD, itime, trun*unit_t0_in_fs, vtk_step, p, vtk_results)

      contains

      subroutine vtk_results(d, r, vtkf)
        use module_vtk
        use module_interaction_specific_types
        implicit none

        type(t_particle_data), intent(in) :: d(:)
        type(t_particle_results), intent(in) :: r(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf

        call vtk_write_particle_data_results(d, r, vtkf)
      end subroutine
    end subroutine write_particles_vtk


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_particles_type(particles, itime_in_, binary, ascii, mpiio)
        use physvars
        use module_checkpoint
        use module_namelist
        use module_prepare
        use module_pepc_types
        implicit none
        include 'mpif.h'
        type(t_particle), allocatable, intent(out) :: particles(:)
        logical, intent(in) :: binary, ascii, mpiio
        integer, intent(in) :: itime_in_
        character(255) :: filename

        if (binary .or. ascii .or. mpiio) then

            !!! read particle data as a binary file
            if (binary) write(*,*) "read_particles(): binary mode unsupported" !call read_particles_binary(my_rank, itime, dp)

            !!! read particle date as a text file
            if (ascii)  write(*,*) "read_particles(): ascii mode unsupported" !call read_particles_ascii(my_rank, itime, dp)

            !!! read particle checkpoint data using mpi-io
            if (mpiio) then
                call read_particles_mpiio(itime_in_, MPI_COMM_WORLD, itime, npart_total, particles, filename)
                call read_frontend_parameters_from_file(filename)
            endif

            call frontend_prepare()

        endif

    end subroutine read_particles_type
      
end module module_diagnostics
