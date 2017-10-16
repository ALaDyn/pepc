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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_kinds
  use module_pepc_types
  use variables

  contains


  subroutine init_after_resume()
      
    use module_pepc
    use module_interaction_specific
    use module_checkpoint
    use module_pepc_types

    implicit none
    include 'mpif.h'
      
    vtk=.false.
    cmd_args = COMMAND_ARGUMENT_COUNT()
    if (cmd_args > 1) then
        call GET_COMMAND_ARGUMENT(1, file_in)
        call GET_COMMAND_ARGUMENT(2, file_out)
    end if
    if (file_out=="vtk") vtk=.true.
    call read_particles_mpiio_from_filename(MPI_COMM_WORLD,step,npart,particles,file_in)
    np=size(particles, kind=kind(np))

    if (vtk) then
        filenameh = trim(file_in)//".h"
        open(123, file=trim(filenameh),action='read')
        read(123,NML=pepcf)
        close(123)
    end if
    call pepc_prepare(3_kind_dim)

  end subroutine

!===============================================================================

    subroutine write_particles_vtk(p)
        use module_vtk
        implicit none

        type(t_particle), allocatable, intent(in) :: p(:)

        integer :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time
        real*8 :: ta, tb

        ta = get_time()

        time = dt * step

        vtk_step = VTK_STEP_NORMAL

        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(np, 0)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", p(:)%x(1), p(:)%x(2), p(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
        call vtk%write_data_array("el_field", p(:)%results%e(1),p(:)%results%e(2), p(:)%results%e(3))
        call vtk%write_data_array("el_pot", p(:)%results%pot)
        call vtk%write_data_array("charge", p(:)%data%q)
        call vtk%write_data_array("mass", p(:)%data%m)
        call vtk%write_data_array("work", p(:)%work)
        call vtk%write_data_array("species", p(:)%data%species)
        call vtk%write_data_array("pelabel", p(:)%label)
        call vtk%write_data_array("local index", [(i,i=1,np)])
        !call vtk%write_data_array("processor", p(:)%pid)
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        !if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

    end subroutine write_particles_vtk

!======================================================================================
  real*8 function get_time()
    implicit none
    include 'mpif.h'
    
    get_time = MPI_WTIME()
    
  end function get_time

end module
