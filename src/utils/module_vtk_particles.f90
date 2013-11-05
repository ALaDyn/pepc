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

module module_vtk_particles

! this is a fix to get rid of I/O (type-bound procedures, really) for OMPSs
#ifndef BASE
      implicit none 
 
      integer, public, parameter :: dummy = 0 
end module module_vtk_particles 
#else

  use module_pepc_types
  use module_vtk
  implicit none
  private

  public :: vtk_write_particles
  public :: vtk_write_particles_coulomb_XYZQVM
  public :: vtk_write_particles_coulomb_XYZQVM_helper

  contains

  subroutine vtk_write_particles(fname, mpi_comm, step, time, vtk_step, p, helper_func)
    implicit none

    include 'mpif.h'
    
    character(*), intent(in) :: fname
    integer(kind_default), intent(in) :: mpi_comm
    integer, intent(in) :: step
    real*8, intent(in) :: time
    integer, intent(in) :: vtk_step
    type(t_particle), intent(in) :: p(:)

    interface
      subroutine helper_func(p, vtkf)
        import t_particle, vtkfile_unstructured_grid
        implicit none

        type(t_particle), intent(in) :: p(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func

    integer(kind_particle) :: i, np
    integer(kind_pe) :: mpi_rank, mpi_size
    integer(kind_default) :: ierr
    type(vtkfile_unstructured_grid) :: vtk
    
    call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    np = size(p)

    call vtk%create_parallel(fname, step, mpi_rank, mpi_size, time, vtk_step)
      call vtk%write_headers(np, 0_kind_particle)
      call vtk%startpoints()
        call vtk%write_data_array("xyz", p(:)%x(1), p(:)%x(2), p(:)%x(3))
      call vtk%finishpoints()
      call vtk%startpointdata()
        call vtk%write_data_array("work", p(:)%work)
        call vtk%write_data_array("pelabel", p(:)%label)
        call vtk%write_data_array("local index", [(i,i=1,np)])
        call vtk%write_data_array("processor", p(:)%pid)
        if (present(helper_func)) then; call helper_func(p, vtk); end if
      call vtk%finishpointdata()
      call vtk%dont_write_cells()
      call vtk%write_final()
    call vtk%close()
  end subroutine vtk_write_particles


  subroutine vtk_write_particles_coulomb_XYZQVM(fname, mpi_comm, step, time, vtk_step, p)
    implicit none
    
    character(*), intent(in) :: fname
    integer(kind_default), intent(in) :: mpi_comm
    integer, intent(in) :: step
    real*8, intent(in) :: time
    integer, intent(in) :: vtk_step
    type(t_particle), intent(in) :: p(:)

    call vtk_write_particles(fname, mpi_comm, step, time, vtk_step, p, &
      vtk_write_particles_coulomb_XYZQVM_helper)
  end subroutine vtk_write_particles_coulomb_XYZQVM


  subroutine vtk_write_particles_coulomb_XYZQVM_helper(p, vtkf)
    implicit none

    type(t_particle), intent(in) :: p(:)
    type(vtkfile_unstructured_grid), intent(inout) :: vtkf

    call vtkf%write_data_array("charge", p(:)%data%q)
    call vtkf%write_data_array("velocity", p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
    call vtkf%write_data_array("mass", p(:)%data%m)
    call vtkf%write_data_array("el_pot", p(:)%results%pot)
    call vtkf%write_data_array("el_field", p(:)%results%e(1), p(:)%results%e(2), p(:)%results%e(3))
  end subroutine vtk_write_particles_coulomb_XYZQVM_helper
end module module_vtk_particles
#endif
