!!!!!!!!!!!!!!!!!!!!
!! Output module
!!
!! Enthaelt Methoden fuer den Output
!!!!!!!!!!!!!!!!!!!!

MODULE output
    use helper

    implicit none

    CONTAINS


!===============================================================================

    subroutine write_particles(p)
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

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt-1) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(np, 0)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", np, p(:)%x(1), p(:)%x(2), p(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", np, p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
        call vtk%write_data_array("el_field", np, p(:)%results%e(1),p(:)%results%e(2), p(:)%results%e(3))
        call vtk%write_data_array("el_pot", np, p(:)%results%pot)
        call vtk%write_data_array("charge", np, p(:)%data%q)
        call vtk%write_data_array("mass", np, p(:)%data%m)
        call vtk%write_data_array("work", np, p(:)%work)
        call vtk%write_data_array("pelabel", np, p(:)%label)
        call vtk%write_data_array("local index", np, [(i,i=1,np)])
        call vtk%write_data_array("processor", np, p(:)%pid)
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

    end subroutine write_particles  

END MODULE
