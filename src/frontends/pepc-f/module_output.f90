!!!!!!!!!!!!!!!!!!!!
!! Output module
!!
!! Enthaelt Methoden fuer den Output
!!!!!!!!!!!!!!!!!!!!

MODULE output
    use helper
    use variables


    implicit none

    CONTAINS

!===============================================================================

    SUBROUTINE timing_output(integrator,particlehandling,pepc_grow,pepc_traverse,pepc_rest,output,filehandle)

        implicit none

        integer,intent(in)      :: filehandle
        real(kind=8),intent(in) :: integrator,particlehandling,pepc_grow,output,pepc_traverse,pepc_rest
        real(kind=8)             :: timestep

        timestep=integrator+particlehandling+pepc_grow+output+pepc_traverse+pepc_rest
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in integrator [s], %            : ", integrator,", ",100.*integrator/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in particlehandling [s], %      : ", particlehandling,", ",100.*particlehandling/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in grow_tree routine [s], %     : ", pepc_grow,", ",100.*pepc_grow/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in traverse_tree routine [s], % : ", pepc_traverse,", ",100.*pepc_traverse/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in other pepc routines [s], %   : ", pepc_rest,", ",100.*pepc_rest/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in output routines [s], %       : ", output,", ",100.*output/timestep," %"
        if(root) write(filehandle,'(a,es16.8)') " == total time in timestep [s]           : ", timestep

    END SUBROUTINE timing_output

!===============================================================================

    SUBROUTINE end_of_ts_output(timestep,filehandle)

        implicit none

        integer,intent(in)      :: timestep,filehandle

        if(root) write(filehandle,'(a,i6)') " == finished computing step              : ",timestep
        if(root) write(filehandle,'(a)') " "

    END SUBROUTINE end_of_ts_output

!===============================================================================

    SUBROUTINE recycling_output(i_hits_r,e_hits_r,i_hits_l,e_hits_l,q_r,q_l,filehandle)

        implicit none

        integer,intent(in)      :: i_hits_r,e_hits_r,i_hits_l,e_hits_l,filehandle
        real(kind=8),intent(in) :: q_r,q_l

        if(root) write(filehandle,'(a,i16)')    " == ion hits on right wall               : ", i_hits_r
        if(root) write(filehandle,'(a,i16)')    " == electron hits on right wall          : ", e_hits_r
        if(root) write(filehandle,'(a,es16.8)') " == incoming charge on right wall        : ", (i_hits_r-e_hits_r)*e*fsup
        if(root) write(filehandle,'(a,es16.8)') " == total charge on right wall           : ", q_r
        if(root) write(filehandle,'(a,i16)')    " == ion hits on left wall                : ", i_hits_l
        if(root) write(filehandle,'(a,i16)')    " == electron hits on left wall           : ", e_hits_l
        if(root) write(filehandle,'(a,es16.8)') " == incoming charge on left wall         : ", (i_hits_l-e_hits_l)*e*fsup
        !if(root) write(filehandle,'(a,es18.6)') " == total charge on left wall            : ", q_l

        if (fixed_npp) then
            if (root) write(filehandle,'(a,l16)')   " == refluxing in this timestep           : ", need_to_reflux
            if (need_to_reflux) then
                if(root) write(filehandle,'(a,i16)')    " == influxed ions                        : ", i_hits_l+i_hits_r+new_i_r_last_ts
                if(root) write(filehandle,'(a,i16)')    " == influxed electrons                   : ", e_hits_l+e_hits_r+new_e_r_last_ts
            else
                if(root) write(filehandle,'(a,i16)')    " == influxed ions                        : ", i_hits_l
                if(root) write(filehandle,'(a,i16)')    " == influxed electrons                   : ", e_hits_l
            end if
        else
            if(root) write(filehandle,'(a,i16)')    " == influxed ions                        : ", i_hits_l+tfpp/2
            if(root) write(filehandle,'(a,i16)')    " == influxed electrons                   : ", e_hits_l+tfpp/2
        end if
        if(root) write(filehandle,'(a,i16)')    " == total number of particles            : ", tnp
        if(root) write(filehandle,'(a,i16)')    " == total number of plasma particles     : ", tnpp

    END SUBROUTINE recycling_output

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
    
        time = dt* step

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
        call vtk%write_data_array("species", np, p(:)%data%species)
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        if(root) write(*,'(a,i6)') " == [write_particles] vtk output at timestep",step

    end subroutine write_particles  

END MODULE
