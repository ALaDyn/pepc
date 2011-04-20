!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_diagnostics
   implicit none
   private



     public write_particles_to_vtk
     public write_total_momentum



   contains

      subroutine write_particles_to_vtk(step)
        use physvars
        use module_vtk
        implicit none
        integer, intent(in) :: step
        character(16) :: fn
        integer :: i
        type(vtkfile_unstructured_grid) :: vtk

        write(fn, '("particles_", I6.6)') step

        call vtk%create_parallel(fn, my_rank, n_cpu)
          call vtk%write_headers(np_local)
            call vtk%startpoints()
              call vtk%write_data_array("xyz", np_local, x, y, z)
            call vtk%finishpoints()
            call vtk%startpointdata()
              call vtk%write_data_array("velocity", np_local, ux, uy, uz)
              call vtk%write_data_array("charge", np_local, q)
              call vtk%write_data_array("mass", np_local, m)
              call vtk%write_data_array("pelabel", np_local, pelabel)
              call vtk%write_data_array("local index", np_local, [(i,i=1,np_local)])
              call vtk%write_data_array("processor", np_local, [(my_rank,i=1,np_local)])
            call vtk%finishpointdata()
          call vtk%write_final()
        call vtk%close()

      end subroutine write_particles_to_vtk


      subroutine write_total_momentum(itime_, trun_, mom)
        use physvars
        implicit none
        include 'mpif.h'

        integer, intent(in) :: itime_
        real, intent(in) :: trun_
        real*8, intent(out) :: mom(4)
        real*8 :: r(4)
        integer :: p, ierr

        mom = 0.

        do p = 1,np_local
          if (q(p) < 0) then
            r = [ux(p), uy(p), uz(p), 0._8]
            r(4) = sqrt(dot_product(r,r))
            mom = mom + r
          endif
        end do

        if (my_rank == 0) then
          call MPI_REDUCE(MPI_IN_PLACE, mom, 4, MPI_REAL8, MPI_SUM,  0, MPI_COMM_WORLD, ierr )
        else
          call MPI_REDUCE(mom, MPI_IN_PLACE, 4, MPI_REAL8, MPI_SUM,  0, MPI_COMM_WORLD, ierr )
        endif

        if (my_rank == 0) then
          if (itime_ <= 1) then
             open(87, FILE='momentum.dat',STATUS='UNKNOWN', POSITION = 'REWIND')
          else
             open(87, FILE='momentum.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
           endif
           write(87,'(i10,5g25.12)') itime_, trun_, mom
           close(87)
        endif
      end subroutine momentum_dump

end module module_diagnostics
