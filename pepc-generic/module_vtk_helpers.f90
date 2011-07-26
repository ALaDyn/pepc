!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates helper functions for outputting different data formats to vtk files
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_vtk_helpers
      implicit none
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public vtk_field_on_grid


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

        subroutine vtk_field_on_grid(filename, step, tsim, vtk_step, globaldims, mydims, xcoords, ycoords, zcoords, &
                          scalarvalues, scalarname, vectorvalues, vectorname, my_rank, num_pe, comm)
          use module_vtk
          implicit none
          character(*), intent(in) :: filename, scalarname, vectorname
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          integer, dimension(2,3), intent(in) :: globaldims, mydims
          real*8, intent(in) :: xcoords(:), ycoords(:), zcoords(:)
          real*8, intent(in) :: scalarvalues(:, :, :), vectorvalues(:,:,:,:)
          integer, intent(in) :: comm, my_rank, num_pe
          integer :: nx, ny, nz

          type(vtkfile_rectilinear_grid) :: vtk

            nx = mydims(2,1) - mydims(1,1) + 1
            ny = mydims(2,2) - mydims(1,2) + 1
            nz = mydims(2,3) - mydims(1,3) + 1

            call vtk%create_parallel(trim(filename), step, my_rank, num_pe, tsim, vtk_step)
              call vtk%set_communicator(comm)
              call vtk%write_headers(globaldims, mydims)
                call vtk%startcoordinates()
                  call vtk%write_data_array("x_coordinate", nx, xcoords)
                  call vtk%write_data_array("y_coordinate", ny, ycoords)
                  call vtk%write_data_array("z_coordinate", nz, zcoords)
                call vtk%finishcoordinates()
                call vtk%startpointdata()
                  call vtk%write_data_array(scalarname, nx, ny, nz, scalarvalues)
                  call vtk%write_data_array(vectorname, nx, ny, nz, vectorvalues(:,:,:,1), vectorvalues(:,:,:,2), vectorvalues(:,:,:,3))
                  ! no point data here
                call vtk%finishpointdata()
                call vtk%startcelldata()
                  ! no cell data here
                call vtk%finishcelldata()
              call vtk%write_final()
            call vtk%close()

        end subroutine



end module module_vtk_helpers
