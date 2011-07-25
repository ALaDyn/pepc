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

        subroutine vtk_field_on_grid(filename, step, tsim, vtk_step, nx, ny, nz, xcoords, ycoords, zcoords, scalarvalues, scalarname, vectorvalues, vectorname)
          use module_vtk
          implicit none
          character(*), intent(in) :: filename, scalarname, vectorname
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          integer, intent(in) :: nx, ny, nz
          real*8, intent(in) :: xcoords(nx), ycoords(ny), zcoords(nz)
          real*8, intent(in) :: scalarvalues(nx, ny, nz), vectorvalues(nx,ny,nz,3)

          type(vtkfile_rectilinear_grid) :: vtk

            call vtk%create(trim(filename), step, tsim, vtk_step)
              call vtk%write_headers([0, 0, 0], [nx, ny, nz])
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
