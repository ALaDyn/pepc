!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_vtk
      implicit none

      character(6), parameter :: subfolder = "./vtk/"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      type vtkfile
        private
          character(40) :: filename
          integer :: filehandle = 97
          integer :: filehandle_par = 98
          logical :: parallel = .false.
          character(12) :: byte_order = "LittleEndian"
          character(3) :: version = "0.1"
          character(5) :: format = "ascii"
          integer :: my_rank
          integer :: num_pe


        contains
          procedure :: create => vtkfile_create ! filename
          procedure :: create_parallel => vtkfile_create_parallel ! filename, mpi_comm, rank, num_pe --> rank 0 writes .pvtX-file
          procedure :: close => vtkfile_close
          procedure :: write_data_array_header => vtkfile_write_data_array_header

          procedure :: write_data_array_Real4_1  => vtkfile_write_data_array_Real4_1
          procedure :: write_data_array_Real4_3  => vtkfile_write_data_array_Real4_3
          procedure :: write_data_array_Real8_1  => vtkfile_write_data_array_Real8_1
          procedure :: write_data_array_Real8_3  => vtkfile_write_data_array_Real8_3

          procedure :: write_data_array_Int4_1  => vtkfile_write_data_array_Int4_1
          procedure :: write_data_array_Int4_3  => vtkfile_write_data_array_Int4_3
          procedure :: write_data_array_Int8_1  => vtkfile_write_data_array_Int8_1
          procedure :: write_data_array_Int8_3  => vtkfile_write_data_array_Int4_3

          generic :: write_data_array => write_data_array_Real4_1, & ! name, one-dim real*4, number of entries
                                            write_data_array_Real4_3,  & ! name, three-dim real*4 as three separate arrays, number of entries
                                            write_data_array_Real8_1,  & ! ...
                                            write_data_array_Real8_3,  &
                                            write_data_array_Int4_1,   &
                                            write_data_array_Int4_3,   &
                                            write_data_array_Int8_1,   &
                                            write_data_array_Int8_3
      end type vtkfile


      type, extends(vtkfile) :: vtkfile_unstructured_grid
        contains
          procedure :: write_headers => vtkfile_unstructured_grid_write_headers ! number of particles, writes anything incl. <Piece>
          procedure :: startpoints => vtkfile_unstructured_grid_startpoints
          procedure :: finishpoints => vtkfile_unstructured_grid_finishpoints
          procedure :: startpointdata => vtkfile_unstructured_grid_startpointdata
          procedure :: finishpointdata => vtkfile_unstructured_grid_finishpointdata
          procedure :: write_final => vtkfile_unstructured_grid_write_final ! writes anything from/incl. <Cells>
      end type vtkfile_unstructured_grid


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains


      subroutine vtkfile_create(vtk, filename_)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: filename_
        call vtk%create_parallel(filename_, 0, 0)
      end subroutine vtkfile_create


      subroutine vtkfile_create_parallel(vtk, filename_, my_rank_, num_pe_)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: filename_
        character(50) :: fn
        character(6) :: tmp
        integer :: my_rank_, num_pe_

        vtk%num_pe   = num_pe_
        vtk%my_rank  = my_rank_
        vtk%filename = filename_

        if (vtk%num_pe>0) then
          vtk%parallel = (vtk%my_rank == 0)
          write(tmp,'(I6.6)') vtk%my_rank
          fn = subfolder//trim(vtk%filename)//"."//tmp//".vtu"
        else
          fn = subfolder//trim(vtk%filename)//".vtu"
        endif

        open(vtk%filehandle, file=fn)

        if (vtk%parallel) then
          open(vtk%filehandle_par, file=subfolder//filename_//".pvtu")
        endif
      end subroutine vtkfile_create_parallel


      subroutine vtkfile_close(vtk)
        implicit none
        class(vtkfile) :: vtk
        close(vtk%filehandle)
        if (vtk%parallel) close(vtk%filehandle_par)
     end subroutine vtkfile_close


     subroutine vtkfile_write_data_array_header(vtk, name, number_of_components, type)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name, type
        integer :: number_of_components

        write(vtk%filehandle, '("<DataArray Name=""",a,""" NumberOfComponents=""", I0, """ type=""", a ,""" format=""", a ,""">")') &
                 name, number_of_components, type, vtk%format
        if (vtk%parallel) then
          write(vtk%filehandle_par, '("<DataArray Name=""",a,""" NumberOfComponents=""", I0, """ type=""", a ,""" format=""", a ,"""/>")') &
                 name, number_of_components, type, vtk%format
        end if
     end subroutine


     subroutine vtkfile_write_data_array_Real4_1(vtk, name, ndata, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*4 :: data(:)
        call vtk%write_data_array_header(name, 1, "Float32")
        do i=1,ndata
          write(vtk%filehandle, '(G14.6)') data(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real4_1


     subroutine vtkfile_write_data_array_Real4_3(vtk, name, ndata, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*4 :: data1(:), data2(:), data3(:)
        call vtk%write_data_array_header(name, 1, "Float32")
        do i=1,ndata
          write(vtk%filehandle, '(3G14.6)') data1(i), data2(i), data3(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real4_3


     subroutine vtkfile_write_data_array_Real8_1(vtk, name, ndata, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*8 :: data(:)
        call vtk%write_data_array_header(name, 1, "Float64")
        do i=1,ndata
          write(vtk%filehandle, '(G14.6)') data(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real8_1


     subroutine vtkfile_write_data_array_Real8_3(vtk, name, ndata, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        real*8 :: data1(:), data2(:), data3(:)
        call vtk%write_data_array_header(name, 3, "Float64")
        do i=1,ndata
          write(vtk%filehandle, '(3G14.6)') data1(i), data2(i), data3(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Real8_3


     subroutine vtkfile_write_data_array_Int4_1(vtk, name, ndata, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*4 :: data(:)
        call vtk%write_data_array_header(name, 1, "Int32")
        do i=1,ndata
          write(vtk%filehandle, '(I20)') data(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int4_1


     subroutine vtkfile_write_data_array_Int4_3(vtk, name, ndata, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*4 :: data1(:), data2(:), data3(:)
        call vtk%write_data_array_header(name, 1, "Int32")
        do i=1,ndata
          write(vtk%filehandle, '(3I20)') data1(i), data2(i), data3(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int4_3


     subroutine vtkfile_write_data_array_Int8_1(vtk, name, ndata, data)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*8 :: data(:)
        call vtk%write_data_array_header(name, 1, "Int64")
        do i=1,ndata
          write(vtk%filehandle, '(I20)') data(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int8_1


     subroutine vtkfile_write_data_array_Int8_3(vtk, name, ndata, data1, data2, data3)
        implicit none
        class(vtkfile) :: vtk
        character(*) :: name
        integer :: ndata, i
        integer*8 :: data1(:), data2(:), data3(:)
        call vtk%write_data_array_header(name, 3, "Int64")
        do i=1,ndata
          write(vtk%filehandle, '(3I20)') data1(i), data2(i), data3(i)
        end do
        write(vtk%filehandle, '("</DataArray>")')
     end subroutine vtkfile_write_data_array_Int8_3


     subroutine vtkfile_unstructured_grid_write_headers(vtk, npart)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        integer :: npart

        write(vtk%filehandle, '("<VTKFile type=""UnstructuredGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, vtk%byte_order
        write(vtk%filehandle, '("<UnstructuredGrid GhostLevel=""0"">")')
        write(vtk%filehandle, '("<Piece NumberOfPoints=""", I0, """ NumberOfCells=""0"">")') npart

        if (vtk%parallel) then
          write(vtk%filehandle_par, '("<VTKFile type=""PUnstructuredGrid"" version=""", a, """ byte_order=""", a, """>")') vtk%version, vtk%byte_order
          write(vtk%filehandle_par, '("<PUnstructuredGrid GhostLevel=""0"">")')
        endif
     end subroutine vtkfile_unstructured_grid_write_headers


     subroutine vtkfile_unstructured_grid_write_final(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk
        integer :: i
        character(6) :: tmp
        character(50) :: fn

          write(vtk%filehandle, '("<Cells>")')
          write(vtk%filehandle, '("<DataArray type=""Int32"" Name=""connectivity"" />")')
          write(vtk%filehandle, '("<DataArray type=""Int32"" Name=""offsets"" />")')
          write(vtk%filehandle, '("<DataArray type=""UInt8"" Name=""types"" />")')
          write(vtk%filehandle, '("</Cells>")')
          write(vtk%filehandle, '("<CellData>")')
          write(vtk%filehandle, '("</CellData>")')
          write(vtk%filehandle, '("</Piece>")')
          write(vtk%filehandle, '("</UnstructuredGrid>")')
          write(vtk%filehandle, '("</VTKFile>")')

          if (vtk%parallel) then
            write(vtk%filehandle_par, '("<PCellData>")')
            write(vtk%filehandle_par, '("</PCellData>")')

            do i = 0,vtk%num_pe-1
              write(tmp,'(I6.6)') i
              fn = trim(vtk%filename)//"."//tmp//".vtu"
              write(vtk%filehandle_par, '("<Piece Source=""", a, """/>")') trim(fn)
            end do

            write(vtk%filehandle_par, '("</PUnstructuredGrid>")')
            write(vtk%filehandle_par, '("</VTKFile>")')
          endif
     end subroutine vtkfile_unstructured_grid_write_final


     subroutine vtkfile_unstructured_grid_startpoints(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("<Points>")')
        if (vtk%parallel) write(vtk%filehandle_par, '("<PPoints>")')
     end subroutine vtkfile_unstructured_grid_startpoints


     subroutine vtkfile_unstructured_grid_finishpoints(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("</Points>")')
        if (vtk%parallel) write(vtk%filehandle_par, '("</PPoints>")')
     end subroutine vtkfile_unstructured_grid_finishpoints


     subroutine vtkfile_unstructured_grid_startpointdata(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("<PointData>")')
        if (vtk%parallel) write(vtk%filehandle_par, '("<PPointData>")')
     end subroutine vtkfile_unstructured_grid_startpointdata


     subroutine vtkfile_unstructured_grid_finishpointdata(vtk)
        implicit none
        class(vtkfile_unstructured_grid) :: vtk

        write(vtk%filehandle, '("</PointData>")')
        if (vtk%parallel) write(vtk%filehandle_par, '("</PPointData>")')
     end subroutine vtkfile_unstructured_grid_finishpointdata

end module module_vtk
