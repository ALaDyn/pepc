!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_fields
      use physvars
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, public, dimension(3) :: field_dump_ncells = [50, 50, 50] !< spatial grid resolution
      real :: threshhold = 0.2
      real*8 :: mincoord(3), maxcoord(3), delta(3)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public field_dump
      public momentum_dump

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable, dimension(:,:,:,:) :: ffieldi, ffielde
      integer, allocatable, dimension(:,:,:) :: npartfi, npartfe

    contains

      subroutine momentum_dump(itime_, trun_)
        use physvars
        implicit none
        include 'mpif.h'

        integer, intent(in) :: itime_
        real, intent(in) :: trun_
        real*8 :: mom(4), r(4)
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



      subroutine field_dump(itime)
        implicit none
        include 'mpif.h'
        integer :: ierr

        integer, intent(in) :: itime
        real*8 :: deltathresh(3)
        character(18) :: filename

        mincoord(1) = minval(x(1:np_local))
        mincoord(2) = minval(y(1:np_local))
        mincoord(3) = minval(z(1:np_local))
        maxcoord(1) = maxval(x(1:np_local))
        maxcoord(2) = maxval(y(1:np_local))
        maxcoord(3) = maxval(z(1:np_local))
        call MPI_ALLREDUCE(MPI_IN_PLACE, mincoord, 3, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(MPI_IN_PLACE, maxcoord, 3, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
        deltathresh = threshhold*(maxcoord-mincoord)
        maxcoord = maxcoord + deltathresh
        mincoord = mincoord - deltathresh
        delta = (maxcoord - mincoord)/(field_dump_ncells-1)

        if (my_rank == 0) then
          write(*,*) "DUMPING VTK FILE"
          write(filename,'("field_",i8.8,".vtk")') itime
         open(87, FILE=filename,STATUS='UNKNOWN')
          call field_write_vtk_header(87)
        end if

        call field_dump_data_v(87, 'Exyz', Ex, Ey, Ez,.true.)
        call field_dump_data_s(87, 'Pot', Pot,.true.)
        call field_dump_data_v(87, 'v', ux, uy, uz,.true.)
        call field_dump_data_s(87, 'q', q,.false.)

        close(87)

      end subroutine field_dump


      subroutine field_dump_data_v(ifile, name, xvals, yvals, zvals, renorm)
        implicit none
        real*8, intent(in), dimension(np_local) :: xvals, yvals, zvals
        character(*), intent(in) :: name
        integer, intent(in) :: ifile
        logical, intent(in) :: renorm

        call allocfields(3)

        call field_gather_local_v(xvals, yvals, zvals)
        call field_reduce_global(renorm)
        if (my_rank == 0) call field_write_to_vtk(ifile, name)

        call deallocfields()

      end subroutine

      subroutine field_dump_data_s(ifile, name, vals, renorm)
        implicit none
        real*8, intent(in), dimension(np_local) :: vals
        character(*), intent(in) :: name
        integer, intent(in) :: ifile
        logical, intent(in) :: renorm

        call allocfields(1)

        call field_gather_local_s(vals)
        call field_reduce_global(renorm)
        if (my_rank == 0) call field_write_to_vtk(ifile, name)

        call deallocfields()

      end subroutine


      subroutine allocfields(fdim)
        implicit none
        integer, intent(in) :: fdim

        allocate(ffieldi(field_dump_ncells(1),field_dump_ncells(2),field_dump_ncells(3), fdim))
        allocate(npartfi(field_dump_ncells(1),field_dump_ncells(2),field_dump_ncells(3)))
        allocate(ffielde(field_dump_ncells(1),field_dump_ncells(2),field_dump_ncells(3), fdim))
        allocate(npartfe(field_dump_ncells(1),field_dump_ncells(2),field_dump_ncells(3)))

        ffieldi = 0.
        npartfi = 0
        ffielde = 0.
        npartfe = 0

      end subroutine

      subroutine deallocfields
        implicit none

        deallocate(ffieldi, npartfi)
        deallocate(ffielde, npartfe)
      end subroutine



      subroutine field_gather_local_v(xvals, yvals, zvals)
        use physvars
        implicit none
        real*8, intent(in), dimension(np_local) :: xvals, yvals, zvals
        integer :: p
        integer :: cell(3)
        real*8 :: cellr(3),coord(3)

        do p = 1,np_local
          coord = [x(p), y(p), z(p)]
          cellr = (coord - mincoord)/delta+1
          cell  = nint(cellr)

          if (q(p) > 0) then
            ffieldi(    cell(1),cell(2), cell(3),:) = ffieldi(    cell(1),cell(2), cell(3),:) + [xvals(p), yvals(p), zvals(p)]
            npartfi(cell(1),cell(2), cell(3)  )     = npartfi(cell(1),cell(2), cell(3)  ) + 1
          else
            ffielde(    cell(1),cell(2), cell(3),:) = ffielde(    cell(1),cell(2), cell(3),:) + [xvals(p), yvals(p), zvals(p)]
            npartfe(cell(1),cell(2), cell(3)  )     = npartfe(cell(1),cell(2), cell(3)  ) + 1
          end if
        end do
      end subroutine


      subroutine field_gather_local_s(vals)
        use physvars
        implicit none
        real*8, intent(in), dimension(np_local) :: vals
        integer :: p
        integer :: cell(3)
        real*8 :: cellr(3),coord(3)

        do p = 1,np_local
          coord = [x(p), y(p), z(p)]
          cellr = (coord - mincoord)/delta+1
          cell  = nint(cellr)

          if (q(p) > 0) then
            ffieldi(    cell(1),cell(2), cell(3),:) = ffieldi(    cell(1),cell(2), cell(3),:) + vals(p)
            npartfi(cell(1),cell(2), cell(3)  )     = npartfi(cell(1),cell(2), cell(3)  ) + 1
          else
            ffielde(    cell(1),cell(2), cell(3),:) = ffielde(    cell(1),cell(2), cell(3),:) + vals(p)
            npartfe(cell(1),cell(2), cell(3)  )     = npartfe(cell(1),cell(2), cell(3)  ) + 1
          end if
        end do
      end subroutine




      subroutine field_reduce_global(renorm)
        implicit none
        include 'mpif.h'
        logical, intent(in) :: renorm
        integer:: fdim
        integer :: numdata, ierr, i

        numdata = product(field_dump_ncells)
        fdim = size(ffieldi(1,1,1,:))

        if (my_rank == 0) then
          call MPI_REDUCE(MPI_IN_PLACE,ffieldi, fdim*numdata, MPI_REAL8,   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(MPI_IN_PLACE,npartfi,      numdata, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(MPI_IN_PLACE,ffielde, fdim*numdata, MPI_REAL8,   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(MPI_IN_PLACE,npartfe,      numdata, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        else
          call MPI_REDUCE(ffieldi,MPI_IN_PLACE, fdim*numdata, MPI_REAL8,   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(npartfi,MPI_IN_PLACE,      numdata, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(ffielde,MPI_IN_PLACE, fdim*numdata, MPI_REAL8,   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          call MPI_REDUCE(npartfe,MPI_IN_PLACE,      numdata, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        endif

        ! normalize by number of contributing particles
        if (renorm .and. (my_rank == 0)) then
          do i=1,fdim
            where (npartfi .ne. 0)
              ffieldi(:,:,:,i) = ffieldi(:,:,:,i) / npartfi
            end where
            where (npartfe .ne. 0)
              ffielde(:,:,:,i) = ffielde(:,:,:,i) / npartfe
            end where
          end do
        end if

      end subroutine



      subroutine field_write_to_vtk(ifile, name)
        implicit none
        integer, intent(in) :: ifile
        character(*), intent(in) :: name
        integer :: fdim, i,j,k,l

        fdim = size(ffieldi(1,1,1,:))

        write (ifile,'(/A,A,A,I3)')'SCALARS ', name, '_ions float ', fdim
        write (ifile,'(A)')'LOOKUP_TABLE default'

        do k = 1,field_dump_ncells(3)
          do j = 1,field_dump_ncells(2)
            do i = 1,field_dump_ncells(1)
              do l=1,fdim
                write (ifile,'((1X,1pE11.4))',advance='no') ffieldi(i,j,k,l)
              end do
              write(ifile,'(A)') ""
            end do
          end do
        end do


        write (ifile,'(/A,A,A,I3)')'SCALARS ', name, '_electrons float ', fdim
        write (ifile,'(A)')'LOOKUP_TABLE default'

        do k = 1,field_dump_ncells(3)
          do j = 1,field_dump_ncells(2)
            do i = 1,field_dump_ncells(1)
              do l=1,fdim
                write (ifile,'((1X,1pE11.4))',advance='no') ffielde(i,j,k,l)
              end do
              write(ifile,'(A)') ""
            end do
          end do
        end do
      end subroutine




      subroutine field_calc_laser
        implicit none

      end subroutine





      subroutine field_write_vtk_header(ifile)
        implicit none
        integer, intent(in) :: ifile
        integer :: i

        write (ifile,'(A26)') '# vtk DataFile Version 3.0'
        write (ifile,'(A)') 'PEPC field output'
        write (ifile,'(A5)')  'ASCII'
        write (ifile,'(A24)') 'DATASET RECTILINEAR_GRID'
        write (ifile,'(A11,1X,3I6)')'DIMENSIONS ', field_dump_ncells

        write (ifile,'(A14,1X,I6,1X,A5)') 'X_COORDINATES ',field_dump_ncells(1), 'float'
        do i=0,field_dump_ncells(1)-1
           write (ifile,'(1X,1pE11.4)',advance="no") mincoord(1) + 1.*i*delta(1)
        end do

        write (ifile,'(/A14,1X,I6,1X,A5)') 'Y_COORDINATES ',field_dump_ncells(2), 'float'
        do i=0,field_dump_ncells(2)-1
           write (ifile,'(1X,1pE11.4)',advance="no") mincoord(2) + 1.*i*delta(2)
        end do

        write (ifile,'(/A14,1X,I6,1X,A5)') 'Z_COORDINATES ',field_dump_ncells(3), 'float'
        do i=0,field_dump_ncells(3)-1
           write (ifile,'(1X,1pE11.4)',advance="no") mincoord(3) + 1.*i*delta(3)
        end do

        write (ifile,'(//A,1X,I7/)')'POINT_DATA', product(field_dump_ncells)

      end subroutine



end module module_fields
