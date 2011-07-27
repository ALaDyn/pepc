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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable, dimension(:,:,:,:) :: ffieldi, ffielde
      integer, allocatable, dimension(:,:,:) :: npartfi, npartfe

    contains

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

        call field_dump_laser()

        return

        if (my_rank == 0) then
          write(*,*) "DUMPING VTK FILE"
          write(filename,'("field_",i8.8,".vtk")') itime
         open(87, FILE=trim(filename),STATUS='UNKNOWN')
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




      subroutine field_dump_laser
        use module_laser
        use module_units
        use module_vtk_helpers
        use module_vtk
        use physvars, only : itime, nt, trun, my_rank, n_cpu, MPI_COMM_PEPC
        implicit none
        include 'mpif.h'
        real*8, dimension(:,:,:,:), allocatable :: efield
        real*8, dimension(:,:,:), allocatable :: pot
        real*8, allocatable, dimension(:) :: xcoords, ycoords, zcoords
        real*8 :: E_pon(3), B_em(3), Phi_pon
        integer :: i, j, k, vtk_step, ierr
        integer :: globaldims(2,3), mydims(2,3)
        integer :: dims(3), coords(3)
        logical :: periods
        integer :: comm_cart

        globaldims(1,:) = 0
        globaldims(2,:) = field_dump_ncells(:)

        dims = 0
        call MPI_DIMS_CREATE(n_cpu, 3, dims, ierr)
        periods = .false.
        call MPI_CART_CREATE(MPI_COMM_PEPC, 3, dims, periods, .false., comm_cart, ierr)
        call MPI_CART_GET(comm_cart, 3, dims, periods, coords, ierr)

        mydims(2,:) = globaldims(2,:) / dims
        mydims(1,:) = coords * mydims(2,:)
        mydims(2,:) = mydims(2,:) + mydims(1,:)
        ! shift domains if modulo(globaldims(2,:), dims !=0
        mydims(1,:) = mydims(1,:) + min(my_rank, modulo(globaldims(2,:), dims))
        mydims(2,:) = mydims(2,:) + min(my_rank, modulo(globaldims(2,:), dims))
        where (coords < modulo(globaldims(2,:), dims))
          mydims(2,:) = mydims(2,:) + 1
        end where

        allocate(efield(mydims(1,1):mydims(2,1), &
                         mydims(1,2):mydims(2,2), &
                         mydims(1,3):mydims(2,3), 3))
        allocate(pot(mydims(1,1):mydims(2,1), &
                      mydims(1,2):mydims(2,2), &
                      mydims(1,3):mydims(2,3)))
        allocate(xcoords(mydims(1,1):mydims(2,1)))
        allocate(ycoords(mydims(1,2):mydims(2,2)))
        allocate(zcoords(mydims(1,3):mydims(2,3)))

        do i=mydims(1,1),mydims(2,1)
          xcoords(i) = (1.*i-0.5)*delta(1) + mincoord(1)
        end do

        do j=mydims(1,2),mydims(2,2)
          ycoords(j) = (1.*j-0.5)*delta(2) + mincoord(2)
        end do

        do k=mydims(1,3),mydims(2,3)
          zcoords(k) = (1.*k-0.5)*delta(3) + mincoord(3)
        end do

        do k=mydims(1,3),mydims(2,3)
          do j=mydims(1,2),mydims(2,2)
            do i=mydims(1,1),mydims(2,1)
              call force_laser_at(xcoords(i), ycoords(j), zcoords(k), 0._8, E_pon, B_em, Phi_pon)
              efield(i, j, k, 1:3) = E_pon
                 pot(i, j, k)      = Phi_pon
            end do
          end do
        end do

       if (itime == 1) then
         vtk_step = VTK_STEP_FIRST
       else if (itime == nt) then
         vtk_step = VTK_STEP_LAST
       else
         vtk_step = VTK_STEP_NORMAL
       endif

       call vtk_field_on_grid("laser", itime, trun*unit_t0_in_fs, vtk_step, &
                    globaldims, mydims, xcoords, ycoords, zcoords, pot, "phipon", efield, "epon", &
                    my_rank, n_cpu, MPI_COMM_PEPC)

       deallocate(efield)
       deallocate(pot)
       deallocate(xcoords)
       deallocate(ycoords)
       deallocate(zcoords)

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
