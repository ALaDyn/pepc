!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_fields
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, public, dimension(3) :: field_dump_min = [ -10., -10., -10.] !< minimum coordinates for fielddump grid
      real*8, public, dimension(3) :: field_dump_max = [ -10., -10., -10.] !< maximum coordinates for fielddump grid
      integer, public, dimension(3) :: field_dump_ncells = [50, 50, 50] !< spatial grid resolution

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

      real*8, allocatable, dimension(:,:,:,:) :: fE
      real*8, allocatable, dimension(:,:,:)   :: fPhi

    contains

      subroutine field_dump(itime)
        implicit none

        integer, intent(in) :: itime

        allocate(fE(field_dump_ncells(1),field_dump_ncells(2),field_dump_ncells(3), 3))
        allocate(fPhi(field_dump_ncells(1),field_dump_ncells(2),field_dump_ncells(3)))

        fE = 0.
        fPhi = 0.

        call field_gather_local
        call field_reduce_global
        call field_calc_laser

        ! dump fE and fPhi to file
        call field_write_vtk

        deallocate(fE, fPhi)

      end subroutine field_dump





      subroutine field_gather_local
        implicit none

      end subroutine




      subroutine field_reduce_global
        implicit none

      end subroutine




      subroutine field_calc_laser
        implicit none

      end subroutine





      subroutine field_write_vtk
        implicit none

!  INTEGER(ap), INTENT(in)::koufem
!       INTEGER(ap):: i
!       CHARACTER file*22
!       WRITE (file,'(A15)')'out/vtk/veloci.'
!
!       IF (kOuFEM.LT.10) THEN
!          WRITE (file(16:16),'(I1)') kOuFEM
!          WRITE (file(17:20),'(A4)') '.vtk'
!       ELSEIF (kOuFEM.LT.100) THEN
!          WRITE (file(16:17),'(I2)') kOuFEM
!          WRITE (file(18:21),'(A4)') '.vtk'
!       ELSE
!          WRITE (file(16:18),'(I3)') kOuFEM
!          WRITE (file(18:22),'(A4)') '.vtk'
!       ENDIF
!
!       print *, file
!
!       OPEN (UNIT=8,FILE=file,STATUS='UNKNOWN')
!       WRITE (8,'(A26)')'# vtk DataFile Version 3.0'
!       WRITE (8,'(A12)')'model R-SWMS'
!       WRITE (8,'(A5)')'ASCII'
!       WRITE (8,'(A24)')'DATASET RECTILINEAR_GRID'
!       WRITE (8,'(A11,1X,3I6)')'DIMENSIONS ', nx, ny, nz
!
!       WRITE (8,'(A14,1X,I6,1X,A5)')'X_COORDINATES ',nx, 'float'
!       DO i=1,nx
!          WRITE (8,'(1X,1pE11.4)',advance="no") xgrid(i)
!       END DO
!
!       WRITE (8,'(/A14,1X,I6,1X,A5)')'Y_COORDINATES ',ny, 'float'
!       DO i=1,ny*nx,nx
!          WRITE (8,'(1pE11.4)',advance="no") ygrid(i)
!       END DO
!
!       WRITE (8,'(/A14,1X,I6,1X,A5)')'Z_COORDINATES ',nz, 'float'
!       DO i=1,size(zgrid),nx*ny
!          WRITE (8,'(1pE11.4)',advance="no") zgrid(i)
!       END DO
!
!       WRITE (8,'(/A15,1X,I7)')'POINT_DATA', nPt
!       WRITE (8,'(A24)')'SCALARS velocity float 3'
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!
!       CALL Veloc(hNew)
!       DO 2 i=1,nPt
!        WRITE (8,'(3(1X,1pE11.4))')Vx(i),Vy(i),Vz(i)
!2     END DO
!
!       WRITE (8,'(A30)')'SCALARS pressurehead float '
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!       WRITE (8,'(1X,1pE11.4)',advance="no") hnew
!
!       WRITE (8,'(A30)')'SCALARS wc float '
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!       WRITE (8,'(1X,1pE11.4)',advance="no") theta
!
!       WRITE (8,'(A30)')'SCALARS conc float '
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!       WRITE (8,'(1X,1pE11.4)',advance="no") conc
!
!       WRITE (8,'(A30)')'SCALARS sink float '
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!       WRITE (8,'(1X,1pE11.4)',advance="no") sink
!
!       WRITE (8,'(A30)')'SCALARS csink float '
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!       WRITE (8,'(1X,1pE11.4)',advance="no") csink
!
!       WRITE (8,'(/A15,1X,I7)')'CELL_DATA', nElm
!       WRITE (8,'(A24)')'SCALARS sinkElm float'
!       WRITE (8,'(A20)')'LOOKUP_TABLE default'
!       WRITE (8,'(1X,1pE11.4)',advance="no") sink_cube
!
!       CLOSE (8)


      end subroutine



end module module_fields
