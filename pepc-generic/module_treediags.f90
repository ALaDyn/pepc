!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Helper functions for checkpointing and restarting purposes
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_treediags
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

      public write_branches_to_vtk
      public write_spacecurve_to_vtk
      public draw_lists
      public draw_map
      public draw_domains
      public draw_tree2d
      public draw2d
      public draw2d_hash

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes the tree branches structure into a vtk file
        !> pepc_fields must have been called with no_dealloc=.true. before
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine write_branches_to_vtk(step, tsim, vtk_step)
          use treevars
          use module_vtk
          use module_spacefilling
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          integer :: key2addr

          integer :: i,j, baddr, bnode
          integer*8 :: bkey
          real*8 :: bcocx(nbranch_sum),bcocy(nbranch_sum),bcocz(nbranch_sum), bsize, bq(nbranch_sum)
          real*8, dimension(nbranch_sum*8) :: bcornersx, bcornersy, bcornersz
          integer, dimension(nbranch_sum*8) :: bcornersidx
          integer, dimension(nbranch_sum) :: bcornersoffsets, bcornerstypes, bowner, blevel
          real*8 :: bx, by, bz

          real, parameter, dimension(3,8) :: box_shift = reshape([ 0., 0., 0., &
                                                                       0., 0., 1., &
                                                                       0., 1., 0., &
                                                                       0., 1., 1., &
                                                                       1., 0., 0., &
                                                                       1., 0., 1., &
                                                                       1., 1., 0., &
                                                                       1., 1., 1. ], shape(box_shift))
          real*8, dimension(3) :: bshift
          type(vtkfile_unstructured_grid) :: vtk

          if (me .ne. 0) return

          if (.not. (allocated(htable) .and. allocated(branch_key))) then
            write(*,*) 'write_branches_to_vtk(): pepc_fields() must have been called with no_dealloc=.true. before'
            return
          endif

          do i = 1,nbranch_sum
            bkey      = branch_key(i)
            baddr     = key2addr(bkey, "write_branches_to_vtk")
            bnode     = htable(baddr)%node
            bowner(i) = htable(baddr)%owner
            blevel(i) = int(log( 1._8*bkey )/log(8._8))
            bsize     = boxsize/2**blevel(i)
            !write(*,'(O10, Z8, I12, I8, I8, 1G12.3, " | ", 3G12.3)') bkey, baddr, bnode, bowner, blevel, bsize, bcoc
            ! prepare voxel data structure
            bcornerstypes(i)   = VTK_VOXEL
            bcornersoffsets(i) = 8*i
            bq(i)              = charge(bnode)
            bcocx(i)           = xcoc(bnode)
            bcocy(i)           = ycoc(bnode)
            bcocz(i)           = zcoc(bnode)
            ! compute real center coordinate
            call key_to_coord(bkey, bx, by, bz)

            do j=1,8
              bcornersidx(8*(i-1)+j) = 8*(i-1)+j - 1
              bshift(1:3) = box_shift(1:3,j) * bsize
              bcornersx(8*(i-1)+j)   = bx + bshift(1)
              bcornersy(8*(i-1)+j)   = by + bshift(2)
              bcornersz(8*(i-1)+j)   = bz + bshift(3)
            end do
          end do


            call vtk%create("branches", step, tsim, vtk_step)
              call vtk%write_headers(nbranch_sum*8, nbranch_sum)
                call vtk%startpoints()
                  call vtk%write_data_array("corners", 8*nbranch_sum, bcornersx, bcornersy, bcornersz)
                call vtk%finishpoints()
                call vtk%startpointdata()
                  ! no point data here
                call vtk%finishpointdata()
                call vtk%startcells()
                  call vtk%write_data_array("connectivity", nbranch_sum*8, bcornersidx)
                  call vtk%write_data_array("offsets", nbranch_sum, bcornersoffsets)
                  call vtk%write_data_array("types", nbranch_sum, bcornerstypes)
                call vtk%finishcells()
                call vtk%startcelldata()
                  call vtk%write_data_array("processor", nbranch_sum, bowner)
                  call vtk%write_data_array("key", nbranch_sum, branch_key)
                  call vtk%write_data_array("level", nbranch_sum, blevel)
                  call vtk%write_data_array("center_of_charge", nbranch_sum, bcocx, bcocy, bcocz)
                  call vtk%write_data_array("total_charge", nbranch_sum, bq)
                call vtk%finishcelldata()
              call vtk%write_final()
            call vtk%close()

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes the space filling curve into a parallel set of vtk files
        !> pepc_fields must have been called with no_dealloc=.true. before
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine write_spacecurve_to_vtk(step, tsim, vtk_step)
          use treevars
          use module_vtk
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(vtkfile_unstructured_grid) :: vtk
          integer :: i

          if (.not. (allocated(x) .and. allocated(y) .and. allocated(z))) then
            if (me .eq. 0) write(*,*) 'write_spacecurve_to_vtk(): pepc_fields() must have been called with no_dealloc=.true. before'
            return
          endif

            call vtk%create_parallel("spacecurve", step, me, num_pe, tsim, vtk_step)
              call vtk%write_headers(npp, 1)
                call vtk%startpoints()
                  call vtk%write_data_array("xyz", npp, x, y, z)
                call vtk%finishpoints()
                call vtk%startpointdata()
                  ! no point data here
                call vtk%finishpointdata()
                call vtk%startcells()
                  call vtk%write_data_array("connectivity", npp, [(i,i=0, npp-1)])
                  call vtk%write_data_array("offsets", npp)
                  call vtk%write_data_array("types", VTK_POLY_LINE)
                call vtk%finishcells()
                call vtk%startcelldata()
                  call vtk%write_data_array("processor", me)
                call vtk%finishcelldata()
              call vtk%write_final()
            call vtk%close()
        end subroutine


end module module_treediags
