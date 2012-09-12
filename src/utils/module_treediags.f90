! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

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
      public write_interaction_partners_to_vtk

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
        subroutine write_interaction_partners_to_vtk(step, label,tsim, vtk_step)
          use treevars
          use module_vtk
          use module_spacefilling
          use module_htable
          use module_interaction_specific
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          integer, intent(in) :: label
          real*8, intent(in) :: tsim

          integer :: i,j, baddr, bnode
          integer*8 :: bkey
          real*8 :: bcocx(no_interaction_partners(label)),bcocy(no_interaction_partners(label)),bcocz(no_interaction_partners(label)), bsize(3), bq(no_interaction_partners(label))
          real*8, dimension(no_interaction_partners(label)*8) :: bcornersx, bcornersy, bcornersz
          integer, dimension(no_interaction_partners(label)*8) :: bcornersidx
          integer, dimension(no_interaction_partners(label)) :: bcornersoffsets, bcornerstypes, bowner, blevel
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
          type(vtkfile_unstructured_grid) :: vtk,vtk2
          character(255) :: dateiname

          if (me .ne. 0) return

          if (.not. (allocated(htable) .and. allocated(branch_key) .and. allocated(tree_nodes))) then
            write(*,*) 'write_branches_to_vtk(): pepc_fields() must have been called with no_dealloc=.true. before'
            return
          endif

          do i = 1,no_interaction_partners(label)
            bkey      = interaction_keylist(label,i)
            baddr     = key2addr(bkey, "write_interaction_partners_to_vtk")
            bnode     = htable(baddr)%node
            bowner(i) = htable(baddr)%owner
            blevel(i) = level_from_key(bkey)
            bsize     = boxsize/2**blevel(i)
            !write(*,'(O10, Z8, I12, I8, I8, 1G12.3, " | ", 3G12.3)') bkey, baddr, bnode, bowner, blevel, bsize, bcoc
            ! prepare voxel data structure
            bcornerstypes(i)   = VTK_VOXEL
            bcornersoffsets(i) = 8*i
            bq(i)              = tree_nodes(bnode)%charge
            bcocx(i)           = tree_nodes(bnode)%coc(1)+interaction_vbox(label,i,1)
            bcocy(i)           = tree_nodes(bnode)%coc(2)+interaction_vbox(label,i,2)
            bcocz(i)           = tree_nodes(bnode)%coc(3)+interaction_vbox(label,i,3)
            ! compute real center coordinate
            call key_to_coord(bkey, bx, by, bz)
            bx=bx+interaction_vbox(label,i,1)
            by=by+interaction_vbox(label,i,2)
            bz=bz+interaction_vbox(label,i,3)

            do j=1,8
              bcornersidx(8*(i-1)+j) = 8*(i-1)+j - 1
              bshift(1:3) = box_shift(1:3,j) * bsize(1:3)
              bcornersx(8*(i-1)+j)   = bx + bshift(1)
              bcornersy(8*(i-1)+j)   = by + bshift(2)
              bcornersz(8*(i-1)+j)   = bz + bshift(3)
            end do
          end do

        write(dateiname,"(a,i3.3)") "int_partner_coc",label
        call vtk2%create(trim(dateiname), step, tsim, vtk_step)
        call vtk2%write_headers(no_interaction_partners(label), 0)
          call vtk2%startpoints()
          call vtk2%write_data_array("xyz", no_interaction_partners(label), bcocx, bcocy, bcocz)
        call vtk2%finishpoints()
        call vtk2%startpointdata()
          call vtk2%write_data_array("level", no_interaction_partners(label), blevel)
          call vtk2%write_data_array("charge", no_interaction_partners(label), bq)
        call vtk2%finishpointdata()
        call vtk2%dont_write_cells()
        call vtk2%write_final()
        call vtk2%close()

            write(dateiname,"(a,i3.3)") "int_partner_box",label
            call vtk%create(trim(dateiname), step, tsim, vtk_step)
              call vtk%write_headers(no_interaction_partners(label)*8, no_interaction_partners(label))
                call vtk%startpoints()
                  call vtk%write_data_array("corners", 8*no_interaction_partners(label), bcornersx, bcornersy, bcornersz)
                call vtk%finishpoints()
                call vtk%startpointdata()
                  ! no point data here
                call vtk%finishpointdata()
                call vtk%startcells()
                  call vtk%write_data_array("connectivity", no_interaction_partners(label)*8, bcornersidx)
                  call vtk%write_data_array("offsets", no_interaction_partners(label), bcornersoffsets)
                  call vtk%write_data_array("types", no_interaction_partners(label), bcornerstypes)
                call vtk%finishcells()
                call vtk%startcelldata()
                  call vtk%write_data_array("processor", no_interaction_partners(label), bowner)
                  call vtk%write_data_array("key", no_interaction_partners(label), branch_key)
                  call vtk%write_data_array("level", no_interaction_partners(label), blevel)
                  call vtk%write_data_array("center_of_charge", no_interaction_partners(label), bcocx, bcocy, bcocz)
                  call vtk%write_data_array("total_charge", no_interaction_partners(label), bq)
                call vtk%finishcelldata()
              call vtk%write_final()
            call vtk%close()

        end subroutine



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
          use module_htable
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim

          integer :: i,j, baddr, bnode
          integer*8 :: bkey
          real*8 :: bcocx(nbranch_sum),bcocy(nbranch_sum),bcocz(nbranch_sum), bsize(3), bq(nbranch_sum)
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

          if (.not. (allocated(htable) .and. allocated(branch_key) .and. allocated(tree_nodes))) then
            write(*,*) 'write_branches_to_vtk(): pepc_fields() must have been called with no_dealloc=.true. before'
            return
          endif

          do i = 1,nbranch_sum
            bkey      = branch_key(i)
            baddr     = key2addr(bkey, "write_branches_to_vtk")
            bnode     = htable(baddr)%node
            bowner(i) = htable(baddr)%owner
            blevel(i) = level_from_key(bkey)
            bsize     = boxsize/2**blevel(i)
            !write(*,'(O10, Z8, I12, I8, I8, 1G12.3, " | ", 3G12.3)') bkey, baddr, bnode, bowner, blevel, bsize, bcoc
            ! prepare voxel data structure
            bcornerstypes(i)   = VTK_VOXEL
            bcornersoffsets(i) = 8*i
            bq(i)              = tree_nodes(bnode)%charge
            bcocx(i)           = tree_nodes(bnode)%coc(1)
            bcocy(i)           = tree_nodes(bnode)%coc(2)
            bcocz(i)           = tree_nodes(bnode)%coc(3)
            ! compute real center coordinate
            call key_to_coord(bkey, bx, by, bz)

            do j=1,8
              bcornersidx(8*(i-1)+j) = 8*(i-1)+j - 1
              bshift(1:3) = box_shift(1:3,j) * bsize(1:3)
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
        subroutine write_spacecurve_to_vtk(step, tsim, vtk_step, particles)
          use treevars, only : npp, me, num_pe
          use module_vtk
          use module_pepc_types
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_particle), intent(in) :: particles(:)
          type(vtkfile_unstructured_grid) :: vtk
          integer :: i

            call vtk%create_parallel("spacecurve", step, me, num_pe, tsim, vtk_step)
              call vtk%write_headers(npp, 1)
                call vtk%startpoints()
                  call vtk%write_data_array("xyz", npp, particles(1:npp)%x(1), particles(1:npp)%x(2), particles(1:npp)%x(3))
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
