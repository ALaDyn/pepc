! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

!>
!> Helper functions for checkpointing and restarting purposes
!>
module module_treediags
      implicit none
      private

      public write_branches_to_vtk
      public write_spacecurve_to_vtk
      public write_interaction_partners_to_vtk

      contains

        !>
        !> Writes a collection of nodes into vtk files, once as boxes, once
        !> as points
        !>
        subroutine write_nodes_to_vtk(step, tsim, vtk_step, num_nodes, &
          node_keys, filename_box, filename_point, node_vbox)
          use treevars
          use module_vtk
          use module_spacefilling
          use module_mirror_boxes
          use module_pepc_types, only: t_tree_node
          use module_tree, only: t_tree, tree_lookup_node_critical
          use module_tree_node
          use module_pepc, only: global_tree
          implicit none

          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          integer, intent(in) :: num_nodes
          integer*8, dimension(:), intent(in) :: node_keys
          character(*), intent(in), optional :: filename_box
          character(*), intent(in), optional :: filename_point
          real*8, dimension(:,:), intent(in), optional :: node_vbox
          
          type(t_tree_node), pointer :: bnode
          type(t_tree), pointer :: t
          integer :: i,j
          integer*8 :: bkey
          real*8 :: bsize(3)
          real*8, dimension(:), allocatable  :: bcocx, bcocy, bcocz, bq, &
                                                bcornersx, bcornersy, bcornersz
          integer, dimension(:), allocatable :: bcornersidx, bcornersoffsets, &
                                                bcornerstypes, bowner, blevel, &
                                                mirror_level
          integer, dimension(:, :), allocatable :: mirror_indices
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
          type(vtkfile_unstructured_grid) :: vtk_box, vtk_point

          ! TODO: generalize!
          t => global_tree

          if (.not. t%comm_env%first) return

          allocate(bcocx(num_nodes), bcocy(num_nodes), bcocz(num_nodes), &
            bq(num_nodes))
          allocate(bcornersx(num_nodes*8), bcornersy(num_nodes*8), &
            bcornersz(num_nodes*8))
          allocate(bcornersidx(num_nodes * 8))
          allocate(bcornersoffsets(num_nodes), bcornerstypes(num_nodes), &
            bowner(num_nodes), blevel(num_nodes), mirror_level(num_nodes))
          allocate(mirror_indices(num_nodes, 3))

          do i = 1, num_nodes
            bkey      = node_keys(i)
            call tree_lookup_node_critical(t, bkey, bnode, 'write_nodes_to_vtk')
            bowner(i) = bnode%owner
            blevel(i) = level_from_key(bkey)
            bsize     = t%bounding_box%boxsize / 2**blevel(i)
            !write(*,'(O10, Z8, I12, I8, I8, 1G12.3, " | ", 3G12.3)') bkey, baddr, bnode, bowner, blevel, bsize, bcoc

            ! prepare voxel data structure
            bcornerstypes(i)   = VTK_VOXEL
            bcornersoffsets(i) = 8*i
            bq(i)              = bnode%interaction_data%charge
            bcocx(i)           = bnode%interaction_data%coc(1)
            bcocy(i)           = bnode%interaction_data%coc(2)
            bcocz(i)           = bnode%interaction_data%coc(3)

            if ( present(node_vbox) ) then
              bcocx(i) = bcocx(i) + node_vbox(i, 1)
              bcocy(i) = bcocy(i) + node_vbox(i, 2)
              bcocz(i) = bcocz(i) + node_vbox(i, 3)

              mirror_indices(i, 1:3) = lattice_indices(node_vbox(i, 1:3))
              mirror_level(i) = maxval(abs(mirror_indices(i, :)))
            end if
            !write(*,*) node_vbox(i,1:3),mirror_indices(i, 1:3),mirror_level(i)

            ! compute real center coordinate
            call key_to_coord(t%bounding_box, bkey, bx, by, bz)
            
            if ( present(node_vbox) ) then
              bx = bx + node_vbox(i, 1)
              by = by + node_vbox(i, 2)
              bz = bz + node_vbox(i, 3)
            end if

            do j=1,8
              bcornersidx(8*(i-1)+j) = 8*(i-1)+j - 1
              bshift(1:3) = box_shift(1:3,j) * bsize(1:3)
              bcornersx(8*(i-1)+j)   = bx + bshift(1)
              bcornersy(8*(i-1)+j)   = by + bshift(2)
              bcornersz(8*(i-1)+j)   = bz + bshift(3)
            end do
          end do

          if ( present(filename_point) ) then
            call vtk_point%create(filename_point, step, tsim, vtk_step)
            call vtk_point%write_headers(num_nodes, 0)
            call vtk_point%startpoints()
              call vtk_point%write_data_array("xyz", num_nodes, bcocx, bcocy, bcocz)
            call vtk_point%finishpoints()
            call vtk_point%startpointdata()
              call vtk_point%write_data_array("level", num_nodes, blevel)
              if ( present(node_vbox) ) then
                call vtk_point%write_data_array("mirror_level", num_nodes, mirror_level)
                call vtk_point%write_data_array("mirror_indices", num_nodes, &
                  mirror_indices(:, 1), mirror_indices(:, 2), mirror_indices(:, 3))
              end if
              call vtk_point%write_data_array("charge", num_nodes, bq)
            call vtk_point%finishpointdata()
            call vtk_point%dont_write_cells()
            call vtk_point%write_final()
            call vtk_point%close()
          end if

          if ( present(filename_box) ) then
            call vtk_box%create(filename_box, step, tsim, vtk_step)
            call vtk_box%write_headers(num_nodes*8, num_nodes)
            call vtk_box%startpoints()
              call vtk_box%write_data_array("corners", 8*num_nodes, bcornersx, bcornersy, bcornersz)
            call vtk_box%finishpoints()
            call vtk_box%startpointdata()
              ! no point data here
            call vtk_box%finishpointdata()
            call vtk_box%startcells()
              call vtk_box%write_data_array("connectivity", num_nodes*8, bcornersidx)
              call vtk_box%write_data_array("offsets", num_nodes, bcornersoffsets)
              call vtk_box%write_data_array("types", num_nodes, bcornerstypes)
            call vtk_box%finishcells()
            call vtk_box%startcelldata()
              call vtk_box%write_data_array("processor", num_nodes, bowner)
              call vtk_box%write_data_array("level", num_nodes, blevel)
              if ( present(node_vbox) ) then
                call vtk_box%write_data_array("mirror_level", num_nodes, mirror_level)
                call vtk_point%write_data_array("mirror_indices", num_nodes, &
                  mirror_indices(:, 1), mirror_indices(:, 2), mirror_indices(:, 3))
              end if
              call vtk_box%write_data_array("center_of_charge", num_nodes, bcocx, bcocy, bcocz)
              call vtk_box%write_data_array("total_charge", num_nodes, bq)
            call vtk_box%finishcelldata()
            call vtk_box%write_final()
            call vtk_box%close()
          end if

          deallocate(bcocx, bcocy, bcocz, bq)
          deallocate(bcornersx, bcornersy, bcornersz)
          deallocate(bcornersidx)
          deallocate(bcornersoffsets, bcornerstypes, bowner, blevel, &
            mirror_level)
          deallocate(mirror_indices)

        end subroutine


        !>
        !> Writes the interaction partners of the particle with the 
        !> specified label into vtk files, once as boxes, once as points
        !>
        subroutine write_interaction_partners_to_vtk(step, label,tsim, vtk_step)
          use treevars
          use module_interaction_specific
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          integer, intent(in) :: label
          real*8, intent(in) :: tsim

          character(255) :: fn_box, fn_point

          if (me .ne. 0) return

          write(fn_point,"(a,i3.3)") "int_partner_coc",label
          write(fn_box  ,"(a,i3.3)") "int_partner_box",label

          call write_nodes_to_vtk(step, tsim, vtk_step, &
            no_interaction_partners(label), interaction_keylist(label, :), &
            node_vbox = interaction_vbox(label, :, :), &
            filename_box = trim(fn_box), filename_point = trim(fn_point))

        end subroutine


        !>
        !> Writes the tree branches structure into a vtk file.
        !>
        subroutine write_branches_to_vtk(step, tsim, vtk_step)
          use treevars
          use module_pepc, only: global_tree
          use module_pepc_types, only: t_tree_node
          use module_tree, only: t_tree, tree_lookup_root, tree_allocated
          use module_debug
          implicit none

          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim

          type(t_tree), pointer :: t
          type(t_tree_node), pointer :: r
          integer*8, allocatable :: branch_keys(:)
          integer*8 :: i

          ! TODO: generalize!
          t => global_tree

          if (.not. t%comm_env%first) return

          if (.not. tree_allocated(t)) then
            write(*,*) 'write_branches_to_vtk(): tree is not allocated, aborting branch output.'
            return
          endif

          allocate(branch_keys(t%nbranch))
          i = 0
          call tree_lookup_root(t, r)
          call collect_branches(r)
          DEBUG_ASSERT(i == t%nbranch)

          call write_nodes_to_vtk(step, tsim, vtk_step, int(t%nbranch), &
            branch_keys, filename_box = "branches")

          deallocate(branch_keys)

          contains

          recursive subroutine collect_branches(n)
            use module_tree_node
            use module_tree, only: tree_node_get_first_child, tree_node_get_next_sibling
            implicit none

            type(t_tree_node), intent(in) :: n

            type(t_tree_node), pointer :: s, ns

            s => null()
            ns => null()

            if (btest(n%flags, TREE_NODE_FLAG_IS_BRANCH_NODE)) then
              i = i + 1
              branch_keys(i) = n%key
            else if (tree_node_get_first_child(t, n, s)) then
              do
                call collect_branches(s)
                if (.not. tree_node_get_next_sibling(t, s, ns)) then
                  exit
                end if
                s => ns
              end do
            end if
          end subroutine collect_branches
        end subroutine

        !>
        !> Writes the space filling curve into a parallel set of vtk files
        !> pepc_fields must have been called with no_dealloc=.true. before
        !>
        subroutine write_spacecurve_to_vtk(step, tsim, vtk_step, particles)
          use module_vtk
          use module_pepc_types, only: t_particle
          use module_pepc, only: global_tree
          use module_tree, only: t_tree
          implicit none

          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_particle), intent(in) :: particles(:)

          type(t_tree), pointer :: t
          type(vtkfile_unstructured_grid) :: vtk
          integer :: i, npp

          t => global_tree
          npp = int(t%npart_me)

            call vtk%create_parallel("spacecurve", step, t%comm_env%rank, t%comm_env%size, tsim, vtk_step)
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
                  call vtk%write_data_array("processor", t%comm_env%rank)
                call vtk%finishcelldata()
              call vtk%write_final()
            call vtk%close()
        end subroutine


end module module_treediags
