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
      public write_leaves_to_vtk
      public write_spacecurve_to_vtk
      public write_interaction_partners_to_vtk

      contains

        !>
        !> Writes a collection of nodes into vtk files as points
        !>
        subroutine write_nodes_to_vtk_as_points(fname, mpi_comm, step, tsim, vtk_step, t, nodes, node_vbox)
          use treevars
          use module_vtk
          use module_spacefilling
          use module_mirror_boxes
          use module_pepc_types, only: t_tree_node
          use module_tree, only: t_tree
          use module_tree_node
          use module_pepc, only: global_tree
          use module_pepc_types, only: kind_node
          implicit none

          include 'mpif.h'

          character(*), intent(in) :: fname
          integer(kind_default) :: mpi_comm
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_tree), intent(in) :: t
          integer(kind_node), dimension(:), intent(in) :: nodes
          real*8, dimension(:,:), intent(in), optional :: node_vbox

          type(t_tree_node), pointer :: bnode
          integer(kind_node) :: i, num_nodes
          integer(kind_default) :: mpi_rank, mpi_size, ierr
          integer :: j
          real*8, dimension(:), allocatable  :: bcocx, bcocy, bcocz
          integer, dimension(:), allocatable :: bowner, blevel, mirror_level
          integer, dimension(:, :), allocatable :: mirror_indices

          type(vtkfile_unstructured_grid) :: vtk

          call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
          call MPI_Comm_size(mpi_comm, mpi_size, ierr)

          num_nodes = size(nodes)

          allocate(bcocx(num_nodes), bcocy(num_nodes), bcocz(num_nodes))
          allocate(bowner(num_nodes), blevel(num_nodes))
          allocate(mirror_level(num_nodes), mirror_indices(num_nodes, 3))

          do i = 1, num_nodes
            bnode => t%nodes(nodes(i))
            bowner(i) = bnode%owner
            blevel(i) = bnode%level

            bcocx(i) = bnode%interaction_data%coc(1)
            bcocy(i) = bnode%interaction_data%coc(2)
            bcocz(i) = bnode%interaction_data%coc(3)

            if (present(node_vbox)) then
              bcocx(i) = bcocx(i) + node_vbox(i, 1)
              bcocy(i) = bcocy(i) + node_vbox(i, 2)
              bcocz(i) = bcocz(i) + node_vbox(i, 3)

              mirror_indices(i, 1:3) = lattice_indices(node_vbox(i, 1:3))
              mirror_level(i) = maxval(abs(mirror_indices(i, :)))
            end if
          end do

          call vtk%create_parallel(fname, step, mpi_rank, mpi_size, tsim, vtk_step)
          call vtk%write_headers(num_nodes, 0_kind_node)
          call vtk%startpoints()
            call vtk%write_data_array("xyz", bcocx, bcocy, bcocz)
          call vtk%finishpoints()
          call vtk%startpointdata()
            call vtk%write_data_array("level", blevel)
            if ( present(node_vbox) ) then
              call vtk%write_data_array("mirror_level", mirror_level)
              call vtk%write_data_array("mirror_indices", &
                mirror_indices(:, 1), mirror_indices(:, 2), mirror_indices(:, 3))
            end if
          call vtk%finishpointdata()
          call vtk%dont_write_cells()
          call vtk%write_final()
          call vtk%close()

          deallocate(bcocx, bcocy, bcocz, bowner, blevel, mirror_level, mirror_indices)
        end subroutine write_nodes_to_vtk_as_points


        !>
        !> Writes a collection of nodes into vtk files as boxes
        !>
        subroutine write_nodes_to_vtk_as_boxes(fname, mpi_comm, step, tsim, vtk_step, t, nodes, node_vbox)
          use treevars
          use module_vtk
          use module_spacefilling
          use module_mirror_boxes
          use module_pepc_types, only: t_tree_node
          use module_tree, only: t_tree
          use module_tree_node
          use module_pepc, only: global_tree
          use module_pepc_types, only: kind_node
          implicit none

          include 'mpif.h'

          character(*), intent(in) :: fname
          integer(kind_default) :: mpi_comm
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_tree), intent(in) :: t
          integer(kind_node), dimension(:), intent(in) :: nodes
          real*8, dimension(:,:), intent(in), optional :: node_vbox

          type(t_tree_node), pointer :: bnode
          integer(kind_node) :: i, num_nodes
          integer(kind_default) :: mpi_rank, mpi_size, ierr
          integer :: j
          integer(kind_key) :: bkey
          real*8 :: bsize(3)
          real*8, dimension(:), allocatable  :: bcornersx, bcornersy, bcornersz
          integer(kind_node), dimension(:), allocatable :: bcornersidx, bcornersoffsets
          integer, dimension(:), allocatable :: bcornerstypes, bowner, blevel, &
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
          type(vtkfile_unstructured_grid) :: vtk

          call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
          call MPI_Comm_size(mpi_comm, mpi_size, ierr)

          num_nodes = size(nodes)

          allocate(bcornersx(num_nodes*8), bcornersy(num_nodes*8), &
            bcornersz(num_nodes*8))
          allocate(bcornersidx(num_nodes * 8))
          allocate(bcornersoffsets(num_nodes), bcornerstypes(num_nodes), &
            bowner(num_nodes), blevel(num_nodes), mirror_level(num_nodes))
          allocate(mirror_indices(num_nodes, 3))

          do i = 1, num_nodes
            bnode => t%nodes(nodes(i))
            bkey      = bnode%key
            bowner(i) = bnode%owner
            blevel(i) = bnode%level
            bsize     = t%bounding_box%boxsize / 2**blevel(i)

            ! prepare voxel data structure
            bcornerstypes(i)   = VTK_VOXEL
            bcornersoffsets(i) = 8*i

            ! compute real center coordinate
            call key_to_coord(t%bounding_box, bkey, bx, by, bz)
            
            if ( present(node_vbox) ) then
              bx = bx + node_vbox(i, 1)
              by = by + node_vbox(i, 2)
              bz = bz + node_vbox(i, 3)

              mirror_indices(i, 1:3) = lattice_indices(node_vbox(i, 1:3))
              mirror_level(i) = maxval(abs(mirror_indices(i, :)))
            end if

            do j=1,8
              bcornersidx(8*(i-1)+j) = 8*(i-1)+j - 1
              bshift(1:3) = box_shift(1:3,j) * bsize(1:3)
              bcornersx(8*(i-1)+j)   = bx + bshift(1)
              bcornersy(8*(i-1)+j)   = by + bshift(2)
              bcornersz(8*(i-1)+j)   = bz + bshift(3)
            end do
          end do

          call vtk%create_parallel(fname, step, mpi_rank, mpi_size, tsim, vtk_step)
          call vtk%write_headers(num_nodes*8, num_nodes)
          call vtk%startpoints()
            call vtk%write_data_array("corners", bcornersx, bcornersy, bcornersz)
          call vtk%finishpoints()
          call vtk%startpointdata()
            ! no point data here
          call vtk%finishpointdata()
          call vtk%startcells()
            call vtk%write_data_array("connectivity", bcornersidx)
            call vtk%write_data_array("offsets", bcornersoffsets)
            call vtk%write_data_array("types", bcornerstypes)
          call vtk%finishcells()
          call vtk%startcelldata()
            call vtk%write_data_array("processor", bowner)
            call vtk%write_data_array("level", blevel)
            if ( present(node_vbox) ) then
              call vtk%write_data_array("mirror_level", mirror_level)
              call vtk%write_data_array("mirror_indices", &
                mirror_indices(:, 1), mirror_indices(:, 2), mirror_indices(:, 3))
            end if
          call vtk%finishcelldata()
          call vtk%write_final()
          call vtk%close()

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
          use module_pepc, only: global_tree
          use module_tree, only: tree_lookup_node_critical
          use treevars
          use module_interaction_specific
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          integer, intent(in) :: label
          real*8, intent(in) :: tsim

          integer(kind_node), allocatable :: i, partner_nodes(:)
          character(255) :: fn_box, fn_point

          write(fn_point,"(a,i3.3)") "int_partner_coc",label
          write(fn_box  ,"(a,i3.3)") "int_partner_box",label

          if (0 == me) then
            allocate(partner_nodes(no_interaction_partners(label)))

            do i = 1, no_interaction_partners(label)
              call tree_lookup_node_critical(global_tree, interaction_keylist(label, i), partner_nodes(i), &
                'write_interaction_partners_to_vtk()')
            end do

            call write_nodes_to_vtk_as_points(fn_point, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes, &
              node_vbox = interaction_vbox(label, :, :))

            call write_nodes_to_vtk_as_boxes(fn_box, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes, &
              node_vbox = interaction_vbox(label, :, :))

            deallocate(partner_nodes)
          else
            allocate(partner_nodes(0))
            call write_nodes_to_vtk_as_points(fn_point, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes)
            call write_nodes_to_vtk_as_boxes(fn_box, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes)
            deallocate(partner_nodes)
          end if
        end subroutine


        !>
        !> Writes the tree leaves into a vtk file.
        !>
        subroutine write_leaves_to_vtk(step, tsim, vtk_step, t_)
          use module_pepc, only: global_tree
          use module_pepc_types, only: t_tree_node, kind_node
          use module_tree, only: t_tree, tree_lookup_root, tree_allocated
          use module_debug
          implicit none

          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_tree), optional, target, intent(in) :: t_

          type(t_tree), pointer :: t
          integer(kind_node) :: r, i
          integer(kind_node), allocatable :: leaves(:)

          if (present(t_)) then
            t => t_
          else
            t => global_tree
          end if

          if (.not. tree_allocated(t)) then
            write(*,*) 'write_leaves_to_vtk(): tree is not allocated, aborting leaves output.'
            return
          end if

          allocate(leaves(t%nleaf_me))
          i = 0
          call tree_lookup_root(t, r)
          call collect_leaves(t, r)
          DEBUG_ASSERT(i == t%nleaf_me)

          call write_nodes_to_vtk_as_boxes("leaves", t%comm_env%comm, step, tsim, vtk_step, t, leaves)

          deallocate(leaves)

          contains

          recursive subroutine collect_leaves(t, nidx)
            use module_tree_node
            use module_tree, only: t_tree
            use module_pepc_types, only: kind_node
            implicit none

            type(t_tree), intent(in) :: t
            integer(kind_node), intent(in) :: nidx

            integer(kind_node) :: s, ns

            s  = NODE_INVALID
            ns = NODE_INVALID

            associate (n => t%nodes(nidx))
              if (btest(n%flags_local1, TREE_NODE_FLAG_LOCAL1_HAS_LOCAL_CONTRIBUTIONS)) then
                if (tree_node_is_leaf(n)) then
                  i = i + 1
                  leaves(i) = nidx
                else if (tree_node_get_first_child(n, s)) &
                  then
                  do
                    call collect_leaves(t, s)
                    if (.not. tree_node_get_next_sibling(t%nodes(s), ns)) then
                      exit
                    end if
                    s = ns
                  end do
                end if
              end if
            end associate
          end subroutine collect_leaves
        end subroutine


        !>
        !> Writes the tree branches structure into a vtk file.
        !>
        subroutine write_branches_to_vtk(step, tsim, vtk_step, t_)
          use module_pepc, only: global_tree
          use module_pepc_types, only: t_tree_node, kind_node
          use module_tree, only: t_tree, tree_lookup_root, tree_allocated
          use module_debug
          implicit none

          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_tree), optional, target, intent(in) :: t_

          type(t_tree), pointer :: t
          integer(kind_node) :: r, i
          integer(kind_node), allocatable :: branch_nodes(:)

          if (present(t_)) then
            t => t_
          else
            t => global_tree
          end if

          if (.not. tree_allocated(t)) then
            write(*,*) 'write_branches_to_vtk(): tree is not allocated, aborting branch output.'
            return
          end if

          if (t%comm_env%first) then
            allocate(branch_nodes(t%nbranch))
            i = 0
            call tree_lookup_root(t, r)
            call collect_branches(t, r)
            DEBUG_ASSERT(i == t%nbranch)

            call write_nodes_to_vtk_as_boxes("branches", t%comm_env%comm, step, tsim, vtk_step, t, branch_nodes)

            deallocate(branch_nodes)
          else
            allocate(branch_nodes(0))
            call write_nodes_to_vtk_as_boxes("branches", t%comm_env%comm, step, tsim, vtk_step, t, branch_nodes)
            deallocate(branch_nodes)
          end if

          contains

          recursive subroutine collect_branches(t, nidx)
            use module_tree_node
            use module_tree, only: t_tree
            use module_pepc_types, only: kind_node
            implicit none

            type(t_tree), intent(in) :: t
            integer(kind_node), intent(in) :: nidx

            integer(kind_node) :: s, ns

            s  = NODE_INVALID
            ns = NODE_INVALID

            associate (n => t%nodes(nidx))
              if (btest(n%flags_global, TREE_NODE_FLAG_GLOBAL_IS_BRANCH_NODE)) then
                i = i + 1
                branch_nodes(i) = nidx
              else if (tree_node_get_first_child(n, s)) then
                do
                  call collect_branches(t, s)
                  if (.not. tree_node_get_next_sibling(t%nodes(s), ns)) then
                    exit
                  end if
                  s = ns
                end do
              end if
            end associate
          end subroutine collect_branches
        end subroutine


        !>
        !> Writes the space filling curve into a parallel set of vtk files
        !> pepc_fields must have been called with no_dealloc=.true. before
        !>
        subroutine write_spacecurve_to_vtk(step, tsim, vtk_step, particles)
          use module_vtk
          use module_pepc_types
          use module_pepc, only: global_tree
          use module_tree, only: t_tree
          implicit none

          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          type(t_particle), intent(in) :: particles(:)

          type(t_tree), pointer :: t
          type(vtkfile_unstructured_grid) :: vtk
          integer(kind_particle) :: i, npp

          t => global_tree
          npp = t%npart_me

            call vtk%create_parallel("spacecurve", step, t%comm_env%rank, t%comm_env%size, tsim, vtk_step)
              call vtk%write_headers(npp, 1_kind_particle)
                call vtk%startpoints()
                  call vtk%write_data_array("xyz", particles(1:npp)%x(1), particles(1:npp)%x(2), particles(1:npp)%x(3))
                call vtk%finishpoints()
                call vtk%startpointdata()
                  ! no point data here
                call vtk%finishpointdata()
                call vtk%startcells()
                  call vtk%write_data_array("connectivity", [(i,i=0, npp-1)])
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
