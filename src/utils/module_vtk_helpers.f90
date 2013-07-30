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
!>  Encapsulates helper functions for outputting different data formats to vtk files
!>
module module_vtk_helpers
  implicit none
  private

  public vtk_write_particles

  public vtk_write_branches
  public vtk_write_leaves
  public vtk_write_interaction_partners

  public vtk_write_spacecurve

  public vtk_write_field_on_grid
  public vtk_write_densities_on_grid

  contains

  !>
  !> Writes particles to a VTK file
  !>
  subroutine vtk_write_particles(fname, mpi_comm, step, time, vtk_step, p, helper_func, coord_scale)
    use module_vtk
    use module_pepc_types
    use module_interaction_specific_types
    implicit none

    include 'mpif.h'
    
    character(*), intent(in) :: fname
    integer(kind_default), intent(in) :: mpi_comm
    integer, intent(in) :: step
    real*8, intent(in) :: time
    integer, intent(in) :: vtk_step
    type(t_particle), intent(in) :: p(:)

    interface
      subroutine helper_func(d, r, vtkf)
        use module_interaction_specific_types, only: t_particle_data, t_particle_results
        use module_vtk, only: vtkfile_unstructured_grid
        implicit none
        type(t_particle_data), intent(in) :: d(:)
        type(t_particle_results), intent(in) :: r(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func
    real*8, intent(in), optional :: coord_scale

    integer(kind_particle) :: i, np
    integer(kind_pe) :: mpi_rank, mpi_size
    integer(kind_default) :: ierr
    type(vtkfile_unstructured_grid) :: vtk
    
    call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    np = size(p)

    call vtk%create_parallel(fname, step, mpi_rank, mpi_size, time, vtk_step)
      call vtk%write_headers(np, np)
      call vtk%startpoints()
        call vtk%write_data_array("xyz", p(:)%x(1), p(:)%x(2), p(:)%x(3), coord_scale)
      call vtk%finishpoints()
      call vtk%startpointdata()
        call vtk%write_data_array("work", p(:)%work)
        call vtk%write_data_array("label", p(:)%label)
        call vtk%write_data_array("local index", [(i,i=1,np)])
        call vtk%write_data_array("processor", int(np, kind = 4), mpi_rank)
        if (present(helper_func)) then
          call helper_func(p(:)%data, p(:)%results, vtk)
        else
          call vtk_write_particle_data_results(p(:)%data, p(:)%results, vtk)
        end if
      call vtk%finishpointdata()
      call vtk%startcells()
        call vtk%write_data_array("connectivity", [(i, i=0,np-1)])
        call vtk%write_data_array("offsets", [(i,i=1,np)])
        call vtk%write_data_array("types", [(VTK_VERTEX,i=1,np)])
      call vtk%finishcells()
      call vtk%startcelldata()
       ! no cell data here as cells correspond to points anyway, in case of problems use PointDataToCellData Filter in Paraview
      call vtk%finishcelldata()
      call vtk%write_final()
    call vtk%close()
  end subroutine vtk_write_particles


  !>
  !> Writes a collection of nodes into vtk files as points
  !>
  subroutine vtk_write_nodes_as_points(fname, mpi_comm, step, tsim, vtk_step, t, nodes, node_vbox, helper_func, coord_scale)
    use treevars
    use module_vtk
    use module_mirror_boxes
    use module_pepc_types
    use module_interaction_specific_types
    use module_tree, only: t_tree
    use module_tree_node
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
    real*8, intent(in), optional :: coord_scale

    interface
      subroutine helper_func(d, vtkf)
        use module_interaction_specific_types, only: t_tree_node_interaction_data
        use module_vtk, only: vtkfile_unstructured_grid
        implicit none
        type(t_tree_node_interaction_data), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func

    type(t_tree_node), pointer :: bnode
    integer(kind_node) :: i, num_nodes
    integer(kind_default) :: mpi_rank, mpi_size, ierr
    real*8, dimension(:), allocatable  :: bcocx, bcocy, bcocz
    integer, dimension(:), allocatable :: bowner, blevel, mirror_level
    integer*8, dimension(:), allocatable :: bleaves, bdescendants
    type(t_tree_node_interaction_data), dimension(:), allocatable :: bdata
    integer, dimension(:, :), allocatable :: mirror_indices

    type(vtkfile_unstructured_grid) :: vtk

    call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    num_nodes = size(nodes)

    allocate(bcocx(num_nodes), bcocy(num_nodes), bcocz(num_nodes))
    allocate(bowner(num_nodes), blevel(num_nodes), bleaves(num_nodes), &
      bdescendants(num_nodes), bdata(num_nodes))
    allocate(mirror_level(num_nodes), mirror_indices(num_nodes, 3))

    do i = 1, num_nodes
      bnode => t%nodes(nodes(i))
      bowner(i)       = bnode%owner
      blevel(i)       = bnode%level
      bleaves(i)      = bnode%leaves
      bdescendants(i) = bnode%descendants
      bdata(i)        = bnode%interaction_data

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
        call vtk%write_data_array("xyz", bcocx, bcocy, bcocz, coord_scale)
      call vtk%finishpoints()
      call vtk%startpointdata()
        call vtk%write_data_array("processor", bowner)
        call vtk%write_data_array("level", blevel)
        call vtk%write_data_array("leaves", bleaves)
        call vtk%write_data_array("descendants", bdescendants)
        if (present(node_vbox)) then
          call vtk%write_data_array("mirror_level", mirror_level)
          call vtk%write_data_array("mirror_indices", &
            mirror_indices(:, 1), mirror_indices(:, 2), mirror_indices(:, 3))
        end if
        if (present(helper_func)) then
          call helper_func(bdata, vtk)
        else
          call vtk_write_node_interaction_data(bdata, vtk)
        end if
      call vtk%finishpointdata()
      call vtk%startcells()
        call vtk%write_data_array("connectivity", [(i, i=0,num_nodes-1)])
        call vtk%write_data_array("offsets", [(i,i=1,num_nodes)])
        call vtk%write_data_array("types", [(VTK_VERTEX,i=1,num_nodes)])
      call vtk%finishcells()
      call vtk%startcelldata()
       ! no cell data here as cells correspond to points anyway, in case of problems use PointDataToCellData Filter in Paraview
      call vtk%finishcelldata()
      call vtk%write_final()
    call vtk%close()

    deallocate(bcocx, bcocy, bcocz, bowner, blevel, mirror_level, mirror_indices)
    deallocate(bleaves, bdescendants, bdata)
  end subroutine vtk_write_nodes_as_points


  !>
  !> Writes a collection of nodes into vtk files as boxes
  !>
  subroutine vtk_write_nodes_as_boxes(fname, mpi_comm, step, tsim, vtk_step, t, nodes, node_vbox, helper_func, coord_scale)
    use treevars
    use module_vtk
    use module_spacefilling
    use module_mirror_boxes
    use module_pepc_types
    use module_interaction_specific_types
    use module_tree, only: t_tree
    use module_tree_node
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
    real*8, intent(in), optional :: coord_scale

    interface
      subroutine helper_func(d, vtkf)
        use module_interaction_specific_types, only: t_tree_node_interaction_data
        use module_vtk, only: vtkfile_unstructured_grid
        implicit none
        type(t_tree_node_interaction_data), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func

    type(t_tree_node), pointer :: bnode
    integer(kind_node) :: i, num_nodes
    integer(kind_default) :: mpi_rank, mpi_size, ierr
    integer :: j
    integer(kind_key) :: bkey
    real*8 :: bsize(3)
    real*8, dimension(:), allocatable  :: bcornersx, bcornersy, bcornersz
    integer(kind_node), dimension(:), allocatable :: bcornersidx, bcornersoffsets
    integer, dimension(:), allocatable :: bowner, blevel, mirror_level
    integer*8, dimension(:), allocatable :: bleaves, bdescendants
    type(t_tree_node_interaction_data), dimension(:), allocatable :: bdata

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
    allocate(bcornersoffsets(num_nodes), bleaves(num_nodes), bdescendants(num_nodes), &
      bowner(num_nodes), blevel(num_nodes), bdata(num_nodes), mirror_level(num_nodes))
    allocate(mirror_indices(num_nodes, 3))

    do i = 1, num_nodes
      bnode => t%nodes(nodes(i))
      bkey            = bnode%key
      bowner(i)       = bnode%owner
      blevel(i)       = bnode%level
      bleaves(i)      = bnode%leaves
      bdescendants(i) = bnode%descendants
      bdata(i)        = bnode%interaction_data
      bsize           = t%bounding_box%boxsize / 2**blevel(i)

      ! prepare voxel data structure
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
        call vtk%write_data_array("corners", bcornersx, bcornersy, bcornersz, coord_scale)
      call vtk%finishpoints()
      call vtk%startpointdata()
        ! no point data here
      call vtk%finishpointdata()
      call vtk%startcells()
        call vtk%write_data_array("connectivity", bcornersidx)
        call vtk%write_data_array("offsets", bcornersoffsets)
        call vtk%write_data_array("types", int(num_nodes, kind = 4), VTK_VOXEL)
      call vtk%finishcells()
      call vtk%startcelldata()
        call vtk%write_data_array("processor", bowner)
        call vtk%write_data_array("level", blevel)
        call vtk%write_data_array("leaves", bleaves)
        call vtk%write_data_array("descendants", bdescendants)
        if ( present(node_vbox) ) then
          call vtk%write_data_array("mirror_level", mirror_level)
          call vtk%write_data_array("mirror_indices", &
            mirror_indices(:, 1), mirror_indices(:, 2), mirror_indices(:, 3))
        end if
        if (present(helper_func)) then
          call helper_func(bdata, vtk)
        else
          call vtk_write_node_interaction_data(bdata, vtk)
        end if
      call vtk%finishcelldata()
      call vtk%write_final()
    call vtk%close()

    deallocate(bcornersx, bcornersy, bcornersz)
    deallocate(bcornersidx)
    deallocate(bcornersoffsets, bowner, blevel, bleaves, bdescendants, mirror_level)
    deallocate(bdata, mirror_indices)
  end subroutine vtk_write_nodes_as_boxes


  !>
  !> Writes the interaction partners of the particle with the 
  !> specified label into vtk files, once as boxes, once as points
  !>
  subroutine vtk_write_interaction_partners(step, label,tsim, vtk_step, interaction_nodelist, no_interaction_partners, interaction_vbox, helper_func)
    use module_pepc, only: global_tree
    use module_pepc_types
    use treevars
    integer, intent(in) :: step
    integer, intent(in) :: vtk_step
    integer, intent(in) :: label
    real*8, intent(in) :: tsim
    
    integer(kind_node), intent(in) :: interaction_nodelist(:,:)
    integer(kind_node), intent(in) :: no_interaction_partners(:)
    real*8, intent(in) :: interaction_vbox(:,:,:)
    
    interface
      subroutine helper_func(d, vtkf)
        use module_interaction_specific_types, only: t_tree_node_interaction_data
        use module_vtk, only: vtkfile_unstructured_grid
        implicit none
        type(t_tree_node_interaction_data), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func

    integer(kind_node), allocatable :: i, partner_nodes(:)
    character(255) :: fn_box, fn_point

    write(fn_point,"(a,i3.3)") "int_partner_coc",label
    write(fn_box  ,"(a,i3.3)") "int_partner_box",label

    if (0 == me) then
      allocate(partner_nodes(no_interaction_partners(label)))

      do i = 1, no_interaction_partners(label)
        partner_nodes(i) = interaction_nodelist(label, i)
      end do

      call vtk_write_nodes_as_points(fn_point, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes, &
        node_vbox = interaction_vbox(label, :, :), helper_func = helper_func)

      call vtk_write_nodes_as_boxes(fn_box, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes, &
        node_vbox = interaction_vbox(label, :, :), helper_func = helper_func)

      deallocate(partner_nodes)
    else
      allocate(partner_nodes(0))
      call vtk_write_nodes_as_points(fn_point, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes)
      call vtk_write_nodes_as_boxes(fn_box, MPI_Comm_lpepc, step, tsim, vtk_step, global_tree, partner_nodes)
      deallocate(partner_nodes)
    end if
  end subroutine


  !>
  !> Writes the tree leaves into a vtk file.
  !>
  subroutine vtk_write_leaves(step, tsim, vtk_step, t_, helper_func, coord_scale)
    use module_pepc, only: global_tree
    use module_pepc_types, only: t_tree_node, kind_node
    use module_tree, only: t_tree, tree_allocated
    use module_debug
    implicit none

    integer, intent(in) :: step
    integer, intent(in) :: vtk_step
    real*8, intent(in) :: tsim
    type(t_tree), optional, target, intent(in) :: t_
    real*8, intent(in), optional :: coord_scale

    interface
      subroutine helper_func(d, vtkf)
        use module_interaction_specific_types, only: t_tree_node_interaction_data
        use module_vtk, only: vtkfile_unstructured_grid
        implicit none
        type(t_tree_node_interaction_data), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func

    type(t_tree), pointer :: t
    integer(kind_node) :: i
    integer(kind_node), allocatable :: leaves(:)

    if (present(t_)) then
      t => t_
    else
      t => global_tree
    end if

    if (.not. tree_allocated(t)) then
      write(*,*) 'vtk_write_leaves(): tree is not allocated, aborting leaves output.'
      return
    end if

    allocate(leaves(t%nleaf_me))
    i = 0
    call collect_leaves(t, t%node_root)
    DEBUG_ASSERT(i == t%nleaf_me)

    call vtk_write_nodes_as_boxes("leaves", t%comm_env%comm, step, tsim, vtk_step, t, leaves, &
      helper_func = helper_func, coord_scale = coord_scale)

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
        if (btest(n%flags_local, TREE_NODE_FLAG_LOCAL_HAS_LOCAL_CONTRIBUTIONS)) then
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
  subroutine vtk_write_branches(step, tsim, vtk_step, t_, helper_func, coord_scale)
    use module_pepc, only: global_tree
    use module_pepc_types, only: t_tree_node, kind_node
    use module_tree, only: t_tree, tree_allocated
    use module_debug
    implicit none

    integer, intent(in) :: step
    integer, intent(in) :: vtk_step
    real*8, intent(in) :: tsim
    type(t_tree), optional, target, intent(in) :: t_
    real*8, intent(in), optional :: coord_scale

    interface
      subroutine helper_func(d, vtkf)
        use module_interaction_specific_types, only: t_tree_node_interaction_data
        use module_vtk, only: vtkfile_unstructured_grid
        implicit none
        type(t_tree_node_interaction_data), intent(in) :: d(:)
        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      end subroutine helper_func
    end interface

    optional :: helper_func

    type(t_tree), pointer :: t
    integer(kind_node) :: i
    integer(kind_node), allocatable :: branch_nodes(:)

    if (present(t_)) then
      t => t_
    else
      t => global_tree
    end if

    if (.not. tree_allocated(t)) then
      write(*,*) 'vtk_write_branches(): tree is not allocated, aborting branch output.'
      return
    end if

    if (t%comm_env%first) then
      allocate(branch_nodes(t%nbranch))
      i = 0
      call collect_branches(t, t%node_root)
      DEBUG_ASSERT(i == t%nbranch)

      call vtk_write_nodes_as_boxes("branches", t%comm_env%comm, step, tsim, vtk_step, t, branch_nodes, &
        helper_func = helper_func, coord_scale = coord_scale)

      deallocate(branch_nodes)
    else
      allocate(branch_nodes(0))
      call vtk_write_nodes_as_boxes("branches", t%comm_env%comm, step, tsim, vtk_step, t, branch_nodes, &
        helper_func = helper_func, coord_scale = coord_scale)
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
  !>
  subroutine vtk_write_spacecurve(step, tsim, vtk_step, particles, coord_scale)
    use treevars, only: MPI_COMM_lpepc
    use module_vtk
    use module_pepc_types
    implicit none

    include 'mpif.h'

    integer, intent(in) :: step
    integer, intent(in) :: vtk_step
    real*8, intent(in) :: tsim
    type(t_particle), intent(in) :: particles(:)
    real*8, intent(in), optional :: coord_scale

    type(vtkfile_unstructured_grid) :: vtk
    integer(kind_particle) :: i, npp
    integer :: mpi_rank, mpi_size, ierr

    call MPI_Comm_rank(MPI_COMM_lpepc, mpi_rank, ierr)
    call MPI_Comm_size(MPI_COMM_lpepc, mpi_size, ierr)
    npp = size(particles, kind = kind(npp))

    call vtk%create_parallel("spacecurve", step, mpi_rank, mpi_size, tsim, vtk_step)
      call vtk%write_headers(npp, 1_kind_particle)
        call vtk%startpoints()
          call vtk%write_data_array("xyz", particles(:)%x(1), particles(:)%x(2), particles(:)%x(3), coord_scale)
        call vtk%finishpoints()
        call vtk%startpointdata()
          ! no point data here
        call vtk%finishpointdata()
        call vtk%startcells()
          call vtk%write_data_array("connectivity", [(i, i = 0, npp - 1)])
          call vtk%write_data_array("offsets", npp)
          call vtk%write_data_array("types", VTK_POLY_LINE)
        call vtk%finishcells()
        call vtk%startcelldata()
          call vtk%write_data_array("processor", mpi_rank)
        call vtk%finishcelldata()
      call vtk%write_final()
    call vtk%close()
  end subroutine


  !>
  !> Writes scalar and vector field values on a regular rectangular grid to VTK.
  !>
  subroutine vtk_write_field_on_grid(filename, step, tsim, vtk_step, globaldims, mydims, xcoords, ycoords, zcoords, &
                    scalarvalues, scalarname, vectorvalues, vectorname, mpi_comm)
    use module_vtk
    implicit none
    include 'mpif.h'
    character(*), intent(in) :: filename, scalarname, vectorname
    integer, intent(in) :: step
    integer, intent(in) :: vtk_step
    real*8, intent(in) :: tsim
    integer, dimension(2,3), intent(in) :: globaldims, mydims
    real*8, intent(in) :: xcoords(:), ycoords(:), zcoords(:)
    real*8, intent(in) :: scalarvalues(:, :, :), vectorvalues(:,:,:,:)
    integer, intent(in) :: mpi_comm

    integer :: mpi_rank, mpi_size, ierr

    type(vtkfile_rectilinear_grid) :: vtk
    
    call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    call vtk%create_parallel(trim(filename), step, mpi_rank, mpi_size, tsim, vtk_step)
      call vtk%set_communicator(mpi_comm)
      call vtk%write_headers(globaldims, mydims)
        call vtk%startcoordinates()
          call vtk%write_data_array("x_coordinate", xcoords)
          call vtk%write_data_array("y_coordinate", ycoords)
          call vtk%write_data_array("z_coordinate", zcoords)
        call vtk%finishcoordinates()
        call vtk%startpointdata()
          call vtk%write_data_array(scalarname, scalarvalues)
          call vtk%write_data_array(vectorname, vectorvalues(:,:,:,1), vectorvalues(:,:,:,2), vectorvalues(:,:,:,3))
          ! no point data here
        call vtk%finishpointdata()
        call vtk%startcelldata()
          ! no cell data here
        call vtk%finishcelldata()
      call vtk%write_final()
    call vtk%close()
  end subroutine


  !>
  !> Writes two scalar fields on a regular rectangular grid to VTK.
  !> dens1 and dens2 are intepreted as cell data instead as point data as above
  !> accordingly, size(densX) = size([xcoords, ycoords, zcoords])-1
  !>
  subroutine vtk_write_densities_on_grid(filename, step, tsim, vtk_step, globaldims, mydims, xcoords, ycoords, zcoords, &
                    dens1, name1, dens2, name2, mpi_comm, coord_scale)
    use module_vtk
    implicit none
    include 'mpif.h'
    character(*), intent(in) :: filename, name1, name2
    integer, intent(in) :: step
    integer, intent(in) :: vtk_step
    real*8, intent(in) :: tsim
    integer, dimension(2,3), intent(in) :: globaldims, mydims
    real*8, intent(in) :: xcoords(:), ycoords(:), zcoords(:)
    real*8, intent(in) :: dens1(:,:,:), dens2(:,:,:)
    integer, intent(in) :: mpi_comm
    real*8, intent(in), optional :: coord_scale
    
    integer :: mpi_rank, mpi_size, ierr

    type(vtkfile_rectilinear_grid) :: vtk
    
    call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
    call MPI_Comm_size(mpi_comm, mpi_size, ierr)

    call vtk%create_parallel(trim(filename), step, mpi_rank, mpi_size, tsim, vtk_step)
      call vtk%set_communicator(mpi_comm)
      call vtk%write_headers(globaldims, mydims)
        call vtk%startcoordinates()
          call vtk%write_data_array("x_coordinate", xcoords, coord_scale)
          call vtk%write_data_array("y_coordinate", ycoords, coord_scale)
          call vtk%write_data_array("z_coordinate", zcoords, coord_scale)
        call vtk%finishcoordinates()
        call vtk%startpointdata()
          ! no point data here
        call vtk%finishpointdata()
        call vtk%startcelldata()
          call vtk%write_data_array(name1, dens1)
          call vtk%write_data_array(name2, dens2)
        call vtk%finishcelldata()
      call vtk%write_final()
    call vtk%close()
  end subroutine

end module module_vtk_helpers
