! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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
!> A simple dual tree traversal
!>
module module_dual_tree_walk
  implicit none

  public :: dual_tree_walk_sow
  public :: dual_tree_walk_reap
  public :: dual_tree_walk_prepare
  public :: dual_tree_walk_finalize
  public :: dual_tree_walk_read_parameters
  public :: dual_tree_walk_write_parameters
  public :: dual_tree_walk_statistics

  contains

  subroutine dual_tree_walk_statistics(u)
    use treevars, only: me
    implicit none

    integer, intent(in) :: u

    if (0 == me) then; write (u, *) "module_walk_dual: no statistics for now."; end if
  end subroutine dual_tree_walk_statistics


  subroutine dual_tree_walk_read_parameters(fh)
    implicit none
    integer, intent(in) :: fh
  end subroutine dual_tree_walk_read_parameters


  subroutine dual_tree_walk_write_parameters(fh)
    implicit none
    integer, intent(in) :: fh
  end subroutine dual_tree_walk_write_parameters


  subroutine dual_tree_walk_prepare()
    implicit none
  end subroutine dual_tree_walk_prepare


  subroutine dual_tree_walk_finalize()
    implicit none
  end subroutine dual_tree_walk_finalize


  subroutine dual_tree_walk_sow(t, lattice_vector)
    use module_pepc_kinds, only: kind_physics, kind_node
    use module_tree, only: t_tree
    use module_pepc_types, only: t_tree_node
    use module_interaction_specific
    use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
      tree_node_get_first_child, tree_node_get_next_sibling, NODE_INVALID
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    real(kind_physics), intent(in) :: lattice_vector(3) !< lattice vector

    call pepc_status('WALK DUAL SOW')

    call sow_aux(t%node_root, t%node_root)

    contains

    recursive subroutine sow_aux(src, dst)
      implicit none
      integer(kind_node), intent(in) :: src, dst

      logical :: src_is_leaf, dst_is_leaf
      real(kind_physics) :: shifted_src_center(3)

      associate (src_node => t%nodes(src), dst_node => t%nodes(dst), &
        rsrc => t%boxdiaglength(t%nodes(src)%level), rdst => t%boxdiaglength(t%nodes(dst)%level))

        src_is_leaf = tree_node_is_leaf(src_node)
        dst_is_leaf = tree_node_is_leaf(dst_node)

        shifted_src_center = src_node%center - lattice_vector

        if ((src_is_leaf .and. dst_is_leaf) &
            .or. dual_mac(shifted_src_center, rsrc, src_node%multipole_moments, dst_node%center, rdst, dst_node%multipole_moments)) then ! TODO: Leafs and well separated nodes together?
          call multipole_to_local(shifted_src_center, src_node%multipole_moments, dst_node%center, dst_node%local_coefficients)
        else if (dst_is_leaf) then
          call split_src(src, dst)
        else if (src_is_leaf) then
          call split_dst(src, dst)
        else if (rsrc > rdst) then
          call split_src(src, dst)
        else
          call split_dst(src, dst)
        end if

      end associate

    end subroutine sow_aux

    recursive subroutine split_src(src, dst)
      use pthreads_stuff, only: pthreads_sched_yield
      use module_tree_communicator, only: tree_node_fetch_children
      use module_atomic_ops, only: atomic_read_barrier
      implicit none
      integer(kind_node), intent(in) :: src, dst

      integer(kind_node) :: ns

      associate (src_node => t%nodes(src))
        ! source node to be split may never be a leaf, leaves cannot be split
        DEBUG_ASSERT(.not. tree_node_is_leaf(src_node))

        if (.not. tree_node_children_available(src_node)) then
          call tree_node_fetch_children(t, src_node, src)
          do ! loop and yield until children have been fetched
            ERROR_ON_FAIL(pthreads_sched_yield())
            if (tree_node_children_available(src_node)) exit
            !!$omp taskyield
          end do
        end if

        call atomic_read_barrier()

        ns = tree_node_get_first_child(src_node)
        do
          call sow_aux(ns, dst)
          ns = tree_node_get_next_sibling(t%nodes(ns))
          if (ns == NODE_INVALID) exit
        end do

      end associate

    end subroutine split_src

    recursive subroutine split_dst(src, dst)
      use module_tree_node, only: tree_node_has_local_contributions
      implicit none
      integer(kind_node), intent(in) :: src, dst

      integer(kind_node) :: ns

      associate (dst_node => t%nodes(dst))
        ! destination node to be split may never be a leaf, leaves cannot be split
        DEBUG_ASSERT(.not. tree_node_is_leaf(dst_node))
        ! destination nodes must only come from the local part of the tree, results are only computed for local leaves
        DEBUG_ASSERT(tree_node_has_local_contributions(dst_node))

        ! since destination nodes are part of the local tree, their children (in fact all descendants) must be available
        DEBUG_ASSERT(tree_node_children_available(dst_node))
        ns = tree_node_get_first_child(dst_node)
        do
          if (tree_node_has_local_contributions(t%nodes(ns))) call sow_aux(src, ns)
          ns = tree_node_get_next_sibling(t%nodes(ns))
          if (ns == NODE_INVALID) exit
        end do

      end associate

    end subroutine split_dst

  end subroutine dual_tree_walk_sow

  subroutine dual_tree_walk_reap(t, p)
    use module_pepc_kinds, only: kind_physics, kind_node, kind_particle
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle, t_tree_node
    use module_interaction_specific
    use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
      tree_node_get_first_child, tree_node_get_next_sibling, tree_node_get_particle, tree_node_has_local_contributions, NODE_INVALID
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)

    call pepc_status('WALK DUAL REAP')

    call reap_aux(t%node_root)

    contains

    recursive subroutine reap_aux(n)
      implicit none
      integer(kind_node), intent(in) :: n

      type(t_tree_node), pointer :: node
      integer(kind_node) :: ns
      integer(kind_particle) :: ps

      associate (node => t%nodes(n))
        ! destination nodes must only come from the local part of the tree, results are only computed for local leaves
        DEBUG_ASSERT(tree_node_has_local_contributions(node))

        if (tree_node_is_leaf(node)) then
          ps = tree_node_get_particle(node)
          call evaluate_at_particle(node%local_coefficients, p(ps)%results)
        else
          ! since destination nodes are part of the local tree, their children (in fact all descendants) must be available
          DEBUG_ASSERT(tree_node_children_available(node))
          ns = tree_node_get_first_child(node)
          do
            if (tree_node_has_local_contributions(t%nodes(ns))) then
              call shift_coefficients_down(node%center, node%local_coefficients, t%nodes(ns)%center, t%nodes(ns)%local_coefficients)
              call reap_aux(ns)
            end if
            ns = tree_node_get_next_sibling(t%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        end if

      end associate

    end subroutine reap_aux

  end subroutine dual_tree_walk_reap

end module module_dual_tree_walk
