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

  public :: dual_tree_walk_run
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


  subroutine dual_tree_walk_run(t, p, vbox)
    use module_pepc_kinds, only: kind_physics, kind_node, kind_particle
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle, t_tree_node
    use module_interaction_specific
    use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
      tree_node_get_first_child, tree_node_get_next_sibling, tree_node_get_particle, NODE_INVALID
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)
    real(kind_physics), intent(in) :: vbox(3) !< lattice vector

    real(kind_physics) :: boxsizelength

    call pepc_status('WALK DUAL')

    boxsizelength = sqrt(dot_product(t%bounding_box%boxsize, t%bounding_box%boxsize))

    call dual_tree_walk_run_sow(t%node_root, t%node_root)

    call dual_tree_walk_run_reap(t%node_root)

    contains

      recursive subroutine dual_tree_walk_run_sow(src, dst)
        implicit none
        integer(kind_node), intent(in) :: src, dst

        type(t_tree_node), pointer :: src_node, dst_node
        logical :: src_is_leaf, dst_is_leaf
        real(kind_physics) :: shifted_src_center(3), rsrc, rdst
        integer(kind_node) :: ns

        src_node => t%nodes(src)
        dst_node => t%nodes(dst)
        src_is_leaf = tree_node_is_leaf(src_node)
        dst_is_leaf = tree_node_is_leaf(dst_node)

        ! TODO: Optimization
        rsrc = 0.5**(src_node%level) * boxsizelength
        rdst = 0.5**(dst_node%level) * boxsizelength

        shifted_src_center = src_node%center - vbox

        if ((src_is_leaf .and. dst_is_leaf) &
            .or. dual_mac(shifted_src_center, rsrc, src_node%multipole_moments, dst_node%center, rdst, dst_node%multipole_moments)) then ! TODO: Leafs and well separated nodes together?
          call multipole_to_local(shifted_src_center, src_node%multipole_moments, dst_node%center, dst_node%local_coefficients)
        else if (dst_is_leaf) then
          ns = tree_node_get_first_child(src_node)
          do
            call dual_tree_walk_run_sow(ns, dst)
            ns = tree_node_get_next_sibling(t%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        else if (src_is_leaf) then
          ns = tree_node_get_first_child(dst_node)
          do
            call dual_tree_walk_run_sow(src, ns)
            ns = tree_node_get_next_sibling(t%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        else if (rsrc > rdst) then 
          ns = tree_node_get_first_child(src_node)
          do
            call dual_tree_walk_run_sow(ns, dst)
            ns = tree_node_get_next_sibling(t%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        else 
          ns = tree_node_get_first_child(dst_node)
          do
            call dual_tree_walk_run_sow(src, ns)
            ns = tree_node_get_next_sibling(t%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        end if
 
      end subroutine dual_tree_walk_run_sow


      recursive subroutine dual_tree_walk_run_reap(n)
        implicit none
        integer(kind_node), intent(inout) :: n

        type(t_tree_node), pointer :: node
        integer(kind_node) :: ns
        integer(kind_particle) :: ps

        node => t%nodes(n)

        if (tree_node_is_leaf(node)) then
          ps = tree_node_get_particle(node)
          call evaluate_at_particle(node%local_coefficients, p(ps)%results)
        else 
          ns = tree_node_get_first_child(node)
          do
            call shift_coefficients_down(node%center, node%local_coefficients, t%nodes(ns)%center, t%nodes(ns)%local_coefficients)
            call dual_tree_walk_run_reap(ns)
            ns = tree_node_get_next_sibling(t%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        end if

      end subroutine dual_tree_walk_run_reap

  end subroutine dual_tree_walk_run

end module module_dual_tree_walk
