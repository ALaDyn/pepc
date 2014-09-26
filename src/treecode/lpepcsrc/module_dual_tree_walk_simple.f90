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
  private

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

    if (0 == me) then; write (u, *) "dual_tree_walk_simple: no statistics for now."; end if
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


  !>
  !> Perform a dual tree traversal of tree `t` through itself. `p` is the list of particles from which `t` has been built.
  !>
  !> This has to be called once for every neighbor box with the associated `lattice_vector`.
  !>
  !> The "sow" part of the algorithm identifies pairs `(src, dst)` of interaction partners in `t` as per the multipole acceptance
  !> criterion and populates the local coefficients of `dst` if `dst` is a tree node or directly computes an interaction result
  !> if `dst` is a particle.
  !>
  subroutine dual_tree_walk_sow(t, p, lattice_vector)
    use module_pepc_kinds, only: kind_physics
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use treevars, only: num_threads
    use pthreads_stuff, only: place_thread, THREAD_TYPE_WORKER, THREAD_TYPE_MAIN
    use omp_lib, only: omp_get_thread_num
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)
    real(kind_physics), intent(in) :: lattice_vector(3)

    call pepc_status('WALK DUAL SOW')

    !$omp parallel default(shared) num_threads(num_threads)
    call place_thread(THREAD_TYPE_WORKER, omp_get_thread_num() + 1)
    !$omp master
    call sow_aux(t, p, lattice_vector, t%node_root, t%node_root)
    !$omp end master
    !$omp end parallel

    call place_thread(THREAD_TYPE_MAIN, 0)
  end subroutine dual_tree_walk_sow


  !>
  !> Consider two nodes `src` and `dst` for interaction.
  !>
  recursive subroutine sow_aux(t, p, lattice_vector, src, dst)
    use module_pepc_kinds, only: kind_physics, kind_node
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use module_interaction_specific
    use module_tree_node, only: tree_node_get_particle, tree_node_is_leaf
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)
    real(kind_physics), intent(in) :: lattice_vector(3)
    integer(kind_node), intent(in) :: src, dst

    logical :: src_is_leaf, dst_is_leaf
    real(kind_physics) :: delta(3), dist2, rsrc, rdst

    associate (src_node => t%nodes(src), dst_node => t%nodes(dst))
      src_is_leaf = tree_node_is_leaf(src_node)
      dst_is_leaf = tree_node_is_leaf(dst_node)

      ! The next two blocks assume that leaves contain a single particle.
      if (src_is_leaf) then
        rsrc = 0.0_kind_physics
      else
        rsrc = t%boxdiaglength(src_node%level)
      end if

      if (dst_is_leaf) then
        rdst = 0.0_kind_physics
      else
        rdst = t%boxdiaglength(dst_node%level)
      end if

      delta = dst_node%center - (src_node%center - lattice_vector)
      dist2 = dot_product(delta, delta)

      if (src_is_leaf .and. dst_is_leaf) then ! A pair of leaves interacts directly via P2P.
        associate(dst_particle => p(tree_node_get_particle(dst_node)))
          dst_particle%work = dst_particle%work + 1.0_8
          if (src /= dst) then
            call calc_force_per_interaction_with_leaf(dst_particle, src_node%multipole_moments, src, delta, dist2, lattice_vector)
          else
            call calc_force_per_interaction_with_self(dst_particle, src_node%multipole_moments, src, delta, dist2, lattice_vector)
          end if
        end associate
      else if (dual_mac(delta, dist2, rsrc, src_node%multipole_moments, rdst, dst_node%multipole_moments)) then ! Approximate...
        if (dst_is_leaf) then ! ...via P2M for destinations that are leaves, or...
          associate(dst_particle => p(tree_node_get_particle(dst_node)))
            dst_particle%work = dst_particle%work + 5.0_8
            call calc_force_per_interaction_with_twig(dst_particle, src_node%multipole_moments, src, delta, dist2, lattice_vector)
          end associate
        else ! ...via M2L for interactions between twigs.
          dst_node%work = dst_node%work + 20.0_8
          call multipole_to_local(delta, dist2, src_node%multipole_moments, dst_node%local_coefficients)
        end if
      else if (dst_is_leaf) then ! Destination cannot be split, so split source.
        call sow_split_src(t, p, lattice_vector, src, dst)
      else if (src_is_leaf) then ! Source cannot be split, so split destination.
        call sow_split_dst(t, p, lattice_vector, src, dst)
      else if (rsrc > rdst) then ! Source is bigger, so split it.
        call sow_split_src(t, p, lattice_vector, src, dst)
      else ! Destination is bigger, so split it.
        call sow_split_dst(t, p, lattice_vector, src, dst)
      end if
    end associate
  end subroutine sow_aux


  !>
  !> `src` and `dst` may not interact. Consider `(c, dst)` for all children `c` of `src` for interaction.
  !>
  recursive subroutine sow_split_src(t, p, lattice_vector, src, dst)
    use module_pepc_kinds, only: kind_physics, kind_node
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use pthreads_stuff, only: pthreads_sched_yield
    use module_tree_communicator, only: tree_node_fetch_children
    use module_atomic_ops, only: atomic_read_barrier
    use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
      tree_node_get_first_child, tree_node_get_next_sibling, NODE_INVALID
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)
    real(kind_physics), intent(in) :: lattice_vector(3)
    integer(kind_node), intent(in) :: src, dst

    integer(kind_node) :: ns

    associate (src_node => t%nodes(src))
      ! Source node to be split may never be a leaf, leaves cannot be split.
      DEBUG_ASSERT(.not. tree_node_is_leaf(src_node))

      if (.not. tree_node_children_available(src_node)) then
        call tree_node_fetch_children(t, src_node, src)
        do ! Loop and yield until children have been fetched.
          ERROR_ON_FAIL(pthreads_sched_yield())
          if (tree_node_children_available(src_node)) exit
          !$omp taskyield
        end do
      end if

      call atomic_read_barrier()

      ns = tree_node_get_first_child(src_node)
      do
        call sow_aux(t, p, lattice_vector, ns, dst)
        ns = tree_node_get_next_sibling(t%nodes(ns))
        if (ns == NODE_INVALID) exit
      end do
    end associate
  end subroutine sow_split_src


  !>
  !> `src` and `dst` may not interact. Consider, in parallel, `(src, c)` for all children `c` of `dst` for interaction.
  !>
  recursive subroutine sow_split_dst(t, p, lattice_vector, src, dst)
    use module_pepc_kinds, only: kind_physics, kind_node, kind_particle
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, tree_node_has_local_contributions, &
      tree_node_get_first_child, tree_node_get_next_sibling, NODE_INVALID
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)
    real(kind_physics), intent(in) :: lattice_vector(3)
    integer(kind_node), intent(in) :: src, dst

    integer(kind_node) :: ns

    associate (dst_node => t%nodes(dst))
      ! Destination node to be split may never be a leaf, leaves cannot be split.
      DEBUG_ASSERT(.not. tree_node_is_leaf(dst_node))
      ! Destination nodes must only come from the local part of the tree, results are only computed for local leaves.
      DEBUG_ASSERT(tree_node_has_local_contributions(dst_node))

      ! Since destination nodes are part of the local tree, their children (in fact all descendants) must be available.
      DEBUG_ASSERT(tree_node_children_available(dst_node))
      ns = tree_node_get_first_child(dst_node)
      do
        if (tree_node_has_local_contributions(t%nodes(ns))) then
          if (t%nodes(ns)%leaves >= 100) then
            !$omp task default(none) firstprivate(ns) shared(t, p, lattice_vector, src) untied
            call sow_aux(t, p, lattice_vector, src, ns)
            !$omp end task
          else
            call sow_aux(t, p, lattice_vector, src, ns)
          end if
        end if
        ns = tree_node_get_next_sibling(t%nodes(ns))
        if (ns == NODE_INVALID) exit
      end do
    end associate
    !$omp taskwait
  end subroutine sow_split_dst


  !>
  !> Perform a top-down traversal of tree `t` to shift local coefficients down towards the leaves and finally evaluate them at the
  !> particles `p`.
  !>
  subroutine dual_tree_walk_reap(t, p)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use treevars, only: num_threads
    use pthreads_stuff, only: place_thread, THREAD_TYPE_WORKER, THREAD_TYPE_MAIN
    use omp_lib, only: omp_get_thread_num
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)

    call pepc_status('WALK DUAL REAP')

    !$omp parallel default(shared) num_threads(num_threads)
    call place_thread(THREAD_TYPE_WORKER, omp_get_thread_num() + 1)
    !$omp master
    call reap_aux(t, p, t%node_root)
    !$omp end master
    !$omp end parallel

    call place_thread(THREAD_TYPE_MAIN, 0)
  end subroutine dual_tree_walk_reap


  !>
  !> Shift local coefficients of node `n` to its children and recurse in parallel.
  !>
  recursive subroutine reap_aux(t, p, n)
    use module_tree, only: t_tree
    use module_pepc_types, only: t_particle
    use module_pepc_kinds, only: kind_physics, kind_node, kind_particle
    use module_interaction_specific
    use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
      tree_node_get_first_child, tree_node_get_next_sibling, tree_node_get_particle, tree_node_has_local_contributions, NODE_INVALID
    use module_debug
    implicit none

    type(t_tree), intent(inout) :: t
    type(t_particle), intent(inout) :: p(:)
    integer(kind_node), intent(in) :: n

    integer(kind_node) :: ns
    real(kind_physics) :: delta(3)

    associate (node => t%nodes(n))
      ! Destination nodes must only come from the local part of the tree, results are only computed for local leaves.
      DEBUG_ASSERT(tree_node_has_local_contributions(node))

      if (tree_node_is_leaf(node)) then
        associate(particle => p(tree_node_get_particle(node)))
          particle%work = particle%work + node%work
          call evaluate_at_particle(node%local_coefficients, particle%results)
        end associate
      else
        ! Since destination nodes are part of the local tree, their children (in fact all descendants) must be available.
        DEBUG_ASSERT(tree_node_children_available(node))
        ns = tree_node_get_first_child(node)
        do
          associate(child_node => t%nodes(ns))
            if (tree_node_has_local_contributions(child_node)) then
              if (child_node%leaves >= 100) then
                !$omp task default(none) firstprivate(ns, child_node) private(delta) shared(t, p, node) untied
                child_node%work = child_node%work + real(child_node%leaves, kind = 8) * node%work / real(node%leaves, kind = 8)
                delta = child_node%center - node%center
                call shift_coefficients_down(delta, node%local_coefficients, child_node%local_coefficients)
                call reap_aux(t, p, ns)
                !$omp end task
              else
                child_node%work = child_node%work + real(child_node%leaves, kind = 8) * node%work / real(node%leaves, kind = 8)
                delta = child_node%center - node%center
                call shift_coefficients_down(delta, node%local_coefficients, child_node%local_coefficients)
                call reap_aux(t, p, ns)
              end if
            end if
          end associate
          ns = tree_node_get_next_sibling(t%nodes(ns))
          if (ns == NODE_INVALID) exit
        end do
      end if
      !$omp taskwait
    end associate
  end subroutine reap_aux
end module module_dual_tree_walk
