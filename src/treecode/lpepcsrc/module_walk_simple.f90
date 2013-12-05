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
!> A simple force calculation / tree traversal
!>
!> The algorithm is single threaded and blocks while waiting for remote nodes to
!> arrive.
!>
module module_walk
  use module_tree, only: t_tree
  use module_pepc_types, only: t_particle
  implicit none
  private

  type(t_tree), pointer :: walk_tree
  type(t_particle), pointer :: walk_particles(:)
  integer :: num_walk_threads

  real*8, public :: interactions_local, mac_evaluations_local

  public :: tree_walk_run
  public :: tree_walk_init
  public :: tree_walk_uninit
  public :: tree_walk_prepare
  public :: tree_walk_finalize
  public :: tree_walk_read_parameters
  public :: tree_walk_write_parameters
  public :: tree_walk_statistics

  contains

  subroutine tree_walk_statistics(u)
    use treevars, only: me
    implicit none
    
    integer, intent(in) :: u

    if (0 == me) then; write (u, *) "module_walk_simple: no statistics for now."; end if
  end subroutine tree_walk_statistics


  subroutine tree_walk_read_parameters(fh)
    implicit none
    integer, intent(in) :: fh
  end subroutine tree_walk_read_parameters


  subroutine tree_walk_write_parameters(fh)
    implicit none
    integer, intent(in) :: fh
  end subroutine tree_walk_write_parameters


  subroutine tree_walk_prepare()
    implicit none
  end subroutine tree_walk_prepare


  subroutine tree_walk_finalize()
    implicit none
  end subroutine tree_walk_finalize


  subroutine tree_walk_init(t, p, num_threads)
    implicit none

    type(t_tree), target, intent(in) :: t
    type(t_particle), target, intent(in) :: p(:)
    integer, intent(in) :: num_threads

    walk_tree => t
    walk_particles => p
    num_walk_threads = num_threads
  end subroutine tree_walk_init


  subroutine tree_walk_uninit(t, p)
    implicit none

    type(t_tree), intent(in) :: t
    type(t_particle), intent(in) :: p(:)
  end subroutine tree_walk_uninit


  subroutine tree_walk_run(vbox)
    use module_pepc_types, only: kind_particle
    use treevars, only: num_threads
    use module_debug
    implicit none

    real*8, intent(in) :: vbox(3) !< lattice vector

    integer(kind_particle) :: i

    call pepc_status('WALK SIMPLE')

    interactions_local = 0.0_8
    mac_evaluations_local = 0.0_8

    !$omp parallel default(shared) private(i) num_threads(num_walk_threads)
    !$omp master
    do i = 1, size(walk_particles, kind = kind_particle)
      !$omp task default(shared) firstprivate(i)
      call tree_walk_single(walk_particles(i), vbox)
      !$omp end task
    end do
    !$omp end master
    !$omp end parallel
  end subroutine tree_walk_run


  subroutine tree_walk_single(p_, vbox)
    use module_pepc_types, only: t_particle, kind_node
    use module_debug
    use treevars, only: nlev
    implicit none

    type(t_particle), intent(inout) :: p_
    real*8, intent(in) :: vbox(3)

    type(t_particle) :: p
    integer(kind_node) :: ni
    integer :: i
    real*8 :: b2(0:nlev), num_int, num_mac

    num_int = 0.0_8
    num_mac = 0.0_8

    b2(0) = maxval(walk_tree%bounding_box%boxsize)**2
    do i = 1, nlev
      b2(i) = b2(i - 1) / 4.0_8
    end do

    p = p_ ! make a local copy of the particle (false sharing and such)

    ni = 0_kind_node
    call tree_walk_single_aux(walk_tree%node_root)
    DEBUG_ASSERT(ni == walk_tree%npart)

    !$omp critical
    interactions_local = interactions_local + num_int
    mac_evaluations_local = mac_evaluations_local + num_mac
    !$omp end critical

    p_ = p

  contains

    recursive subroutine tree_walk_single_aux(n)
      use module_pepc_types, only: t_tree_node
      use module_atomic_ops, only: atomic_read_barrier
      use module_interaction_specific
      use module_tree_communicator, only: tree_node_fetch_children
      use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
        tree_node_get_first_child, tree_node_get_next_sibling, NODE_INVALID
      use pthreads_stuff, only: pthreads_sched_yield
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
      use module_mirror_boxes, only: spatial_interaction_cutoff
      #endif
      implicit none

      integer(kind_node), intent(in) :: n

      type(t_tree_node), pointer :: node
      integer(kind_node) :: ns
      real*8 :: d2, d(3)
      logical :: is_leaf

      node => walk_tree%nodes(n)
      is_leaf = tree_node_is_leaf(node)

      d = (p%x - vbox) - node%interaction_data%coc
      d2 = dot_product(d, d)

      if (is_leaf) then
        ni = ni + 1
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (any(abs(d) >= spatial_interaction_cutoff)) return
        #endif
        num_int = num_int + 1.0_8
        if (d2 > 0.0_8) then ! not self
          call calc_force_per_interaction_with_leaf(p, node%interaction_data, n, d, d2, vbox)
        else ! self
          call calc_force_per_interaction_with_self(p, node%interaction_data, n, d, d2, vbox)
        end if
      else ! not a leaf, evaluate MAC
        num_mac = num_mac + 1.0_8
        if (mac(IF_MAC_NEEDS_PARTICLE(p) node%interaction_data, d2, b2(node%level))) then ! MAC OK: interact
          ni = ni + node%leaves
          #ifndef NO_SPATIAL_INTERACTION_CUTOFF
          if (any(abs(d) >= spatial_interaction_cutoff)) return
          #endif
          num_int = num_int + 1.0_8
          call calc_force_per_interaction_with_twig(p, node%interaction_data, n, d, d2, vbox)
        else ! MAC fails: resolve
          if (.not. tree_node_children_available(node)) then
            call tree_node_fetch_children(walk_tree, node, n)
            do ! loop and yield until children have been fetched
              ERROR_ON_FAIL(pthreads_sched_yield())
              if (tree_node_children_available(node)) exit
              !$omp taskyield
            end do
          end if

          call atomic_read_barrier()

          ns = tree_node_get_first_child(node)
          DEBUG_ASSERT_MSG(ns /= NODE_INVALID, *, "walk_simple: unexpectedly, this twig had no children")
          do
            call tree_walk_single_aux(ns)
            ns = tree_node_get_next_sibling(walk_tree%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        end if
      end if
    end subroutine tree_walk_single_aux
  end subroutine tree_walk_single
end module module_walk
