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
  implicit none
  private

  real*8, public :: interactions_local, mac_evaluations_local

  public :: tree_walk
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


  subroutine tree_walk(t, ps, twalk, twalk_loc, vbox)
    use module_pepc_types, only: t_particle, kind_node
    use module_tree, only: t_tree
    use treevars, only: nlev
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), target, intent(inout) :: t !< tree to be traversed
    type(t_particle), target, intent(inout) :: ps(:) !< list of particles
    real*8, intent(inout) :: twalk !< time until completion
    real*8, intent(inout) :: twalk_loc !< also time until completion
    real*8, intent(in) :: vbox(3) !< lattice vector

    type(t_particle), pointer :: p
    integer :: i
    integer(kind_node) :: ni
    logical :: in_central_box
    real*8 :: b2(0:nlev)

    call pepc_status('WALK SIMPLE')

    twalk = - MPI_WTIME()
    interactions_local = 0.0_8
    mac_evaluations_local = 0.0_8
    in_central_box = dot_product(vbox, vbox) == 0.0_8

    b2(0) = maxval(t%bounding_box%boxsize)**2
    do i = 1, nlev
      b2(i) = b2(i - 1) / 4.0_8
    end do

    do i = lbound(ps, 1), ubound(ps, 1)
      p => ps(i)
      ni = 0_kind_node
      call tree_walk_single(t%node_root)
      DEBUG_ASSERT(ni == t%npart)
    end do

    twalk = MPI_WTIME() + twalk
    twalk_loc = twalk

    contains

    recursive subroutine tree_walk_single(n)
      use module_atomic_ops, only: atomic_read_barrier
      use module_spacefilling, only: is_ancestor_of_particle
      use module_interaction_specific, only: mac, calc_force_per_interaction
      use module_tree_communicator, only: tree_node_fetch_children
      use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
        tree_node_get_first_child, tree_node_get_next_sibling
      use pthreads_stuff, only: pthreads_sched_yield
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
      use module_mirror_boxes, only: spatial_interaction_cutoff
      #endif
      implicit none

      integer(kind_node), intent(in) :: n

      integer(kind_node) :: s, ns
      real*8 :: d2, d(3)
      logical :: is_leaf

      is_leaf = tree_node_is_leaf(t%nodes(n))

      d = (p%x - vbox) - t%nodes(n)%interaction_data%coc
      d2 = dot_product(d, d)

      if (is_leaf) then
        if (d2 > 0.0_8) then ! not self, interact
          goto 1
        else ! self, count as interaction partner, otherwise ignore
          goto 2
        end if
      else ! not a leaf, evaluate MAC
        mac_evaluations_local = mac_evaluations_local + 1.0_8
        if (mac(p, t%nodes(n)%interaction_data, d2, b2(t%nodes(n)%level))) then ! MAC OK: interact
          goto 1
        else ! MAC fails: resolve
          goto 3
        end if
      end if

      DEBUG_ASSERT_MSG(.false., *, "The block of ifs above should be exhaustive!")
      return

      ! interact
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
1     if (all(abs(d) < spatial_interaction_cutoff)) then
      #endif
        call calc_force_per_interaction(p, t%nodes(n)%interaction_data, n, d, d2, vbox, is_leaf)
        interactions_local = interactions_local + 1.0_8
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
      end if
      #endif
      ! count partner
2     ni = ni + t%nodes(n)%leaves
      return

      ! resolve
3     if (.not. tree_node_children_available(t%nodes(n))) then
        call tree_node_fetch_children(t, t%nodes(n), n)
        do
          ERROR_ON_FAIL(pthreads_sched_yield())
          if (tree_node_children_available(t%nodes(n))) then
            call atomic_read_barrier()
            exit
          end if
        end do
      end if

      if (tree_node_get_first_child(t%nodes(n), ns)) then
        do
          s = ns
          call tree_walk_single(s)
          if (.not. tree_node_get_next_sibling(t%nodes(s), ns)) then
            exit
          end if
        end do
      else
        DEBUG_ERROR(*, "walk_simple: unexpectedly, this twig had no children")
      end if
      return
    end subroutine tree_walk_single
  end subroutine tree_walk
end module module_walk
