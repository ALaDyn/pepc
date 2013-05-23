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
    use module_interaction_specific_types
    use module_pepc_types, only: t_particle, kind_node
    use module_tree, only: t_tree, tree_lookup_root
    use treevars, only: nlev
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_tree), target, intent(inout) :: t !< tree to be traversed
    type(t_particle), target, intent(inout) :: ps(:) !< list of particles
    real*8, intent(inout) :: twalk !< time until completion
    real*8, intent(inout) :: twalk_loc !< also time until completion
    real*8, intent(in) :: vbox(3) !< lattice vector

    integer(kind_node) :: r
    type(t_particle_thread) :: p
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

    call tree_lookup_root(t, r)

    do i = lbound(ps, 1), ubound(ps, 1)
       p%x        = ps(i)%x 
       p%work     = ps(i)%work 
       p%key      = ps(i)%key 
       p%key_leaf = ps(i)%key_leaf 
       p%label    = ps(i)%label 
       p%pid      = ps(i)%pid 
       p%data     = ps(i)%data 
       p%results  = ps(i)%results 
       p%queued   = -1 
       p%my_idx   = i 
       ni = 0_kind_node
       call tree_walk_single(r)
       DEBUG_ASSERT(ni == t%npart)
    end do

    twalk = MPI_WTIME() + twalk
    twalk_loc = twalk

    contains

    recursive subroutine tree_walk_single(n)
      use module_atomic_ops, only: atomic_read_barrier
      use module_spacefilling, only: is_ancestor_of_particle
      use module_interaction_specific
      use module_tree_communicator, only: tree_node_fetch_children
      use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
        tree_node_get_first_child, tree_node_get_next_sibling
      use pthreads_stuff, only: pthreads_sched_yield
      use module_mirror_boxes, only: spatial_interaction_cutoff
      implicit none

      integer(kind_node), intent(in) :: n

      integer(kind_node) :: s, ns
      real*8 :: d2, d(3)
      logical :: is_leaf, is_related

      is_leaf = tree_node_is_leaf(t%nodes(n))
      is_related = in_central_box .and. is_ancestor_of_particle(p%key, t%nodes(n)%key, t%nodes(n)%level)

      if (.not. (is_leaf .or. is_related)) then ! A twig that is not an ancestor
        d = (p%x - vbox) - t%nodes(n)%interaction_data%coc
        d2 = dot_product(d, d)
        mac_evaluations_local = mac_evaluations_local + 1.0_8
        if (mac(p, t%nodes(n)%interaction_data, d2, b2(t%nodes(n)%level))) then ! MAC OK: interact
          go to 1 ! interact
        else ! MAC fails: resolve
          go to 3 ! resolve
        end if
      else if (is_leaf .and. (.not. is_related)) then ! Always interact with leaves
        d = (p%x - vbox) - t%nodes(n)%interaction_data%coc
        d2 = dot_product(d, d)
        go to 1 ! interact
      else if ((.not. is_leaf) .and. is_related) then ! This is a parent: resolve
        go to 3 ! resolve
      else if (is_leaf .and. is_related) then ! self
        go to 2 ! ignore, but count
      end if

      if(p%queued .gt. 0) call compute_iact_list(p)
      if(p%queued_l .gt. 0) call compute_iact_list_leaf(p)

      DEBUG_ASSERT_MSG(.false., *, "The block of ifs above should be exhaustive!")
      return

      ! interact
1     if (all(abs(d) < spatial_interaction_cutoff)) then
        call calc_force_per_interaction(p, t%nodes(n)%interaction_data, t%nodes(n)%key, d, d2, vbox, is_leaf)
        interactions_local = interactions_local + 1.0_8
      end if
      ! count partner
2     ni = ni + t%nodes(n)%leaves
      return

      ! resolve
3     if (.not. tree_node_children_available(t%nodes(n))) then
        call tree_node_fetch_children(t, t%nodes(n))
        do
          DEBUG_ERROR_ON_FAIL(pthreads_sched_yield())
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
