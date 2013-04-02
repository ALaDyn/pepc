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

  public :: tree_walk
  public :: tree_walk_prepare
  public :: tree_walk_finalize
  public :: tree_walk_read_parameters
  public :: tree_walk_write_parameters
  public :: tree_walk_statistics

  contains

  subroutine tree_walk_statistics(t, fh, perform_output)
    use module_tree, only: t_tree
    implicit none
    
    type(t_tree), intent(in) :: t
    integer, intent(in) :: fh
    logical, intent(in) :: perform_output

    if (perform_output) then
      write (20, *) "module_walk_simple: no statistics available."
    end if
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
    use module_pepc_types, only: t_particle, t_tree_node
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

    type(t_tree_node), pointer :: r
    type(t_particle), pointer :: p
    integer :: i
    integer*8 :: ni
    logical :: in_central_box
    real*8 :: b2(0:nlev)

    call pepc_status('WALK SIMPLE')

    twalk = - MPI_WTIME()
    in_central_box = dot_product(vbox, vbox) == 0.0_8

    b2(0) = maxval(t%bounding_box%boxsize)**2
    do i = 1, nlev
      b2(i) = b2(i - 1) / 4.0_8
    end do

    call tree_lookup_root(t, r)

    do i = lbound(ps, 1), ubound(ps, 1)
      p => ps(i)
      ni = 0_8
      call tree_walk_single(r)
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
      use module_tree, only: tree_node_get_first_child, tree_node_get_next_sibling
      use module_tree_node, only: tree_node_children_available, tree_node_is_leaf
      implicit none

      type(t_tree_node), intent(inout) :: n

      type(t_tree_node), pointer :: s, ns
      real*8 :: d2, d(3)
      logical :: is_leaf, is_same_particle, is_ancestor

      d = (p%x - vbox) - n%interaction_data%coc
      d2 = dot_product(d, d)

      is_leaf = tree_node_is_leaf(n)
      is_ancestor = in_central_box .and. is_ancestor_of_particle(p%key, n%key, n%level)
      is_same_particle = is_ancestor .and. is_leaf
      is_ancestor = is_ancestor .and. (.not. is_leaf)

      if (is_same_particle) then ! ignore same particle
        ni = ni + 1
      else if ((.not. is_ancestor) .and. (is_leaf .or. mac(p, n%interaction_data, d2, b2(n%level)))) then ! mac ok -> interact
        call calc_force_per_interaction(p, n%interaction_data, n%key, d, d2, vbox, is_leaf)
        ni = ni + n%leaves
      else ! mac fails -> resolve
        if (.not. tree_node_children_available(n)) then
          call tree_node_fetch_children(t, n)
          do
            if (tree_node_children_available(n)) then
              call atomic_read_barrier()
              exit
            end if
          end do
        end if

        if (tree_node_get_first_child(t, n, ns)) then
          do
            s => ns
            call tree_walk_single(s)
            if (.not. tree_node_get_next_sibling(t, s, ns)) then
              exit
            end if
          end do
        else
          DEBUG_ERROR(*, "walk_simple: unexpectedly, this twig had no children")
        end if
      end if
    end subroutine tree_walk_single
  end subroutine tree_walk
end module module_walk
