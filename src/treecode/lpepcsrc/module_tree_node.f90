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
!>  Encapsulates functions for accessing, manipulating, and verifying hash table data
!>
module module_tree_node
    use module_pepc_types
    implicit none
    private

    ! bits in flags to be set when children are requested, the request has been sent, and they have arrived
    integer, public, parameter :: TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE       = 0 !< bit is used in flags_local to denote that children information for the node is available in the local hashtable
    integer, public, parameter :: TREE_NODE_FLAG_LOCAL_HAS_LOCAL_CONTRIBUTIONS  = 2 !< bit is set in flags_local for all nodes that contain some local nodes beneath them
    integer, public, parameter :: TREE_NODE_FLAG_LOCAL_HAS_REMOTE_CONTRIBUTIONS = 3 !< bit is set in flags_local for all nodes that contain some remote nodes beneath them
    integer, public, parameter :: TREE_NODE_FLAG_GLOBAL_IS_BRANCH_NODE           = 0 !< bit is set in flags_global for all branch nodes (set in tree_exchange)
    integer, public, parameter :: TREE_NODE_FLAG_GLOBAL_IS_FILL_NODE             = 1 !< bit is set in flags_global for all nodes that are above (towards root) branch nodes

    integer(kind_node), public, parameter :: NODE_INVALID = -1

    public tree_node_get_first_child
    public tree_node_get_next_sibling
    public tree_node_is_leaf
    public tree_node_is_root
    public tree_node_children_available
    public tree_node_get_num_children
    public tree_node_get_num_preceding_children
    public tree_node_get_childkeys
    public tree_node_has_child
    public tree_node_pack
    public tree_node_unpack

    contains

    !>
    !> Returns the first child of node `p`.
    !>
    !> If `p` has children that are locally available, returns `.true.` and `fc`
    !> points to the child.
    !> Otherwise, `.false.` is returned and `fc` points to `null()`.
    !>
    function tree_node_get_first_child(p, fc)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_spacefilling, only: child_key_from_parent_key
      use module_debug
      implicit none

      logical :: tree_node_get_first_child
      type(t_tree_node), intent(in) :: p
      integer(kind_node), intent(out) :: fc

      if (tree_node_children_available(p)) then
        fc = p%first_child
      else
        fc = NODE_INVALID
      endif
      tree_node_get_first_child = (fc .ne. NODE_INVALID)
    end function tree_node_get_first_child


    !>
    !> Returns the next sibling of node `p`.
    !>
    !> If there is a next sibling to `n` returns `.true.` and `s` points
    !> to the sibling node.
    !> Otherwise `.false.` is returned and `s` points to `null()`.
    !>
    !> In this context, "next" is defined by the ordering of the node keys.
    !>
    function tree_node_get_next_sibling(n, s)
      use module_pepc_types, only: t_tree_node, kind_node
      use module_debug
      implicit none

      logical :: tree_node_get_next_sibling
      type(t_tree_node), intent(in) :: n
      integer(kind_node), intent(out) :: s

      s = n%next_sibling
      tree_node_get_next_sibling = (s .ne. NODE_INVALID)
    end function tree_node_get_next_sibling


    !>
    !> checks whether `n` is a leaf or twig node
    !>
    function tree_node_is_leaf(n)
      implicit none
      type(t_tree_node), intent(in) :: n
      logical :: tree_node_is_leaf

      tree_node_is_leaf = (0 == n%childcode)
    end function tree_node_is_leaf


    !>
    !> returns the number of children between (including) first_child 
    !> and (including) child `c` below `n`
    !>
    function tree_node_get_num_preceding_children(n, c)
      implicit none
      type(t_tree_node), intent(in) :: n
      integer(kind_byte), intent(in) :: c
      integer(kind_byte) :: tree_node_get_num_preceding_children
      
      tree_node_get_num_preceding_children = popcnt(ishft(n%childcode, bit_size(n%childcode)-(c+1)))-1
    end function
      

    !>
    !> checks whether `n` is a root node
    !>
    function tree_node_is_root(n)
      implicit none
      type(t_tree_node), intent(in) :: n

      logical :: tree_node_is_root

      tree_node_is_root = n%level == 0
    end function tree_node_is_root


    !>
    !> checks whether the children of `n` are locally available or have
    !> to be requested from remote ranks
    !>
    function tree_node_children_available(n)
      implicit none
      type(t_tree_node), intent(in) :: n
      logical :: tree_node_children_available

      tree_node_children_available = btest(n%flags_local, TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE)
    end function tree_node_children_available


    !>
    !> Returns `.true.` if tree node `n` has a child with number `i`.
    !>
    !> Child number as defined by `child_number_from_key()`.
    !>
    function tree_node_has_child(n, i)
      implicit none

      logical :: tree_node_has_child
      type(t_tree_node), intent(in) :: n
      integer, intent(in) :: i

      tree_node_has_child = btest(n%childcode, i)
    end function tree_node_has_child


    !>
    !> returns the number of children of node `n`
    !>
    function tree_node_get_num_children(n) result(res)
      implicit none

      integer :: res
      type(t_tree_node), intent(in) :: n

      res = popcnt(n%childcode)
    end function


    !>
    !> returns the keys of all children, that are attached to the
    !> node `n`
    !>
    subroutine tree_node_get_childkeys(n, childnum, childkeys)
      use treevars, only: idim
      use module_spacefilling, only: shift_key_by_level
      implicit none
      type(t_tree_node), intent(in) :: n
      integer, intent(out) :: childnum
      integer(kind_key), dimension(:), intent(out) :: childkeys

      integer   :: i
      integer(kind_key) :: keyhead

      keyhead   = shift_key_by_level(n%key, 1_kind_level)
      childnum = 0

      do i = 0, 2**idim - 1
        if (tree_node_has_child(n, i)) then
          childnum            = childnum + 1
          childkeys(childnum) = ior(keyhead, 1_kind_key*i)
        end if
      end do
    end subroutine tree_node_get_childkeys


    subroutine tree_node_pack(n, p)
      use module_pepc_types, only: t_tree_node_package
      implicit none

      type(t_tree_node), intent(in) :: n
      type(t_tree_node_package), intent(out) :: p

      p%key              = n%key
      p%childcode        = n%childcode
      p%flags_global     = n%flags_global
      p%level            = n%level
      p%owner            = n%owner
      p%leaves           = n%leaves
      p%descendants      = n%descendants
      p%parent           = NODE_INVALID ! this is to be filled by answer_request()
      p%first_child      = n%first_child
      p%interaction_data = n%interaction_data
      
    end subroutine tree_node_pack


    subroutine tree_node_unpack(p, n)
      use module_pepc_types, only: t_tree_node_package
      implicit none

      type(t_tree_node_package), intent(in) :: p
      type(t_tree_node), intent(out) :: n

      n%key              = p%key
      n%childcode        = p%childcode
      n%flags_global     = p%flags_global
      n%level            = p%level
      n%flags_local      = 0_kind_byte
      n%owner            = p%owner
      n%leaves           = p%leaves
      n%descendants      = p%descendants
      n%parent           = p%parent
      n%first_child      = p%first_child
      n%next_sibling     = NODE_INVALID
      n%request_sent     = .false.
      n%interaction_data = p%interaction_data
      
    end subroutine tree_node_unpack


end module module_tree_node
