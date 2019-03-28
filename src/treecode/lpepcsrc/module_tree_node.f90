! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2017 Juelich Supercomputing Centre, 
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
    use module_pepc_kinds
    use iso_fortran_env
    implicit none
    private

    ! bits in flags to be set when children are requested, the request has been sent, and they have arrived
    integer, public, parameter :: TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE       = 0 !< bit is used in flags_local to denote that children information for the node is available in the local hashtable
    integer, public, parameter :: TREE_NODE_FLAG_LOCAL_REQUEST_SENT             = 1 !< bit is set in flags_local if request for child nodes has actually been sent
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
    public tree_node_pack
    public tree_node_unpack

    !! CRC stuff for debugging
    logical :: crc_tabled = .false.
    integer(int32) :: crc_table(0:255)

    contains

    !>
    !> Returns the first child of node `p`.
    !>
    !> If `p` has children that are locally available, returns `.true.` and `fc`
    !> points to the child.
    !> Otherwise, `.false.` is returned and `fc` points to `null()`.
    !>
    function tree_node_get_first_child(p) result(fc)
      use module_pepc_types, only: t_tree_node
      use module_atomic_ops
      implicit none

      integer(kind_node) :: fc
      type(t_tree_node), intent(in) :: p

      if (tree_node_children_available(p)) then
        call atomic_read_barrier()
        fc = p%first_child
      else
        fc = NODE_INVALID
      end if
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
    function tree_node_get_next_sibling(n) result(s)
      use module_pepc_types, only: t_tree_node
      implicit none

      integer(kind_node) :: s
      type(t_tree_node), intent(in) :: n

      s = n%next_sibling
    end function tree_node_get_next_sibling


    !>
    !> checks whether `n` is a leaf or twig node
    !>
    function tree_node_is_leaf(n)
      use module_pepc_types, only: t_tree_node
      implicit none
      type(t_tree_node), intent(in) :: n
      logical :: tree_node_is_leaf

      tree_node_is_leaf = (0 == n%descendants)
    end function tree_node_is_leaf


    !>
    !> checks whether `n` is a root node
    !>
    function tree_node_is_root(n)
      use module_pepc_types, only: t_tree_node
      implicit none
      type(t_tree_node), intent(in) :: n

      logical :: tree_node_is_root

      tree_node_is_root = (0 == n%level)
    end function tree_node_is_root


    !>
    !> checks whether the children of `n` are locally available or have
    !> to be requested from remote ranks
    !>
    function tree_node_children_available(n)
      use module_pepc_types, only: t_tree_node
      implicit none
      type(t_tree_node), intent(in) :: n
      logical :: tree_node_children_available

      tree_node_children_available = btest(n%flags_local, TREE_NODE_FLAG_LOCAL_CHILDREN_AVAILABLE)
    end function tree_node_children_available


    subroutine tree_node_pack(n, p)
      use module_pepc_types, only: t_tree_node_package, t_tree_node
      implicit none

      type(t_tree_node), intent(in) :: n
      type(t_tree_node_package), intent(out) :: p
      character(len=10) :: fst
      character(len=20) :: scnd
      character(len=168) :: cmplt
      integer(int32) :: crc1, crc2
      integer(int8) :: dummy(1:4)

      p%key              = n%key
      p%flags_global     = n%flags_global
      p%level            = n%level
      p%owner            = n%owner
      p%leaves           = n%leaves
      p%descendants      = n%descendants
      p%parent           = NODE_INVALID ! this is to be filled by answer_request()
      p%first_child      = n%first_child
      p%interaction_data = n%interaction_data

      cmplt = transfer(p, cmplt)
      fst = transfer(cmplt(1:10), fst)
      scnd = transfer(cmplt(13:32), scnd)
      call update_crc(fst, crc1)
      call update_crc(scnd, crc2)
      dummy = transfer(crc1, dummy)
      p%dummy(1) = dummy(1)
      dummy = transfer(crc2, dummy)
      p%dummy(2) = dummy(2)

    end subroutine tree_node_pack


    subroutine tree_node_unpack(p, n)
      use module_pepc_types, only: t_tree_node_package, t_tree_node
      implicit none

      type(t_tree_node_package), intent(in) :: p
      type(t_tree_node), intent(out) :: n
      character(len=10) :: fst
      character(len=20) :: scnd
      character(len=168) :: cmplt
      integer(int32) :: crc1, crc2
      integer(int8) :: dummy(1:4)

      n%key              = p%key
      n%flags_global     = p%flags_global
      n%level            = p%level
      n%flags_local      = 0_kind_byte
      n%owner            = p%owner
      n%leaves           = p%leaves
      n%descendants      = p%descendants
      n%parent           = p%parent
      n%first_child      = p%first_child
      n%next_sibling     = NODE_INVALID
      n%request_posted   = .false.
      n%interaction_data = p%interaction_data

      cmplt = transfer(p, cmplt)
      fst = transfer(cmplt(1:10), fst)
      scnd = transfer(cmplt(13:32), scnd)
      call update_crc(fst, crc1)
      call update_crc(scnd, crc2)
      dummy = transfer(crc1, dummy)
      if (dummy(1) /= p%dummy(1)) then
         write(*,*) 'CRC error1'
         write(*,*) dummy(1), p%dummy(1)
         write(*,*) p
      end if
      dummy = transfer(crc2, dummy)
      if (dummy(2) /= p%dummy(2)) then
         write(*,*) 'CRC error2'
         write(*,*) dummy(2), p%dummy(2)
         write(*,*) p
      end if

    end subroutine tree_node_unpack

    !!! adding CRC routine for debugging purposes
    subroutine update_crc(a, crc)
       integer :: n, i
       character(*) :: a
       integer(int32) :: crc

       if (.not. crc_tabled) then
          crc_tabled = .true.
          call init_table
       end if

       crc = 0
       crc = not(crc)
       n = len(a)
       do i = 1, n
          crc = ieor(shiftr(crc, 8), crc_table(iand(ieor(crc, iachar(a(i:i))), 255)))
       end do
       crc = not(crc)
    end subroutine

    subroutine init_table
       integer :: i, j
       integer(int32) :: k

       do i = 0, 255
          k = i
          do j = 1, 8
             if (btest(k, 0)) then
                k = ieor(shiftr(k, 1), -306674912)
             else
                k = shiftr(k, 1)
             end if
          end do
          crc_table(i) = k
       end do
    end subroutine

end module module_tree_node
