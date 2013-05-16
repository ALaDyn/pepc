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
!> Encapsulates functions for accessing, manipulating, and verifying hash
!> table data. This implements Coalesced Hashing [ch].
!>
!> [ch]: http://en.wikipedia.org/wiki/Coalesced_hashing
!>
module module_htable
    use module_pepc_types
    use module_tree_node, only: NODE_INVALID
    implicit none
    private

    type, public :: t_htable
      integer(kind_node)                 :: maxentries    !< max number of entries
      integer(kind_node)                 :: nentries      !< number of entries present
      type(t_tree_node), public, pointer :: values(:)     !< array of entry values
    end type t_htable

    public htable_create
    public htable_allocated
    public htable_maxentries
    public htable_add
    public htable_clear
    public htable_destroy

    contains


    !>
    !> allocates a hash table `t` of size `n` and initializes the 
    !> collision resolution lists
    !>
    subroutine htable_create(t, n)
      use module_debug
      implicit none

      type(t_htable), intent(inout) :: t
      integer(kind_node) :: n

      DEBUG_ASSERT(.not. htable_allocated(t))

      t%maxentries = max(n, 2_kind_node**15)
      t%nentries   = 0_kind_node

      allocate( t%values(1:t%maxentries) )

      call htable_clear(t)
    end subroutine htable_create


    !>
    !> empties the hash table `t`
    !>
    subroutine htable_clear(t)
      use module_debug
      implicit none
      type(t_htable), intent(inout) :: t

      t%nentries = 0_kind_node

      DEBUG_ASSERT(htable_allocated(t))
    end subroutine


    !>
    !> returns .true. if memory has been allocated for hash table `t`
    !>
    pure function htable_allocated(t)
      implicit none
      type(t_htable), intent(in) :: t

      logical :: htable_allocated

      htable_allocated = associated(t%values)
    end function htable_allocated


    !>
    !> returns the maximum number of entries that will fit in hash table `t`
    !>
    pure function htable_maxentries(t)
      implicit none
      type(t_htable), intent(in) :: t

      integer(kind_node) :: htable_maxentries

      htable_maxentries = t%maxentries
    end function htable_maxentries


    !>
    !> deallocates the hash table `t`
    !>
    subroutine htable_destroy(t)
      use module_debug
      implicit none
      type(t_htable), intent(inout) :: t

      DEBUG_ASSERT(htable_allocated(t))

      deallocate(t%values)
      t%maxentries = 0_kind_node
      t%nentries   = 0_kind_node
    end subroutine htable_destroy


    !>
    !> Add an entry (`k`, `v`) to the hash table `t`.
    !> If present, `entry_pointer` points to the inserted value `v`.
    !> Resolve collision if necessary
    !> `htable_add == .true.` if everything went fine
    !> `htable_add == .false.` if the key already exists in `t` and
    !> `entry_pointer` points to the current value, but the entry 
    !> itself is _not_ modified
    !>
    function htable_add(t, k, v, entry_pointer)
      use treevars
      use module_debug
      implicit none

      type(t_htable), intent(inout) :: t
      integer(kind_key), intent(in) :: k
      type(t_tree_node), intent(in) :: v
      integer(kind_node), intent(out), optional :: entry_pointer
      logical :: htable_add

      integer(kind_node) :: hashaddr

      DEBUG_ASSERT(htable_allocated(t))
      if (t%nentries >= t%maxentries) then
        DEBUG_ERROR('("Tree arrays full. # Entries: ", I0,"/",I0)', t%nentries, t%maxentries)
      end if

      t%nentries = t%nentries + 1_kind_node
      t%values(t%nentries) = v
      htable_add = .true.

      if (present(entry_pointer)) then
        entry_pointer = t%nentries
      end if
    end function htable_add

end module module_htable
