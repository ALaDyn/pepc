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
    use module_pepc_types, only: t_tree_node
    implicit none
    private

    integer*8, public, parameter :: HTABLE_KEY_INVALID = -1_8
    integer*8, public, parameter :: HTABLE_KEY_EMPTY   =  0_8
    integer*8, parameter :: HTABLE_FREE_LO = 1024_8 !< min address allowed for resolving collisions (from 4th level up)

    type :: t_htable_bucket
      integer*8                  :: key = HTABLE_KEY_EMPTY  !< the entry key
      integer*8                  :: link = -1_8             !< link for collision resolution
    end type t_htable_bucket
    type (t_htable_bucket), parameter :: HTABLE_EMPTY_BUCKET = t_htable_bucket(HTABLE_KEY_EMPTY, -1_8) !< constant for empty hashentry

    type, public :: t_htable
      private
      integer*8                          :: hashconst     !< bitmask for hashfunction
      integer*8                          :: maxentries    !< max number of entries
      integer*8                          :: nentries      !< number of entries present
      integer*8                          :: sum_unused    !< number of free addresses
      type(t_htable_bucket), allocatable :: buckets(:)    !< hashtable buckets
      integer*8, allocatable             :: free_addr(:)  !< list of free hash table addresses
      integer*8, allocatable             :: point_free(:) !< pointer to free address index
      type(t_tree_node), pointer         :: values(:)     !< array of entry values
    end type t_htable

    type, public :: t_htable_iterator
      private
      type(t_htable), pointer :: t => null()
      integer*8               :: i =  0_8
    end type t_htable_iterator

    public htable_create
    public htable_allocated
    public htable_entries
    public htable_maxentries
    public htable_add
    public htable_contains
    public htable_lookup
    public htable_lookup_critical
    public htable_remove_keys
    public htable_clear
    public htable_destroy
    public htable_check
    public htable_dump

    contains


    !>
    !> allocates a hash table `t` of size `n` and initializes the 
    !> collision resolution lists
    !>
    subroutine htable_create(t, n)
      implicit none

      type(t_htable), intent(inout) :: t
      integer*8 :: n

      t%maxentries = max(n, 2_8**15)
      t%nentries   = 0_8
      t%hashconst  = 2_8**int(log(1. * t%maxentries) / log(2.)) - 1_8

      allocate( &
        t%buckets(0:t%maxentries - 1), &
        t%free_addr(0:t%maxentries - 1), &
        t%point_free(0:t%maxentries - 1), &
        t%values(0:t%maxentries - 1) &
      )

      call htable_clear(t)
    end subroutine htable_create


    !>
    !> checks whether the given htable-entry is valid and not empty
    !>
    elemental function htable_entry_is_valid(e)
      implicit none
      type(t_htable_bucket), intent(in) :: e
      logical:: htable_entry_is_valid

      htable_entry_is_valid = (e%key .ne. HTABLE_KEY_EMPTY) .and. (e%key .ne. HTABLE_KEY_INVALID)
    end function


    !>
    !> empties the hash table `t`
    !>
    subroutine htable_clear(t)
      implicit none
      type(t_htable), intent(inout) :: t

      t%nentries = 0_8
      t%buckets  = HTABLE_EMPTY_BUCKET

      call htable_prepare_address_list(t)
    end subroutine


    !>
    !> returns .true. if memory has been allocated for hash table `t`
    !>
    pure function htable_allocated(t)
      implicit none
      type(t_htable), intent(in) :: t

      logical :: htable_allocated

      htable_allocated = allocated(t%buckets)
    end function htable_allocated


    !>
    !> returns the number of entries contained in `t`
    !>
    pure function htable_entries(t)
      implicit none
      type(t_htable), intent(in) :: t

      integer*8 :: htable_entries

      htable_entries = t%nentries
    end function htable_entries


    !>
    !> returns the maximum number of entries that will fit in hash table `t`
    !>
    pure function htable_maxentries(t)
      implicit none
      type(t_htable), intent(in) :: t

      integer*8 :: htable_maxentries

      htable_maxentries = t%maxentries
    end function htable_maxentries


    !>
    !> deallocates the hash table `t`
    !>
    subroutine htable_destroy(t)
      implicit none
      type(t_htable), intent(inout) :: t

      deallocate(t%buckets, t%free_addr, t%point_free, t%values)
      t%maxentries = 0_8
      t%nentries   = 0_8
      t%hashconst  = 0_8
    end subroutine htable_destroy


    !>
    !> prepares the collision resolution list by traversing the full
    !> hash table
    !>
    subroutine htable_prepare_address_list(t)
      implicit none
      type(t_htable), intent(inout) :: t
      integer*8 :: i

      ! build list of free addresses for faster collision resolution on insertion into htable
      t%sum_unused = 0_8
      do i = lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
        if (t%buckets(i)%key == HTABLE_KEY_EMPTY .and. i >= HTABLE_FREE_LO) then
          t%free_addr(t%sum_unused) = i                ! Free address list for resolving collisions
          t%point_free(i)           = t%sum_unused     ! Index
          t%sum_unused              = t%sum_unused + 1_8
        else
          t%point_free(i)           = -1_8
        end if
      end do
    end subroutine htable_prepare_address_list


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
      integer*8, intent(in) :: k
      type(t_tree_node), intent(in) :: v
      type(t_tree_node), intent(out), pointer, optional :: entry_pointer
      logical :: htable_add

      integer*8 :: hashaddr

      if (t%nentries >= t%maxentries) then
        DEBUG_ERROR('("Tree arrays full. # Entries: ", I0,"/",I0)', t%nentries, t%maxentries)
      end if

      if (.not. testaddr(t, k, hashaddr)) then
        ! this key does not exist in the htable 
        htable_add = .true.

        if (hashaddr .eq. -1_8) then
          ! the first entry is already empty
          hashaddr = iand(k, t%hashconst)

          if (t%point_free(hashaddr) /= -1_8) then ! Check if new address in collision res. list
            if (t%sum_unused <= 0_8) then
              DEBUG_ERROR(*, "Hash table collision resolution list exhausted. ", t%nentries, " / ", t%maxentries, " entries")
            end if
            t%sum_unused                            = t%sum_unused - 1_8
            t%free_addr(t%point_free(hashaddr))     = t%free_addr(t%sum_unused) ! Replace free address with last on list
            t%point_free(t%free_addr(t%sum_unused)) = t%point_free(hashaddr)    ! Reset pointer
            t%point_free(hashaddr)                  = -1_8
          end if
        else
          ! we are at the end of a linked list --> create new entry
          if (t%sum_unused <= 0_8) then
            DEBUG_ERROR(*, "Hash table collision resolution list exhausted. ", t%nentries, " / ", t%maxentries, " entries")
          end if
          t%sum_unused             = t%sum_unused - 1_8
          t%buckets(hashaddr)%link = t%free_addr(t%sum_unused)
          hashaddr                 = t%buckets( hashaddr )%link
          t%point_free(hashaddr)   = -1_8
        end if

        ! check if new entry is really empty
        if (t%buckets(hashaddr)%key /= HTABLE_KEY_EMPTY) then
          write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
          write (*,*) 'PE ',me,' key ',k,' entry',hashaddr,' unused ',t%sum_unused
          write (*,*) "desired entry:     ", v
          call debug_mpi_abort()
        end if

        t%buckets(hashaddr)%key = k
        t%nentries = t%nentries + 1_8
        t%values(hashaddr) = v
      else
        ! this key does already exists in the htable - as 'hashaddr' we return its current address
        htable_add = .false.
      end if

      if (present(entry_pointer)) then
        entry_pointer => t%values(hashaddr)
      end if

    end function htable_add


    !>
    !> returns `.true.` if hash table `t` contains a value for key `k`,
    !> `.false.` otherwise
    !>
    function htable_contains(t, k)
      implicit none

      type(t_htable), intent(in) :: t
      integer*8, intent(in) :: k

      logical :: htable_contains

      htable_contains = testaddr(t, k)

    end function htable_contains


    !>
    !> returns `.true.` if hash table `t` contains a value for key `k`
    !> and makes `v` point to the corresponding value, `.false.` and 
    !> `null()` otherwise
    !>
    function htable_lookup(t, k, v)
      implicit none

      type(t_htable), intent(in) :: t
      integer*8, intent(in) :: k
      type(t_tree_node), pointer, intent(out) :: v

      logical :: htable_lookup
      integer*8 :: addr

      htable_lookup = testaddr(t, k, addr)
      if (htable_lookup) then
        v => t%values(addr)
      else
        v => null()
      end if
    end function htable_lookup


    !>
    !> makes `v` point to the entry corresponding to `k` if there is
    !> one, otherwise halts the whole program
    !>
    !> @exception if key does not exist, the whole program is aborted
    !>
    subroutine htable_lookup_critical(t, k, v, caller)
      use treevars, only: me
      use module_debug
      implicit none

      type(t_htable), intent(in) :: t
      integer*8, intent(in) :: k
      type(t_tree_node), pointer, intent(out) :: v
      character(LEN=*), intent(in) :: caller

      integer*8 :: addr

      if (.not. testaddr(t, k, addr)) then
        ! could not find key
        DEBUG_WARNING_ALL('("Key not resolved in htable_lookup_critical at ",a)', caller)
        DEBUG_WARNING_ALL('("Bad address, check hash table and key list for PE", I7)', me)
        DEBUG_WARNING_ALL('("key                  (oct) = ", o22)', k)
        DEBUG_WARNING_ALL('("initial address      (dez) = ", i22)', iand( k, t%hashconst))
        DEBUG_WARNING_ALL('("   last address      (dez) = ", i22)', addr)
        if (.not. (addr == -1)) then
          DEBUG_WARNING_ALL('("htable(lastaddr)%key (oct) = ", o22)', t%buckets(addr)%key)
        end if
        DEBUG_WARNING_ALL('("# const              (dez) = ", i22)', t%hashconst)
        DEBUG_WARNING_ALL('("     maxentries      (dez) = ", i22)', t%maxentries)
        call htable_dump(t)
        ! TODO: this goes in module_tree, htables are not distributed!
        call debug_mpi_abort()
      end if

      v => t%values(addr)

    end subroutine htable_lookup_critical


    !>
    !> Removes the entry for `key` from `t`, i.e., `htable_contains(t, key) == .false.`
    !> afterwards
    !>
    !> \note Entries are actually only invalidated. This means that no space in
    !> the hash table is freed and `htable_entries(t)` remains unchanged.
    !>
    subroutine htable_remove_key(t, key)
      implicit none

      type(t_htable), intent(inout) :: t
      integer*8, intent(in) :: key

      call htable_remove_keys(t, (/ key /), 1)
    end subroutine htable_remove_key


    !>
    !> Removes the entry for `keys` from `t`, i.e., `htable_contains(t, k) == .false.`
    !> for any `k` in `keys` afterwards
    !>
    !> \note Entries are actually only invalidated. This means that no space in
    !> the hash table is freed and `htable_entries(t)` remains unchanged.
    !>
    subroutine htable_remove_keys(t, keys, num_keys)
      implicit none

      type(t_htable), intent(inout) :: t
      integer*8, intent(in) :: keys(num_keys)
      integer, intent(in) :: num_keys

      integer :: i
      integer*8 :: addr

      do i = 1, num_keys
        if (testaddr(t, keys(i), addr)) then
          t%buckets( addr )%key = HTABLE_KEY_INVALID
        end if
      end do

    end subroutine


    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] t hash table
    !> @param[in] keyin inverse mapping candidate.
    !> @return `.true.` if the key has been found, `.false.` otherwise
    !> @param[out] addr address if candidate exists, address of last entry in linked list otherwise, -1 if already the first lookup failed
    !>
    function testaddr(t, keyin, addr)
      use treevars
      implicit none

      type(t_htable), intent(in) :: t
      integer*8, intent(in)  :: keyin
      integer*8, intent(out), optional :: addr
      logical :: testaddr

      integer :: ires
      integer*8 :: nextaddr, lastaddr

      nextaddr = iand( keyin, t%hashconst) ! cell address hash function

      if (t%buckets( nextaddr )%key == HTABLE_KEY_EMPTY) then ! already the first entry is empty
        testaddr                = .false.  ! key does not exist in htable
        if (present(addr)) addr = -1_8     ! we return -1
        return
      end if

      ires     =  1 ! counter for number of htable lookups

      do while ( t%buckets( nextaddr )%key .ne. keyin )
        lastaddr = nextaddr
        nextaddr = t%buckets( nextaddr )%link    ! look at next linked entry
        ires     = ires + 1

        if (   (nextaddr == -1_8) & ! reached end of linked list without finding the key --> node is not in htable or htable is corrupt
          .or. (ires >= t%maxentries) ) & ! we probed at as many positions as the htable has entries --> circular linked list or htable corrupt
          then
            testaddr                = .false.  ! key does not exist in htable
            if (present(addr)) addr = lastaddr ! we return last entry in the linked list
            return
        end if

      end do

      testaddr                = .true.   ! key has been found in htable
      if (present(addr)) addr = nextaddr ! we return its address

    end function testaddr


    !>
    !> returns an iterator for traversing the entries in hash table `t` linearly
    !>
    function htable_iterator(t)
      implicit none

      type(t_htable), target, intent(in) :: t

      type(t_htable_iterator) :: htable_iterator

      htable_iterator = t_htable_iterator(t, lbound(t%buckets, dim = 1))
    end function htable_iterator


    !>
    !> Returns the next entry in the hash table associated with iterator `it` as 
    !> a pair of key and tree node in `k` and `v` and returns `.true.`.
    !> Once the iterator is exhausted, `.false.` is returned.
    !>
    !> \note If the associated hash table is modified during traversal,
    !> behaviour is undefined.
    !>
    function htable_iterator_next(it, k, v)
      implicit none

      type(t_htable_iterator), intent(inout) :: it
      integer*8, intent(out) :: k
      type(t_tree_node), pointer, intent(out) :: v

      logical htable_iterator_next

      do while (it%i <= ubound(it%t%buckets, dim = 1))
        if (htable_entry_is_valid(it%t%buckets(it%i))) then 
          htable_iterator_next = .true.
          k    =  it%t%buckets(it%i)%key
          v    => it%t%values(it%i)
          it%i =  it%i + 1_8
          return
        end if

        it%i = it%i + 1_8
      end do

      htable_iterator_next = .false.
      k =  HTABLE_KEY_EMPTY
      v => null()
    end function htable_iterator_next


    !>
    !> performs a sanity check of the internals of hash table `t`. if an error is
    !> encountered, the whole hash table is dumped to disk.
    !>
    subroutine htable_check(t, caller)
      use module_debug
      use module_tree_node
      implicit none

      type(t_htable), intent(in) :: t
      character(*), intent(in) :: caller

      integer*8 :: i, nentries_check, sum_unused_check
      logical :: error

      call pepc_status('CHECK TABLE')

      error = .false.
      nentries_check   = count(t%buckets(:)%key /= HTABLE_KEY_EMPTY)
      sum_unused_check = t%maxentries - HTABLE_FREE_LO - &
        count(t%buckets(HTABLE_FREE_LO:)%key /= HTABLE_KEY_EMPTY)
      
      do i = lbound(t%point_free, dim = 1), ubound(t%point_free, dim = 1)
        if (t%point_free(i) /= -1_8) then
          if (t%point_free(i) < lbound(t%free_addr, dim = 1) .or. &
              t%point_free(i) > ubound(t%free_addr, dim = 1)) then
            DEBUG_WARNING("(2a,/,a,i0,a,i0,a,i0)", "htable_check() called from ", caller, " point_free(", i, ") out of bounds of free_addr(", lbound(t%free_addr, dim = 1), ":", ubound(t%free_addr, dim = 1),")")
            error = .true.
          else if (t%free_addr(t%point_free(i)) /= i) then
            DEBUG_WARNING("(2a,/,a,i0,a,i0)", "htable_check() called from ", caller, " free_addr(point_free(", i, ")) is ", t%free_addr(t%point_free(i)))
            error = .true.
          end if
        end if
      end do

      if (nentries_check /= t%nentries) then
        DEBUG_WARNING("(2a,/,a,i0,a,i0)", "htable_check() called from ", caller, "nentries is ", t%nentries, " should be ", nentries_check)
        error = .true.
      end if

      if (sum_unused_check /= t%sum_unused) then
        DEBUG_WARNING("(2a,/,a,i0,a,i0)", "htable_check() called from ", caller, "sum_unused is ", t%sum_unused, " should be ", sum_unused_check)
        error = .true.
      end if

      if (error) then
        call htable_dump(t)
      end if

    end subroutine htable_check


    !>
    !> Do some quick checks on the tree structure
    !>
    ! TODO: Check nleaf, ntwig somewhere else.
    !subroutine check_table(t, callpoint)
    !  use treevars
    !  use module_debug
    !  use module_tree_node
    !  implicit none

    !  type(t_htable), intent(in) :: t
    !  character(*), intent(in) :: callpoint

    !  type(t_htable_iterator) :: it
    !  integer*8 :: key
    !  type(t_tree_node), pointer :: node
    !  integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check
    !  logical :: error

    !  call pepc_status('CHECK TABLE')

    !  error = .false.

    !  nleaf_check = 0
    !  nleaf_me_check = 0
    !  ntwig_check = 0
    !  ntwig_me_check = 0
    !  it = htable_iterator(t)
    !  do while (htable_iterator_next(it, key, node))
    !    if (tree_node_is_leaf(node)) then
    !      nleaf_check = nleaf_check + 1
    !      if (node%owner == me) then
    !        nleaf_me_check = nleaf_me_check + 1
    !      end if
    !    else
    !      ntwig_check = ntwig_check + 1
    !      if (node%owner == me) then
    !        ntwig_me_check = ntwig_me_check + 1
    !      end if
    !    end if
    !  end do

    !  if (nleaf /= nleaf_check) then
    !    DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
    !    '# leaves in table = ',nleaf_check,' vs ',nleaf,' accumulated',
    !    'Fixing and continuing for now..')
    !    !     nleaf = nleaf_check
    !    error = .true.
    !  end if

    !  if (ntwig /= ntwig_check) then
    !    DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
    !    '# twigs in table = ',ntwig_check,' vs ',ntwig,' accumulated',
    !    'Fixing and continuing for now..')
    !    !     ntwig = ntwig_check
    !    error = .true.
    !  end if

    !  if (nleaf_me /= nleaf_me_check) then
    !    DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
    !    '# own leaves in table = ',nleaf_me_check,' vs ',nleaf_me,' accumulated',
    !    'Fixing and continuing for now..')
    !    nleaf_me = nleaf_me_check
    !    error = .true.
    !  end if
    !  if (ntwig_me /= ntwig_me_check) then
    !    DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
    !    '# own twigs in table = ',ntwig_me_check,' vs ',ntwig_me,' accumulated',
    !    'Fixing and continuing for now..')

    !    ntwig_me = ntwig_me_check
    !    error = .true.
    !  end if

    !  if (error) then
    !    !call diagnose_tree(t)
    !  end if

    !end subroutine check_table


    !>
    !> Print tree structure from hash table to ipefile
    !>
    subroutine htable_dump(t, particles)
        use treevars
        use module_pepc_types
        use module_spacefilling
        use module_utils
        use module_debug, only : debug_ipefile_open, debug_ipefile_close, debug_ipefile, pepc_status
        use module_tree_node
        implicit none

        type(t_htable), intent(in) :: t
        type(t_particle), optional, intent(in) :: particles(:)

        character(1) :: collision
        integer*8 :: i

        call pepc_status('DIAGNOSE')
        call debug_ipefile_open()

        ! output hash table

        write(debug_ipefile,'(/a)') 'Hash table'

        write(debug_ipefile,'(153x,a35)') &
                    "IS_FILL_NODE              ", &
                    "|IS_BRANCH_NODE           ", &
                    "||HAS_REMOTE_CONTRIBUTIONS", &
                    "|||HAS_LOCAL_CONTRIBUTIONS", &
                    "||||REQUEST_SENT          ", &
                    "|||||CHILDREN_AVAILABLE   ", &
                    "||||||REQUEST_POSTED      "

        write(debug_ipefile,'(4(x,a10),3(x,a22),x,a14,x,a10,4x,a5,a30,/,162("-"),7("V")," 76543210")') &
                     'entry_10', &
                     'entry_8', &
                     'owner', &
                     'level', &
                     'key_8', &
                     'key_10', &
                     'parent_8', &
                     'collision link', &
                     'leaves', &
                     'flags', &
                     '||||||| childcod'

        ! write(debug_ipefile,'(154x,a)') " 3      .   2    .     1  .        "
        ! write(debug_ipefile,'(154x,a)') "10987654.32109876.54321098.76543210"

        do i = lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
          if (htable_entry_is_valid(t%buckets(i))) then
            ! flag  collisions
            if (t%buckets(i)%link /= -1_8 ) then
              collision="C"
            else
              collision=" "
            end if

            write (debug_ipefile,'(x,i10,x,o10,2(x,i10),x,o22,x,i22,x,o22,x,a1,x,i12,x,i10,4x,3(b8.8,"."),b8.8)') &
                    i, &
                    i, &
                    t%values(i)%owner, &
                    level_from_key(t%buckets(i)%key), &
                    t%buckets(i)%key, &
                    t%buckets(i)%key, &
                    parent_key_from_key(t%buckets(i)%key), &
                    collision, &
                    t%buckets(i)%link, &
                    t%values(i)%leaves, &
                    ishft(iand(t%values(i)%flags, Z'FF000000'), -24), &
                    ishft(iand(t%values(i)%flags, Z'00FF0000'), -16), &
                    ishft(iand(t%values(i)%flags, Z'0000FF00'), -08), &
                    ishft(iand(t%values(i)%flags, Z'000000FF'), -00)
          end if
        end do

        write (debug_ipefile,'(///a)') 'Tree structure'

        write(debug_ipefile,'(//a/,x,a,/,179("-"))') 'Twigs from hash-table', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

        do i = lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
          if (htable_entry_is_valid(t%buckets(i)) .and. .not. tree_node_is_leaf(t%values(i))) then
            write(debug_ipefile,*) t%values(i)%interaction_data
          end if
        end do

        write(debug_ipefile,'(//a/,x,a,/,179("-"))') 'Leaves from hash-table', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

        do i = lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
          if (htable_entry_is_valid(t%buckets(i)) .and. tree_node_is_leaf(t%values(i))) then
            write(debug_ipefile,*) t%values(i)%interaction_data
          end if
        end do

        if (present(particles)) then
          ! local particles
          write(debug_ipefile,'(//a/,x,a10,x,a,/,189("-"))') 'Local particles', 'index', 'data (see module_module_pepc_types::t_particle for meaning of the columns)'

          do i = 1_8, ubound(particles, 1)
            write(debug_ipefile,'(x,i10,x)',advance='no') i
            write(debug_ipefile,*) particles(i)
          end do
        end if

        call debug_ipefile_close()

    end subroutine htable_dump


end module module_htable
