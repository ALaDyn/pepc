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

    integer*8, public, parameter :: HTABLE_KEY_INVALID = -1
    integer*8, public, parameter :: HTABLE_KEY_EMPTY   =  0
    integer, parameter :: HTABLE_FREE_LO = 1024 !< min address allowed for resolving collisions (from 4th level up)

    type :: t_htable_bucket
      integer*8                  :: key = HTABLE_KEY_EMPTY  !< the entry key
      type(t_tree_node), pointer :: value_pointer => null() !< index into the values array
      integer                    :: link = -1               !< link for collision resolution
    end type t_htable_bucket
    type (t_htable_bucket), parameter :: HTABLE_EMPTY_BUCKET = t_htable_bucket(HTABLE_KEY_EMPTY, null(), -1) !< constant for empty hashentry

    type, public :: t_htable
      private
      integer*8                          :: hashconst     !< bitmask for hashfunction
      integer                            :: maxvalues     !< max number of entries
      integer                            :: nvalues       !< number of entries present
      integer                            :: iused         !< counter for collision resolution array free_addr()
      integer                            :: sum_unused    !< number of free addresses
      type(t_htable_bucket), allocatable :: buckets(:)    !< hashtable buckets
      integer, allocatable               :: free_addr(:)  !< list of free hash table addresses
      integer, allocatable               :: point_free(:) !< pointer to free address index
      type(t_tree_node), pointer         :: values(:)     !< array of entry values
    end type t_htable

    type(t_htable), save, public :: global_htable !< global singleton hash table for gradual transition

    public htable_init
    public htable_allocated
    public htable_add
    public htable_contains
    public htable_lookup
    public htable_lookup_critical
    public htable_remove_keys
    public htable_clear
    public htable_destroy
    public check_table
    public diagnose_tree

    contains


    !>
    !> allocates a hash table `t` of size `n` and initializes the 
    !> collision resolution lists
    !>
    subroutine htable_init(t, n)
      implicit none

      type(t_htable), intent(inout) :: t
      integer :: n

      t%maxvalues = n
      t%nvalues   = 0
      t%hashconst = 2**int(max(log(1. * n) / log(2.), 15.)) - 1

      allocate(t%buckets(0:n), t%free_addr(0:n), t%point_free(0:n), t%values(n+1))

      call htable_clear(t)
    end subroutine htable_init


    !>
    !> checks whether the given htable-entry is valid and not empty
    !>
    ! TODO: use this function everywhere consistently
    function htable_entry_is_valid(e)
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

      t%nvalues = 0
      t%buckets = HTABLE_EMPTY_BUCKET ! TODO: need list of 'live' adresses to speed this up
                                      ! possible solution: use a "bitmap" of occupied addresses in htable
                                      ! then, only this bitmap has to be cleared upon startup
                                      ! and every test for occupancy of a htable entry is done in this bitmap

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
    !> returns the maximum number of entries that will fit in hash table `t`
    !>
    pure function htable_maxentries(t)
      implicit none
      type(t_htable), intent(in) :: t

      integer :: htable_maxentries

      htable_maxentries = t%maxvalues
    end function htable_maxentries


    !>
    !> deallocates the hash table `t`
    !>
    subroutine htable_destroy(t)
      implicit none
      type(t_htable), intent(inout) :: t

      deallocate(t%buckets, t%free_addr, t%point_free, t%values)
      t%maxvalues = 0
      t%nvalues   = 0
      t%hashconst = 0
    end subroutine htable_destroy


    !>
    !> prepares the collision resolution list by traversing the full
    !> hash table
    !>
    subroutine htable_prepare_address_list(t)
        implicit none
        type(t_htable), intent(inout) :: t
        integer :: i

        ! build list of free addresses for faster collision resolution on insertion into htable
        t%sum_unused = 0
        t%iused      = 1   ! reset used-address counter
        do i=0, t%maxvalues
           if ((.not. associated(t%buckets(i)%value_pointer)) .and. t%buckets(i)%key /=-1 .and. i > HTABLE_FREE_LO) then
              t%sum_unused              = t%sum_unused + 1
              t%free_addr(t%sum_unused) = i                ! Free address list for resolving collisions
              t%point_free(i)           = t%sum_unused     ! Index
           else
              t%point_free(i)           = 0
           endif
        enddo

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

        integer :: hashaddr

        if (t%nvalues >= t%maxvalues) then
          DEBUG_ERROR('("Tree arrays full. # Entries: ", I0,"/",I0)', t%nvalues, t%maxvalues)
        end if

        if (.not. testaddr(t, k, hashaddr)) then
          ! this key does not exist in the htable 
          htable_add = .true.

          if (hashaddr .eq. -1) then
            ! the first entry is already empty
            hashaddr = int(IAND(k, t%hashconst))

            if (t%point_free(hashaddr) /= 0) then     ! Check if new address in collision res. list
                t%free_addr( t%point_free(hashaddr) ) = t%free_addr(t%sum_unused)  ! Replace free address with last on list
                t%point_free(t%free_addr(t%sum_unused)) = t%point_free(hashaddr)   ! Reset pointer
                t%point_free(hashaddr)              = 0
                t%sum_unused = t%sum_unused - 1
            endif
          else
            ! we are at the end of a linked list --> create new entry
            t%buckets( hashaddr )%link = t%free_addr(t%iused)
            hashaddr                   = t%buckets( hashaddr )%link
            t%iused                    = t%iused + 1
          endif

          ! check if new entry is really empty
          if (associated(t%buckets(hashaddr)%value_pointer) .or. (t%buckets( hashaddr )%key/=0)) then
            write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
            write (*,*) 'PE ',me,' key ',k,' entry',hashaddr,' used ',t%iused,'/',t%sum_unused
            !write (*,*) "htable(hashaddr):  ", t%buckets(hashaddr)
            write (*,*) "desired entry:     ", v
            call debug_mpi_abort()
          endif

          t%buckets(hashaddr)%key = k
          t%nvalues = t%nvalues + 1
          t%values(t%nvalues) = v
          t%buckets(hashaddr)%value_pointer => t%values(t%nvalues)
        else
          ! this key does already exists in the htable - as 'hashaddr' we return its current address
          htable_add = .false.
        endif

        if (present(entry_pointer)) then
          entry_pointer => t%buckets(hashaddr)%value_pointer
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
      integer :: dummy

      htable_contains = testaddr(t, k, dummy)

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
      integer :: addr

      htable_lookup = testaddr(t, k, addr)
      if (htable_lookup) then
        v => t%buckets(addr)%value_pointer
      else
        v => null()
      end if
    end function htable_lookup


    !>
    !> makes `v` point to the entry corresponding to `k` if there is
    !> one, otherwise halts the whole program
    !>
    subroutine htable_lookup_critical(t, k, v, caller)
      implicit none

      type(t_htable), intent(in) :: t
      integer*8, intent(in) :: k
      type(t_tree_node), pointer, intent(out) :: v
      character(LEN=*), intent(in) :: caller

      integer :: addr

      addr = key2addr(t, k, caller)
      v => t%buckets(addr)%value_pointer
    end subroutine htable_lookup_critical


    !>
    !> Invalidates entries in the hash table `t`, i.e. sets their key to KEY_INVALID
    !> TODO: this function does not free the nodelist-storage and does
    !> not fix any connections inside the hash table, esp. concerning the
    !> children-available-flag (it even does not care for them)
    !> additionally, it does not modify `nleaf` or `ntwig` which would be
    !> necessary to survive `check_table()` if the entry really was removed
    !>
    subroutine htable_remove_keys(t, keys, num_keys)
      implicit none

      type(t_htable), intent(inout) :: t
      integer*8, intent(in) :: keys(num_keys)
      integer, intent(in) :: num_keys

      integer :: i

      do i=1,num_keys
        t%buckets(  key2addr(t, keys(i), 'htable_remove_keys')  )%key = HTABLE_KEY_INVALID
      end do

    end subroutine

    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] t hash table
    !> @param[in] keyin inverse mapping candidate.
    !> @param[in] cmark a description.
    !> @param[out] key2addr the adress if the key exists
    !> @exception if key does not exist, the whole program is aborted
    !>
    function key2addr(t, keyin,cmark)

        use treevars
        use module_debug
        implicit none

        type(t_htable), intent(in) :: t
        integer*8, intent(in)  :: keyin
        character(LEN=*) :: cmark
        integer :: key2addr

        if (.not. testaddr(t, keyin, key2addr)) then
          ! could not find key
          DEBUG_WARNING_ALL('("Key not resolved in KEY2ADDR at ",a)', cmark)
          DEBUG_WARNING_ALL('("Bad address, check #-table and key list for PE", I7)', me)
          DEBUG_WARNING_ALL('("key                  (oct) = ", o22)', keyin)
          DEBUG_WARNING_ALL('("initial address      (dez) = ", i22)', int(IAND( keyin, t%hashconst)))
          DEBUG_WARNING_ALL('("   last address      (dez) = ", i22)', key2addr)
          if (.not. (key2addr == -1)) then
            DEBUG_WARNING_ALL('("htable(lastaddr)%key (oct) = ", o22)', t%buckets(key2addr)%key)
          end if
          DEBUG_WARNING_ALL('("# const              (dez) = ", i22)', t%hashconst)
          DEBUG_WARNING_ALL('("     maxvalues       (dez) = ", i22)', t%maxvalues)
          call diagnose_tree(t)
          call debug_mpi_abort()
        endif

    end function key2addr


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
      integer, intent(out), optional :: addr
      logical :: testaddr

      integer :: ires
      integer :: nextaddr, lastaddr

      nextaddr = int(IAND( keyin, t%hashconst))     ! cell address hash function

      if (t%buckets( nextaddr )%key == HTABLE_KEY_EMPTY) then ! already the first entry is empty
        testaddr                = .false.  ! key does not exist in htable
        if (present(addr)) addr = -1       ! we return -1
        return
      endif

      ires     =  1 ! counter for number of htable lookups

      do while ( t%buckets( nextaddr )%key .ne. keyin )
        lastaddr = nextaddr
        nextaddr = t%buckets( nextaddr )%link    ! look at next linked entry
        ires     = ires + 1

        if (   (nextaddr == -1) & ! reached end of linked list without finding the key --> node is not in htable or htable is corrupt
          .or. (ires >= t%maxvalues) ) & ! we probed at as many positions as the htable has entries --> circular linked list or htable corrupt
          then
            testaddr                = .false.  ! key does not exist in htable
            if (present(addr)) addr = lastaddr ! we return last entry in the linked list
            return
        endif

      end do

      testaddr                = .true.   ! key has been found in htable
      if (present(addr)) addr = nextaddr ! we return its address

    end function testaddr


    !>
    !> Do some quick checks on the tree structure
    !>
    subroutine check_table(t, callpoint)
        use treevars
        use module_debug
        use module_tree_node
        implicit none

        type(t_htable), intent(in) :: t
        character(*), intent(in) :: callpoint

        integer :: i, nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check
        logical :: error

        call pepc_status('CHECK TABLE')

        error = .false.

        nleaf_check = 0
        nleaf_me_check = 0
        ntwig_check = 0
        ntwig_me_check = 0
        do i = lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1) 
          if (htable_entry_is_valid(t%buckets(i))) then
            if (tree_node_is_leaf(t%buckets(i)%value_pointer)) then
              nleaf_check = nleaf_check + 1
              if (t%buckets(i)%value_pointer%owner == me) then
                nleaf_me_check = nleaf_me_check + 1
              end if
            else
              ntwig_check = ntwig_check + 1
              if (t%buckets(i)%value_pointer%owner == me) then
                ntwig_me_check = ntwig_me_check + 1
              end if
            end if
          end if
        end do

        if (nleaf /= nleaf_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# leaves in table = ',nleaf_check,' vs ',nleaf,' accumulated',
                                                       'Fixing and continuing for now..')
        !     nleaf = nleaf_check
            error = .true.
        endif

        if (ntwig /= ntwig_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# twigs in table = ',ntwig_check,' vs ',ntwig,' accumulated',
                                                       'Fixing and continuing for now..')
        !     ntwig = ntwig_check
            error = .true.
        endif

        if (nleaf_me /= nleaf_me_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# own leaves in table = ',nleaf_me_check,' vs ',nleaf_me,' accumulated',
                                                       'Fixing and continuing for now..')
            nleaf_me = nleaf_me_check
            error = .true.
        endif
        if (ntwig_me /= ntwig_me_check) then
            DEBUG_WARNING('(3a,i0,/,a,i0,a,i0,a,/,a)', 'Table check called ',callpoint,' by PE',me,
                                                       '# own twigs in table = ',ntwig_me_check,' vs ',ntwig_me,' accumulated',
                                                       'Fixing and continuing for now..')

            ntwig_me = ntwig_me_check
            error = .true.
        endif

        if (error) then
          call diagnose_tree(t)
        endif

    end subroutine check_table


    !>
    !> Print tree structure from hash table to ipefile
    !>
    subroutine diagnose_tree(t, particles)
        use treevars
        use module_pepc_types
        use module_spacefilling
        use module_utils
        use module_debug, only : debug_ipefile_open, debug_ipefile_close, debug_ipefile, pepc_status
        implicit none

        type(t_htable), intent(in) :: t
        type(t_particle), optional, intent(in) :: particles(1:npp)

        integer*8, dimension(ntwig) :: node_twig      ! twig-nodes
        integer*8, dimension(nleaf) :: node_leaf      ! leaf-nodes
        character(1) :: collision
        integer :: i

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

        do i=lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
            if (htable_entry_is_valid(t%buckets(i))) then
              ! flag  collisions
              if (t%buckets(i)%link/= -1 ) then
                collision="C"
              else
                collision=" "
              endif

              write (debug_ipefile,'(x,i10,x,o10,2(x,i10),x,o22,x,i22,x,o22,x,a1,x,i12,x,i10,4x,3(b8.8,"."),b8.8)') &
                      i,                   &
                      i,                   &
                      t%buckets(i)%value_pointer%owner,     &
                      level_from_key(t%buckets(i)%key), &
                      t%buckets(i)%key,       &
                      t%buckets(i)%key,       &
                      parent_key_from_key(t%buckets(i)%key), &
                      collision,           &
                      t%buckets(i)%link,      &
                      t%buckets(i)%value_pointer%leaves,    &
                      ishft(iand(t%buckets(i)%value_pointer%info_field, Z'FF000000'), -24), &
                      ishft(iand(t%buckets(i)%value_pointer%info_field, Z'00FF0000'), -16), &
                      ishft(iand(t%buckets(i)%value_pointer%info_field, Z'0000FF00'), -08), &
                      ishft(iand(t%buckets(i)%value_pointer%info_field, Z'000000FF'), -00)
           end if
        end do

        write (debug_ipefile,'(///a)') 'Tree structure'

        write(debug_ipefile,'(//a/,x,a,/,179("-"))') 'Twigs from hash-table', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

        do i=lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
          if (htable_entry_is_valid(t%buckets(i))) then
            write(debug_ipefile,*) t%buckets(i)%value_pointer%interaction_data
          end if
        end do

        write(debug_ipefile,'(//a/,x,a,/,179("-"))') 'Leaves from hash-table', 'data (see module_interaction_specific::t_tree_node_interaction_data for meaning of the columns)'

        do i=lbound(t%buckets, dim = 1), ubound(t%buckets, dim = 1)
          if (htable_entry_is_valid(t%buckets(i))) then
            write(debug_ipefile,*) t%buckets(i)%value_pointer%interaction_data
          end if
        end do

        if (present(particles)) then
          ! local particles
          write(debug_ipefile,'(//a/,x,a10,x,a,/,189("-"))') 'Local particles', 'index', 'data (see module_module_pepc_types::t_particle for meaning of the columns)'

          do i=1,npp
            write(debug_ipefile,'(x,i10,x)',advance='no') i
            write(debug_ipefile,*) particles(i)
          end do
        endif

        call debug_ipefile_close()

    end subroutine diagnose_tree


end module module_htable
