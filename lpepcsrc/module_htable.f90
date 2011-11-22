!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates functions for accessing, manipulating, and verifying hash table data
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_htable
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Hash table datatype - 36 bytes per entry
    type :: t_hash ! TODO: exchange order of node and key
        integer   :: node          !< Address of particle/pseudoparticle data
        integer*8 :: key           !< Key
        integer   :: link          !< Pointer to next empty address in table in case of collision
        integer   :: leaves        !< # leaves contained within twig (=1 for leaf, npart for root)
        integer   :: childcode     !< Byte code indicating position of children (twig node); particle label (leaf node)
        integer   :: owner         !< Node owner (for branches)
    end type t_hash

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type (t_hash), public, target, allocatable :: htable(:) !< hash table
    integer*8,     public ::  hashconst  !< hashing constants
    integer*8,     public ::  hashchild = b'111' !< bits that contain the child index in a key

    ! bits in childcode to be set when children are requested, the request has been sent, and they have arrived
    integer, public, parameter :: CHILDCODE_BIT_REQUEST_POSTED     =  8 !< this bit is used inside the childcode to denote that a request for children information is already in the request queue
    integer, public, parameter :: CHILDCODE_BIT_CHILDREN_AVAILABLE =  9 !< this bit is used inside the childcode to denote that children information for the node is available in the local hashtable
    integer, public, parameter :: CHILDCODE_BIT_REQUEST_SENT       = 10 !< this bit is used inside the childcode to denote that children information has already been requested from the owner
    integer, public, parameter :: CHILDCODE_CHILDBYTE            = b'11111111' !< bits that contain the children information for this node

    integer, public :: maxaddress                    !< max address allowed in #table

    ! TODO: make the following private
    integer, public, allocatable :: free_addr(:)    !< List of free #table addresses (for HASHENTRY routine)
    integer, public, allocatable :: point_free(:)   !< Pointer to free address index
    integer, parameter :: free_lo = 1024             !< min address allowed for resolving collisions (from 4th level up)
    integer :: iused                                  !< counter for collision resolution array free_addr()
    integer :: sum_unused                             !< # free addresses


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public children_available
    public get_next_node
    public get_childkeys
    public make_hashentry
    public key2addr
    public testaddr
    public htable_clear
    public htable_clear_and_insert_root
    public htable_prepare_address_list
    public check_table
    public diagnose_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, private, parameter :: start_child_idx = 0 !< index of first child to be used in traversal - do not change, currently not completely implemented
    type (t_hash), private, parameter :: HASHENTRY_EMPTY = t_hash(0,0_8,-1,0,0,0) !< constant for empty hashentry

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks whether the given htable-entry is a leaf or a twig
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TODO: use this function everywhere consistently
    function is_leaf(hashaddr)
      implicit none
      integer, intent(in) :: hashaddr
      logical:: is_leaf
      ! TODO: this way of identifying leaves/twigs is not good -- use number of leaves or (even better) a flag in the childcode, instead
      is_leaf = (htable(hashaddr)%node > 0)
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> empties the htable
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_clear()
        implicit none

        htable = HASHENTRY_EMPTY ! TODO: need list of 'live' adresses to speed this up

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> empties the htable, inserts the root node and prepares the
    !> collision resolution list
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_clear_and_insert_root()
        use treevars, only : ntwig, me
        implicit none

        call htable_clear()

        ntwig = 1

        htable(1) = t_hash(-1, 1_8, -1, 0, IBSET(0, CHILDCODE_BIT_CHILDREN_AVAILABLE), me)

        call htable_prepare_address_list()

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> prepares the collision resolution list by traversing the full htable
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_prepare_address_list()
        implicit none
        integer :: i

        ! build list of free addresses for faster collision resolution on insertion into htable ! TODO: move free_addr and related fields to module_htable and add appropriate access routines
        sum_unused = 0
        iused      = 1   ! reset used-address counter
        do i=0, maxaddress
           if (htable(i)%node == 0 .and. htable(i)%key /=-1 .and. i > free_lo) then
              sum_unused            = sum_unused + 1
              free_addr(sum_unused) = i            ! Free address list for resolving collisions
              point_free(i)         = sum_unused   ! Index
           else
              point_free(i)         = 0
           endif
        enddo

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks whether children for certain htable-address are locally available
    !> or have to be requested
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function children_available(addr)
      implicit none
      logical :: children_available
      integer, intent(in) :: addr

      children_available = btest(htable( addr )%childcode, CHILDCODE_BIT_CHILDREN_AVAILABLE)

    end function


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Function to return key of next node in local tree-walk,
    !> i.e. search for next sibling, uncle, great-uncle etc
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_next_node(keyin)

        use treevars

        implicit none
        integer*8 :: get_next_node
        integer*8, intent(in) :: keyin

        integer*8 :: search_key, parent_key
        integer :: parent_addr
        integer :: parent_child_byte, search_child_idx

        search_key = keyin

        ! search for next sibling, uncle, great-uncle etc
        do while (search_key > 1) ! loop through parent nodes up to root
            parent_key        = ishft(search_key,-3)
            parent_addr       = key2addr( parent_key ,"next_node(), get parent_addr" )
            parent_child_byte = ibits( htable( parent_addr ) % childcode, 0, 8)

            search_child_idx  = int(ibits( search_key, 0, 3), kind(search_child_idx) ) ! lower three bits of key

            do ! loop over all siblings
                search_child_idx   = modulo(search_child_idx + 1, 8) ! get next sibling, possibly starting again from first one

                ! if sibling-loop wrapped and reached starting point again --> go up one level
                if ( search_child_idx == start_child_idx ) then
                    search_key = parent_key      ! go up one level
                    exit
                endif

                ! if sibling exists: next_node has been found
                if ( btest(parent_child_byte, search_child_idx) ) then
                    get_next_node = ior(int(ishft(parent_key, 3), kind(search_child_idx)), search_child_idx) ! assemble next_node out of parent-key and new sibling-index
                    return
                endif
            end do
        end do

        get_next_node  = 1 ! nothing has been found, i.e. top-right corner reached: set pointer=root

    end function get_next_node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> returns the keys of all children, that are attached to the
    !> node at a particular htable address
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_childkeys(addr, childnum, childkeys)
        use treevars
        implicit none
        integer, intent(in) :: addr
        integer, intent(out) :: childnum
        integer*8, dimension(:), intent(out) :: childkeys
        integer :: i

        integer*8 :: keyhead
        integer :: childcode

        keyhead   = ishft(htable(addr)%key, 3)
        childcode = htable(addr)%childcode
        childnum = 0

        do i=0,7
          if (btest(childcode, i)) then
            childnum            = childnum + 1
            childkeys(childnum) = ior(keyhead, 1_8*i)
          end if
        end do

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Make entry in hash-table - returns address 'newentry'
    !>  Resolve collision if necessary
    !>  ierror == 0 if anything went fine
    !>  ierror == 1 if key already exists in htable, newentry is set to
    !>               return current address, but the htable-entry itself is *not* modified
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> TODO: maybe use t_hash as parameter type instead of different parameters?
    subroutine make_hashentry( keyin, nodein, leavesin, codein, ownerin, newentry, ierror)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in) :: keyin
        integer, intent(in) :: nodein, leavesin, codein, ownerin      ! call input parameters
        integer, intent(out) :: newentry  ! address in # table returned to calling routine
        integer, intent(out) :: ierror

        integer :: ierr

        if (.not. testaddr(keyin, newentry)) then
          ! this key does not exist in the htable 
          ierror = 0

          if (newentry .eq. -1) then
            ! the first entry is already empty
            newentry = int(IAND( keyin, hashconst))

            if (point_free(newentry) /= 0) then     ! Check if new address in collision res. list
                free_addr( point_free(newentry) ) = free_addr(sum_unused)  ! Replace free address with last on list
                point_free(free_addr(sum_unused)) = point_free(newentry)  ! Reset pointer
                point_free(newentry)              = 0
                sum_unused = sum_unused - 1
            endif
          else
            ! we are at the end of a linked list --> create new entry
            htable( newentry )%link = free_addr(iused)
            newentry                = htable( newentry )%link
            iused                   = iused + 1
          endif

          ! check if new entry is really empty
          if ((htable(newentry)%node /= 0 ) .or. (htable( newentry )%key/=0)) then
            write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
            write (*,*) 'PE ',me,' key ',keyin,' entry',newentry,' used ',iused,'/',sum_unused
            write (*,*) "htable(newentry):  ", htable(newentry)
            write (*,*) "desired entry:     ", t_hash( nodein, keyin, -1, leavesin, codein, ownerin )
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          endif

          htable( newentry ) = t_hash( nodein, keyin, htable( newentry )%link, leavesin, codein, ownerin )
        else
          ! this key does already exists in the htable - as 'newentry' we return its current address
          ierror = 1
        endif

    end subroutine make_hashentry



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] keyin inverse mapping candidate.
    !> @param[in] cmark a description.
    !> @param[out] key2addr the adress if the key exists
    !> @exception if key does not exist, the whole program is aborted
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function key2addr(keyin,cmark)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in)  :: keyin
        integer :: ierr
        character(LEN=*) :: cmark
        integer :: key2addr

        if (.not. testaddr(keyin, key2addr)) then
          ! could not find key
          write(*,'("Key not resolved in KEY2ADDR at ",a)') cmark
          write(*,'("Bad address, check #-table and key list for PE", I7)') me
          write(*,'("key (octal)           = ", o22)') keyin
          write(*,'("initial address (dez) = ", i22)') int(IAND( keyin, hashconst))
          write(*,'("   last address (dez) = ", i22)') key2addr
          write(*,'("# const         (dez) = ", i22)') hashconst
          write(*,'("maxaddress      (dez) = ", i22)') maxaddress
          call diagnose_tree()
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

    end function key2addr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] keyin inverse mapping candidate.
    !> @return .true. if the key has been found, .false. otherwise
    !> @param[out] addr address if candidate exists, address of last entry in linked list otherwise, -1 if already the first lookup failed
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function testaddr(keyin, addr)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in)  :: keyin
        integer :: ires
        logical :: testaddr
        integer, intent(out), optional :: addr
        integer :: nextaddr, lastaddr

        nextaddr = int(IAND( keyin, hashconst))     ! cell address hash function

        if (htable( nextaddr )%key == 0) then ! already the first entry is empty
           testaddr                = .false.  ! key does not exist in htable
           if (present(addr)) addr = -1       ! we return -1
           return
        endif

        ires     =  1 ! counter for number of htable lookups

        do while ( htable( nextaddr )%key .ne. keyin )
          lastaddr = nextaddr
          nextaddr = htable( nextaddr )%link    ! look at next linked entry
          ires     = ires + 1

          if (   (nextaddr == -1) & ! reached end of linked list without finding the key --> node is not in htable or htable is corrupt
            .or. (ires >= maxaddress) ) & ! we probed at as many positions as the htable has entries --> circular linked list or htable corrupt
            then
              testaddr                = .false.  ! key does not exist in htable
              if (present(addr)) addr = lastaddr ! we return last entry in the linked list
              return
          endif

        end do

        testaddr                = .true.   ! key has been found in htable
        if (present(addr)) addr = nextaddr ! we return its address

    end function testaddr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Do some quick checks on the tree structure
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_table(callpoint)
        use treevars
        implicit none
        character(*) :: callpoint
        integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check

        ntwig_check = count(mask =  htable%node <0 )
        nleaf_check = count(mask =  htable%node >0 )
        nleaf_me_check = count(mask = htable%owner==me .and. htable%node >0 )
        ntwig_me_check = count(mask = htable%owner==me .and. htable%node <0 )

        if (nleaf /= nleaf_check) then
            write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
            write(*,*) '# leaves in table = ',nleaf_check,'vs ',nleaf,'accumulated'
            write(*,*) 'Fixing and continuing for now..'
        !     nleaf = nleaf_check
        endif

        if (ntwig /= ntwig_check) then
            write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
            write(*,*) ' # twigs in table = ',ntwig_check,'vs ',ntwig,'accumulated'
            write(*,*) 'Fixing and continuing for now..'
        !     ntwig = ntwig_check
        endif

        if (nleaf_me /= nleaf_me_check) then
            write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
            write(*,*) ' # own leaves in table = ',nleaf_me_check,'vs ',nleaf_me,'accumulated'
            write(*,*) 'Fixing and continuing for now..'
            nleaf_me = nleaf_me_check
        endif
        if (ntwig_me /= ntwig_me_check) then
            write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
            write(*,*) ' # own twigs in table = ',ntwig_me_check,'vs ',ntwig_me,'accumulated'
            write(*,*) 'Fixing and continuing for now..'
            ntwig_me = ntwig_me_check
        endif

    end subroutine check_table



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Check integrity of tree structure from hash table
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine diagnose_tree
        use treevars
        use module_spacefilling
        implicit none

        integer*8 :: key_twig(ntwig), key_leaf(nleaf)
        integer, dimension(ntwig) :: child_twig, addr_twig, ind_twig      ! twig-nodes
        real*8 :: rcoc2(1:ntwig)
        integer, dimension(nleaf) :: plist_leaf, ind_leaf, owner_leaf       ! leaf-nodes

        character(1) :: csnap, collision

        integer :: i, isnap

        save isnap
        data isnap/1/

        csnap=achar(mod(isnap,10)+48)

        call status('DIAGNOSE')

        ! output hash table

        write (ipefile,'(/a/8x,a,a)') 'Hash table ', &
        'entry,    owner    node,            key_8     key_10        parent       next   ', &
        '    link   #leaves  childcode  collision', &
        '----------------------------------------------------------------------------------------------- '

        ! flag  collisions

        do i=0,maxaddress
            collision=" "
            if (htable(i)%node/=0 .and. htable(i)%link/= -1 ) collision="C"

            if (htable(i)%node /= 0) write (ipefile,'(3i10,o22,i10,o22,i8,i10,z4,4x,a1)') &
            i,htable(i)%owner,htable(i)%node,htable(i)%key,htable(i)%key,ishft( htable(i)%key,-3 ), &
            htable(i)%link,htable(i)%leaves,htable(i)%childcode,collision
        end do


        ! get keys of twig nodes from hash table
        key_twig(1:ntwig)  = pack(htable(0:maxaddress)%key,mask=htable(0:maxaddress)%node<0)
        ! get levels of twigs
        addr_twig(1:ntwig) = (/( key2addr( key_twig(i),'DIAGNOSE_TREE' ),i=1,ntwig)/)   !  Table address
        child_twig(1:ntwig) = (/( htable( key2addr( key_twig(i),'DIAGNOSE_TREE' ) )%childcode,i=1,ntwig )/)   !  Children byte-code
        ind_twig(1:ntwig) = (/( htable( key2addr( key_twig(i),'DIAGNOSE_TREE' ) )%node,i=1,ntwig )/)   !  Twig node pointers

        do i=1,ntwig
            rcoc2(i) = dot_product(tree_nodes(ind_twig(i))%coc, tree_nodes(ind_twig(i))%coc)
        end do

        write (ipefile,'(///a)') 'Tree structure'

        !  write (ipefile,'(/a/a/(3i5,2i10,2i8,b11,i2,i8,i10,9(1pe15.4)))') 'Twigs from hash-table:', &
        write (ipefile,'(/a/a,a/(3i5,2o15,2i8,o15,i8,14(1pe30.19)))') 'Twigs from hash-table:', &
        '    i  level  owner        key     parent-key       #    node  code      1st child #leaves ', &
        ' abs_charge    charge   xcoc   ycoc   zcoc   xdip   ydip   zdip   xxquad   yyquad   zzquad   xyquad   yzquad   zxquad', &
        (i,level_from_key(key_twig(i)), &              !  index, level
        htable( addr_twig(i) )%owner, &                            ! Owner-PE of node
        key_twig(i),ishft( key_twig(i),-3 ), &                             ! key, parent key
        addr_twig(i), ind_twig(i), &    ! Table address and node number
        child_twig(i), &                         ! Children byte-code
        htable( addr_twig(i) )%leaves, &                           ! # leaves contained in branch
        tree_nodes(ind_twig(i))%abs_charge, &    ! Twig absolute charge
        tree_nodes(ind_twig(i))%charge, &    ! Twig  charge
        tree_nodes(ind_twig(i))%coc(1), & ! Centre of charge
        tree_nodes(ind_twig(i))%coc(2), &
        tree_nodes(ind_twig(i))%coc(3), &
        tree_nodes(ind_twig(i))%dip(1), &
        tree_nodes(ind_twig(i))%dip(2), &
        tree_nodes(ind_twig(i))%dip(3), &
        tree_nodes(ind_twig(i))%quad(1), &
        tree_nodes(ind_twig(i))%quad(2), &
        tree_nodes(ind_twig(i))%quad(3), &
        tree_nodes(ind_twig(i))%xyquad, &
        tree_nodes(ind_twig(i))%yzquad, &
        tree_nodes(ind_twig(i))%zxquad, &
        i=1,ntwig)


        ! get keys of local leaf nodes from hash table
        key_leaf(1:nleaf_me) = pack(htable%key,mask=(htable%node>0 .and. htable%owner == me))
        ind_leaf(1:nleaf_me) = pack(htable%node,mask=(htable%node>0 .and. htable%owner == me))         ! particle/leaf index
        plist_leaf(1:nleaf_me) = pack(htable%childcode,mask=(htable%node>0 .and. htable%owner == me))   ! particle label
        owner_leaf(1:nleaf_me) = pack(htable%owner,mask=(htable%node>0 .and. htable%owner == me))   ! who owns leaf node


        write (ipefile,'(/a/3a5,2a10,2a15,a25,4a11/(3i5,2i10,2o15,o25,4f30.19))') 'Local leaves from hash-table:', &
        'i','owner','plab','i-leaf','lev','key','parent','pkey','x','y','z','q', &
        (i,owner_leaf(i),plist_leaf(i),ind_leaf(i),level_from_key(key_leaf(i)),key_leaf(i), &
        ishft( key_leaf(i),-3 ), &      ! parent
        particles(ind_leaf(i))%key, &  ! particle key
        particles(ind_leaf(i))%x, particles(ind_leaf(i))%data%q, &
        i=1,nleaf_me)

        ! get keys of NON-local leaf nodes from hash table
        key_leaf(1:nleaf-nleaf_me) = pack(htable%key,mask=(htable%node>0 .and. htable%owner /= me))
        ind_leaf(1:nleaf-nleaf_me) = pack(htable%node,mask=(htable%node>0 .and. htable%owner /= me))         ! leaf index
        plist_leaf(1:nleaf-nleaf_me) = pack(htable%childcode,mask=(htable%node>0 .and. htable%owner /= me))   ! global particle label
        owner_leaf(1:nleaf-nleaf_me) = pack(htable%owner,mask=(htable%node>0 .and. htable%owner /= me))   ! who owns leaf node


        write (ipefile,'(//a/a/(4i5,2o15,i5,2f11.4,f6.1,f11.4))') 'Non-local leaves from hash-table:', &
        '    i   owner    i-leaf    lev    key    parent  plabel  xcoc  ycoc  charge      ', &
        (i,owner_leaf(i),ind_leaf(i),level_from_key(key_leaf(i)),key_leaf(i), &
          ishft( key_leaf(i),-3 ), &      ! parent
          plist_leaf(i), & ! global particle label
          tree_nodes(ind_twig(i))%coc(1),&
          tree_nodes(ind_twig(i))%coc(2),&
          tree_nodes(ind_twig(i))%charge,&
          tree_nodes(ind_twig(i))%dip(1), &
        i=1,nleaf-nleaf_me)

    end subroutine diagnose_tree


end module module_htable
