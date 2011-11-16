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
    type :: t_hash
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
    integer, public, parameter :: CHILDCODE_NODE_TOUCHED           = 11 !< this bit is used inside the childcode to notify of nodes, that already contain valid multipole information and may not be set to zero in tree_global
    integer, public, parameter :: CHILDCODE_CHILDBYTE            = b'11111111' !< bits that contain the children information for this node

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
    !> empties the htable
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine htable_clear()
        implicit none

        htable = HASHENTRY_EMPTY ! TODO: need list of 'live' adresses to speed this up

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
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine make_hashentry( keyin, nodein, leavesin, codein, ownerin, newentry, ierror)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in) :: keyin
        integer, intent(in) :: nodein, leavesin, codein, ownerin      ! call input parameters
        integer, intent(out) :: newentry  ! address in # table returned to calling routine
        integer :: ierr
        integer :: link_addr, cell_addr, ierror, free, addr_count
        logical :: resolved

        ierror = 0
        addr_count = 1
        free = free_addr(iused)

        ! perform address hashing on key
        cell_addr = int(IAND( keyin, hashconst))

        if ( htable( cell_addr )%node == 0 .and. htable( cell_addr )%key /=-1 ) then       ! Is entry empty?
            newentry = cell_addr                 ! Yes, so create new entry:
            htable( cell_addr )%node = nodein          !   local pointer
            htable( cell_addr )%key =  keyin            !   key
            htable( cell_addr )%leaves = leavesin       !   # contained nodes
            htable( cell_addr )%childcode = codein       !  child byte-code or particle #
            htable( cell_addr )%owner = ownerin       ! PE owning branch node

            if (point_free(cell_addr) /= 0) then     ! Check if new address in collision res. list
                free_addr( point_free(cell_addr) ) = free_addr(sum_unused)  ! Replace free address with last on list
                point_free(free_addr(sum_unused)) = point_free(cell_addr)   ! Reset pointer
                point_free(cell_addr) = 0
                sum_unused = sum_unused - 1
            endif

        else if ( htable( cell_addr )%node /= 0 .AND. htable(cell_addr)%key == keyin ) then
            ! Entry exists and keys match
            ! => local node or already inserted
            ierror = 1
            newentry = cell_addr

        else            ! Entry exists and keys do not match: COLLISION

            if ( htable( cell_addr )%link == -1) then     ! Entry was occupied without existing link

                newentry =  free  ! Pick next free # address from list of unused table positions (computed at end of treebuild)

                if (htable(free)%node /= 0 ) then
                    write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
                    write (*,*) 'PE ',me,' key ',keyin,' entry',newentry,' used ',iused,'/',sum_unused
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                endif
                htable( free )%node = nodein
                htable( free )%key = keyin
                htable( free )%childcode = codein
                htable( free )%leaves = leavesin
                htable( free )%owner = ownerin
                htable( free )%link = -1
                htable( cell_addr )%link = free     ! Create link from 1st duplicate entry to new entry
                iused = iused + 1        ! increment free address counter


            else if ( htable( cell_addr )%link /= -1 ) then     ! Address occupied with link already

                link_addr = cell_addr                          ! Start of chain
                resolved = .false.                                     ! Resolve flag

                do while (  .not. resolved .and. addr_count < maxaddress)       ! Find end of chain
                    link_addr = htable(link_addr)%link

                    if ( htable( link_addr )%key == keyin ) then
                        ! Occupied with same key -> local or boundary node, so skip
                        resolved = .true.
                        ierror = 1
                        newentry = link_addr

                    else if ( htable(link_addr)%node == 0 .and. htable (link_addr)%link == -1 ) then
                        ! Found end of chain: entry was occupied by boundary node, so reuse entry
                        ! link from previous chain entry still valid
                        newentry = link_addr
                        htable( link_addr )%node = nodein
                        htable( link_addr )%key = keyin
                        htable( link_addr )%childcode = codein
                        htable( link_addr )%leaves = leavesin
                        htable( link_addr )%owner = ownerin

                        if (point_free(link_addr) /= 0) then     ! Check if new address in collision res. list
                            free_addr( point_free(link_addr) ) = free_addr(sum_unused)  ! Replace free address with last on list
                            point_free(free_addr(sum_unused)) = point_free(link_addr)   ! Reset pointer
                            point_free(link_addr) = 0
                            sum_unused = sum_unused - 1
                        endif
                        resolved = .true.

                    else if ( htable(link_addr)%node /= 0 .and. htable (link_addr)%link == -1 ) then
                        ! Found end of chain: entry occupied with no link
                        newentry = free

                        if (htable(free)%node /= 0 ) then
                            write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
                            write (*,*) 'PE ',me,' key ',keyin,' entry',newentry,' used ',iused,'/',sum_unused
                            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                        endif

                        htable( free )%node = nodein
                        htable( free )%key = keyin
                        htable( free )%childcode = codein
                        htable( free )%leaves = leavesin
                        htable( free )%owner = ownerin
                        htable( free )%link = -1
                        htable( link_addr )%link = free     ! Create link from 1st duplicate entry to new entry
                        iused = iused + 1                  ! Increment free_address counter
                        resolved = .true.
                    else
                        addr_count = addr_count + 1
                       ! not yet resolved - go to next link in chain
                    endif
                end do

                if (addr_count >= maxaddress) then
                    write (ipefile,'(a5,o20,a)') 'Key ',keyin,' not resolved in MAKE_HASHENTRY'
                    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                endif
            else
                write (ipefile,'(a5,o20,a)') 'Key ',keyin,' not resolved in MAKE_HASHENTRY'
                ierror = 2
            end if

        end if


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

        key2addr = testaddr(keyin) ! cell address hash function

        if (key2addr == -1) then
          ! could not find key
          write(*,'("Key not resolved in KEY2ADDR at ",a)') cmark
          write(*,'("Bad address, check #-table and key list for PE", I7)') me
          write(*,'("key             = ", o22)') keyin
          write(*,'("initial address = ", i22)') int(IAND( keyin, hashconst))
          write(*,'("# const         = ", i22)') hashconst
          write(*,'("maxaddress      = ", i22)') maxaddress
          call diagnose_tree()
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

    end function key2addr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] keyin inverse mapping candidate.
    !> @return address if candidate exists, else -1
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function testaddr(keyin)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in)  :: keyin
        integer :: ires
        integer :: testaddr

        testaddr = int(IAND( keyin, hashconst))     ! cell address hash function
        ires     = 1 ! counter for number of htable lookups

        do while ( htable( testaddr )%key .ne. keyin )
          testaddr = htable( testaddr )%link    ! look at next linked entry

          ires = ires + 1

          if (   (testaddr == -1) & ! reached end of linked list without finding the key --> node is not in htable or htable is corrupt
            .or. (ires     >= maxaddress) ) & ! we probed at as many positions as the htable has entries --> circular linked list or htable corrupt
            then
              testaddr = -1
              return
          endif

        end do

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

        if (me==0) write (*,*) 'DIAGNOSE TREE'
        write (ipefile,*) 'DIAGNOSE TREE'

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
