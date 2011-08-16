!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates functions for accessing, manipulating, and verifying hash table data
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_htable
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Hash table datatype - 36 bytes per entry
    type :: hash
        integer   :: node          !< Address of particle/pseudoparticle data
        integer*8 :: key           !< Key
        integer   :: link          !< Pointer to next empty address in table in case of collision
        integer   :: leaves        !< # leaves contained within twig (=1 for leaf, npart for root)
        integer   :: childcode     !< Byte code indicating position of children (twig node); particle label (leaf node)
        integer*8 :: next          !< Pointer to next key to examine in tree-walk
        integer   :: owner         !< Node owner (for branches)
    end type hash

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type (hash), public, allocatable :: htable(:) !< hash table
    integer*8,   public ::  hashconst  !< hashing constants
    integer*8,   public ::  hashchild=7_8 !< bits that contain the child index in a key

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public get_next_node
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
    type (hash), private, parameter :: HASHENTRY_EMPTY = hash(0,0_8,-1,0,0,0_8,0) !< constant for empty hashentry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

    subroutine htable_clear()
        implicit none

        htable = HASHENTRY_EMPTY ! TODO: need list of 'live' adresses to speed this up

    end subroutine


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
                    get_next_node = ior(ishft(parent_key, 3), search_child_idx) ! assemble next_node out of parent-key and new sibling-index
                    return
                endif
            end do
        end do

        get_next_node  = 1 ! nothing has been found, i.e. top-right corner reached: set pointer=root

    end function get_next_node



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

        tablehigh = max(tablehigh,cell_addr)                 ! Track highest address

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
            tablehigh = max(tablehigh,free)                 ! Track highest address

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
        integer :: cell_addr, link_addr, ires,i, ierr
        logical :: resolved
        character(LEN=*) :: cmark
        integer :: key2addr

        cell_addr = int(IAND( keyin, hashconst))     ! cell address hash function

        if ( htable( cell_addr )%key == keyin ) then
            key2addr = cell_addr       ! Keys match -> found entry
            return
        else
            resolved = .false.
            link_addr = cell_addr
            ires = 0

            do while (  .not. resolved .and. ires <= maxaddress )       ! Repeat until keys match or run out of links
                link_addr = htable(link_addr)%link    ! Next linked entry
                if (link_addr == -1 ) then
                    write (*,'(a,a20)') 'Key not resolved in KEY2ADDR at ',cmark
                    write (*,*) 'check #-table and key list for PE ',me
                    write(*,*) 'Bad address'
                    write(*,'(a15,o22)') 'Key = ',keyin
                    write(*,*) 'Initial address =',cell_addr
                    write(*,*) '# const =',hashconst
                    write(*,*) 'ires =',ires
                    exit
                endif
                ires = ires + 1
                if ( htable( link_addr )%key == keyin ) then
                    key2addr = link_addr      ! Keys match -> found entry
                    return
                endif
            end do
            ! Not resolved - something wrong: invalid key or #-table wrong
            write (*,'(a5,o24,a2,i20,a1,a12,i15)') 'Key #: ',keyin,' (',keyin,')',' Address: ',cell_addr
            write (ipefile,*) 'Keys in table:'
            do i=0,maxaddress
                if (htable(i)%key/=0) write(ipefile,'(i8,o25,i10)') i,htable(i)%key,htable(i)%link
            end do

            !     call diagnose_tree
            close(75)
            !     call closefiles
            !     pause
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            stop
        endif
    end function key2addr


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculates inverse mapping from the hash-key to the hash-address.
    !>
    !> @param[in] keyin inverse mapping candidate.
    !> @return true if candidate exists
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function testaddr(keyin)

        use treevars
        implicit none
        include 'mpif.h'

        integer*8, intent(in)  :: keyin

        integer :: cell_addr, link_addr, ires
        logical :: resolved
        logical :: testaddr

        ! cell address hash function
        cell_addr = INT(IAND(keyin,hashconst))

        ! Keys match -> found entry
        if ( htable( cell_addr )%key == keyin ) then
            testaddr = .true.
            return
        else
            resolved = .false.
            link_addr = cell_addr
            ires = 0

            do while (.not.resolved .and. ires <= maxaddress )       ! Repeat until keys match or run out of links
                link_addr = htable(link_addr)%link    ! Next linked entry

                ! no more links, thus no entry for keyin
                if (link_addr == -1 ) then
                    testaddr = .false.
                    return
                endif

                ires = ires + 1

                ! Keys match -> found entry
                if ( htable( link_addr )%key == keyin ) then
                    testaddr = .true.
                    return
                endif
            end do
        endif
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
            if (htable(i)%node /= 0 .and. htable(i)%next >=0) write (ipefile,'(3i10,o22,i10,2o22,i8,i10,z4,4x,a1)') &
            me,htable(i)%owner,htable(i)%node,htable(i)%key,htable(i)%key,ishft( htable(i)%key,-3 ), htable(i)%next, &
            htable(i)%link,htable(i)%leaves,htable(i)%childcode,collision
            if (htable(i)%node /= 0 .and. htable(i)%next <0) write (ipefile,'(3i10,2o15,i15,i15,i5,z6,4x,a1)') &
            me,htable(i)%owner,htable(i)%node,htable(i)%key,ishft( htable(i)%key,-3 ), htable(i)%next, &
            htable(i)%link,htable(i)%leaves,htable(i)%childcode,collision
        end do


        ! get keys of twig nodes from hash table
        key_twig(1:ntwig)  = pack(htable(0:maxaddress)%key,mask=htable(0:maxaddress)%node<0)
        ! get levels of twigs
        addr_twig(1:ntwig) = (/( key2addr( key_twig(i),'DIAGNOSE_TREE' ),i=1,ntwig)/)   !  Table address
        child_twig(1:ntwig) = (/( htable( key2addr( key_twig(i),'DIAGNOSE_TREE' ) )%childcode,i=1,ntwig )/)   !  Children byte-code
        ind_twig(1:ntwig) = (/( htable( key2addr( key_twig(i),'DIAGNOSE_TREE' ) )%node,i=1,ntwig )/)   !  Twig node pointers

        do i=1,ntwig
            rcoc2(i) = xcoc(ind_twig(i))**2+ycoc(ind_twig(i))**2 + zcoc(ind_twig(i))**2
        end do

        write (ipefile,'(///a)') 'Tree structure'

        !  write (ipefile,'(/a/a/(3i5,2i10,2i8,b11,i2,i8,i10,9(1pe15.4)))') 'Twigs from hash-table:', &
        write (ipefile,'(/a/a,a/(3i5,2o15,2i8,z6,o15,i8,15(1pe30.19)))') 'Twigs from hash-table:', &
        '    i  level  owner        key     parent-key       #    node  code      1st child #leaves ', &
        ' abs_charge    charge   xcoc   ycoc   zcoc   xdip   ydip   zdip   sqrtbla   xxquad   yyquad   zzquad   xyquad   yzquad   zxquad', &
        (i,node_level(ind_twig(i)), &              !  index, level
        htable( key2addr( key_twig(i),'DIAGNOSE_TREE' ) )%owner, &                            ! Owner-PE of node
        key_twig(i),ishft( key_twig(i),-3 ), &                             ! key, parent key
        addr_twig(i), ind_twig(i), &    ! Table address and node number
        child_twig(i), &                         ! Children byte-code
        first_child( ind_twig(i) ), &            ! key of 1st child
        htable( addr_twig(i) )%leaves, &                           ! # leaves contained in branch
        abs_charge(ind_twig(i)), &    ! Twig absolute charge
        charge(ind_twig(i)), &    ! Twig  charge
        xcoc(ind_twig(i)), & ! Centre of charge
        ycoc(ind_twig(i)), &
        zcoc(ind_twig(i)), &
        xdip(ind_twig(i)), &
        ydip(ind_twig(i)), &
        zdip(ind_twig(i)), &
        sqrt(size_node(ind_twig(i))/htable(addr_twig(i))%leaves-rcoc2(i)), &
        xxquad(ind_twig(i)), &
        yyquad(ind_twig(i)), &
        zzquad(ind_twig(i)), &
        xyquad(ind_twig(i)), &
        yzquad(ind_twig(i)), &
        zxquad(ind_twig(i)), &
        i=1,ntwig)


        ! get keys of local leaf nodes from hash table
        key_leaf(1:nleaf_me) = pack(htable%key,mask=(htable%node>0 .and. htable%owner == me))
        ind_leaf(1:nleaf_me) = pack(htable%node,mask=(htable%node>0 .and. htable%owner == me))         ! particle/leaf index
        plist_leaf(1:nleaf_me) = pack(htable%childcode,mask=(htable%node>0 .and. htable%owner == me))   ! particle label
        owner_leaf(1:nleaf_me) = pack(htable%owner,mask=(htable%node>0 .and. htable%owner == me))   ! who owns leaf node


        write (ipefile,'(/a/3a5,2a10,2a15,a25,4a11/(3i5,2i10,2o15,o25,4f30.19))') 'Local leaves from hash-table:', &
        'i','owner','plab','i-leaf','lev','key','parent','pkey','x','y','z','q', &
        (i,owner_leaf(i),plist_leaf(i),ind_leaf(i),node_level(ind_leaf(i)),key_leaf(i), &
        ishft( key_leaf(i),-3 ), &      ! parent
        pekey(ind_leaf(i)), &  ! particle key
        x(ind_leaf(i)),y(ind_leaf(i)),z(ind_leaf(i)), q(ind_leaf(i)), &
        i=1,nleaf_me)

        ! get keys of NON-local leaf nodes from hash table
        key_leaf(1:nleaf-nleaf_me) = pack(htable%key,mask=(htable%node>0 .and. htable%owner /= me))
        ind_leaf(1:nleaf-nleaf_me) = pack(htable%node,mask=(htable%node>0 .and. htable%owner /= me))         ! leaf index
        plist_leaf(1:nleaf-nleaf_me) = pack(htable%childcode,mask=(htable%node>0 .and. htable%owner /= me))   ! global particle label
        owner_leaf(1:nleaf-nleaf_me) = pack(htable%owner,mask=(htable%node>0 .and. htable%owner /= me))   ! who owns leaf node


        write (ipefile,'(//a/a/(4i5,2o15,i5,2f11.4,f6.1,f11.4))') 'Non-local leaves from hash-table:', &
        '    i   owner    i-leaf    lev    key    parent  plabel  xcoc  ycoc  charge      ', &
        (i,owner_leaf(i),ind_leaf(i),node_level(ind_leaf(i)),key_leaf(i), &
        ishft( key_leaf(i),-3 ), &      ! parent
        plist_leaf(i), & ! global particle label
        xcoc(ind_leaf(i)),ycoc(ind_leaf(i)), charge(ind_leaf(i)), xdip(ind_leaf(i)), &
        i=1,nleaf-nleaf_me)

    end subroutine diagnose_tree


end module module_htable
