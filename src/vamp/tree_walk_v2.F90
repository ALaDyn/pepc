! ===========================================
!
!           TREE_WALK
!
!   Perform tree walk for all local particles
!
!  Algorithm follows Warren & Salmon's 'latency-hiding' concept,
!  but retains list-based tree-walk from vectorised code.
!
!  Structure:
!
!    while (any particle still active) do
!
!      while (any particle not finished walk) do
!
!         if (MAC OK) 
!            put node on interaction list
!         else if (MAC not OK and local)
!            resolve node
!         else 
!            put node on defer list
!         endif
!
!      end while
!      
!      process requests for non-local nodes
!      insert new children into local tree (htable)
!      copy defer list -> new walk list
!      count remaining active particles
!
!    end while
!
! ===========================================

subroutine tree_walk(pshort,npshort)

  use treevars
  use utils
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none

  integer, intent(in) :: npshort
  integer, intent(in) :: pshort(nshortm)

  integer :: npackm   ! Max # children shipped
  integer :: nchild_shipm

! Key arrays (64-bit)

  integer*8,  dimension(nshortm) :: walk_key, add_key, walk_next 
  integer*8, dimension(size_tree*2)  :: request_key, ask_key, process_key
  integer*8, dimension(nppm,0:num_pe-1) ::  ship_key
  integer*8, dimension(8) :: sub_key, key_child, next_child
  integer*8, dimension(nintmax,nshortm) :: defer_list, walk_list
  integer*8, dimension(size_tree) :: last_child                      ! List of 'last' children fetched from remote PEs
  
  integer, dimension(nshortm) :: plist

  integer, dimension(nshortm) ::  walk_addr, walk_node 
  integer, dimension(nshortm) :: nlocal, ndefer,  nwalk, defer_ctr


  integer, dimension(size_tree) ::  process_addr, request_owner
  integer, dimension(5*nppm) :: childbyte
  integer, dimension(8) :: addr_child, node_child, byte_child, leaves_child

  real, dimension(8) :: xcoc_child, ycoc_child, zcoc_child

  real, dimension(0:nlev) :: boxlength
  logical, dimension(nshortm) :: finished, ignore, mac_ok, requested

  integer, dimension(num_pe) ::   nactives, nreq
  integer, dimension(0:num_pe-1) :: ntoship, nrequested, istart

! Key working vars
  integer*8 :: node_key, kchild, kparent, search_key,  nxchild

  integer :: i, j, k, ic, ipe, iwait, inner_pass          ! loop counters
  integer :: nnew, nshare, nplace, newrequest, nreqs, i1, i2, send_count, send_prop_count, ioff, rec_count, ipack, source_id
  integer :: nchild, newleaf, newtwig, nactive, maxactive, ntraversals

  integer ::  bchild, nodchild, lchild, hashaddr, nlast_child
  integer :: max_nplace, max_pack
  real :: sbox, theta2, dx, dy, dz, s2

! stuff for tree-patch after traversals complete
  integer ::  node_addr, parent_addr, parent_node, child_byte
  integer :: jmatch(1)
  logical :: resolved, keymatch(8), emulate_blocking=.false.

  integer :: key2addr        ! Mapping function to get hash table address from key

!VAMPINST subroutine_start
       CALL VTENTER(IF_tree_walk,VTNOSCL,VTIERR)
!      write(*,*) 'VT: tree_walk S>',VTIERR,
!     *    IF_tree_walk,ICLASSH
!
  npackm = size_tree
  nchild_shipm = size_tree
  walk_debug = .true.

  if (walk_debug) write(ipefile,'(/a)') 'TREE WALK'

  sbox = boxsize


  theta2 = theta**2               ! Clumping parameter**2 for MAC

  boxlength(0:nlev) = (/ (sbox/2**i, i=0,nlev ) /)  ! Preprocessed box sizes for each level

  walk_key(1:npshort) = 1                    ! initial walk list starts at root
  nlist = npshort                 ! Inner loop list length
  plist(1:npshort) = (/ (i,i=1,npshort) /)       ! initial local particle indices 
  nterm = 0
  nactive = npshort
!  Find global max active particles - necessary if some PEs enter walk on dummy pass

  call mpi_allreduce( nactive, maxactive, one, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr )
 

  ntraversals = 0
  nlast_child = 0

  do while (maxactive > 0)        ! Outer loop over 'active' traversals

     ntraversals = ntraversals + 1  ! # Tree-walks
     finished(1:nlist) = .false.  ! label all particles not finished
     ndefer(1:npshort) = 0          ! Zero # deferrals (absolute particle number)

     nshare = 0                   ! Total number of requests for nonlocal keys
     nplace = 0                   ! Total # non-local children to be fetched
     newleaf = 0
     newtwig = 0                  ! local bookkeeping
     inner_pass = 0


     do while (nlist > 0 )        ! Inner loop - tree walk

        inner_pass = inner_pass+1
        walk_addr(1:nlist) = (/( key2addr( walk_key(j) ),j=1,nlist)/)      ! get htable addresses
        walk_node(1:nlist) = htable( walk_addr(1:nlist) )%node             ! Walk node index - points to multipole moments  
        walk_next(1:nlist) = htable( walk_addr(1:nlist) )%next             ! Next node pointer 

        childbyte(1:nlist) = htable( walk_addr(1:nlist) )%childcode        ! Children byte-code

        ! set ignore flag if leaf node corresponds to particle itself
        ignore(1:nlist) = (/ ( pshort(plist(i)) == htable( walk_addr(i) )%node,i=1,nlist) /)

        ! children of local/non-local parents already fetched: HERE flag


        do i=1,nlist
           dx = x( pshort(plist(i)) ) - xcoc( walk_node(i) )      ! Separations
           dy = y( pshort(plist(i)) ) - ycoc( walk_node(i) )
           dz = z( pshort(plist(i)) ) - zcoc( walk_node(i) )

           s2 = boxlength( node_level(walk_node(i)) )**2
           mac_ok(i) = ( s2 < theta2*(dx**2+dy**2+dz**2 ) )             ! Preprocess MAC
        end do
        add_key(1:nlist) = walk_key(1:nlist)                                ! Remember current key


        ! Possible courses of action:


        do i=1,nlist


           ! 1) MAC test OK, so put cell on interaction list and find next node for tree walk
           if ( mac_ok(i) .or. (walk_node(i) >0 .and. .not.ignore(i) ) ) then
              walk_key(i) = walk_next(i)
              nterm(plist(i)) = nterm(plist(i)) + 1
              intlist( nterm( plist(i)), plist(i) ) = add_key(i)      ! Augment interaction list
              nodelist( nterm( plist(i)), plist(i) ) = walk_node(i)   ! Node number


              ! 2) MAC fails at node for which children present, so resolve cell & put 1st child on walk_list
           else  if ( .not.mac_ok(i) .and. walk_node(i) < 0 .and. btest(childbyte(i),9) ) then
              ! if local put 1st child node on walk_list
              walk_key(i) = first_child( walk_node(i) )


              ! 3) MAC fails at node for which children _absent_, so put node on REQUEST list (flag with add=2)
           else if ( .not.mac_ok(i) .and. walk_node(i) < 0 .and.  .not. btest(childbyte(i),9) ) then
              walk_key(i) = walk_next(i)  ! Continue with walk for now
              ndefer(plist(i)) = ndefer(plist(i)) + 1
              defer_list( ndefer( plist(i)), plist(i) ) = add_key(i)  ! Deferred list of nodes to search, pending request
              ! for data from nonlocal PEs
              if (.not. BTEST( htable(walk_addr(i))%childcode, 8 ) ) then  ! Check if node already requested
                 nshare = nshare + 1
                 nplace = nplace + n_children(walk_node(i))     ! Total # children to be fetched from remote PEs
                 request_key(nshare) = add_key(i)       ! New request key
                 request_owner(nshare) = htable( walk_addr(i) )%owner           ! Owner 
                 htable(walk_addr(i))%childcode   =  IBSET( htable(walk_addr(i))%childcode, 8 ) ! Set requested flag

              endif

              ! 1) particle and leaf node identical, so skip

           else
              walk_key(i) = walk_next(i)

           endif
        end do



        do i=1,nlist
           if (walk_key(i) == -1 )  then         ! Trap condition for last of nonlocal children - skip to next deferred node
              defer_ctr(plist(i)) = defer_ctr(plist(i)) + 1
              walk_key(i) = walk_list (defer_ctr(plist(i)),plist(i) ) ! Select next deferred node from walk list for particle plist(i)
           endif
        end do

        where ( walk_key(1:nlist) == 1 ) 
           finished(1:nlist) = .true.    ! Flag particles whose walks are complete
           elsewhere
           finished(1:nlist) = .false.  ! label all particles not finished

        end where

        nnew = count( mask = .not.finished(1:nlist) )                       ! Count remaining particles

        plist(1:nnew) =  pack( plist(1:nlist), mask = .not.finished(1:nlist) )    ! Compress particle index list
        walk_key(1:nnew) =  pack( walk_key(1:nlist), mask = .not.finished(1:nlist) )       ! Compress walk lists etc.

        nlist = nnew
     end do

     if (walk_debug) then
        write(ipefile,'(a/(o15,i7))') 'Shared request list: ',(request_key(i),htable( key2addr( request_key(i) ) )%owner,i=1,nshare)
!        if (me==2) then
!           write(*,'(a/(o15,i7))') 'Shared request list: ',(request_key(i),htable( key2addr( request_key(i) ) )%owner,i=1,nshare)
!        endif
     endif
     ! At this point all particles have completed their walks with the locally available tree data.
     ! Each PE has list of nonlocal nodes which have been requested during these walks - request_key(1:nshare)
     ! Need to ship these requests to the PEs on which the child-data is actually held (cf particle swapping).

     ! First find out how many requests are to be sent to each PE.

     ntoship(0:num_pe-1) = (/ (count( mask = request_owner(1:nshare) == ipe ), ipe=0,num_pe-1) /)

     call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! Wait for other PEs to catch up


     ! Send off the request numbers to the appropriate destination & initiate receives
     do ipe=0,num_pe-1
        call MPI_IRECV( nreq(ipe+1), one, MPI_INTEGER, MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, pe_handle(ipe+1), ierr)
        call MPI_ISEND( ntoship(ipe), one, MPI_INTEGER, ipe, tag1, MPI_COMM_WORLD, pe_handle(num_pe+ipe+1), ierr)
     end do

     ! Wait for incoming requests and sort according to rank
     nbuf = 2*num_pe

     do i=1,2*num_pe
        call MPI_WAITANY( nbuf, pe_handle, index, status, ierr)  ! Wait for one of receives or sends to complete
        if (index .le. num_pe) then
           ! receive complete
           source_id = status(MPI_SOURCE)
           nrequested(source_id) = nreq(index)
        else
           ! send complete
        endif
     end do



  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first




     ! Can now ship the keys in the request lists to the PEs where the data sits:
     !   the recipients now know how many requests to expect.

     send_count=0  ! Handle counter

     do ipe = 0,num_pe-1 
        if ( ipe /= me .and. ntoship(ipe) > 0 ) then   ! avoid shipping to oneself and shipping nothing
           ! Extract sub-list according to location - don't overwrite buffer!
           ship_key(1:ntoship(ipe),ipe) = pack(request_key(1:nshare), mask = request_owner(1:nshare) == ipe ) 
	   send_count=send_count+1
           nbuf = ntoship(ipe)
           destpe = ipe
           call MPI_ISEND(ship_key(1:ntoship(ipe),ipe), nbuf, MPI_INTEGER8, destpe, tag1, &
                MPI_COMM_WORLD, send_key_handle(send_count), ierr ) ! Ship to data location
           call MPI_REQUEST_FREE( send_key_handle(send_count), ierr)
        endif

     end do

!     do i=1,send_count
 !       call MPI_WAIT( send_key_handle(send_count), ierr)
 !    end do

     ! Initiate receives for request keys: could reuse shared request list
     istart(me) = 0    ! Starting point of key list for given rank
     i1=1
     rec_count = 0
     ipack = 0

     do ipe = 0,num_pe-1
        if ( ipe /= me .and. nrequested(ipe) > 0 ) then   ! skip self or if nothing requested
           istart(ipe) = i1
           i2 = i1+nrequested(ipe)-1
           rec_count = rec_count + 1  ! receive counter
           nbuf =  nrequested(ipe)
           destpe = ipe
           call MPI_IRECV( process_key(i1:i2), nbuf, MPI_INTEGER8, destpe, tag1, &
                           MPI_COMM_WORLD, recv_key_handle(rec_count), ierr)
           i1 = i1+nrequested(ipe)
        endif
     end do



     ! Now have complete list of requests from all PEs in rank order. 
     ! Fetch child data from local htable and send back to requesting PE (recipient already knows how many children to expect)

     send_prop_count = 0

     do iwait = 1,rec_count
        nbuf = rec_count
        call MPI_WAITANY( nbuf, recv_key_handle, index, status, ierr)  ! Wait for one of receives to complete
        ipe = status(MPI_SOURCE)
        nreqs = nrequested(ipe)  ! # parent keys
        i1 = istart(ipe)
        i2 = istart(ipe) + nrequested(ipe) - 1
        if (walk_debug) then
           write(ipefile,'(a,3i7/(i7))') '# Received request from PE: ',ipe,nrequested(ipe),istart(ipe),process_key(i1:i2) 
        endif
        process_addr(1:nreqs) = (/( key2addr( process_key(j) ),j=i1,i2)/)    ! get htable addresses
        childbyte(1:nreqs) = htable( process_addr(1:nreqs) )%childcode        !  Children byte-code


        ! For each key in the request list, fetch and package tree info for children

        do i=1,nreqs
           ipack = i1+i-1   ! Request number: this needs to be unique for each SEND, otherwise buffer pack_child(..,ipack)
           ! may get overwritten before send actually completes
           nchild = SUM( (/ (ibits(childbyte(i),j,1),j=0,2**idim-1) /) )                     ! Get # children
           sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(childbyte(i),j),j=0,7) /) )      ! Extract sub key from byte code
           key_child(1:nchild) = IOR( ishft( process_key(ipack),idim ), sub_key(1:nchild) ) ! Construct keys of children
           addr_child(1:nchild) = (/( key2addr( key_child(j) ),j=1,nchild)/)                 ! Table address of children
           node_child(1:nchild) = htable( addr_child(1:nchild) )%node                        ! Child node index  
           byte_child(1:nchild) = IAND( htable( addr_child(1:nchild) )%childcode,255 )        ! Catch lowest 8 bits of childbyte - filter off requested and here flags 
           leaves_child(1:nchild) = htable( addr_child(1:nchild) )%leaves                    ! # contained leaves
           next_child(1:nchild-1) = htable( addr_child(1:nchild-1) )%next                    ! # next-node pointer
           next_child(nchild) = -1                             ! Last child gets pointed back to _parent_ for non-local nodes
           ! This is used to distinguish particles' walks during 'defer' phase

           ! Package children properties into user-defined multipole array for shipping
           do ic = 1,nchild
              send_prop_count = send_prop_count+1
              pack_child(ic,ipack) = multipole ( key_child(ic), &
               				 	 byte_child(ic), &
               					 leaves_child(ic), &
              	   			 	 next_child(ic), &
              				 	 charge( node_child(ic) ), &
              					 abs_charge( node_child(ic) ), &
             				 	 xcoc( node_child(ic)), &
              				 	 ycoc( node_child(ic)), &
              					 zcoc( node_child(ic)), &
              					 xdip( node_child(ic)), &
              					 ydip( node_child(ic)), &
              					 zdip( node_child(ic)), &
              				 	 xxquad( node_child(ic)), &
              				 	 yyquad( node_child(ic)), &
              					 zzquad( node_child(ic)), &
              					 xyquad( node_child(ic)), &
              					 yzquad( node_child(ic)), &
              				 	 zxquad( node_child(ic)) )

              ! Ship child data back to PE that requested it
              destpe = ipe

              call MPI_ISEND( pack_child(ic,ipack), 1, MPI_type_multipole, destpe, tag1,&
                   		MPI_COMM_WORLD, send_child_handle(send_prop_count), ierr )
	      call MPI_REQUEST_FREE(send_child_handle(send_prop_count), ierr) 

            end do

	end do
        if (walk_debug) write(ipefile,*) 'Total children shipped to processor ',ipe,' from ',me,' was', send_prop_count 

     end do


     ! Initiate receives for returning child nodes

     do i=1,nplace
        call MPI_IRECV( get_child(i), 1, MPI_type_multipole, MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, recv_child_handle(i), ierr)
     end do

!     if (emulate_blocking) then
!        call MPI_WAITALL( nplace, recv_handle, stat_pe, ierr )  ! Simulate blocking recv: wait for ALL requested data to arrive before processing
!     endif

     ! Wait for data to arrive

     do i=1,nplace
        nbuf = nplace
        call MPI_WAITANY( nbuf, recv_child_handle, index, status, ierr)  ! Wait for one of receives to complete
        source_id = status(MPI_SOURCE)
        kchild = get_child(index)%key
        kparent = ishft( kchild,-idim )
        bchild = get_child(index)%byte
        lchild = get_child(index)%leaves
        nxchild = get_child(index)%next

        if (lchild ==1 ) then
           newleaf = newleaf + 1
           nleaf = nleaf + 1
           nodchild = nleaf
           n_children(nodchild) = 0
           first_child(nodchild) = kchild

        else if (lchild > 1) then
           newtwig = newtwig + 1
           ntwig = ntwig + 1
           nodchild = -ntwig
           nchild = SUM( (/ (ibits(bchild,j,1),j=0,2**idim-1) /) )   ! Get # children
           n_children( nodchild ) = nchild       
           sub_key(1:nchild) = pack( bitarr(0:7), mask=(/ (btest(bchild,j),j=0,7) /) )  ! Extract child sub-keys from byte code
           first_child( nodchild ) = IOR( ishft( kchild,idim), sub_key(1) )              ! Construct key of 1st (grand)child

        else
           write(ipefile,'(a,o15,a,i7)') '# leaves <=0 for received child node ',kchild,' from PE ',source_id
        endif

        if (walk_debug) write(ipefile,'(a,o15,a,i7,a,o13)') &
	     'Child data arrived:',kchild,' from ',source_id,' requested for key ',kparent

        ! Insert new node into local #-table

        call make_hashentry( kchild, nodchild, lchild, bchild, source_id, hashaddr, ierr )

        htable(hashaddr)%next = nxchild           ! Fill in special next-node pointer for non-local children
        htable( key2addr( kparent) )%childcode = IBSET(  htable( key2addr( kparent) )%childcode, 9) ! Set children_HERE flag for parent node

        node_level( nodchild ) = log(1.*kchild)/log(2.**idim)  ! get level from keys and prestore as node property

        ! Physical properties

        charge( nodchild ) = get_child(index)%q
        abs_charge( nodchild ) = get_child(index)%absq
        xcoc( nodchild ) = get_child(index)%xcoc
        ycoc( nodchild ) = get_child(index)%ycoc
        zcoc( nodchild ) = get_child(index)%zcoc
        xdip( nodchild ) = get_child(index)%xdip
        ydip( nodchild ) = get_child(index)%ydip
        zdip( nodchild ) = get_child(index)%zdip
        xxquad( nodchild ) = get_child(index)%xxquad
        yyquad( nodchild ) = get_child(index)%yyquad
        zzquad( nodchild ) = get_child(index)%zzquad
        xyquad( nodchild ) = get_child(index)%xyquad
        yzquad( nodchild ) = get_child(index)%yzquad
        zxquad( nodchild ) = get_child(index)%zxquad

 ! Put last child onto list for post-traversal processing
	if (nxchild == -1) then
	   nlast_child = nlast_child + 1
	   last_child(nlast_child) = kchild
	endif

     end do


     ! Copy defer lists to new walk lists for next tree-walk iteration

     do i=1,npshort
        nwalk(i) = ndefer(i)
        if (nwalk(i) /= 0) then
           nlist = nlist + 1
           walk_list( 1:nwalk(i), i ) = defer_list( 1:nwalk(i), i )
           walk_list( nwalk(i) + 1, i ) = 1        ! Last node root for correct 'next-node'
           walk_key(nlist) = defer_list(1, i)       ! Start node for next walk
           plist(nlist) = i                        ! Particle index
           defer_ctr(i) = 1                        ! Deferral counter (1 ... nwalk)
        endif
     enddo

     nactive = count( mask = nwalk(1:npshort) /= 0 )     ! Count remaining 'active' particles - those still with deferred nodes to search

     call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! Wait for other PEs to catch up

     ! Broadcast # remaining particles to other PEs

     call MPI_ALLGATHER( nactive, one, MPI_INTEGER, nactives, one, MPI_INTEGER, MPI_COMM_WORLD, ierr )

     maxactive = maxval(nactives)

! Determine global max
     call MPI_ALLREDUCE( ipack, max_pack, one, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
     call MPI_ALLREDUCE( nplace, max_nplace, one, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )

     if (walk_debug .and. mod(itime,idump)==0 .and. ntraversals==1) then
        write (ipefile,*) 'Traversal # ',ntraversals,' :'
        write (ipefile,*) 'max children shipped: ',max_pack,' /',npackm
        write (ipefile,*) 'max non-local children fetched: ',max_nplace,' /',nchild_shipm
     endif
     if (walk_debug) then
        write (ipefile,'(3(/a30,i6)/a/(2i5))') 'New twigs: ',newtwig, &
             'New leaves:',newleaf, &
             'New list length: ',nlist, &
             '# remaining active particles on each PE: ',(i,nactives(i+1),i=0,num_pe-1)
        write(ipefile,'(a/(2i5))') 'New shortlist: ',(plist(i),pshort(plist(i)),i=1,nlist)

     endif
  end do



  !  Determine 'next_node' pointers for 'last child' list & update hash table (tree).
  !  This  ensures that traversals in next pass treat already-fetched nodes as local,
  !  avoiding deferral list completely.

 ! ** TODO:  Put this in separate routine FIND_NEXTNODE

  do i = 1,nlast_child

     search_key = last_child(i)                   
     node_key = search_key                     ! keep key, address of node 
     node_addr = key2addr(search_key)
     resolved = .false.

     !   Search for next sibling, uncle, great-uncle etc

     do while (.not. resolved .and. search_key > 1)
        kparent =  ishft(search_key,-idim)                 ! parent
	parent_addr = key2addr( kparent)
        parent_node = htable( parent_addr )%node   ! parent node pointer

        child_byte = htable( parent_addr )%childcode                           !  Children byte-code
        nchild = SUM( (/ (ibits(child_byte,j,1),j=0,2**idim-1) /) )                   ! # children = sum of bits in byte-code
        sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(child_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
        key_child(1:nchild) = IOR( ishft(kparent,idim), sub_key(1:nchild) )         ! Construct keys of children

        keymatch=.false.
        keymatch(1:nchild) = (/ (key_child(j) == search_key,j=1,nchild) /)
        jmatch = pack(bitarr, mask = keymatch ) + 1                                  ! Pick out position of current search key in family

        if (jmatch(1) < nchild ) then                                                ! if search_key has 'elder' sibling then
           htable( node_addr )%next  = key_child(jmatch(1)+1)                        ! store next_node as sibling of parent/grandparent
           resolved = .true.
        else
           search_key = ishft(search_key, -idim)                                     ! Go up one level 
        endif
     end do

     if (.not. resolved .and. search_key == 1)  htable( node_addr )%next  = 1        ! Top-right corner reached: set pointer=root

  end do

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: tree_walk S<',VTIERR,ICLASSH
!
end subroutine tree_walk
