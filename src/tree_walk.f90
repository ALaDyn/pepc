! ===========================================
!
!           TREE_WALK
!
!   $Revision$
!
!   Perform tree walk for all local particles
!
!  Algorithm follows Warren & Salmon's 'latency-hiding' concept,
!  retaining list-based tree-walk from vectorised code by Pfalzner & Gibbon.
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

subroutine tree_walk(pshort,npshort,pass,theta,itime)

  use treevars
  use utils
  !VAMPINST include
  !  INCLUDE 'VTcommon.h'
  !  INTEGER VTIERR
  !

  implicit none

  real, intent(in) :: theta
  integer, intent(in) :: npshort,itime
  integer, intent(in) :: pshort(npshort)

  integer :: npackm   ! Max # children shipped
  integer :: nchild_shipm

  ! Key arrays (64-bit)

  integer*8,  dimension(npshort) :: walk_key, walk_last 
  integer*8, dimension(size_tree)  :: request_key, ask_key, process_key
  integer*8, dimension(size_tree/3,0:num_pe-1) ::  ship_key
  integer*8, dimension(8) :: sub_key, key_child, next_child
  integer*8, dimension(nintmax,npshort) :: defer_list, walk_list
  integer*8, dimension(size_tree) :: last_child   ! List of 'last' children fetched from remote PEs

  integer, dimension(npshort) :: plist

  integer ::  walk_addr, walk_node, entry_next 
  integer, dimension(npshort) :: nlocal, ndefer,  nwalk, defer_ctr


  integer, dimension(size_tree) ::  process_addr, request_owner
  integer, dimension(size_tree) :: childbyte
  integer, dimension(8) :: addr_child, node_child, byte_child, leaves_child

  real, dimension(8) :: xcoc_child, ycoc_child, zcoc_child

  real, dimension(0:nlev) :: boxlength
  logical, dimension(npshort) :: finished, requested
  integer :: hops(21) ! array to control max # iterations in single traversal 
  integer, dimension(num_pe) ::   nactives
  integer, dimension(0:num_pe-1) :: ntoship, &              ! # keys needed
       nrequested, &           ! # keys requested from elsewhere
       istart, ic_start, &     ! # fenceposts
       nplace,&                ! # children (new entries) to place in table
       nchild_ship       ! # children shipped to others

  ! Key working vars
  integer*8 :: node_key, add_key, walk_next, kchild, kparent, search_key,  nxchild

  integer :: i, j, k, ic, ipe, iwait, inner_pass, nhops          ! loop counters
  integer :: p,pass
  integer :: nnew, nshare, newrequest, nreqs, i1, i2, ioff, ipack
  integer :: nchild, newleaf, newtwig, nactive, maxactive, ntraversals, own
  integer :: ic1, ic2, ihand, nchild_ship_tot, nplace_max
  integer, save ::  sum_nhops, sum_inner_pass, sum_nhops_old=0, sum_inner_old=0
  integer :: request_count, fetch_pe_count, send_prop_count  ! buffer counters

  integer ::  bchild, nodchild, lchild, hashaddr, nlast_child, cbyte
  integer :: max_nplace, max_pack
  real :: sbox, theta2, theta2_ion, dx, dy, dz, s2

  ! stuff for tree-patch after traversals complete
  integer ::  node_addr, parent_addr, parent_node, child_byte
  integer :: jmatch(1)
  logical :: resolved, keymatch(8), emulate_blocking=.false.
  logical :: ignore, mac_ok

  integer :: nrest, ndef

  integer :: key2addr        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk


  !
  npackm = size_tree
  nchild_shipm = size_tree
  !walk_debug = .false.
  ! ipefile = 6
  if (walk_debug .or. walk_summary) write(ipefile,'(/2(a,i6))') '*** TREE WALK for timestep ',itime,' pass ',pass
  if (me.eq.0 .and. walk_summary) write(*,'(/2(a,i6))') 'TREE WALK for timestep ',itime,' pass ',pass

  sbox = boxsize


  theta2 = theta**2               ! Clumping parameter**2 for MAC
!  theta2_ion = min(1.0,2*theta2)  ! Ion MAC 50% larger than electron MAC
  theta2_ion=theta2
  boxlength(0:nlev) = (/ (sbox/2**i, i=0,nlev ) /)  ! Preprocessed box sizes for each level

  walk_key(1:npshort) = 1                    ! initial walk list starts at root
  nwalk(1:npshort) = 0   ! # keys on deferred list
  defer_ctr(1:npshort) = 1   ! Deferral counter
  nlist = npshort                 ! Inner loop list length
  plist(1:npshort) = (/ (i,i=1,npshort) /)       ! initial local particle indices 
  nterm = 0
  nactive = npshort

  !  Find global max active particles - necessary if some PEs enter walk on dummy pass

  !  call MPI_ALLREDUCE( nactive, maxactive, one, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLGATHER( nactive, one, MPI_INTEGER, nactives, one, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  maxactive = maxval(nactives)

  ntraversals = 0
  nlast_child = 0
  nchild_ship = 0
  sum_nhops = 0
  sum_inner_pass = 0
  if (sum_nhops_old < 0) then
     hops(1) = sum_nhops_old*0.9 
     hops(2:4) = hops(1:3)*0.2
     hops(5:21) = 2**30
  else
     hops = 2**30    ! First time round set infinite - finish all traversals
     !     hops(1) = 500
  endif

  do while (maxactive > 0)        ! Outer loop over 'active' traversals


     ntraversals = ntraversals + 1  ! # Tree-walks
     if (walk_debug) write(ipefile,'(a,i6)') 'Start of traversal ',ntraversals
     !     if (walk_debug) write(*,'(a,i6)') 'Start of traversal ',ntraversals

     finished(1:nlist) = .false.  ! label all particles not finished
     ndefer(1:npshort) = 0          ! Zero # deferrals (absolute particle number)

     nshare = 0                   ! Total number of requests for nonlocal keys
     nplace = 0                   ! Total # non-local children to be fetched
     newleaf = 0
     newtwig = 0                  ! local bookkeeping
     inner_pass = 0
     nhops = 0

     !VAMPINST subroutine_start
     !    CALL VTENTER(IF_tree_walk,VTNOSCL,VTIERR)
     !      write(*,*) 'VT: tree_walk S>',VTIERR,
     !     *    IF_tree_walk,ICLASSH

     do while (nlist>0 .and. nhops <= hops(ntraversals) )        ! Inner loop - single hop in tree walk
        inner_pass = inner_pass+1        ! statistics
        nhops = nhops + nlist
        sum_nhops = sum_nhops + nlist
        sum_inner_pass = sum_inner_pass + 1

        do i=1,nlist
           p = plist(i)  ! particle index
           walk_addr =  key2addr( walk_key(i) )     ! get htable address
           walk_node = htable( walk_addr )%node             ! Walk node index - points to multipole moments  
           walk_next = htable( walk_addr )%next             ! Next node pointer 

           cbyte = htable( walk_addr )%childcode        ! Children byte-code

           ! set ignore flag if leaf node corresponds to particle itself (number in pshort)
           ignore =  ( pshort(p) == htable( walk_addr )%node )

           ! children of local/non-local parents already fetched: HERE flag


           dx = x( pshort(p) ) - xcoc( walk_node )      ! Separations
           dy = y( pshort(p) ) - ycoc( walk_node )
           dz = z( pshort(p) ) - zcoc( walk_node )

           s2 = boxlength( node_level(walk_node) )**2
 
           mac_ok = ( s2 < theta2*(dx**2+dy**2+dz**2 ) )             ! Preprocess MAC
 
           add_key = walk_key(i)                                ! Remember current key

           ! Possible courses of action:


           ! 1) MAC test OK, so put cell on interaction list and find next node for tree walk
           if ( mac_ok .or. (walk_node >0 .and. .not.ignore ) ) then
              walk_key(i) = walk_next
	      entry_next = nterm(p) + 1
              intlist( entry_next, p ) = add_key      ! Augment interaction list
              nodelist( entry_next, p ) = walk_node   ! Node number
              nterm(p) = entry_next


              ! 2) MAC fails at node for which children present, so resolve cell & put 1st child on walk_list
           else  if ( .not.mac_ok .and. walk_node < 0 .and. btest(cbyte,9) ) then
              ! if local put 1st child node on walk_list
              walk_key(i) = first_child( walk_node )


              ! 3) MAC fails at node for which children _absent_, so put node on REQUEST list (flag with add=2)
           else if ( .not.mac_ok .and. walk_node < 0 .and.  .not. btest(cbyte,9) ) then
              walk_key(i) = walk_next  ! Continue with walk for now
              ndefer(p) = ndefer(p) + 1
              defer_list( ndefer( p), p ) = add_key  ! Deferred list of nodes to search, pending request
              ! for data from nonlocal PEs
              if (.not. BTEST( htable(walk_addr)%childcode, 8 ) ) then  ! Check if node already requested
                 nshare = nshare + 1

                 request_key(nshare) = add_key       ! New request key
                 own = htable( walk_addr )%owner           ! Owner 
                 request_owner(nshare) = own
                 nplace(own) = nplace(own) + n_children(walk_node)     ! Total # children to be fetched from remote PEs
                 htable(walk_addr)%childcode   =  IBSET( htable(walk_addr)%childcode, 8 ) ! Set requested flag

              endif

              ! 1) particle and leaf node identical, so skip

           else
              walk_key(i) = walk_next

           endif


           ! Trap condition for last of nonlocal children - fetch next deferred node from walk_list
           if (walk_key(i) == -1 )  then       
              defer_ctr(p) = defer_ctr(p) + 1
              walk_key(i) = walk_list (defer_ctr(p),p ) ! Select next deferred node from walk list for particle plist(i)
           endif

           walk_last(p) = walk_key(i)  ! Store last reached in traversal


           ! Check for completed traversals
           if ( walk_key(i) == 1 )  then  ! Reached root  
              finished(i) = .true.    ! Flag particles whose local walks are complete
              defer_ctr(p) = defer_ctr(p) + 1

           else 
              finished(i) = .false. 
           endif
        end do

!        nnew = count( mask = .not.finished(1:nlist) )            ! Count remaining particles

!        plist(1:nnew) =  pack( plist(1:nlist), mask = .not.finished(1:nlist) )    ! Compress particle index list
!        walk_key(1:nnew) =  pack( walk_key(1:nlist), mask = .not.finished(1:nlist) )       ! Compress walk lists etc.
	nnew=0
	do i=1,nlist
	  if (.not.finished(i)) then
	    nnew=nnew+1
	    plist(nnew) = plist(i)
	    walk_key(nnew) = walk_key(i)
	  endif
        end do
  
        nlist = nnew

     end do   ! END DO_WHILE

     ! For remaining unfinished particles, need to copy rest of walk_list (still not inspected)
     ! back onto defer list.

     do i=1,npshort
        ndef = ndefer(i)  ! # new deferrals added to defer_list during traversal
        nrest = nwalk(i)-defer_ctr(i) + 1  ! # unprocessed deferrals from previous fetch
        if (nrest>0) then  
           ! Augment defer list with rest of walk_list
           defer_list(ndef+1:ndef+nrest,i) = walk_list(defer_ctr(i):nwalk(i),i)
           ndefer(i) = ndefer(i) + nrest
        endif
        if (walk_last(i) /= 1) then
           defer_list(ndef+nrest+1,i) = walk_last(i)  ! Tack on key where local walk was interrupted
           ndefer(i) = ndefer(i) + 1
        endif

     end do



     !VAMPINST subroutine_end
     !CALL VTLEAVE(ICLASSH,VTIERR)
     !      write(*,*) 'VT: tree_walk S<',VTIERR,ICLASSH
     !

     !CALL VTENTER(IF_tree_nlswap,VTNOSCL,VTIERR)
     !      write(*,*) 'VT: tree_walk S>',VTIERR,
     !     *    IF_tree_walk,ICLASSH
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


     ! Exchange numbers of keys to be shipped and requested

     call MPI_ALLTOALL( ntoship, 1, MPI_INTEGER, nrequested, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

     if (walk_debug) then
        write (ipefile,*) 
        write (ipefile,*) 'PE ',me,': Keys to request: ',ntoship(0:num_pe-1),' to process: ',nrequested(0:num_pe-1)
     endif


     ! Initiate receives for incoming request keys: could reuse shared request list
     istart(me) = 0    ! Starting point of key list for given rank
     i1=1
     request_count = 0
     ipack = 0

     do ipe = 0,num_pe-1
        if ( ipe /= me .and. nrequested(ipe) > 0 ) then   ! skip self or if nothing requested
           istart(ipe) = i1
           request_count = request_count + 1  ! receive counter
           call MPI_IRECV( process_key(i1), nrequested(ipe), MPI_INTEGER8, ipe, tag1, &
                MPI_COMM_WORLD, recv_key_handle(request_count), ierr)
           i1 = i1+nrequested(ipe)
        endif
     end do


     ! Can now ship the keys in the request lists to the PEs where the data sits:
     !   the recipients already know how many requests to expect.

     i1=1
     fetch_pe_count = 0

     do ipe = 0,num_pe-1 

        if ( ipe /= me .and. ntoship(ipe) > 0 ) then   ! avoid shipping to oneself and shipping nothing

           ic_start(ipe) = i1        ! fencepost
           fetch_pe_count = fetch_pe_count + 1  ! receive counter

           ! First initiate receives for returning child info
           call MPI_IRECV( get_child(i1), nplace(ipe), MPI_type_multipole, ipe, tag1, &
                MPI_COMM_WORLD, recv_child_handle(fetch_pe_count), ierr)
           i1 = i1+nplace(ipe)

           ! Extract sub-list of keys to request according to location - don't overwrite buffer!

           ship_key(1:ntoship(ipe),ipe) = pack(request_key(1:nshare), mask = request_owner(1:nshare) == ipe )

           if (emulate_blocking) then
              call MPI_SEND(ship_key(1,ipe), ntoship(ipe), MPI_INTEGER8, ipe, tag1, &
                   MPI_COMM_WORLD,  ierr ) ! Ship to data location
           else
              call MPI_ISEND(ship_key(1,ipe), ntoship(ipe), MPI_INTEGER8, ipe, tag1, &
                   MPI_COMM_WORLD, send_key_handle(fetch_pe_count), ierr ) ! Ship to data location
              call MPI_REQUEST_FREE( send_key_handle(fetch_pe_count), ierr)
           endif

        endif
     end do

     ! Now have complete list of requests from all PEs in rank order. 
     ! Fetch child data from local htable and send back to requesting PE (recipient already knows how many children to expect)


     ! Wait for child request and ship

     ic1 = 1
     send_prop_count = 0


     do iwait = 1,request_count
        call MPI_WAITANY( request_count, recv_key_handle, ihand, status, ierr)  ! Wait for one of receives to complete
        ipe = status(MPI_SOURCE)    ! which PE sent it?
        nreqs = nrequested(ipe)  ! # parent keys
        i1 = istart(ipe)
        i2 = istart(ipe) + nrequested(ipe) - 1
        if (walk_debug) then
           write(ipefile,'(a,i4,a,3i7/(i7))') 'PE ',me,' received request from: ',ipe,nrequested(ipe),istart(ipe),process_key(i1:i2) 
        endif
        process_addr(1:nreqs) = (/( key2addr( process_key(j) ),j=i1,i2)/)    ! get htable addresses
        childbyte(1:nreqs) = htable( process_addr(1:nreqs) )%childcode        !  Children byte-code
        nchild_ship(ipe) = 0


        ! For each key in the request list, fetch and package tree info for children

        do i=1,nreqs
           ipack = i1+i-1   ! Request number: this needs to be unique for each SEND, otherwise buffer pack_child(..,ipack)
           ! may get overwritten before send actually completes
           nchild = SUM( (/ (ibits(childbyte(i),j,1),j=0,2**idim-1) /) )                     ! Get # children
           nchild_ship(ipe) = nchild_ship(ipe) + nchild  ! Total # children to be shipped

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
              pack_child(send_prop_count) = multipole ( key_child(ic), &
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

           end do

           !  Keep record of requested keys
!           requested_keys( nreqs_total(ipe)+1:nreqs_total(ipe)+nchild,ipe ) = key_child(1:nchild)
           nreqs_total(ipe) = nreqs_total(ipe) + nchild  ! Record cumulative total of # children requested 

	end do

        ! Ship child data back to PE that requested it

        if (emulate_blocking) then
           call MPI_SEND( pack_child(ic1), nchild_ship, MPI_type_multipole, ipe, tag1, MPI_COMM_WORLD, ierr )
        else
           call MPI_ISEND( pack_child(ic1), nchild_ship(ipe), MPI_type_multipole, ipe, tag1,&
                MPI_COMM_WORLD, send_child_handle(iwait), ierr )
           call MPI_REQUEST_FREE(send_child_handle(iwait), ierr) 
        endif

        ic1 = ic1 + nchild_ship(ipe)  ! increment start position

        if (walk_debug) write(ipefile,*) 'Total children shipped to processor ',ipe,' from ',me,' was', nchild_ship(ipe) 

     end do

     nchild_ship_tot = ic1-1  ! Total # children shipped to all non-local procs


     ! Wait for data to arrive

     do i=1, fetch_pe_count  ! loop over # PEs originally sent requests

        call MPI_WAITANY( fetch_pe_count, recv_child_handle, ihand, status, ierr)  ! Wait for one of receives to complete
        ipe = status(MPI_SOURCE)   ! which PE?
        ic1 = ic_start(ipe)         ! fenceposts
        ic2 = ic1 + nplace(ipe) -1
        do ic = ic1, ic2
           kchild = get_child(ic)%key
           kparent = ishft( kchild,-idim )
           bchild = get_child(ic)%byte
           lchild = get_child(ic)%leaves
           nxchild = get_child(ic)%next

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
              write(ipefile,'(a,o15,a,i7)') '# leaves <=0 for received child node ',kchild,' from PE ',ipe
           endif

           if (walk_debug) write(ipefile,'(a,o15,a,i7,a,o13)') &
                'Child data arrived:',kchild,' from ',ipe,' requested for key ',kparent

           ! Insert new node into local #-table

           call make_hashentry( kchild, nodchild, lchild, bchild, ipe, hashaddr, ierr )

           htable(hashaddr)%next = nxchild           ! Fill in special next-node pointer for non-local children
           htable( key2addr( kparent) )%childcode = IBSET(  htable( key2addr( kparent) )%childcode, 9) ! Set children_HERE flag for parent node

           node_level( nodchild ) = log(1.*kchild)/log(2.**idim)  ! get level from keys and prestore as node property

           ! Physical properties

           charge( nodchild ) = get_child(ic)%q
           abs_charge( nodchild ) = get_child(ic)%absq
           xcoc( nodchild ) = get_child(ic)%xcoc
           ycoc( nodchild ) = get_child(ic)%ycoc
           zcoc( nodchild ) = get_child(ic)%zcoc
           xdip( nodchild ) = get_child(ic)%xdip
           ydip( nodchild ) = get_child(ic)%ydip
           zdip( nodchild ) = get_child(ic)%zdip
           xxquad( nodchild ) = get_child(ic)%xxquad
           yyquad( nodchild ) = get_child(ic)%yyquad
           zzquad( nodchild ) = get_child(ic)%zzquad
           xyquad( nodchild ) = get_child(ic)%xyquad
           yzquad( nodchild ) = get_child(ic)%yzquad
           zxquad( nodchild ) = get_child(ic)%zxquad

           ! Put last child onto list for post-traversal processing
           if (nxchild == -1) then
              nlast_child = nlast_child + 1
              last_child(nlast_child) = kchild
           endif

           !  Add child key to list of fetched nodes
           nfetch_total(ipe) = nfetch_total(ipe) + 1
 !          fetched_keys( nfetch_total(ipe),ipe ) = kchild

        end do
     end do


     ! Copy defer lists to new walk lists for next tree-walk iteration
     nlist = 0
     do i=1,npshort
        nwalk(i) = ndefer(i)   ! # deferrals still to process
        if (nwalk(i) /= 0) then
           nlist = nlist + 1
           walk_list( 1:nwalk(i), i ) = defer_list( 1:nwalk(i), i )  ! Walk list to inspect
           walk_list( nwalk(i)+1, i ) = 1        ! Last node root for correct 'next-node'
           walk_key(nlist) = walk_list(1, i)       ! Start node for next walk
           plist(nlist) = i                        ! Particle index
           defer_ctr(i) = 1                        ! Deferral counter (1 ... nwalk)
        endif
     enddo

     nactive = count( mask = nwalk(1:npshort) /= 0 )     ! Count remaining 'active' particles - those still with deferred nodes to search

     call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! Wait for other PEs to catch up

     ! Broadcast # remaining particles to other PEs

     call MPI_ALLGATHER( nactive, one, MPI_INTEGER, nactives, one, MPI_INTEGER, MPI_COMM_WORLD, ierr )

     maxactive = maxval(nactives)
     nplace_max = maxval(nplace)
     nchild_ship_tot = maxval(nchild_ship)

     ! cumulative totals

     !     nfetch_total = nfetch_total + nplace


     if (walk_summary ) then
        write (ipefile,'(/a,i8,a2)') 'Summary for traversal # ',ntraversals,' :'
        ! Determine global max
        call MPI_ALLREDUCE( nchild_ship_tot, max_pack, one, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
        call MPI_ALLREDUCE( nplace_max, max_nplace, one, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
        write (ipefile,*) ' # inner loop iterations: ', inner_pass,', sum: ',sum_inner_pass,' previous ave: ',sum_inner_old
        write (ipefile,*) ' # tree hops in inner loop: ',nhops,', sum: ',sum_nhops,' previous: ',sum_nhops_old
        write (ipefile,*) ' # local children shipped: ',nchild_ship,', global max:',max_pack
        write (ipefile,*) ' cumulative total shipped: ',nreqs_total
        write (ipefile,*) ' # non-local children fetched: ',nplace,', global max:',max_nplace
        write (ipefile,*) ' cumulative total fetched: ',nfetch_total

        write (ipefile,*) 'array limit',npackm

        write (ipefile,'(3(/a30,i6)/a/(2i5))') &
             'New twigs: ',newtwig, &
             'New leaves:',newleaf, &
             'New list length: ',nlist, &
             '# remaining active particles on each PE: ',SUM(nactives),MAXVAL(nactives)
!             '# remaining active particles on each PE: ',(i,nactives(i+1),i=0,num_pe-1)
        !        write(ipefile,'(a/(2i5))') 'New shortlist: ',(plist(i),pshort(plist(i)),i=1,nlist)

     endif

     !     CALL VTLEAVE(ICLASSH,VTIERR)

  end do



  !  Determine 'next_node' pointers for 'last child' list & update hash table (tree).
  !  This  ensures that traversals in next pass treat already-fetched nodes as local,
  !  avoiding deferral list completely.


  do i = 1,nlast_child
     search_key = last_child(i)                   
     node_addr = key2addr(search_key)
     htable( node_addr )%next = next_node(search_key)  !   Get next sibling, uncle, great-uncle in local tree
  end do

  sum_inner_old = sum_inner_pass
  call MPI_ALLREDUCE( sum_nhops, sum_nhops_old, one, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( sum_inner_pass, sum_inner_old, one, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  sum_nhops_old = sum_nhops_old/num_pe
  sum_inner_old = sum_inner_old/num_pe

end subroutine tree_walk
