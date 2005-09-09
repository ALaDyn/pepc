! ===========================================
!
!           TREE_WALK
!
!   $ Revision: 1.19 $
!
!   Perform tree walk for all local particles
!
!  Fully collective version: multipole info exchanged en masse for each pass (particle block)
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

subroutine tree_walkc(pshort,npshort, pass,theta,itime,mac,twalk,tfetch)

  use treevars
  use tree_utils
  use utils
  implicit none
  include 'mpif.h'

  real, intent(in) :: theta
  integer, intent(in) :: npshort,itime
  integer, intent(in) :: pshort(npshort)
  integer, intent(in) :: mac
  ! integer, intent(out) :: nodelist(nintm, npshort)
  integer :: npackm   ! Max # children shipped
  integer :: nchild_shipm
  real :: twalk, tfetch, tw1, tw2, tc1, tf1, tf2

  ! Key arrays (64-bit)

  integer*8,  dimension(npshort) :: walk_key, walk_last 
  integer*8, dimension(maxaddress)  :: fetch_key, work_key, request_key
  integer*8, dimension(8) :: sub_key, key_child, child_sub, next_child
  integer*8, dimension(nintmax,npshort) :: defer_list, walk_list
  integer*8, dimension(maxaddress) :: last_child   ! List of 'last' children fetched from remote PEs

  integer, dimension(npshort) :: plist

  integer ::  walk_addr, walk_node, entry_next 
  integer, dimension(npshort) :: nlocal, ndefer,  nwalk, defer_ctr


  integer, dimension(maxaddress) ::  request_owner, fetch_owner, sort_owner, indx
  integer, dimension(8) :: addr_child, node_child, leaves_child, byte_child

  real, dimension(8) :: xcoc_child, ycoc_child, zcoc_child

  real, dimension(0:nlev) :: boxlength2
  logical, dimension(npshort) :: finished, requested
  integer :: hops(21) ! array to control max # iterations in single traversal 
  integer, dimension(num_pe) ::   nactives
  integer, dimension(0:num_pe-1) :: nfetches, &              ! # remote keys to fetch
       nrequests, &           ! # keys requested from elsewhere
       sstrides,  rstrides,  &  ! fence posts for fetch/request lists
       istart, ic_start, &     ! # fenceposts
       nplace,&                ! # children (new entries) to place in table
       nchild_ship       ! # children shipped to others

  ! Key working vars
  integer*8 :: node_key, add_key, walk_next, kchild, kparent, search_key,  nxchild
  integer*8 :: recv_key, recv_parent, recv_next
  integer :: recv_leaves,recv_byte, recv_count

  integer*8 :: ship_key, ship_parent, ship_next
  integer :: ship_leaves, ship_byte, ship_count, ship_node, addr_parent, ship_address

  integer :: i, j, k, ic, ipe, iwait, inner_pass, nhops          ! loop counters
  integer :: p,pass, ipost
  integer :: nnew, newrequest, nreqs, i1, i2, ioff, ipack
  integer :: nchild, newleaf, newtwig, nactive, maxactive, ntraversals, own
  integer :: ic1, ic2, ihand
  integer ::  nchild_ship_sum, nplace_sum  ! Local # total ships and fetches
  integer, save ::  sum_nhops, sum_inner_pass, sum_nhops_old=0, sum_inner_old=0
  integer :: nfetch_sum, nreqs_sum

  integer ::  bchild, nodchild, lchild, hashaddr, nlast_child, cbyte
  integer :: max_nplace, max_pack
  real :: sbox, theta2, theta2_ion, dx, dy, dz, s2, dist2

  ! stuff for tree-patch after traversals complete
  integer ::  node_addr, parent_addr, parent_node, child_byte
  integer :: jmatch(1)
  logical :: resolved, keymatch(8), emulate_blocking=.false.
  logical :: ignore, mac_ok

  integer :: nrest, ndef
  integer :: ierr, nbuf, status(MPI_STATUS_SIZE)
  integer :: tag1=40, iofile

  integer :: key2addr        ! Mapping function to get hash table address from key
  integer :: key2addr_db        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk


  !
  twalk=0.
  tfetch=0.

  npackm = maxaddress
  nchild_shipm = maxaddress
!  walk_debug = .true.
  ! ipefile = 6
  if (walk_debug .or. walk_summary) write(ipefile,'(/2(a,i6))') '*** TREE WALK for timestep ',itime,' pass ',pass
  if (me.eq.0 .and. walk_summary) write(*,'(2(a,i6))') 'LPEPC | TREE WALK for timestep ',itime,' pass ',pass

  sbox = boxsize


  theta2 = theta**2               ! Clumping parameter**2 for MAC
  !  theta2_ion = min(1.0,2*theta2)  ! Ion MAC 50% larger than electron MAC
  theta2_ion=theta2
  boxlength2(0)=sbox**2
  do i=1,nlev
     boxlength2(i) =  boxlength2(i-1)/4.  ! Preprocessed box sizes for each level
  end do

  walk_key(1:npshort) = 1                    ! initial walk list starts at root
  nwalk(1:npshort) = 0   ! # keys on deferred list
  defer_ctr(1:npshort) = 1   ! Deferral counter
  nlist = npshort                 ! Inner loop list length
  plist(1:npshort) = (/ (i,i=1,npshort) /)       ! initial local particle indices 
  nterm = 0
  nactive = npshort

  !  Find global max active particles - necessary if some PEs enter walk on dummy pass

!  call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! Wait for other PEs to catch up

  call MPI_ALLGATHER( nactive, 1, MPI_INTEGER, nactives, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  maxactive = maxval(nactives)

  ntraversals = 0
  nlast_child = 0
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

     !POMP$ INST BEGIN(walk_local)

     call cputime(tw1)
     ntraversals = ntraversals + 1  ! # Tree-walks
     if (walk_debug) write(ipefile,'(/a,i6)') 'WALK-CO: Start of traversal ',ntraversals
     !     if (walk_debug) write(*,'(a,i6)') 'Start of traversal ',ntraversals

     finished(1:nlist) = .false.  ! label all particles not finished
     ndefer(1:npshort) = 0          ! Zero # deferrals (absolute particle number)

     nfetch_sum = 0                   ! Total number of requests for nonlocal keys
     nplace = 0                   ! Total # non-local children to be fetched
     newleaf = 0
     newtwig = 0                  ! local bookkeeping
     inner_pass = 0
     nhops = 0


     do while (nlist>0 .and. nhops <= hops(ntraversals) )        ! Inner loop - single hop in tree walk
        inner_pass = inner_pass+1        ! statistics
        nhops = nhops + nlist
        sum_nhops = sum_nhops + nlist
        sum_inner_pass = sum_inner_pass + 1

        do i=1,nlist
           p = plist(i)  ! particle index
           walk_addr =  key2addr_db( walk_key(i),'WALK: local ' )     ! get htable address
           walk_node = htable( walk_addr )%node             ! Walk node index - points to multipole moments  
           walk_next = htable( walk_addr )%next             ! Next node pointer 

           cbyte = htable( walk_addr )%childcode        ! Children byte-code


           ! children of local/non-local parents already fetched: HERE flag


           dx = x( pshort(p) ) - xcoc( walk_node )      ! Separations
           dy = y( pshort(p) ) - ycoc( walk_node )
           dz = z( pshort(p) ) - zcoc( walk_node )

           s2 = boxlength2( node_level(walk_node) )
           dist2 = theta2*(dx**2+dy**2+dz**2)
           mac_ok = ( s2 < dist2 )             ! Preprocess MAC

           ! set ignore flag if leaf node corresponds to particle itself (number in pshort)
           ignore =  ( pshort(p) == htable( walk_addr )%node )

           ! Wakefield QSA mac condition: prevent forward transmission of pw info
           if (mac==5) ignore = (ignore .or. dx<0) 

           add_key = walk_key(i)                                ! Remember current key

           ! Possible courses of action:


           ! 1) MAC test OK, so put cell on interaction list and find next node for tree walk

           if ( mac_ok .or. (walk_node >0 .and. .not.ignore ) ) then
              walk_key(i) = walk_next
	      entry_next = nterm(p) + 1
              intlist( entry_next, p ) = add_key      ! Augment interaction list - only need keys for diagnosis
              nodelist( entry_next, p ) = walk_node   ! Node number for sum_force
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
                 nfetch_sum = nfetch_sum + 1

                 fetch_key(nfetch_sum) = add_key       ! New request key
                 own = htable( walk_addr )%owner           ! Owner 
                 fetch_owner(nfetch_sum) = own
                 nplace(own) = nplace(own) + n_children(walk_node)     ! Total # children to be fetched from remote PEs
                 ! and to be inserted into local tree
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

     call cputime(tw2)
     twalk=twalk+tw2-tw1
     !POMP$ INST END(walk_local)


     if (walk_debug) then
        write(ipefile,'(a/(o15,i7))') 'Shared defer list to be fetched: ',(fetch_key(i),htable( key2addr( fetch_key(i) ) )%owner,i=1,nfetch_sum)
     endif


     ! At this point all particles have completed their walks with the locally available tree data.
     ! Each PE has list of nonlocal nodes which have been requested during these walks - fetch_key(1:nfetch_sum)
     ! Need to ship these requests to the PEs on which the child-data is actually held (cf particle swapping).

     ! First find out how many requests are to be sent to each PE 
 
! Create dummy keys to prevent zero length buffer in all2allv
   do i=1,num_pe
     nfetch_sum=nfetch_sum+1 
     fetch_key(nfetch_sum)=0
     fetch_owner(nfetch_sum) = i-1
     nplace(i-1) = nplace(i-1)+1 
  enddo

     nfetches(0:num_pe-1) = (/ (count( mask = fetch_owner(1:nfetch_sum) == ipe ), ipe=0,num_pe-1) /)

     !POMP$ INST BEGIN(exchange)

     call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! Wait for other PEs to catch up


     ! First exchange numbers of keys to be shipped and requested

     call MPI_ALLTOALL( nfetches, 1, MPI_INTEGER, nrequests, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )


     ! Derive strides needed for all2all

     sstrides = (/ 0, nfetches(0), ( SUM( nfetches(0:i-1) ),i=2,num_pe-1 ) /)
     rstrides = (/ 0, nrequests(0), ( SUM( nrequests(0:i-1) ),i=2,num_pe-1 ) /)




     ! Sort fetch_keys according to owner
     if (nfetch_sum>0) call indexsort(fetch_owner, indx, nfetch_sum, maxaddress)

     do i=1,nfetch_sum
        work_key(i)=fetch_key(indx(i))
        sort_owner(i)=fetch_owner(indx(i))
     end do

     if (walk_debug) then
        write (ipefile,*) '# keys to fetch: ',nfetches(0:num_pe-1)
        write (ipefile,*) ' fetch strides: ',sstrides(0:num_pe-1)
!        write(ipefile,'(a/(i5,i20))') 'unsorted: ',(fetch_owner(i),fetch_key(i),i=1,nfetch_sum) 
	write(ipefile,'(a/(i5,o15))') 'sorted: ',(sort_owner(i),work_key(i),i=1,nfetch_sum) 
     endif

     ! Exchange key lists

     call MPI_ALLTOALLV( work_key,   nfetches, sstrides, MPI_INTEGER8, &
          request_key, nrequests, rstrides, MPI_INTEGER8, &
          MPI_COMM_WORLD, ierr)

     nreqs_sum = SUM(nrequests)

     ! Recover ownership from # requests - should be neater way of doing this
     ipe=0
     ipost=nrequests(0)
     do i=1,nreqs_sum
        if (i>ipost) then
           ipe=ipe+1
!           if (ipe==me) ipe=ipe+1  ! Skip self
           ipost=ipost+nrequests(ipe)
        endif
        request_owner(i)=ipe
     end do

     if (walk_debug) then

        write (ipefile,*) ' # requested: ',nrequests(0:num_pe-1)
        write (ipefile,*) ' req strides: ',rstrides(0:num_pe-1)
        write(ipefile,'(a/(2i5,o15))') 'requests ',(i,request_owner(i),request_key(i),i=1,nreqs_sum) 
     endif

     ! Now have complete list of requests from all PEs in rank order. 
     ! -- ready for all-to-all multipole swap


!call cleanup
     ! Prepare send buffer: pack multipole info together.

     ship_count=0
     nchild_ship = 0              ! Total # local children to be shipped

     do i=1,nreqs_sum

        if (request_key(i)==0) then
! Setup dummy node
	  nchild=1
	  key_child(1)=0_8
	  next_child(1)=0_8
	  leaves_child(1)=1
	  byte_child(1)=1
	  node_child(1)=1  ! Root node
	  ipe=request_owner(i)  ! Owner 	
          nchild_ship(ipe) = nchild_ship(ipe) + nchild  ! Total # children to be shipped back to this PE
	else

        ! For each key in the request list, fetch and package tree info for children
        ship_address = key2addr_db(request_key(i),'WALK-COL: preship')  ! # address
        ship_node = htable(ship_address)%node
        ship_byte = htable( ship_address )%childcode   ! child byte code
        ship_leaves = htable( ship_address )%leaves                    ! # contained leaves
        ipe = request_owner(i)  ! who requested key

        nchild = SUM( (/ (ibits(ship_byte,j,1),j=0,7) /) )            ! Get # children from 1st 8 bits
        nchild_ship(ipe) = nchild_ship(ipe) + nchild  ! Total # children to be shipped back to this PE

        sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(ship_byte,j),j=0,7) /) )      ! Extract sub key from byte code
        key_child(1:nchild) = IOR( ishft( request_key(i),3 ), sub_key(1:nchild) ) ! Construct keys of children
        addr_child(1:nchild) = (/( key2addr_db( key_child(j),'WALK-COL: ship child ' ),j=1,nchild)/)                 ! Table address of children
        node_child(1:nchild) = htable( addr_child(1:nchild) )%node                        ! Child node index  
        byte_child(1:nchild) = IAND( htable( addr_child(1:nchild) )%childcode,255 )        ! Catch lowest 8 bits of childbyte - filter off requested and here flags 
        leaves_child(1:nchild) = htable( addr_child(1:nchild) )%leaves                    ! # contained leaves
        next_child(1:nchild-1) = htable( addr_child(1:nchild-1) )%next                    ! # next-node pointer
        next_child(nchild) = -1                             ! Last child gets pointed back to _parent_ for non-local nodes
        ! This is used to distinguish particles' walks during 'defer' phase
       endif

        ! Package children properties into user-defined multipole array for shipping

        do ic = 1,nchild
           ship_count = ship_count+1
           pack_child(ship_count) = multipole ( key_child(ic), &
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
                zxquad( node_child(ic)), &
                jx( node_child(ic)), &
                jy( node_child(ic)), &
                jz( node_child(ic)), &
                magmx( node_child(ic)), &
                magmy( node_child(ic)), &
                magmz( node_child(ic)) )

        end do

        !  Keep record of requested child keys - exclude dummy
	if (key_child(1)<>0) then
          nreqs_total(ipe) = nreqs_total(ipe) + nchild  ! Cummulative total of # children shipped
	  requested_keys(sum_ships+1:sum_ships+nchild) = key_child(1:nchild)  ! Keep record of multipole ships
	  requested_owner(sum_ships+1:sum_ships+nchild) = ipe  ! where requests came from
    	  sum_ships = sum_ships + nchild
	endif

        if (walk_debug) then
           write(ipefile,'(a,i4,i7/(o12,i15))') 'Keys requested from: ',ipe,nchild,(key_child(j),next_child(j),j=1,nchild) 
        endif

     end do

     ! Determine strides for # child nodes
     ! First exchange numbers of children to be shipped back


     ! Derive strides needed for all2all
     ! nplace(ipe) = # new nodes to be received
     ! nchild_ship(ipe) = Total # children shipped back
 
     sstrides = (/ 0, nchild_ship(0), ( SUM( nchild_ship(0:i-1) ),i=2,num_pe-1 ) /)
     rstrides = (/ 0, nplace(0), ( SUM( nplace(0:i-1) ),i=2,num_pe-1 ) /)

     if (walk_debug) then
        write (ipefile,*) '# children to ship back: ',nchild_ship(0:num_pe-1)
        write (ipefile,*) ' ship strides: ',sstrides(0:num_pe-1)
        write (ipefile,*) '# children to insert:  ',nplace(0:num_pe-1)
     endif

     nplace_sum = SUM(nplace)
     nchild_ship_sum = SUM(nchild_ship)

     ! Ship multipole data of children
! Problem here if nplace_sum=0 - send dummy to  self

!     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
     call MPI_ALLTOALLV( pack_child,   nchild_ship, sstrides, MPI_TYPE_MULTIPOLE, &
          get_child, nplace, rstrides, MPI_TYPE_MULTIPOLE, &
          MPI_COMM_WORLD, ierr)

! if (ntraversals==3) call cleanup
     ! Make new hash entries with newly-fetched nodes
     newleaf = 0
     newtwig = 0

     ! Recover ownership from # requests - should be neater way of doing this
     ipe=0
     ipost=nplace(0)
     do i=1,nplace_sum
        if (i>ipost) then
           ipe=ipe+1
 !          if (ipe==me) ipe=ipe+1  ! Skip self
           ipost=ipost+nplace(ipe)
        endif
        fetch_owner(i)=ipe
     end do

     do i=1, nplace_sum
        ipe = fetch_owner(i)  
        recv_key = get_child(i)%key
	if (recv_key<>0) then   ! Skip dummies
        recv_parent = ishft( recv_key,-3 )
        recv_byte = get_child(i)%byte
        recv_leaves = get_child(i)%leaves
        recv_next = get_child(i)%next

        if (recv_leaves ==1 ) then
           newleaf = newleaf + 1
           nleaf = nleaf + 1
           nodchild = nleaf
           n_children(nodchild) = 0
           first_child(nodchild) = recv_key

        else if (recv_leaves > 1) then
           newtwig = newtwig + 1
           ntwig = ntwig + 1
           nodchild =  -ntwig
           nchild = SUM( (/ (ibits(recv_byte,j,1),j=0,7) /) )   ! Get # children
           n_children( nodchild ) = nchild       
           sub_key(1:nchild) = pack( bitarr(0:7), mask=(/ (btest(recv_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
           first_child( nodchild ) = IOR( ishft( recv_key,3), sub_key(1) )              ! Construct key of 1st (grand)child

        else
           write(*,'(a,i5,a,i5)') 'LPEPC | WALK-CO on ',me,' Bad twig received from PE ',ipe
           write(*,'(a,o25)') '# leaves <=0 for received child node ',recv_key
           write(*,'(a,2i10)') 'recv count, i ',recv_count,i
           write(*,'(a,i10)') '# leaves ',recv_leaves
           write(*,*) get_child(recv_count)
           !           write(*,'(a,i5)') 'LPEPC | ... key will be ignored' 
           call cleanup  ! Abort
        endif

        ! Insert new node into local #-table

        call make_hashentry( recv_key, nodchild, recv_leaves, recv_byte, ipe, hashaddr, ierr )

        htable(hashaddr)%next = recv_next           ! Fill in special next-node pointer for non-local children
        node_addr =  key2addr_db( recv_parent,'WALK-CO: MNHE ')
        htable( node_addr )%childcode = IBSET(  htable( node_addr )%childcode, 9) ! Set children_HERE flag for parent node

        node_level( nodchild ) = log(1.*recv_key)/log(8.)  ! get level from keys and prestore as node property

        ! Store physical properties

        charge( nodchild ) = get_child(i)%q
        abs_charge( nodchild ) = get_child(i)%absq
        xcoc( nodchild ) = get_child(i)%xcoc
        ycoc( nodchild ) = get_child(i)%ycoc
        zcoc( nodchild ) = get_child(i)%zcoc
        xdip( nodchild ) = get_child(i)%xdip
        ydip( nodchild ) = get_child(i)%ydip
        zdip( nodchild ) = get_child(i)%zdip

        xxquad( nodchild ) = get_child(i)%xxquad
        yyquad( nodchild ) = get_child(i)%yyquad
        zzquad( nodchild ) = get_child(i)%zzquad
        xyquad( nodchild ) = get_child(i)%xyquad
        yzquad( nodchild ) = get_child(i)%yzquad
        zxquad( nodchild ) = get_child(i)%zxquad

        jx( nodchild ) = get_child(i)%jx
        jy( nodchild ) = get_child(i)%jy
        jz( nodchild ) = get_child(i)%jz

        magmx( nodchild ) = get_child(i)%magmx
        magmy( nodchild ) = get_child(i)%magmy
        magmz( nodchild ) = get_child(i)%magmz

        ! Put last child onto list for post-traversal processing
        if (recv_next == -1) then
           nlast_child = nlast_child + 1
           last_child(nlast_child) = recv_key
        endif

        ! Bookkeeping
	nfetch_total(ipe)=nfetch_total(ipe)+1
	sum_fetches=sum_fetches+1  ! Total fetches for iteration
        fetched_keys(sum_fetches) = recv_key
        if (walk_debug) write(ipefile,'(a6,i6,o15,a6,i7,a9,o15,a6,i15)') &
             'Walk: ',sum_fetches,recv_key,' from ',ipe,' parent ',recv_parent,' next ',recv_next
	endif
     end do

! Correct sums for dummies
    nplace_sum=nplace_sum-num_pe
    nchild_ship_sum=nchild_ship_sum-num_pe

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


     ! Broadcast # remaining particles to other PEs

     call MPI_ALLGATHER( nactive, 1, MPI_INTEGER, nactives, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

     maxactive = maxval(nactives)

!     sum_fetches = sum_fetches + nplace_sum  ! Total # fetches/iteration
!     sum_ships = sum_ships + nchild_ship_sum ! Total # shipments/iteration

     if (walk_summary ) then
        write (ipefile,'(/a,i8,a2)') 'LPEPC | Summary for traversal # ',ntraversals,' :'
        ! Determine global max
        call MPI_ALLREDUCE( nplace_sum, max_nplace, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
        write (ipefile,'(a40,i8,a7,i8,a12,i8)') ' # inner loop iterations: ', inner_pass,', sum: ',sum_inner_pass,' previous ave: ',sum_inner_old
        write (ipefile,'(a40,i8,a7,i8,a12,i8)') ' # tree hops in inner loop: ',nhops,', sum: ',sum_nhops,' previous: ',sum_nhops_old
        write (ipefile,*) ' # local children shipped:     ',nchild_ship,', traversal sum:',nchild_ship_sum
        write (ipefile,*) ' # non-local children fetched: ',nplace,', traversal sum:',nplace_sum
        write (ipefile,*) ' cumulative # requested keys:  ',nreqs_total
        write (ipefile,*) ' cumulative # fetched keys:    ',nfetch_total

        write (ipefile,*) 'array limit',size_fetch

        write (ipefile,'(3(/a30,i6)/a/(2i5))') &
             'New twigs: ',newtwig, &
             'New leaves:',newleaf, &
             'New list length: ',nlist, &
             '# remaining active particles on each PE: ',SUM(nactives),MAXVAL(nactives)
        !           '# remaining active particles on each CPU: ',(i,nactives(i+1),i=0,num_pe-1)
        !        write(ipefile,'(a/(2i5))') 'New shortlist: ',(plist(i),pshort(plist(i)),i=1,nlist)

     endif

     ! Array bound checks
     if (nleaf>.95*maxaddress/2) then
        write (6,*) 'LPEPC | WARNING: tree arrays >90% full on CPU ',me
        write (6,*) 'LPEPC | nleaf = ',nleaf,' / ',maxaddress/2
        call cleanup
     endif
     if (ntwig>.95*maxaddress/2) then
        write (6,*) 'LPEPC | WARNING: tree arrays >90% full on CPU ',me
        write (6,*) 'LPEPC | ntwig = ',ntwig,' / ',maxaddress/2
        call cleanup
     endif

     call cputime(tc1)
     tfetch=tfetch+tc1-tw2  ! timing for 2nd half of walk
     !POMP$ INST END(exchange)
  end do



  !  Determine 'next_node' pointers for 'last child' list & update hash table (tree).
  !  This  ensures that traversals in next pass treat already-fetched nodes as local,
  !  avoiding deferral list completely.
!   write(*,*) 'nlast_child = ',nlast_child

  do i = 1,nlast_child
     search_key = last_child(i)                   
     node_addr = key2addr_db(search_key,'WALK: NN search ')
     htable( node_addr )%next = next_node(search_key)  !   Get next sibling, uncle, great-uncle in local tree
!  if (walk_debug) then
!	write(ipefile,'(a,o15,a4,o15)') 'Changed next node for ',search_key,' to ',htable(node_addr)%next
!  endif
  end do

  sum_inner_old = sum_inner_pass
!  call MPI_ALLREDUCE( sum_nhops, sum_nhops_old, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
!  call MPI_ALLREDUCE( sum_inner_pass, sum_inner_old, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  sum_nhops_old = sum_nhops_old/num_pe
  sum_inner_old = sum_inner_old/num_pe
  maxtraverse = max(maxtraverse,ntraversals)

end subroutine tree_walkc
