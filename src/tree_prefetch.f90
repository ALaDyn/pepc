! ===========================================
!
!           TREE_PREFETCH
!
!   $Revision$
!
!  Exchange non-local multipole information used at last timestep
!  as 'prefetch' for tree walk.
!
!
! ===========================================

subroutine tree_prefetch(itime)

  use treevars
  use utils

  implicit none

  integer, intent(in) :: itime
  ! Key arrays (64-bit)

  integer*8, dimension(size_tree/2,0:num_pe-1) :: remove_keys, nofetch_keys ! List of deleted keys
  integer*8, dimension(size_tree) :: req_parent, req_compress, fetch_parent, fetch_comp, absent   ! work arrays

  integer*8, dimension(8) :: sub_key, key_child, child_sub, child_key, next_child, siblings

  integer, dimension(0:num_pe-1) :: ntoship, &              ! # keys needed
       nrequested, &           ! # keys requested from elsewhere
       nremove, &              ! # keys deleted from request list
       nadd, &                 ! # new ones added
       nnofetch, &             ! # keys to delete from fetch list
       sstrides, &              ! fence posts for remove lists
       rstrides, &              ! fence posts for remove lists
       istart, ic_start, &     ! # fenceposts
       nplace,&                ! # children (new entries) to place in table
       nchild_ship       ! # children shipped to others

  integer*8, dimension(size_tree) :: last_child   ! List of 'last' children fetched from remote PEs

  logical, dimension(size_tree) :: key_present, duplicate

  ! Key working vars
  integer :: ship_byte, ship_leaves, ship_node, ship_address
  integer*8 :: ship_next, ship_key
  integer :: recv_byte, recv_leaves, recv_node, recv_address
  integer*8 :: recv_next, recv_key, recv_parent
  integer*8 :: kparent, kin, search_key
  integer :: nodchild, nlast_child, newtwig, newleaf, hashaddr, nuniq, npar, nabsent
  integer :: cchild, nchild, node_addr, addr_parent, child_byte
  integer :: i, j, k, ic, ipe, iwait, inner_pass, nhops, nreqs_new, nreqs_old, nfetch_new         ! loop counters
  integer :: size_remove, nreq_max, nfetch_max, timestamp, send_prop_count, recv_count, nnot_local
  character*1 :: ctick
  character(30) :: cfile, ccol1, ccol2, ccol0

  ! external functions
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk
  logical :: key_local   ! Tests whether key present in local # table

  !
  if (prefetch_debug) write(ipefile,'(/a,i6)') 'TREE PREFETCH for timestep ',itime
  ! if (walk_debug) write(*,'(/a,i6)') 'TREE PREFETCH for timestep ',itime
  timestamp=itime



  ! First examine fetched-key lists to check whether parents exist locally
  do ipe = 0,num_pe-1
     if (prefetch_debug) write(ipefile,*) 'PASS 1: Fetches from PE ',ipe
     npar = nfetch_total(ipe)
     fetch_comp(1:npar) = fetched_keys(1:npar,ipe)  ! Work array
     nfetch_new = npar  ! New list length

     do while (npar>0)  ! Loop recursively back up tree until all absent ancestors identified

        if (prefetch_debug) call sort(fetch_comp(1:npar))   ! Sort key list for debugging
        fetch_parent(1:npar) = ishft( fetch_comp(1:npar),-idim)  ! Parent keys of fetch-list
        key_present(1:npar) = (/ (key_local( fetch_parent(i) ),i=1,npar) /) ! Check whether parent node present locally
        nnot_local = count(mask = .not.key_present(1:npar))     ! Count absentees
        if (prefetch_debug) write(ipefile,'(a,i6/(2o15,l3))') &
	'Old fetch list   Parents    Present ',npar,(fetch_comp(i),fetch_parent(i),key_present(i),i=1,npar)

        absent(1:nnot_local) = &
             pack( fetch_parent(1:npar), mask = .not.key_present(1:npar) ) ! Make list of absent parents

        call sort(absent(1:nnot_local))   ! Sort absentee list

        absent(nnot_local+1) = 0 ! dummy at end
        duplicate(1:nnot_local) = (/ (absent(i) /= absent(i+1), i=1,nnot_local) /)  ! Identify unique keys     
        nabsent= count(mask = duplicate(1:nnot_local))          ! Count them
        absent(1:nabsent) = pack(absent(1:nnot_local), mask = duplicate(1:nnot_local))   ! Store compressed list
        fetch_comp(1:nabsent) = absent(1:nabsent) ! buffer for next level up
        npar = nabsent
        fetched_keys(nfetch_new+1:nfetch_new+nabsent,ipe) = absent(1:nabsent)  ! Add absent keys to list
        nfetch_new = nfetch_new + nabsent

     end do


!     if (prefetch_debug) write(ipefile,'(a,i6/(o15))') 'New fetch list: ',nfetch_new,(fetched_keys(i,ipe),i=1,nfetch_new)

     nfetch_total(ipe) = nfetch_new   ! New # keys to fetch

  end do

  !  Get fence posts for key swap

  sstrides = (/ (i*size_tree,i=0,num_pe-1) /)
  rstrides = sstrides

  !  Send new fetch-list totals to destination PEs
  call MPI_ALLTOALL( nfetch_total, 1, MPI_INTEGER, nreqs_total, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)

  !  Do key-swap to update request lists
  call MPI_ALLTOALLV( fetched_keys,   nfetch_total, sstrides, MPI_INTEGER8, &
       requested_keys, nreqs_total, rstrides, MPI_INTEGER8, &
       MPI_COMM_WORLD, ierr)



  ! Now compare request lists with actual local tree nodes:
  !  - flag and remove keys which no longer exist

  do ipe = 0,num_pe-1
     nremove(ipe) = 0
     if (prefetch_debug) write(ipefile,*) 'PASS 2: Requests from PE ',ipe
     do i=1,nreqs_total(ipe)
        key_present(i) = key_local(requested_keys(i,ipe))
        if (key_present(i)) then
           ctick = ''
        else
           ! tag key as non existent and add to remove list
           nremove(ipe) = nremove(ipe) + 1
           remove_keys(nremove(ipe),ipe) = requested_keys(i,ipe)
           ctick = 'N'
        endif
        if (prefetch_debug)  write(ipefile,'(o15,2x,a1)') requested_keys(i,ipe),ctick
     end do
     write(ipefile,*)

     ! remove non-existent keys
     ! could do something more sophisticated here, like suggesting alternative key
     ! but for now let tree_walk routine deal with missing keys

     nreqs_new = count(mask = key_present(1:nreqs_total(ipe)))     ! Count them
     req_compress(1:nreqs_new) = &
          pack(requested_keys(1:nreqs_total(ipe),ipe), mask = key_present(1:nreqs_total(ipe)))        ! Compress list
     nreqs_total(ipe) = nreqs_new
     if (prefetch_debug) write(ipefile,'(a,i8/(o15))') 'Keys removed: ', nremove(ipe),(remove_keys(i,ipe),i=1,nremove(ipe))

     ! TODO: Check that sibling sets complete so that parent's 'HERE' flag can be set in childcode.
     ! Need this so that tree_walk MAC options still function as-is

     ! First make list of parent keys from compressed list

     call sort(req_compress(1:nreqs_new))   ! Sort key list
     req_parent(1:nreqs_new) = ishft( req_compress(1:nreqs_new),-idim)

     if (prefetch_debug) write(ipefile,'(a,i6/(2o15))') & 
	'Compressed list   Parents ',nreqs_new,(req_compress(i),req_parent(i),i=1,nreqs_new)
     req_parent(nreqs_new+1) = 0
     duplicate(1:nreqs_new) = (/ (req_parent(i) /= req_parent(i+1), i=1,nreqs_new) /)  ! Identify unique keys     
     nuniq= count(mask = duplicate(1:nreqs_new))            ! Count them
     req_parent(1:nuniq) = pack(req_parent(1:nreqs_new), mask = duplicate(1:nreqs_new))   ! Compress list
     if (prefetch_debug) write(ipefile,'(a,i6/(o15))') 'Unique parents ',nuniq,(req_parent(i),i=1,nuniq)

     ! Make new list of requested children from parent list
     nreqs_old = nreqs_new
     nreqs_new = 0

     do i=1,nuniq
        cchild = htable( key2addr( req_parent(i) ) )%childcode   !  Children byte-code
        nchild = SUM( (/ (ibits(cchild,j,1),j=0,7) /) ) ! # children = sum of bits in byte-code
        sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(cchild,j),j=0,7) /) )  ! Extract child sub-keys from byte code           
        requested_keys(nreqs_new+1:nreqs_new+nchild,ipe) = IOR( ishft(req_parent(i),idim ), sub_key(1:nchild) ) ! New siblings of original requested key
        nreqs_new = nreqs_new + nchild
     end do
     nreqs_total(ipe) = nreqs_new
     nadd(ipe) = nreqs_new-nreqs_old
  end do


  !  Send new request-list totals to requesting PEs
  call MPI_ALLTOALL( nreqs_total, 1, MPI_INTEGER, nfetch_total, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)


  call MPI_ALLTOALLV( requested_keys,   nreqs_total, sstrides, MPI_INTEGER8, &
       fetched_keys, nfetch_total, rstrides, MPI_INTEGER8, &
       MPI_COMM_WORLD, ierr)

  ! Final sweep of fetched-key lists to eliminate nodes which might have been 
  ! created locally during tree_fill.

  do ipe = 0,num_pe-1
     if (ipe /= me) then
        if (prefetch_debug) write(ipefile,*) 'PASS 3: Fetches from PE ',ipe
        npar = nfetch_total(ipe)
        fetch_comp(1:npar) = fetched_keys(1:npar,ipe)  ! Work array
        key_present(1:npar) = (/ (key_local( fetch_comp(i) ),i=1,npar) /) ! Check whether node already present locally
        nnot_local = count(mask = .not.key_present(1:npar))     ! Count absentees
        if (prefetch_debug) write(ipefile,'(a,i6/(o15,l3))') 'Fetch list  Present ',npar,(fetch_comp(i),key_present(i),i=1,npar)
        nfetch_new = nnot_local
        fetched_keys(1:nfetch_new,ipe) = &
             pack( fetch_comp(1:npar), mask = .not.key_present(1:npar) ) ! Make list of absent keys
        nfetch_total(ipe) = nfetch_new   ! New # keys to fetch
     endif
  end do


  !  Send new fetch-list totals to destination PEs
  call MPI_ALLTOALL( nfetch_total, 1, MPI_INTEGER, nreqs_total, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)

  !  Do key-swap to update request lists
  call MPI_ALLTOALLV( fetched_keys,   nfetch_total, sstrides, MPI_INTEGER8, &
       requested_keys, nreqs_total, rstrides, MPI_INTEGER8, &
       MPI_COMM_WORLD, ierr)


  if (walk_debug) then
     nreq_max=maxval(nreqs_total)
     nfetch_max=maxval(nfetch_total)
     ! Insert dummy values
     do ipe=0,num_pe-1
        requested_keys(nreqs_total(ipe)+1:nreq_max,ipe) = 0
        fetched_keys(nfetch_total(ipe)+1:nfetch_max,ipe) = 0
     end do

     ! formatting for fetch/request lists
     if (num_pe<10) then
        ccol0='(//a,i7/3('//achar(mod(num_pe,10)+48)//'i15/),a/)'
        ccol1='(//a,i7/'//achar(mod(num_pe,10)+48)//'i15/a/)'
        ccol2='('//achar(mod(num_pe,10)+48)//'o15)'
     else
        ccol0='(//a,i7/3('//achar(num_pe/10+48)//achar(mod(num_pe,10)+48)//'i15/),a/)'
        ccol1='(//a,i7/'//achar(num_pe/10+48)//achar(mod(num_pe,10)+48)//'i15/a/)'
        ccol2='('//achar(num_pe/10+48)//achar(mod(num_pe,10)+48)//'o15)'
     endif

     write(ipefile,ccol0) 'Keys requested/removed/added at timestep ', &
	timestamp,nreqs_total(0:num_pe-1),nremove(0:num_pe-1),nadd(0:num_pe-1), &
        '------------------------------------------------------------'
     do i=1,nreq_max
        write(ipefile,ccol2) (requested_keys(i,ipe),ipe=0,num_pe-1)
     end do

     write(ipefile,ccol1) 'Keys fetched at timestep ',timestamp,nfetch_total(0:num_pe-1), &
          '------------------------------------------------------------'
     do i=1,nfetch_max
        write(ipefile,ccol2) (fetched_keys(i,ipe),ipe=0,num_pe-1)
     end do
  endif

 !  return
  ! Now ready for all-to-all multipole swap

  ! Prepare send buffer: pack multipole info together as in tree_walk.
  send_prop_count = 0
  nlast_child = 0

  do ipe=0,num_pe-1
     !     sstrides(ipe) = send_prop_count  ! prestore fencepost
     do i=1,nreqs_total(ipe)
        send_prop_count = send_prop_count + 1
        ship_key = requested_keys(i,ipe)
        ship_address = key2addr(ship_key)  ! # address
        ship_node = htable(ship_address)%node
        ship_byte = IAND( htable( ship_address )%childcode,255 ) ! Catch lowest 8 bits of childbyte - filter off requested and here flags 
        ship_leaves = htable( ship_address )%leaves                    ! # contained leaves

        !  Need to reset ship_next=-1 if node is last of children
        kparent = ishft(ship_key ,-idim )
        addr_parent = key2addr(kparent)
        child_byte = htable( addr_parent )%childcode            !  Children byte-code
        nchild = SUM( (/ (ibits(child_byte,j,1),j=0,7) /) )       ! # children = sum of bits in byte-code
        child_sub(1:nchild) = pack( bitarr, mask=(/ (btest(child_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
        child_key(1:nchild) = IOR( ishft(kparent,idim), child_sub(1:nchild) )         ! Construct keys of children

        if (ship_key == child_key(nchild) ) then   
           ship_next = -1      ! last of siblings or only child
        else
           ship_next = htable( ship_address )%next   ! Node has elder siblings
        endif



        pack_child(send_prop_count) =  multipole ( ship_key, &
             ship_byte, &
             ship_leaves, &
             ship_next, &
             charge( ship_node ), &
             abs_charge( ship_node ), &
             xcoc( ship_node), &
             ycoc( ship_node), &
             zcoc( ship_node), &
             xdip( ship_node), &
             ydip( ship_node), &
             zdip( ship_node), &
             xxquad( ship_node), &
             yyquad( ship_node), &
             zzquad( ship_node), &
             xyquad( ship_node), &
             yzquad( ship_node), &
             zxquad( ship_node), &
             jx( ship_node), &
             jy( ship_node), &
             jz( ship_node), &
             magmx( ship_node), &
             magmy( ship_node), &
             magmz( ship_node) )
     end do
  end do

  sstrides = (/ 0, nreqs_total(0), ( SUM( nreqs_total(0:i-1) ),i=2,num_pe-1 ) /)

  rstrides = (/ 0, nfetch_total(0), ( SUM( nfetch_total(0:i-1) ),i=2,num_pe-1 ) /)

  ! write (ipefile,'((2i8))') (sstrides(i), rstrides(i), i=0,num_pe-1) 

  ! Ship multipole data

  call MPI_ALLTOALLV( pack_child,   nreqs_total, sstrides, MPI_TYPE_MULTIPOLE, &
       get_child, nfetch_total, rstrides, MPI_TYPE_MULTIPOLE, &
       MPI_COMM_WORLD, ierr)

  ! Make new hash entries with newly-fetched nodes
  recv_count = 0
  nlast_child = 0
  newleaf = 0
  newtwig = 0

  do ipe=0,num_pe-1
     do i=1, nfetch_total(ipe)
        recv_count=recv_count+1
        recv_key = get_child(recv_count)%key
        recv_parent = ishft( recv_key,-idim )
        recv_byte = get_child(recv_count)%byte
        recv_leaves = get_child(recv_count)%leaves
        recv_next = get_child(recv_count)%next

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
           first_child( nodchild ) = IOR( ishft( recv_key,idim), sub_key(1) )              ! Construct key of 1st (grand)child

        else
           write(ipefile,'(a,o15,a,i7)') '# leaves <=0 for received child node ',recv_key,' from PE ',ipe
        endif


        ! Insert new node into local #-table

        call make_hashentry( recv_key, nodchild, recv_leaves, recv_byte, ipe, hashaddr, ierr )

        htable(hashaddr)%next = recv_next           ! Fill in special next-node pointer for non-local children
        htable( key2addr( recv_parent) )%childcode = IBSET(  htable( key2addr( recv_parent) )%childcode, 9) ! Set children_HERE flag for parent node

        node_level( nodchild ) = log(1.*recv_key)/log(2.**idim)  ! get level from keys and prestore as node property

        ! Physical properties

        charge( nodchild ) = get_child(recv_count)%q
        abs_charge( nodchild ) = get_child(recv_count)%absq
        xcoc( nodchild ) = get_child(recv_count)%xcoc
        ycoc( nodchild ) = get_child(recv_count)%ycoc
        zcoc( nodchild ) = get_child(recv_count)%zcoc
        xdip( nodchild ) = get_child(recv_count)%xdip
        ydip( nodchild ) = get_child(recv_count)%ydip
        zdip( nodchild ) = get_child(recv_count)%zdip

        xxquad( nodchild ) = get_child(recv_count)%xxquad
        yyquad( nodchild ) = get_child(recv_count)%yyquad
        zzquad( nodchild ) = get_child(recv_count)%zzquad
        xyquad( nodchild ) = get_child(recv_count)%xyquad
        yzquad( nodchild ) = get_child(recv_count)%yzquad
        zxquad( nodchild ) = get_child(recv_count)%zxquad

        jx( nodchild ) = get_child(recv_count)%jx
        jy( nodchild ) = get_child(recv_count)%jy
        jz( nodchild ) = get_child(recv_count)%jz

        magmx( nodchild ) = get_child(recv_count)%magmx
        magmy( nodchild ) = get_child(recv_count)%magmy
        magmz( nodchild ) = get_child(recv_count)%magmz

 ! Put last child onto list for post-traversal processing
        if (recv_next == -1) then
           nlast_child = nlast_child + 1
           last_child(nlast_child) = recv_key
        endif

        if (prefetch_debug) write(ipefile,'(a,o15,a,i7,a,o13)') &
             'Prefetch: ',recv_key,' from ',ipe,' parent key ',recv_parent

     end do
  end do

  !  Determine 'next_node' pointers for 'last child' list in local tree
  ! - these differ from pointers in NL tree from which nodes were sent
  !  This  ensures that traversals in tree_walk treat pre-fetched nodes as local,
  !  avoiding deferral list completely.

  do i = 1,nlast_child
     search_key = last_child(i)                   
     node_addr = key2addr(search_key)
     htable( node_addr )%next = next_node(search_key)  !   Get next sibling, uncle, great-uncle in local tree
  end do


end subroutine tree_prefetch

