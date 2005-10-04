! ===========================================
!
!           TREE_UPDATE
!
!   $Revision$
!
!  Exchange non-local multipole information used at last timestep
!  as 'prefetch' for tree walk in freeze mode.
!
!
! ===========================================

subroutine tree_update(itime)

  use treevars
  use tree_utils
  implicit none
  include 'mpif.h'

  integer, intent(in) :: itime
  integer*8, dimension(8) :: sub_key, key_child, child_sub, child_key, next_child, siblings
  
  integer, dimension(0:num_pe-1) ::  & 
       sstrides, &              ! fence posts for remove lists
       rstrides, &              ! fence posts for remove lists
       istart, ic_start, &     ! # fenceposts
       nplace,&                ! # children (new entries) to place in table
       nchild_ship       ! # children shipped to others
  integer :: ierr
  integer*8, dimension(size_fetch) :: work_key, last_child   ! List of 'last' children fetched from remote PEs
  integer*8, dimension(size_fetch) :: sort_reqs, sort_fetch
  integer, dimension(size_fetch) :: indxr, indxf, sort_reqowner, sort_fetchowner
  logical, dimension(size_fetch) :: key_present, duplicate

  ! Key working vars
  integer :: ship_byte, ship_leaves, ship_node, ship_address
  integer*8 :: ship_next, ship_key
  integer :: recv_byte, recv_leaves, recv_node, recv_address
  integer*8 :: recv_next, recv_key, recv_parent
  integer*8 :: kparent, kin, search_key
  integer :: nodchild, nlast_child, newtwig, newleaf, hashaddr, nuniq, npar, nabsent
  integer :: cchild, nchild, node_addr, addr_parent, child_byte
  integer :: i, j, k, ic, ipe, iwait, inner_pass, nhops, nreqs_new, nreqs_old, nfetch_new         ! loop counters
  integer :: iofile, sum_reqs, sum_fetch
  integer :: size_remove, nreq_max, nfetch_max,  timestamp, send_prop_count, recv_count, nnot_local
  character*1 :: ctick
  character(30) :: cfile, ccol1, ccol2, ccol0
  integer, save :: sumfetches

  ! external functions
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer :: key2addr_db        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk
  logical :: key_local   ! Tests whether key present in local # table
  logical :: update_debug=.false.

  iofile = ipefile
!  iofile = 6
  !
  if (update_debug) write(iofile,'(/a,i6)') 'TREE NONLOCAL UPDATE for timestep ',itime

! TODO prepare buffers as in walk: requested_keys 1D; needs sorting first

! nfetch_total(ipe) contains the # key fetches
! nreqs_total(ipe) contains the # key requests



     call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! Wait for other PEs to catch up


     ! Sort fetched and requested_keys according to owner

     if (sum_ships>0) call indexsort(requested_owner, indxr, sum_ships, maxaddress)
     if (sum_fetches>0) call indexsort(fetched_owner, indxf, sum_fetches, maxaddress)

     do i=1,sum_ships
        sort_reqs(i)=requested_keys(indxr(i))
        sort_reqowner(i)=requested_owner(indxr(i))
     end do

     do i=1,sum_fetches
        sort_fetch(i)=fetched_keys(indxf(i))
        sort_fetchowner(i)=fetched_owner(indxf(i))
     end do

     ! Derive strides needed for all2all

     sstrides = (/ 0, nfetch_total(0), ( SUM( nfetch_total(0:i-1) ),i=2,num_pe-1 ) /)
     rstrides = (/ 0, nreqs_total(0), ( SUM( nreqs_total(0:i-1) ),i=2,num_pe-1 ) /)

     if (update_debug) then
        write (ipefile,*) '# keys to ship: ',sum_ships,nreqs_total(0:num_pe-1)
        write (ipefile,*) ' ship strides: ',rstrides(0:num_pe-1)
!        write(ipefile,'(a/(i5,i20))') 'unsorted: ',(requested_owner(i),requested_keys(i),i=1,sum_ships)
        write(ipefile,'(a/(i5,o15))') 'sorted: ',(sort_reqowner(i),sort_reqs(i),i=1,sum_ships)
        write (ipefile,*) ' # fetches: ',sum_fetches,nfetch_total(0:num_pe-1)
        write (ipefile,*) ' fetch strides: ',sstrides(0:num_pe-1)
        write(ipefile,'(a/(2i5,o15))') 'fetches ',(i,sort_fetchowner(i),sort_fetch(i),i=1,sum_fetches)
     endif

!if (itime==4) call cleanup
     ! Now have complete list of requests from all PEs in rank order.
     ! -- ready for all-to-all multipole swap


  ! Prepare send buffer with sorted shipping list: pack multipole info together as in tree_walk.

  do i=1,sum_ships
        ship_key = sort_reqs(i)
        ship_address = key2addr_db(ship_key,'UPDATE: pack1 ')  ! # address
        ship_node = htable(ship_address)%node
        ship_byte = IAND( htable( ship_address )%childcode,255 ) ! Catch lowest 8 bits of childbyte - filter off requested and here flags 
        ship_leaves = htable( ship_address )%leaves                    ! # contained leaves

        !  Need to reset ship_next=-1 if node is last of children
        kparent = ishft(ship_key ,-3 )
        addr_parent = key2addr_db(kparent,'UPDATE: pack2 ')
        child_byte = htable( addr_parent )%childcode            !  Children byte-code
        nchild = SUM( (/ (ibits(child_byte,j,1),j=0,7) /) )       ! # children = sum of bits in byte-code
        child_sub(1:nchild) = pack( bitarr, mask=(/ (btest(child_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
        child_key(1:nchild) = IOR( ishft(kparent,3), child_sub(1:nchild) )         ! Construct keys of children

        if (ship_key == child_key(nchild) ) then   
           ship_next = -1      ! last of siblings or only child
        else
           ship_next = htable( ship_address )%next   ! Node has elder siblings
        endif

        pack_child(i) =  multipole ( ship_key, &
             ship_byte, &
             ship_leaves, &
	     me, &
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


  ! Ship multipole data
  call MPI_ALLTOALLV( pack_child,  nreqs_total, rstrides, MPI_TYPE_MULTIPOLE, &
       get_child, nfetch_total, sstrides, MPI_TYPE_MULTIPOLE, &
       MPI_COMM_WORLD, ierr)

  ! Update hash entries with refreshed multipole data 

     do i=1, sum_fetches
	ipe=sort_fetchowner(i)
        recv_key = get_child(i)%key
        recv_parent = ishft( recv_key,-3 )
        recv_byte = get_child(i)%byte
        recv_leaves = get_child(i)%leaves
        recv_next = get_child(i)%next

        ! Check new node exists in local #-table

!        call make_hashentry( recv_key, nodchild, recv_leaves, recv_byte, ipe, hashaddr, ierr )
!	key_present = key_local(recv_key)

     	node_addr = key2addr_db(recv_key,'UPDATE: MNHE ')
	nodchild = htable(node_addr)%node

        ! Physical properties

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


        if (update_debug) write(iofile,'(a,o15,a,i7,a,o13)') &
             'Prefetch: ',recv_key,' from ',ipe,' parent key ',recv_parent

     end do

end subroutine tree_update

