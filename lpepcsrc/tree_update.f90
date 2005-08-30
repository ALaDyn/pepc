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
  integer :: ierr
  integer*8, dimension(size_fetch) :: last_child   ! List of 'last' children fetched from remote PEs

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
  integer :: iofile
  integer :: size_remove, nreq_max, nfetch_max,  timestamp, send_prop_count, recv_count, nnot_local
  character*1 :: ctick
  character(30) :: cfile, ccol1, ccol2, ccol0
  integer, save :: sumfetches

  ! external functions
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer :: key2addr_db        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk
  logical :: key_local   ! Tests whether key present in local # table
!  logical :: update_debug=.true.

  iofile = ipefile
!  iofile = 6
  !
  if (prefetch_debug) write(iofile,'(/a,i6)') 'TREE NONLOCAL UPDATE for timestep ',itime



  ! Prepare send buffer: pack multipole info together as in tree_walk.
  send_prop_count = 0
  nlast_child = 0

  do ipe=0,num_pe-1
     !     sstrides(ipe) = send_prop_count  ! prestore fencepost
     do i=1,nreqs_total(ipe)
        send_prop_count = send_prop_count + 1
        ship_key = requested_keys(i,ipe)
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

  ! write (iofile,'((2i8))') (sstrides(i), rstrides(i), i=0,num_pe-1) 

  ! Ship multipole data
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
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
        recv_parent = ishft( recv_key,-3 )
        recv_byte = get_child(recv_count)%byte
        recv_leaves = get_child(recv_count)%leaves
        recv_next = get_child(recv_count)%next

        ! Check new node exists in local #-table

!        call make_hashentry( recv_key, nodchild, recv_leaves, recv_byte, ipe, hashaddr, ierr )
!	key_present = key_local(recv_key)

     	node_addr = key2addr_db(recv_key,'UPDATE: MNHE ')
	nodchild = htable(node_addr)%node

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


        if (prefetch_debug) write(iofile,'(a,o15,a,i7,a,o13)') &
             'Prefetch: ',recv_key,' from ',ipe,' parent key ',recv_parent

     end do
  end do


! Get total # multipole ships from prefetch 
 sumfetches = MAXVAL(nfetch_total)
 call MPI_ALLREDUCE( sumfetches, max_prefetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr ) 

end subroutine tree_update

