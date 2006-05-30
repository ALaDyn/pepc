!
! ===========================================
!
!           TREE_FILL
!
!   Fill in coarsest levels of tree to complete local data structure 

!   
!
! ===========================================

subroutine tree_fill

  use treevars
  use tree_utils

  implicit none
  include 'mpif.h'

!  integer, parameter :: size_t=1000

 integer*8, dimension(maxaddress) :: sub_key, parent_key
 integer*8, dimension(8) :: child_key, child_sub

  integer, dimension(nbranch_max) :: branch_level
  integer, dimension(maxaddress) :: twig_addr, twig_code,  cell_addr, tree_node, parent_addr
  integer, dimension(8) :: child_addr !  children nodes
  logical :: duplicate(maxaddress), resolved, keymatch(8)

  integer*8 :: node_key, search_key, parent,  child_top
  integer :: maxlevel, ilevel, nsub,i,j,k, nparent, nuniq, child_byte, child_bit, nchild, link_addr, hashaddr
  integer :: maxtwig, maxleaf, nleaf_check, ntwig_check
  integer ::  node_addr, jmatch(1),  parent_node, parent_level, nodtwig
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk
  integer :: ierr

  tree_debug=.true.
  if (tree_debug) write(ipefile,'(/a)') 'TREE FILL'
  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | FILL'

  if (tree_debug) call check_table('after make_branches ')

  ! get levels of branch nodes
  maxlevel=0
  do i=1,nbranch_sum
     branch_level(i) = log( 1.*branch_key(i) )/log(8.)
     maxlevel = max( maxlevel, branch_level(i) )        ! Find maximum level
  end do

  nparent = 0
  parent_key(1) = 0
!  if (tree_debug.and.me>250) then
!     write (*,'(/a,i5/(i5,o16,i3))') 'Global branch list on ',me,&
!          (i,branch_key(i), branch_level(i),i=1,nbranch_sum)
!  endif

  do ilevel = maxlevel,1,-1                                            ! Start at finest branch level
     nsub = count( mask=branch_level(1:nbranch_sum) == ilevel )                       ! Count # branches at this level
!     if (tree_debug.and.me>250) then
!        write(*,*) 'proc ',me,' level= ',ilevel,' nsub= ',nsub,' nbranches=',nbranch_sum
!     endif
     sub_key(1:nsub) =  pack(branch_key(1:nbranch_sum), mask = branch_level(1:nbranch_sum) == ilevel)        ! Pick out branches at current level

     sub_key(nsub+1:nsub+nparent) = parent_key(1:nparent)              ! Augment list with parent keys checked at previous level
     nsub = nsub + nparent
     call sort(sub_key(1:nsub))                                        ! Sort keys

     sub_key(nsub+1) = 0
     duplicate(1:nsub) = (/ (sub_key(i) /= sub_key(i+1), i=1,nsub) /)  ! Identify unique keys     
     nuniq= count(mask = duplicate(1:nsub))                            ! Count them
     sub_key(1:nuniq) = pack( sub_key(1:nsub), mask = duplicate(1:nsub) )        ! Compress list


     do i=1,nuniq

        parent_key(i) = ISHFT( sub_key(i),-3 )             ! Parent key

        child_bit = IAND( sub_key(i), hashchild)                    ! extract child bit from key: which child is it?
        child_byte = 0
        child_byte = IBSET(child_byte, child_bit)    ! convert child bit to byte code
        nodtwig = -ntwig -1     ! predicted twig # 

        call make_hashentry( parent_key(i), nodtwig, 0, child_byte, me, hashaddr, ierr )

        if ( ierr == 1 ) then     
           ! keys match, so node already exists locally
           ! Set child-bit in existing parent byte-code
           hashaddr = key2addr( parent_key(i),'FILL: sweep1'  )
           htable( hashaddr )%childcode = IBSET( htable( hashaddr )%childcode, child_bit )

        else if (ierr == 0 ) then
           ntwig = ntwig + 1
           ntwig_me = ntwig_me+1               ! # local twigs
	   twig_key(ntwig_me) = htable( hashaddr)%key  ! add to list of local twigs

        else
           write (ipefile,*) 'Key number ',i,' not resolved'

        endif

     end do
     nparent = nuniq
  end do


  if (tree_debug) call check_table('End of tree-fill    ')

  !  Go through twig nodes and fix # leaves in #table to include non-local branch nodes
  !  - put this stuff in tree_properties for correct array sizes

  nnodes = ntwig + nleaf

  nnodes_pw = nnodes   ! Store # nodes prior to tree-walk
  ntwig_pw = ntwig     ! - need for re-calculating moments using same tree-structure
  nleaf_pw = nleaf



  do i=1,ntwig_me
    hashaddr = key2addr( twig_key(i),'FILL: twigs' )                         !  Table address
    htable( hashaddr )%leaves = 0                                                ! Reset # leaves to zero for recount including non-local branches
    htable( hashaddr )%childcode =  IBSET( htable( hashaddr )%childcode,9 ) ! Set children_HERE flag for all local twig nodes
  end do


  treekey(1:ntwig) = pack(htable%key,mask = htable%node < 0)                                ! list of all twig keys excluding root
  if (me==1756) then 
     do i=1,nnodes
	write(*,'(i8,o20)') i,treekey(i)
    end do
  endif

  call sort(treekey(1:ntwig))                                                               ! Sort keys
  treekey(ntwig+1:ntwig+nleaf) = pack(htable%key,mask = htable%node > 0)                    ! add list of leaf keys

  if (me==1756) then 
	write(*,*) 'After sort'
     do i=1,nnodes
	write(*,'(i8,o20)') i,treekey(i)
    end do
  endif
call MPI_FINALIZE(ierr)
  stop
  tree_node(1) = -1  ! root node #
  cell_addr(1) = key2addr(1_8,'FILL: root')
 

  do i=2,nnodes
    hashaddr = key2addr( treekey(i),'FILL: nodes' )
    cell_addr(i) = hashaddr 
    tree_node(i) =  htable( hashaddr )%node                     ! node property pointers
    parent_key(i) = ishft(treekey(i),-3 )                    ! Parent keys, skipping root
    parent_addr(i) =  key2addr( parent_key(i),'FILL: node par' )                 ! parents' #table addresses
  end do


!  Sweep back up, and augment leaf count of parent nodes

  do i=nnodes,2,-1
     if (htable( parent_addr(i) )%owner == me ) then
        htable( parent_addr(i) )%leaves = htable( parent_addr(i) )%leaves + htable(cell_addr(i) )%leaves 
     endif
  end do

  node_level( tree_node(1:nnodes) ) = log(1.*treekey(1:nnodes))/log(8.)  ! get levels from keys and prestore as node property
  node_level(0) = 0

  ! Check tree integrity: Root node should now contain all particles!
  if (htable(1)%leaves /= npart) then
     write(*,*) 'Problem with tree on PE ',me
     write(*,*) 'Leaf checksum (',htable(1)%leaves,')  does not match # particles (',npart,')'
	call diagnose_tree
  endif


  !  Determine 'next_node' pointers for tree walk & store in hash table.
  !  Preprocess child info and store as node properties
  !  - already have (approx) sorted key list for twig nodes from leaf count above.

  do i = nnodes,2,-1
     search_key = treekey(i)                   
     node_addr = cell_addr(i)
     htable( node_addr )%next = next_node(search_key)  !   Get next sibling, uncle, great-uncle in local tree
  end do


  ! Fill in 1st child, # children in twig properties
  do i = 1,ntwig
     child_byte = htable( cell_addr(i) )%childcode                           !  Children byte-code
     nchild = SUM( (/ (ibits(child_byte,j,1),j=0,7) /) )                   ! # children = sum of bits in byte-code
     child_sub(1:nchild) = pack( bitarr, mask=(/ (btest(child_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
    child_top = ishft(treekey(i),3) 

    do j=1,nchild
       child_key(j) = IOR( child_top, child_sub(j) )         ! Construct keys of children
    end do

     first_child( tree_node(i) ) = child_key(1)   ! Store 1st child as twig-node property - used in tree_walk
     n_children( tree_node(i) ) = nchild             ! Store # children   "    "
  end do

  !  Dummy values for leaves
  first_child(1:nleaf) = treekey(ntwig+1:ntwig+nleaf) 
  n_children(1:nleaf) = 0

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize



end subroutine tree_fill

