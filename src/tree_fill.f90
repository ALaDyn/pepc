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
  use utils

  implicit none

 integer*8, dimension(size_tree) :: sub_key, parent_key
 integer*8, dimension(8) :: child_key, child_sub

  integer, dimension(size_tree) :: branch_level
  integer, dimension(size_tree) :: twig_addr, twig_code,  cell_addr, tree_node, parent_addr
  integer, dimension(8) :: child_addr !  children nodes

  integer*8 :: node_key, search_key, parent,  child_top
  integer :: maxlevel, ilevel, nsub,i,j,k, nparent, nuniq, child_byte, child_bit, nchild, link_addr, hashaddr
  integer :: maxtwig, maxleaf, nleaf_check, ntwig_check
  integer ::  node_addr, jmatch(1),  parent_node, parent_level, nodtwig
  logical :: duplicate(size_tree), resolved, keymatch(8)
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer*8 :: next_node   ! Function to get next node key for local tree walk
  integer :: ierr

  if (tree_debug) write(ipefile,'(/a)') 'TREE FILL'

  if (tree_debug) call check_table('after make_branches ')

  ! get levels of branch nodes

  branch_level(1:nbranch_sum) = log( 1.*branch_key(1:nbranch_sum) )/log(2.**idim)
  maxlevel = maxval( branch_level(1:nbranch_sum) )        ! Find maximum level
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

     parent_key(1:nuniq) = ISHFT( sub_key(1:nuniq),-idim )             ! Parent keys

     do i=1,nuniq

        child_bit = IAND( sub_key(i), hashchild)                    ! extract child bit from key: which child is it?
        child_byte = 0
        child_byte = IBSET(child_byte, child_bit)    ! convert child bit to byte code
        nodtwig = -ntwig -1     ! predicted twig # 

        call make_hashentry( parent_key(i), nodtwig, 0, child_byte, me, hashaddr, ierr )

        if ( ierr == 1 ) then     
           ! keys match, so node already exists locally
           ! Set child-bit in existing parent byte-code
           hashaddr = key2addr( parent_key(i) )
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

!  treekey(1:ntwig_me) = pack(htable%key,mask = ( htable%owner == me .and. htable%node <0))  ! list of locally owned twigs
!  treekey(1:ntwig_me) = twig_key(1:ntwig_me)   ! list of locally owned twigs
!  cell_addr(1:ntwig_me) = (/( key2addr( treekey(i) ),i=1,ntwig_me) /)                       !  Table address
!  cell_addr(1:ntwig_me) = (/( key2addr( treekey(i) ),i=1,ntwig_me) /)                       !  Table address
!  htable( cell_addr(1:ntwig_me) )%leaves = 0                                                ! Reset # leaves to zero for recount including non-local branches
!  htable( cell_addr(1:ntwig_me) )%childcode =  IBSET( htable( cell_addr(1:ntwig_me) )%childcode,9 ) ! Set children_HERE flag for all local twig nodes

  do i=1,ntwig_me
    hashaddr = key2addr( twig_key(i) )                         !  Table address
    htable( hashaddr )%leaves = 0                                                ! Reset # leaves to zero for recount including non-local branches
    htable( hashaddr )%childcode =  IBSET( htable( hashaddr )%childcode,9 ) ! Set children_HERE flag for all local twig nodes
  end do


  treekey(1:ntwig) = pack(htable%key,mask = htable%node < 0)                                ! list of all twig keys excluding root
  call sort(treekey(1:ntwig))                                                               ! Sort keys

  treekey(ntwig+1:ntwig+nleaf) = pack(htable%key,mask = htable%node > 0)                    ! add list of leaf keys

!  cell_addr(1:nnodes) = (/( key2addr( treekey(i) ),i=1,nnodes) /)                           ! table address
!  tree_node(1:nnodes) = (/ (htable( cell_addr(i) )%node, i=1,nnodes) /)                     ! node property pointers

!  parent_key(2:nnodes) = ishft(treekey(2:nnodes),-idim )                                    ! Parent keys, skipping root
!  parent_addr(2:nnodes) = (/( key2addr( parent_key(i) ),i=2,nnodes) /)                      ! parents' #table addresses

  tree_node(1) = -1  ! root node #
  cell_addr(1) = key2addr(1_8)

  do i=2,nnodes
    hashaddr = key2addr( treekey(i) )
    cell_addr(i) = hashaddr 
    tree_node(i) =  htable( hashaddr )%node                     ! node property pointers
    parent_key(i) = ishft(treekey(i),-idim )                    ! Parent keys, skipping root
    parent_addr(i) =  key2addr( parent_key(i) )                 ! parents' #table addresses
  end do


!  Sweep back up, and augment leaf count of parent nodes

  do i=nnodes,2,-1
     if (htable( parent_addr(i) )%owner == me ) then
        htable( parent_addr(i) )%leaves = htable( parent_addr(i) )%leaves + htable(cell_addr(i) )%leaves 
     endif
  end do

  node_level( tree_node(1:nnodes) ) = log(1.*treekey(1:nnodes))/log(2.**idim)  ! get levels from keys and prestore as node property
  node_level(0) = 0

  ! Check tree integrity: Root node should now contain all particles!
  if (htable(1)%leaves /= npart) then
     write(*,*) 'Problem with tree on PE ',me
     write(*,*) 'Leaf checksum (',htable(1)%leaves,')  does not match # particles (',npart,')'
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
     nchild = SUM( (/ (ibits(child_byte,j,1),j=0,2**idim-1) /) )                   ! # children = sum of bits in byte-code
     child_sub(1:nchild) = pack( bitarr, mask=(/ (btest(child_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
    child_top = ishft(treekey(i),idim) 
    child_key(1:nchild) = IOR( child_top, child_sub(1:nchild) )         ! Construct keys of children

     first_child( tree_node(i) ) = child_key(1)   ! Store 1st child as twig-node property - used in tree_walk
     n_children( tree_node(i) ) = nchild             ! Store # children   "    "
  end do

  !  Dummy values for leaves
  first_child(1:nleaf) = treekey(ntwig+1:ntwig+nleaf) 
  n_children(1:nleaf) = 0

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize



end subroutine tree_fill

