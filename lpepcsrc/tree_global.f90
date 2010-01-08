subroutine tree_global

  use treevars
  use timings
  use tree_utils
  implicit none
  include 'mpif.h'

  real*8 :: ts1b=0., ts1e=0., ta1b=0., ta1e=0.
  real*8 :: xss, yss, zss
  integer :: i,j, ierr, maxlevel, ilevel, nparent, nsub, nuniq, child_byte, child_bit, nodtwig, hashaddr, node_addr, nchild
  integer*8 :: search_key, child_top

  integer*8, dimension(8) :: child_key, child_sub
  
  integer, dimension(nbranch_sum) :: branch_level, branch_addr, branch_node 
  integer*8, dimension(maxaddress) :: sub_key, parent_key
  integer, dimension(maxaddress) ::  tree_node, cell_addr, parent_addr, parent_node
  logical :: duplicate(maxaddress)

  integer, external :: key2addr        ! Mapping function to get hash table address from key
  integer*8, external :: next_node   ! Function to get next node key for local tree walk

  character(30) :: cfile

  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  if (tree_debug) write(ipefile,'(a)') 'TREE GLOBAL'
  if (me==0 .and. tree_debug) then
	write(*,'(a)') 'LPEPC | GLOBAL'
  endif

  if (tree_debug .and. proc_debug.eq.-1) then 
	call check_table('after make_branches ')
  else if (tree_debug .and. proc_debug==me) then
	call check_table('after make_branches ')
  endif
  
! get levels of branch nodes
  maxlevel=0
  do i=1,nbranch_sum
     branch_level(i) = log( 1.*branch_key(i) )/log(8.)
     maxlevel = max( maxlevel, branch_level(i) )        ! Find maximum level
  end do  

  nparent = 0
  parent_key(1) = 0

  do ilevel = maxlevel,1,-1                                            ! Start at finest branch level
     nsub = count( mask=branch_level(1:nbranch_sum) == ilevel )        ! Count # branches at this level
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
           hashaddr = key2addr( parent_key(i),'GLOBAL: sweep1'  )
           htable( hashaddr )%childcode = IBSET( htable( hashaddr )%childcode, child_bit )
        else if (ierr == 0 ) then
           ntwig = ntwig + 1
           ntwig_me = ntwig_me+1               ! # local twigs
	   twig_key(ntwig_me) = htable( hashaddr)%key  ! add to list of local twigs
        else
           write (ipefile,*) 'Key number ',i,' not resolved'
           call MPI_ABORT(MPI_COMM_WORLD,ierr)
           stop       
        endif

        branch_addr(i) = key2addr( sub_key(i),'PROPERTIES: fill' )   !  branches' #table addresses
        branch_node(i) = htable( branch_addr(i) )%node
        parent_addr(i) = key2addr( parent_key(i),'PROPERTIES:fill' )   ! parents' #table addresses
        parent_node(i) = htable( parent_addr(i) )%node          ! parents' node numbers
        
     end do

     ! Compute parent properties from children
     do i=nuniq,1,-1
        abs_charge( parent_node(i) ) = abs_charge( parent_node(i) ) + abs_charge( branch_node(i) )           ! Sum |q|
        charge( parent_node(i) ) = charge( parent_node(i) ) + charge( branch_node(i) )                       ! Sum q
     end do

     ! parent charges should be complete before computing coq's

     do i=nuniq,1,-1
        ! Centres of charge
        xcoc( parent_node(i) ) = xcoc(parent_node(i)) + (xcoc( branch_node(i) ) * abs_charge( branch_node(i) ) ) &
             / abs_charge( parent_node(i)) ! coq
        ycoc( parent_node(i) ) = ycoc(parent_node(i)) + (ycoc( branch_node(i) ) * abs_charge( branch_node(i) ) ) &
             / abs_charge( parent_node(i)) ! coq
        zcoc( parent_node(i) ) = zcoc(parent_node(i)) + (zcoc( branch_node(i) ) * abs_charge( branch_node(i) ) ) &
             / abs_charge( parent_node(i)) ! coq
     end do

     do i=nuniq,1,-1
        ! Shifts and multipole moments
        xss = xcoc( parent_node(i) ) - xshift( branch_node(i)  )     ! Shift vector for current child node
        yss = ycoc( parent_node(i) ) - yshift( branch_node(i) )
        zss = zcoc( parent_node(i) ) - zshift( branch_node(i) )

        xshift( parent_node(i) ) = xcoc( parent_node(i) ) ! Shift variable for next level up
        yshift( parent_node(i) ) = ycoc( parent_node(i) ) 
        zshift( parent_node(i) ) = zcoc( parent_node(i) ) 

        ! dipole moment
        xdip( parent_node(i) ) = xdip( parent_node(i) ) + xdip( branch_node(i) ) - charge( branch_node(i) )*xss 
        ydip( parent_node(i) ) = ydip( parent_node(i) ) + ydip( branch_node(i) ) - charge( branch_node(i) )*yss 
        zdip( parent_node(i) ) = zdip( parent_node(i) ) + zdip( branch_node(i) ) - charge( branch_node(i) )*zss 

        ! quadrupole moment
        xxquad( parent_node(i) ) = xxquad( parent_node(i) ) +  xxquad( branch_node(i) ) - 2*xdip( branch_node(i) )*xss &
             + charge( branch_node(i) )*xss**2 
        yyquad( parent_node(i) ) = yyquad( parent_node(i) ) +  yyquad( branch_node(i) ) - 2*ydip( branch_node(i) )*yss &
             + charge( branch_node(i) )*yss**2
        zzquad( parent_node(i) ) = zzquad( parent_node(i) ) +  zzquad( branch_node(i) ) - 2*zdip( branch_node(i) )*zss &
             + charge( branch_node(i) )*zss**2 
        xyquad( parent_node(i) ) = xyquad( parent_node(i) ) +  xyquad( branch_node(i) ) - xdip( branch_node(i) )*yss &
             - ydip( branch_node(i) )*xss + charge( branch_node(i) )*xss*yss
        yzquad( parent_node(i) ) = yzquad( parent_node(i) ) +  yzquad( branch_node(i) ) - ydip( branch_node(i) )*zss &
             - zdip( branch_node(i) )*yss + charge( branch_node(i) )*yss*zss
        zxquad( parent_node(i) ) = zxquad( parent_node(i) ) +  zxquad( branch_node(i) ) - zdip( branch_node(i) )*xss &
             - xdip( branch_node(i) )*zss + charge( branch_node(i) )*zss*xss

        magmx( parent_node(i) ) = magmx( parent_node(i) ) + magmx( branch_node(i) ) - 0.5*yss*jz( branch_node(i) ) &
		 - 0.5*zss*jy( branch_node(i) ) 
        magmy( parent_node(i) ) = magmy( parent_node(i) ) + magmy( branch_node(i) ) - 0.5*zss*jx( branch_node(i) ) &
		 - 0.5*xss*jz( branch_node(i) ) 
        magmz( parent_node(i) ) = magmz( parent_node(i) ) + magmz( branch_node(i) ) - 0.5*xss*jy( branch_node(i) ) &
		 - 0.5*yss*jx( branch_node(i) ) 

        jx( parent_node(i) ) = jx( parent_node(i) ) + jx( branch_node(i) )
        jy( parent_node(i) ) = jy( parent_node(i) ) + jy( branch_node(i) )
        jz( parent_node(i) ) = jz( parent_node(i) ) + jz( branch_node(i) )

! Multipole extent
     	size_node( parent_node(i) ) = size_node( parent_node(i) ) + size_node(branch_node(i))

     end do

     nparent = nuniq
  end do

  ! Rezero dipole and quadrupole sums of all local leaf nodes
  do i=1,nleaf
     xdip(i) = 0.
     ydip(i) = 0.
     zdip(i) = 0.
     xxquad(i) = 0.
     yyquad(i) = 0.
     zzquad(i) = 0.
     xyquad(i) = 0.
     yzquad(i) = 0.
     zxquad(i) = 0.
     magmx(i) = 0.
     magmy(i) = 0.
     magmz(i) = 0.
  end do

  if (tree_debug) call check_table('End of local fill    ')

  ta1e = MPI_WTIME()
  t_fill_local = ta1e-ta1b  
  ta1b = MPI_WTIME()  

  !  Go through twig nodes and fix # leaves in #table to include non-local branch nodes
  nnodes = ntwig + nleaf

  do i=1,ntwig_me
    hashaddr = key2addr( twig_key(i),'FILL: twigs' )                          ! Table address
    htable( hashaddr )%leaves = 0                                             ! Reset # leaves to zero for recount including non-local branches
    htable( hashaddr )%childcode =  IBSET( htable( hashaddr )%childcode,9 )   ! Set children_HERE flag for all local twig nodes
  end do

  treekey(1:ntwig) = pack(htable%key,mask = htable%node < 0)                                ! list of all twig keys excluding root

  call sort(treekey(1:ntwig))                                                               ! Sort keys
  treekey(ntwig+1:ntwig+nleaf) = pack(htable%key,mask = htable%node > 0)                    ! add list of leaf keys

  tree_node(1) = -1  ! root node #
  cell_addr(1) = key2addr(1_8,'FILL: root')
 
  do i=2,nnodes
    hashaddr = key2addr( treekey(i),'FILL: nodes' )
    cell_addr(i) = hashaddr 
    tree_node(i) =  htable( hashaddr )%node                       ! node property pointers
    parent_key(i) = ishft(treekey(i),-3 )                         ! Parent keys, skipping root
    parent_addr(i) =  key2addr( parent_key(i),'FILL: node par' )  ! parents' #table addresses
  end do

  !  Sweep back up, and augment leaf count of parent node
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
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop           
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

  ta1e = MPI_WTIME()
  t_fill_global = ta1e-ta1b  
  ts1e = MPI_WTIME()  
  t_global = ts1e - ts1b
  

end subroutine tree_global
