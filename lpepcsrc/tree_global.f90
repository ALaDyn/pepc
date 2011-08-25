subroutine tree_global

  use treevars
  use timings
  use tree_utils
  use module_htable
  use module_spacefilling
  implicit none
  include 'mpif.h'

  real*8 :: xss, yss, zss
  integer :: i, ierr, maxlevel, ilevel, nparent, nsub, nuniq, child_byte, child_bit, nodtwig, hashaddr, nchild

  integer, dimension(nbranch_sum) :: branch_level, branch_addr, branch_node 
  integer*8, dimension(maxaddress) :: sub_key, parent_key
  integer, allocatable :: tree_node(:), cell_addr(:), parent_addr(:)
  integer, dimension(maxaddress) ::  parent_node
  logical :: duplicate(maxaddress)

  call timer_start(t_global)
  call timer_start(t_fill_local)

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
     branch_level(i) = level_from_key(branch_key(i))
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
        child_bit = int(IAND( sub_key(i), hashchild))                    ! extract child bit from key: which child is it?
        child_byte = 0
        child_byte = IBSET(child_byte, child_bit)    ! convert child bit to byte code
        nodtwig = -ntwig -1     ! predicted twig # 

        call make_hashentry( parent_key(i), nodtwig, 0, child_byte, me, hashaddr, ierr )

        if ( ierr == 1 ) then     
           ! keys match, so node already exists locally
           ! Set child-bit in existing parent byte-code
           hashaddr = key2addr( parent_key(i),'GLOBAL: sweep1'  )
           htable( hashaddr )%childcode = IBSET( htable( hashaddr )%childcode, child_bit )
           nodtwig = htable( hashaddr )%node 
        else if (ierr == 0 ) then
           ntwig = ntwig + 1
           ntwig_me = ntwig_me+1               ! # local twigs
     	   twig_key(ntwig_me) = htable( hashaddr)%key  ! add to list of local twigs
        else
           write (ipefile,*) 'Key number ',i,' not resolved'
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           stop
        endif

        ! Set mm-arrays to zero (initially) for all twigs that are parents of branches and have not been initially sett in tree_local
        if (.not. BTEST( htable(hashaddr)%childcode, CHILDCODE_NODE_TOUCHED )) then
           ! zero multipole information for any entries that have not been inside the tree before
           associate(twig=>tree_nodes(nodtwig))
               twig%abs_charge = 0.  !
               twig%charge     = 0.
               twig%xcoc       = 0.
               twig%ycoc       = 0.
               twig%zcoc       = 0.
               twig%xdip       = 0.
               twig%ydip       = 0.
               twig%zdip       = 0.
               twig%xxquad     = 0.
               twig%yyquad     = 0.
               twig%zzquad     = 0.
               twig%xyquad     = 0.
               twig%yzquad     = 0.
               twig%zxquad     = 0.
               twig%xshift     = 0.
               twig%yshift     = 0.
               twig%zshift     = 0.
               htable(hashaddr)%childcode = IBSET(htable(hashaddr)%childcode,CHILDCODE_NODE_TOUCHED) ! I will now touch this again
               twig%byte       = htable(hashaddr)%childcode ! TODO: maybe inconsistent with htable data
           end associate
        endif

        branch_addr(i) = key2addr( sub_key(i),'PROPERTIES: fill' )   !  branches` #table addresses
        branch_node(i) = htable( branch_addr(i) )%node
        parent_node(i) = nodtwig                                     ! parents` node numbers
        
     end do

     ! Compute parent properties from children
     do i=nuniq,1,-1
        associate(parent=>tree_nodes(parent_node(i)), branch=>tree_nodes(branch_node(i)))
            parent%abs_charge = parent%abs_charge + branch%abs_charge           ! Sum |q|
            parent%charge     = parent%charge     + branch%charge               ! Sum q
        end associate
     end do

     ! parent charges should be complete before computing coq`s

     do i=nuniq,1,-1
        associate(parent=>tree_nodes(parent_node(i)), branch=>tree_nodes(branch_node(i)))
            ! Centres of charge
            parent%xcoc = parent%xcoc + (branch%xcoc * branch%abs_charge )  / parent%abs_charge ! coq
            parent%ycoc = parent%ycoc + (branch%ycoc * branch%abs_charge )  / parent%abs_charge ! coq
            parent%zcoc = parent%zcoc + (branch%zcoc * branch%abs_charge )  / parent%abs_charge ! coq
        end associate
     end do

     do i=nuniq,1,-1
        associate(parent=>tree_nodes(parent_node(i)), branch=>tree_nodes(branch_node(i)))
            ! Shifts and multipole moments
            xss = parent%xcoc - branch%xshift  ! Shift vector for current child node
            yss = parent%ycoc - branch%yshift
            zss = parent%zcoc - branch%zshift

            parent%xshift = parent%xcoc ! Shift variable for next level up
            parent%yshift = parent%ycoc
            parent%zshift = parent%zcoc

            ! dipole moment
            parent%xdip = parent%xdip + branch%xdip - branch%charge*xss
            parent%ydip = parent%ydip + branch%ydip - branch%charge*yss
            parent%zdip = parent%zdip + branch%zdip - branch%charge*zss

            ! quadrupole moment
            parent%xxquad = parent%xxquad +  branch%xxquad - 2*branch%xdip*xss + branch%charge*xss**2
            parent%yyquad = parent%yyquad +  branch%yyquad - 2*branch%ydip*yss + branch%charge*yss**2
            parent%zzquad = parent%zzquad +  branch%zzquad - 2*branch%zdip*zss + branch%charge*zss**2

            parent%xyquad = parent%xyquad +  branch%xyquad - branch%xdip*yss - branch%ydip*xss + branch%charge*xss*yss
            parent%yzquad = parent%yzquad +  branch%yzquad - branch%ydip*zss - branch%zdip*yss + branch%charge*yss*zss
            parent%zxquad = parent%zxquad +  branch%zxquad - branch%zdip*xss - branch%xdip*zss + branch%charge*zss*xss
        end associate
     end do

     nparent = nuniq
  end do

  ! Rezero dipole and quadrupole sums of all local leaf nodes
  do i=1,nleaf
    associate(leaf=>tree_nodes(i))
      leaf%xdip   = 0.
      leaf%ydip   = 0.
      leaf%zdip   = 0.
      leaf%xxquad = 0.
      leaf%yyquad = 0.
      leaf%zzquad = 0.
      leaf%xyquad = 0.
      leaf%yzquad = 0.
      leaf%zxquad = 0.
    end associate
  end do

  if (tree_debug) call check_table('End of local fill    ')

  call timer_stop(t_fill_local)
  call timer_start(t_fill_global)

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

  allocate(tree_node(nnodes),cell_addr(nnodes),parent_addr(nnodes))

  tree_node(1) = -1  ! root node #
  cell_addr(1) = key2addr(1_8,'FILL: root')

  do i=2,nnodes
    hashaddr = key2addr( treekey(i),'FILL: nodes' )
    cell_addr(i) = hashaddr 
    tree_node(i) =  htable( hashaddr )%node                       ! node property pointers
    parent_key(i) = ishft(treekey(i),-3 )                         ! Parent keys, skipping root
    parent_addr(i) =  key2addr( parent_key(i),'FILL: node par' )  ! parents` #table addresses
  end do

  !  Sweep back up, and augment leaf count of parent node
  do i=nnodes,2,-1
     if (htable( parent_addr(i) )%owner == me ) then
        htable( parent_addr(i) )%leaves = htable( parent_addr(i) )%leaves + htable(cell_addr(i) )%leaves 
     endif
  end do

  do i=1,nnodes
    node_level( tree_node(i) ) = level_from_key(treekey(i))  ! get levels from keys and prestore as node property
  end do
  node_level(0) = 0

  ! Check tree integrity: Root node should now contain all particles!
  if (htable(1)%leaves /= npart) then
     write(*,*) 'Problem with tree on PE ',me
     write(*,*) 'Leaf checksum (',htable(1)%leaves,')  does not match # particles (',npart,')'
     call diagnose_tree
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     stop           
  endif

  deallocate(tree_node,cell_addr,parent_addr)

  call timer_stop(t_fill_global)
  call timer_stop(t_global)

end subroutine tree_global
