
! ======================
!
!   DIAGNOSE_TREE 
!
!   Check integrity of tree structure from hash table
!   
!
! ======================

subroutine diagnose_tree


  use treevars

  implicit none
  integer :: key2addr        ! Mapping function to get hash table address from key

  integer*8 :: key_twig(ntwig), key_leaf(nleaf) 
  integer, dimension(size_tree) :: node_list, owner_list, addr_list
  integer, dimension(ntwig) :: level_twig, nchild_twig, child_twig, addr_twig, ind_twig      ! twig-nodes
  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf, owner_leaf       ! leaf-nodes

  character(30) :: cfile
  character(1) :: csnap, collision
  integer ic(5),jc(5),lc(5)

  integer :: i, ip, j, ilev, isnap, ibt
  real :: s, xt, yt

  save isnap
  data isnap/1/

  csnap=achar(mod(isnap,10)+48)




  ! output hash table

  write (ipefile,'(/a/8x,a/8x,a)') 'Final hash table ', &
       'entry,    owner    node,            key        parent       next       link   # leaves  childcode  collision', &
       '----------------------------------------------------------------------------------------------- '

  ! flag  collisions

  do i=0,maxaddress
     collision=" "
     if (htable(i)%node/=0 .and. htable(i)%link/= -1 ) collision="C"
     if (htable(i)%node /= 0 .and. htable(i)%next >=0) write (ipefile,'(3i10,3o15,i8,i5,z4,4x,a1)') &
	  i,htable(i)%owner,htable(i)%node,htable(i)%key,ishft( htable(i)%key,-idim ), htable(i)%next, &
          htable(i)%link,htable(i)%leaves,htable(i)%childcode,collision
     if (htable(i)%node /= 0 .and. htable(i)%next <0) write (ipefile,'(3i10,2o15,i15,i15,i5,z4,4x,a1)') &
	  i,htable(i)%owner,htable(i)%node,htable(i)%key,ishft( htable(i)%key,-idim ), htable(i)%next, &
          htable(i)%link,htable(i)%leaves,htable(i)%childcode,collision
  end do



 ! get keys of twig nodes from hash table
  key_twig(1:ntwig)  = pack(htable(0:maxaddress)%key,mask=htable(0:maxaddress)%node<0)
  
 ! get levels of twigs
  addr_twig(1:ntwig) = (/( key2addr( key_twig(i) ),i=1,ntwig)/)   !  Table address
  child_twig(1:ntwig) = (/( htable( key2addr( key_twig(i) ) )%childcode,i=1,ntwig )/)   !  Children byte-code
  ind_twig(1:ntwig) = (/( htable( key2addr( key_twig(i) ) )%node,i=1,ntwig )/)   !  Twig node pointers



  write (ipefile,'(///a)') 'Tree structure'

!  write (ipefile,'(/a/a/(3i5,2i10,2i8,b11,i2,i8,i10,9(1pe15.4)))') 'Twigs from hash-table:', &
  write (ipefile,'(/a/a/(3i5,2o15,2i8,z4,i2,o15,i5,10(1pe15.4)))') 'Twigs from hash-table:', &
       '    i  level  owner   key    parent-key    #     node     code    #c  1st child    #leaves ', &
       (i,node_level(ind_twig(i)), &              !  index, level
         htable( key2addr( key_twig(i) ) )%owner, &                            ! Owner-PE of node
         key_twig(i),ishft( key_twig(i),-idim ), &                             ! key, parent key
         addr_twig(i), ind_twig(i), &    ! Table address and node number
         child_twig(i), &                         ! Children byte-code 
         n_children( ind_twig(i) ), &
         first_child( ind_twig(i) ), &            ! key of 1st child
         htable( addr_twig(i) )%leaves, &                           ! # leaves contained in branch 
         abs_charge(ind_twig(i)), &    ! Twig absolute charge
         charge(ind_twig(i)), &    ! Twig  charge
         xcoc(ind_twig(i)), & ! Centre of charge
         ycoc(ind_twig(i)), &  
         xdip(ind_twig(i)), &  
         ydip(ind_twig(i)), &  
         jx(ind_twig(i)), &  
         jy(ind_twig(i)), &  
         magmx(ind_twig(i)), &  
         magmy(ind_twig(i)), &  
!         xxquad(ind_twig(i)), &  
!         yyquad(ind_twig(i)), &  
!         xyquad(ind_twig(i)), &  
         i=1,ntwig)     
       

 ! get keys of local leaf nodes from hash table
  key_leaf(1:nleaf_me) = pack(htable%key,mask=(htable%node>0 .and. htable%owner == me))
  ind_leaf(1:nleaf_me) = pack(htable%node,mask=(htable%node>0 .and. htable%owner == me))         ! particle/leaf index
  plist_leaf(1:nleaf_me) = pack(htable%childcode,mask=(htable%node>0 .and. htable%owner == me))   ! particle label
  owner_leaf(1:nleaf_me) = pack(htable%owner,mask=(htable%node>0 .and. htable%owner == me))   ! who owns leaf node


  write (ipefile,'(/a/3a5,2a10,2a15,a25,4a11/(3i5,2i10,2o15,o25,2f11.4,2f11.4))') 'Local leaves from hash-table:', &
       'i','owner','plab','i-leaf','lev','key','parent','pkey','x','y','q','jx', &
       (i,owner_leaf(i),plist_leaf(i),ind_leaf(i),node_level(ind_leaf(i)),key_leaf(i), &
        ishft( key_leaf(i),-idim ), &      ! parent
        pekey(ind_leaf(i)), &  ! particle key
        x(ind_leaf(i)),y(ind_leaf(i)), q(ind_leaf(i)), jx(ind_leaf(i)), &
	i=1,nleaf_me)

 ! get keys of NON-local leaf nodes from hash table
  key_leaf(1:nleaf-nleaf_me) = pack(htable%key,mask=(htable%node>0 .and. htable%owner /= me))
  ind_leaf(1:nleaf-nleaf_me) = pack(htable%node,mask=(htable%node>0 .and. htable%owner /= me))         ! leaf index
  plist_leaf(1:nleaf-nleaf_me) = pack(htable%childcode,mask=(htable%node>0 .and. htable%owner /= me))   ! global particle label
  owner_leaf(1:nleaf-nleaf_me) = pack(htable%owner,mask=(htable%node>0 .and. htable%owner /= me))   ! who owns leaf node


  write (ipefile,'(//a/a/(4i5,2o15,i5,2f11.4,f6.1,f11.4))') 'Non-local leaves from hash-table:', &
       '    i   owner    i-leaf    lev    key    parent  plabel  xcoc  ycoc  charge      ', &
       (i,owner_leaf(i),ind_leaf(i),node_level(ind_leaf(i)),key_leaf(i), &
        ishft( key_leaf(i),-idim ), &      ! parent
        plist_leaf(i), & ! global particle label
        xcoc(ind_leaf(i)),ycoc(ind_leaf(i)), charge(ind_leaf(i)), xdip(ind_leaf(i)), &
	i=1,nleaf-nleaf_me)

! Interaction lists
  
  write(ipefile,'(//a)') 'Interaction lists'
  do i=1,npp
     write(ipefile,'(//a,i5,a,i5)') 'Particle ',pelabel(i),' # terms: ',nterm(i)
     write(ipefile,'(a/(3i7))') 'List: ',(intlist(j,i),htable( key2addr( intlist(j,i) ) )%owner &
          ,htable( key2addr( intlist(j,i) ) )%node,j=1,nterm(i))
  end do




end subroutine diagnose_tree
