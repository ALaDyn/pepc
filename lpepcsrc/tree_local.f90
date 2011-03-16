subroutine tree_local

  use treevars
  use timings
  use tree_utils
  use module_math_tools

  implicit none
  include 'mpif.h'

  real*8 :: gamma, vx, vy, vz
  integer :: i, k, j, i1, i2, ierr, level_match, level_diff, ibit, iend, level, nbound, newleaf, nres, ncoll, &
       link_addr, ipoint, parent_addr, childbyte, treelevel, ncheck, nsubset, newsub, cchild, nchild, addr_leaf, &
       p_leaf, ntwig_domain, node_leaf, addr_twig, addr_child, nsearch
  integer*8 ::  key_lo, cell1, cell2, parent_key
  logical :: resolved

  integer*8, dimension(8) :: sub_key, key_child   ! Child partial key
  integer, dimension(8) :: node_child   
  real*8, dimension(8) :: xs, ys, zs   ! multipole shift vector

  integer*8, dimension(nppm+2) :: local_key
  integer, dimension(nppm+2)  :: local_plist, local_owner, local_ind
  integer*8, allocatable :: subcell(:)
  integer, allocatable :: cell_addr(:)
  integer*8, dimension(maxaddress) :: res_key
  integer, dimension(maxaddress) :: newentry, res_addr, res_node, res_child, res_owner
  integer*8, dimension(nbranch_max) :: search_key, resolve_key

  ! stuff for getting a virtual domain
  integer*8 :: right_limit_me, right_limit
  integer*8 :: left_limit_me, left_limit
  integer*8 :: left_virt_limit, right_virt_limit
  integer*8 :: left_cell, right_cell
  integer*8 :: search_key_mod

  ! stuff for estimation
  integer*8 :: branch_level(1:nlev)
  integer*8 :: branch_level_D1(1:nlev) ! start at level 1
  integer*8 :: branch_level_D2(1:nlev) 
  integer*8 :: D1, D2               ! sub-domains
  integer*8 :: L                    ! inner limit
  integer*8 :: branch_max_local     ! estimation for local branches
  integer*8 :: branch_max_global    ! estimation for global branches
  integer*8 :: ilevel, pos

  integer,external :: key2addr        ! Mapping function to get hash table address from key

!!! --------------- TREE BUILD ---------------

 
  call timer_start(t_local)
  call timer_start(t_build)
  call timer_start(t_build_neigh)

  call OutputMemUsage(18, "[tree_local]", memory_debug .and. (me==0), 59)

  if (tree_debug) write(ipefile,'(/a)') 'TREE LOCAL'
  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | LOCAL BUILD'

  ! zero table: need list of 'live' addresses to speed up
  htable = hash(0,0_8,-1,0,0,0_8,0)

  do i=1,npp+2
     local_plist(i) = pelabel(i)       ! Particle (global) label for tracking purposes 
     local_key(i) = pekey(i)           ! Particle key
     local_owner(i) = pepid(i)
  end do
  i=0                     ! List of free (unused #table-addresses) 

  !  Check whether boundary particles are close enough to warrant inclusion
  !  in tree construction.
  
  nlist = npp

  ! RH neighbour PE 
  if ( me /= num_pe-1 ) then
     ! First find level shared by last particle pair in list
     key_lo = ieor( pekey(npp), pekey(npp-1) )   ! Picks out 1st position where keys differ
     level_diff =  int(log(1.*key_lo)/log(8.))
     level_match = max( nlev - level_diff, 1 )    ! Excludes low-level end-cells for discontinuous domains
     ibit = nlev-level_match               ! bit shift factor 

     cell1 = ishft( pekey(npp), -3*ibit )    ! Subcell number of RH boundary particle
     cell2 = ishft( pekey(npp+1), -3*ibit )    ! Subcell number of LH boundary particle of neighbouring PE
     if (cell1 == cell2) then
        nlist = nlist+1                    ! If keys match at this level then include boundary particle in local tree
        local_plist(nlist) = pelabel(npp+1)
        local_key(nlist) = pekey(npp+1)
        local_owner(nlist) = me+1
        if (build_debug) write (ipefile,'(a,i5,a)') 'boundary particle from PE',me+1,' included in list'
     endif

  endif

  ! 2nd set: LH neighbour PE
  if ( me /= 0 ) then
     iend = npp+2
     if ( me == num_pe-1 ) iend = npp+1         ! End node only has one boundary particle
     ! First find level shared by first particle pair in list
     key_lo = ieor( pekey(1), pekey(2)  )   ! Picks out lower order bits where keys differ
     level_diff =  int(log(1.0*key_lo)/log(8.))

     level_match = max( nlev - level_diff, 1 )    ! Excludes low-level end-cells for discontinuous domains
     ibit = nlev-level_match               ! bit shift factor 

     cell1 = ishft( pekey(1), -3*ibit )    ! Subcell number of LH boundary particle
     cell2 = ishft( pekey(iend), -3*ibit )    ! Subcell number of RH boundary particle of neighbouring PE
     if (cell1 == cell2) then
        nlist = nlist+1         ! If keys match at this level then include boundary particle in local tree
        local_plist(nlist) = pelabel(iend)
        local_key(nlist) = pekey(iend) 
        local_owner(nlist) = me-1
        if (build_debug) write (ipefile,'(a,i5,a)') 'boundary particle from PE',me-1,' included in list'
     endif
  endif

  call timer_stop(t_build_neigh)
  call timer_start(t_build_part)

  allocate(subcell(1:nlist),cell_addr(1:nlist))

  level = 0                        ! start at root
  local_ind(1:nlist) = (/(k,k=1,nlist)/) ! Sorted local particle index / leaf number: points to particle properties

  !  set up root node

  htable(1)%node   =  -1              !  node #
  htable(1)%owner  =  me              ! Owner
  htable(1)%key    =   1_8            !  key
  htable(1)%link   =  -1              !  collision link
  htable(1)%leaves = npp              ! root contains all leaves, excluding boundary particles
  htable(1)%next   =   1              ! root points to itself as next to abort tree walk even after directly interacting with root node
  ntwig = 1
  nleaf = 0
  tablehigh = 0


  do while ( nlist > 0 )              ! While any particle not finished:

     level = level + 1                ! Next sublevel

     if (level>nlev) then
       write(*,*) 'Problem with tree on PE ',me,' - no more levels '
       write(*,'(a/(i8,o30))') 'Remaining keys: ',(local_plist(i),local_key(i),i=1,nlist)
       call MPI_ABORT(MPI_COMM_WORLD,ierr)
       stop
     endif

     ibit = nlev - level               ! bit shift factor (0=highest leaf node, nlev-1 = root)
     newleaf=0                         ! # new leaves at this level
     nbound=0                          ! # boundary particles located (max 2)

     ! Determine subcell # from key
     ! At a given level, these will be unique   
     do i=1,nlist
        subcell(i) = ishft( local_key(i), -3_8*ibit )    
        cell_addr(i) = int(IAND( subcell(i), hashconst))         ! cell address hash function
     end do

     if (build_debug) then
        write (ipefile,'(/a//a,i5,a5,i5,a5,z20/)') '---------------','Starting level ',level, &
             ' ibit',ibit,' iplace',iplace
        write (ipefile,*) 'i,   p,     key,           owner     cell_key,    cell_addr'
        write (ipefile,'(2i4,z20,i3,o20,i6)') &
             (i,local_plist(i),local_key(i),local_owner(i),subcell(i),cell_addr(i),i=1,nlist) 
     endif

     nres =0          ! # resolved entries
     ncoll = 0

     !  Make set of entries at current level and tag collisions
     do i = 1,nlist
        tablehigh = max(tablehigh,cell_addr(i))                 ! Track highest address
        if ( htable( cell_addr(i) )%node == 0 .and. htable(cell_addr(i))%key /= -1 ) then       ! Is entry empty?
           nres = nres + 1
           newentry(nres) = cell_addr(i)                  ! Yes, so create new entry:
           htable(cell_addr(i))%node = i                  !   local pointer
           htable(cell_addr(i))%key = subcell(i)          !   key
           htable(cell_addr(i))%leaves = 1                !   # children
           htable(cell_addr(i))%owner = local_owner(i)    ! Set owner of node equal to owner of particle

        else if ( htable( cell_addr(i) )%node /= 0 .AND. &   ! Entry exists and keys match
             htable(cell_addr(i))%key == subcell(i) ) then   ! => twig cell  

           htable(cell_addr(i))%node = i               ! Overwrite previous (leaf) entry
           htable(cell_addr(i))%leaves =  htable(cell_addr(i))%leaves + 1
           htable(cell_addr(i))%owner = me             ! Set owner of node equal to local PE

        else                                    ! Entry exists and keys do not match
           ncoll = ncoll + 1                    ! Increment collision count
           res_addr(ncoll) = cell_addr(i)       ! Reduced list of addresses 
           res_key(ncoll) = subcell(i)          !     ..    .. of keys
           res_node(ncoll) = i                  !     ..    .. of pointers
           res_child(ncoll) =1
           res_owner(ncoll) =local_owner(i)
        end if
     end do

     ! Make reduced list of empty addresses to speed up collision resolution
     i2 = maxaddress  ! count back from end of table
     do i=1,ncoll
        i1 = i2
        do while (htable(i1)%node /= 0 .or. htable(i1)%key==-1 ) 
           i1=i1-1  ! skip existing or dummy entries
        end do
        free_addr(i) = all_addr(i1)
        i2 = i1-1
     end do

     if (build_debug) &
          write(ipefile,*) '# collisions ',ncoll,' free addresses: ',(free_addr(i),i=1,ncoll)

     do i = 1,ncoll

        if ( htable( res_addr(i) )%link == -1) then     ! Entry occupied without link

           nres = nres + 1                             ! create entry at empty address
           newentry(nres) = free_addr(i)                 

           htable( free_addr(i) )%node = res_node(i)                     
           htable( free_addr(i) )%key = res_key(i)
           htable( free_addr(i) )%leaves = res_child(i)
           htable( free_addr(i) )%owner = res_owner(i)
           htable( free_addr(i) )%link = -1
           htable( res_addr(i) )%link = free_addr(i)     ! Create link from 1st duplicate entry to new entry

           ! need to check here for colliding twig nodes

        else if ( htable( res_addr(i) )%link /= -1 ) then     ! Occupied with link already

           link_addr = res_addr(i)                            ! Start of chain
           resolved = .false.                                 ! Resolve flag

           do while ( .not. resolved )       ! Find end of chain  
              link_addr = htable(link_addr)%link

              if ( htable( link_addr )%key == res_key(i) ) then
                 ! Occupied with same key -> twig node
                 htable( link_addr )%leaves =  htable( link_addr )%leaves + 1
                 resolved = .true.
              else if ( htable (link_addr)%link == -1 ) then
                 nres = nres + 1                             ! create entry at empty address
                 newentry(nres) = free_addr(i)
                 htable( free_addr(i) )%node = res_node(i)                     
                 htable( free_addr(i) )%key = res_key(i)
                 htable( free_addr(i) )%leaves = res_child(i)
                 htable( free_addr(i) )%owner = res_owner(i)
                 htable( free_addr(i) )%link = -1                
                 htable( link_addr )%link = free_addr(i)      ! Create link from last duplicate entry to displaced entry
                 resolved = .true.
              else
                 ! not yet resolved
              endif
           end do

        else
           write (ipefile,*) 'Key number ',i,' not resolved'
           call MPI_ABORT(MPI_COMM_WORLD,ierr)
           stop           
        end if
        tablehigh = max(tablehigh,free_addr(i))                 ! Track highest address

     end do

     ! Go through new entries and sort them into twigs/leaves
     do i = 1,nres

        if ( htable( newentry(i) )%leaves == 1 .and. htable(newentry(i))%owner == me) then
           ! create new leaf
           newleaf = newleaf + 1 
           nleaf = nleaf + 1 
           ipoint =  htable( newentry(i) )%node                  ! local pointer
           htable( newentry(i) )%leaves = 1                      ! contained leaves
           htable( newentry(i) )%childcode = pelabel(local_ind(ipoint))  ! store label in #-table
           htable( newentry(i) )%node = local_ind(ipoint)        ! store index in #-table
           htable( newentry(i) )%owner = me                      ! Set owner
           local_plist( ipoint ) = 0                             ! label as done
           leaf_key(nleaf) = htable(newentry(i))%key	         ! Add to key list
 
        else if ( htable( newentry(i) )%leaves == 1 .and. htable(newentry(i))%owner /= me) then
           ! unwanted leaf generated by boundary particle, so remove entry
           ! -> this prevents duplicate branch nodes being generated in make_branches
           nbound = nbound + 1
           ipoint =  htable( newentry(i) )%node         ! local pointer
           local_plist( ipoint ) = 0                    ! label as done - removes particle from list
           htable( newentry(i) )%node = 0               ! remove node from #table
           htable( newentry(i) )%key = -1_8             ! but retain %link to dummy entry 
                                                        ! in case it is in the middle of a chain
           htable( newentry(i) )%leaves = 0
           htable( newentry(i) )%childcode = 0 
        else if ( htable( newentry(i) )%leaves > 1 ) then
           ! twig
           twig_key(ntwig) = htable(newentry(i))%key	         ! Add to key list
           ntwig = ntwig + 1
           htable( newentry(i) )%node = -ntwig
           htable( newentry(i) )%owner = me                       ! Set owner
           htable( newentry(i) )%childcode = 0  ! Zero children byte-code in twig nodes
        else
           write (ipefile,*) 'Problem with flagging'
           call MPI_ABORT(MPI_COMM_WORLD,ierr)
           stop           
        endif

     end do

     ! Make new lists from unfinished particles.

     k=0
     do i=1,nlist
     	if (.not. (local_plist(i)==0)) then
	      k=k+1
          local_key(k) = local_key(i)
          local_ind(k) = local_ind(i)
	      local_plist(k) = local_plist(i)
	      local_owner(k) = local_owner(i)
	    endif
     end do
     
     nlist = nlist - newleaf - nbound

     if (build_debug) then
        !  Go through list of unresolved entries
        write (ipefile,*) 'So far: ',nleaf,' leaves, ',newleaf, 'new on level,',ntwig,' twigs, ',ncoll,' collisions'
        write (ipefile,*) 'Particles left: ',nlist,', new list:',local_plist(1:nlist)
     endif

  end do

  deallocate(subcell,cell_addr)

  call timer_stop(t_build_part)
  call timer_start(t_build_byte)

  nnodes = nleaf + ntwig  ! total number of local tree nodes

  treekey(1:nleaf) = leaf_key(1:nleaf)
  twig_key(ntwig) = 1  ! add root for later
  treekey(nleaf+1:nnodes) = twig_key(1:ntwig)  ! add twigs to list 
  nsubset=0

  do i=1,nnodes-1
     childbyte = int(IAND( treekey(i), hashchild))    ! extract last 3 bits from key
     parent_key = ishft( treekey(i),-3 )      ! parent key
     parent_addr = key2addr(parent_key,'BUILD: childbyte')
     ! Construct children byte-code (8 settable bits indicating which children nodes present)
     htable(parent_addr)%childcode = ibset( htable(parent_addr)%childcode, childbyte )
     treelevel  = int(log(1.*treekey(i))/log(8.))     ! node levels
     if (treelevel==1) then
        nsubset = nsubset+1  ! # nodes at level 1
        search_key(nsubset) = treekey(i)   ! Subset of nodes at same level
     endif
  end do
  treelevel  = int(log(1.*treekey(nnodes))/log(8.))     ! node levels
  if (treelevel==1) then
     nsubset = nsubset+1  ! # nodes at level 1
     search_key(nsubset) = treekey(nnodes)   ! Subset of nodes at same level
  endif

  sum_unused = 0
  iused = 1   ! reset used-address counter
  do i=0, maxaddress
     if (htable(i)%node == 0 .and. htable(i)%key /=-1 .and. i> free_lo) then
        sum_unused = sum_unused+1
        free_addr(sum_unused) = i            ! Free address list for resolving collisions
        point_free(i) = sum_unused           ! Index
     else
        point_free(i) = 0
     endif
  enddo
  
  call timer_stop(t_build_byte)
  call timer_stop(t_build)

!!! --------------- TREE BRANCHES (local part) ---------------

  call timer_start(t_branches_find)
  
  nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
  ntwig_me = ntwig
  if (tree_debug .and. (proc_debug==me .or.proc_debug==-1)) call check_table('after treebuild     ')
  
  ! get local key limits
  left_limit_me=pekey(1)
  right_limit_me=pekey(npp)

  ! get key limits for neighbor PE`s
  ! and build virtual limits, so that a minimum set a branch nodes comes arround
  ! boundary PE`s can access their boundary space fully only need one virtual limit
  if(me.eq.0)then
     right_limit=pekey(npp+1)
     right_virt_limit = bpi(right_limit_me,right_limit)
     left_virt_limit=2_8**(nlev*3)
  else if(me.eq.(num_pe-1))then
     left_limit=pekey(npp+1)
     left_virt_limit  = bpi(left_limit,left_limit_me)
     right_virt_limit=2_8**(3*nlev+1)-1
  else
     left_limit=pekey(npp+2)
     right_limit=pekey(npp+1)
     left_virt_limit  = bpi(left_limit,left_limit_me)
     right_virt_limit = bpi(right_limit_me,right_limit)
  end if

  ! First find highest power in the Virtual Domain
  L = bpi(left_virt_limit,right_virt_limit)
  
  ! divide in two sub-domains
  D1 = L-left_virt_limit
  D2 = right_virt_limit-L+1
  
  ! get estimation number of branches at all levels
  do ilevel=1,nlev
     pos=3*(nlev-ilevel)
     branch_level_D1(ilevel)=ibits(D1,pos,3_8)
     branch_level_D2(ilevel)=ibits(D2,pos,3_8)
     branch_level(ilevel)=branch_level_D1(ilevel)+branch_level_D2(ilevel)
  end do

  ! estimate local number
  branch_max_local = SUM(branch_level(1:nlev))

  call MPI_REDUCE(branch_max_local, branch_max_global, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if(tree_debug .and. me.eq.0) write(*,*) "global branch estimation: ", branch_max_global

  ! adapt left Virtual Limit 
  if(me.ne.0)then
     left_virt_limit=left_virt_limit-1
  end if

  ! Determine minimum set of branch nodes making up local domain
  ncheck  = 0
  nbranch = 0
  newsub  = 0
  level   = 1
  nsubset = 0


  ! Artificial Startup: mimic first level cells + placeholder bit
  ! TODO: works only for AND-W&S-hash table
  do i=0,7
     if(htable((8+i))%key.eq.(8+i))then ! entry exists
        nsubset = nsubset+1          ! # nodes at level 1
        search_key(nsubset) = (8+i)  ! Subset of nodes at same level
     end if
  end do

  ! while any particle is not in a branch
  do while ( ncheck < nleaf )

     if (nsubset > 0 ) then
        
        ! calculate parent cell for this level from virtual limits 
        left_cell  = ibits(left_virt_limit, 3_8*(nlev-level),3_8*level)
        right_cell = ibits(right_virt_limit,3_8*(nlev-level),3_8*level)

        ! for all nodes at this level
        do i=1,nsubset
           
           ! Important: discount placeholder-bit from search_key-entries
           search_key_mod=search_key(i)-2**(3*level)
           
           ! check if key between limits
           ! Important: discount placeholder-bit from search_key-entries
           if ( me.eq.0 .and. ((search_key_mod < right_cell) .or. &
		( htable( key2addr( search_key(i),'BRANCHES: search' ) )%node > 0 ))) then
              ! twig node between limits or
              ! node is a leaf
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i),'BRANCHES: ncheck' ) )%leaves  ! Augment checksum

           elseif ( me.eq.num_pe-1 .and. ((search_key_mod > left_cell) .or. &
		( htable( key2addr( search_key(i),'BRANCHES: search' ) )%node > 0 ) )) then
              ! twig node between limits or
              ! node is a leaf
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i),'BRANCHES: ncheck' ) )%leaves  ! Augment checksum

           elseif (  ((search_key_mod > left_cell) .and. (search_key_mod < right_cell)) .or. &
		( htable( key2addr( search_key(i),'BRANCHES: search' ) )%node > 0 )) then
              ! twig node between limits or
              ! node is a leaf
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i),'BRANCHES: ncheck' ) )%leaves  ! Augment checksum

           else 
              ! end node: check for complete twigs; otherwise subdivide
              
              ! get children byte-code from hash-table
              cchild = htable( key2addr( search_key(i),'BRANCHES: cchild' ) )%childcode   

              ! get number of children (sum of 1-bits in childcode)
              nchild = SUM( (/ (ibits(cchild,j,1),j=0,7) /) ) 
              
              ! extract sub key from childcode
              sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(cchild,j),j=0,7) /) )  

              ! put children in list for next level branching
              do j=1,nchild
                 
                 ! Construct keys of children
                 resolve_key(newsub+j) = IOR( ishft( search_key(i),3 ), sub_key(j) )
                 
              end do
              newsub = newsub + nchild

           endif
        end do

        if (branch_debug) then
           write (ipefile,'(/a,i7,a,i7/a/(i5,o16,i6,z5,i8))') 'Branches at level:',level, ' Checksum: ',ncheck, &
                '    i      key         node     code     #leaves', &
                (i,search_key(i), &
                htable( key2addr( search_key(i),'BRANCHES: debug' ) )%node, &      ! Node #
                htable( key2addr( search_key(i),'BRANCHES: debug' ) )%childcode, & ! Children byte-code 
                htable( key2addr( search_key(i),'BRANCHES:debug' ) )%leaves, &     ! # leaves contained in branch 
                i=1,nsubset)
        endif

     endif
     
     ! determine branches in next level
     level = level + 1
     
     ! refresh search list for next level with children of twigs which are not a branch at this level
     search_key(1:newsub) = resolve_key(1:newsub) 
     nsubset = newsub
     newsub = 0
  end do

  if (branch_debug) write (ipefile,'(/a/(i6,o16))') 'Domain branch list:',(i,pebranch(i),i=1,nbranch)
  
  if (ncheck .ne. nleaf) then
     write(*,*) 'Checksum ',ncheck,' /= # leaves on PE ',me
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop               
  endif 

  if (tree_debug .and. (proc_debug==me .or.proc_debug==-1)) call check_table('after local branches     ')
  
  call timer_stop(t_branches_find)


!!! --------------- TREE PROPERTIES (local part) ---------------


  call timer_start(t_props_leafs)

  !  Start with *local* leaf properties
  do i=1, nleaf_me
     addr_leaf = key2addr( leaf_key(i),'PROPERTIES: local' )   !  Table address
     p_leaf = htable( addr_leaf )%node   !  Local particle index  - points to properties on PE
     htable(addr_leaf)%childcode = IBSET( htable(addr_leaf)%childcode, CHILDCODE_NODE_TOUCHED ) ! I have touched this node, do not zero its properties (in tree_global)
     node_leaf = p_leaf   !  Leaf node index is identical to particle index for *local* leaves 

     xcoc( node_leaf ) = x( p_leaf )         ! Centre of charge
     ycoc( node_leaf ) = y( p_leaf )
     zcoc( node_leaf ) = z( p_leaf )

     charge( node_leaf ) = q( p_leaf )       ! Charge
     abs_charge( node_leaf ) = abs( q(p_leaf) )   ! Absolute charge (needed for c.o.c)

     xdip( node_leaf ) = x(p_leaf) * q(p_leaf)   ! Dipole moment
     ydip( node_leaf ) = y(p_leaf) * q(p_leaf) 
     zdip( node_leaf ) = z(p_leaf) * q(p_leaf) 

     xxquad( node_leaf ) = q(p_leaf) * x(p_leaf)**2   ! Quadrupole moment
     yyquad( node_leaf ) = q(p_leaf) * y(p_leaf)**2   
     zzquad( node_leaf ) = q(p_leaf) * z(p_leaf)**2   
     xyquad( node_leaf ) = q(p_leaf) * x(p_leaf) * y(p_leaf) 
     yzquad( node_leaf ) = q(p_leaf) * y(p_leaf) * z(p_leaf)  
     zxquad( node_leaf ) = q(p_leaf) * z(p_leaf) * x(p_leaf)   

     !  get true velocities from momenta
     gamma = sqrt(1.0+ux(p_leaf)**2+uy(p_leaf)**2 + uz(p_leaf)**2)
     vx = ux(p_leaf)/gamma
     vy = uy(p_leaf)/gamma
     vz = uz(p_leaf)/gamma
     jx( node_leaf ) = vx * q(p_leaf)   ! Drift moment
     jy( node_leaf ) = vy * q(p_leaf) 
     jz( node_leaf ) = vz * q(p_leaf) 

     magmx( node_leaf ) = 0.5*q(p_leaf) * ( y(p_leaf)*vz - z(p_leaf)*vy)  ! Magnetic dipole moment
     magmy( node_leaf ) = 0.5*q(p_leaf) * ( z(p_leaf)*vx - x(p_leaf)*vz)  
     magmz( node_leaf ) = 0.5*q(p_leaf) * ( x(p_leaf)*vy - y(p_leaf)*vx)  

     !  Zero shift vector
     xshift( node_leaf ) = 0.
     yshift( node_leaf ) = 0.
     zshift( node_leaf ) = 0.

     size_node( node_leaf ) = x(p_leaf)**2 + y(p_leaf)**2 + z(p_leaf)**2

  end do

  call timer_stop(t_props_leafs)
  call timer_start(t_props_twigs)

  !  Accumulate local twig properties
  !  Make list of twig nodes contained within local branch list (pebranch).  Recursive search adapted from make_branches.

  search_key(1:nbranch) = pebranch(1:nbranch)            ! start with branch list
  ncheck = 0       ! Checksum
  newsub = 0
  ntwig_domain = 0  ! # twigs contained in branch list
  nsearch = nbranch
  res_key = 0

  do while ( ncheck < nleaf_me )     ! Repeat until all local leaves found

     do i=1,nsearch
        if (  htable( key2addr( search_key(i),'PROPERTIES: twigs' ) )%node > 0 ) then
           !  leaf,  so skip and increment checksum
           ncheck = ncheck +  1  

        else 
           ! twig: add to list and  subdivide
           ntwig_domain = ntwig_domain + 1
           res_key(ntwig_domain) = search_key(i)

           cchild = htable( key2addr( search_key(i),'PROPERTIES: tw2' ) )%childcode   !  Children byte-code
           nchild = SUM( (/ (ibits(cchild,j,1),j=0,7) /) ) ! # children = sum of bits in byte-code
           sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(cchild,j),j=0,7) /) )  ! Extract child sub-keys from byte code

           resolve_key(newsub+1:newsub+nchild) = IOR( ishft( search_key(i),3 ), sub_key(1:nchild) ) ! Add keys of children to new search list
           newsub = newsub + nchild

        endif
     end do

     search_key(1:newsub) = resolve_key(1:newsub)        ! Put children into search list
     nsearch = newsub  ! # new nodes to search
     newsub = 0
  end do

  call sort(res_key(1:ntwig_domain))
  
  if (props_debug) then
     write (ipefile,*) '# Twigs contained in local branch list: ',ntwig_domain
     write (ipefile,*) 'Found ',ncheck,' out of ',nleaf,' leaves'
  endif

  do i=1,ntwig_domain
     addr_twig = key2addr( res_key(i),'PROPERTIES: domain' )   !  Table address
     res_node(i) = htable( addr_twig )%node   !  Twig node index  
     res_child(i) = htable( addr_twig )%childcode   !  Twig children byte-code 
     htable(addr_twig)%childcode = IBSET( htable(addr_twig)%childcode, CHILDCODE_NODE_TOUCHED ) ! I have touched this node, do not zero its properties (in tree_global)
  end do  

  ! Go up through tree, starting at deepest level (largest key first)
  ! and accumulate multipole moments onto twig nodes
  do i = ntwig_domain,1,-1
     nchild = SUM( (/ (ibits(res_child(i),j,1),j=0,7) /) )                 ! Get # children
     sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(res_child(i),j),j=0,7) /) )  ! Extract sub key from byte code

     do j=1,nchild
        key_child(j) = IOR( ishft( res_key(i),3 ), sub_key(j) )      ! Construct keys of children
        addr_child = key2addr( key_child(j),'PROPERTIES: domain2' )             ! Table address of children
        node_child(j) = htable( addr_child )%node                     ! Child node index  
     end do

     charge( res_node(i) ) = SUM( charge( node_child(1:nchild) ) )                 ! Sum charge of child nodes and place
     ! result with parent node
     abs_charge( res_node(i) ) = SUM( abs_charge( node_child(1:nchild) ) )                 ! Sum |q|

     ! Centres of charge
     xcoc( res_node(i) ) = 0.
     ycoc( res_node(i) ) = 0.
     zcoc( res_node(i) ) = 0.

     do j=1,nchild
        xcoc( res_node(i) ) = xcoc( res_node(i) ) + ( xcoc(node_child(j)) * abs_charge(node_child(j))) &
             / abs_charge( res_node(i) )
        ycoc( res_node(i) ) = ycoc( res_node(i) ) + ( ycoc(node_child(j)) * abs_charge(node_child(j))) &
             / abs_charge( res_node(i) )
        zcoc( res_node(i) ) = zcoc( res_node(i) ) + ( zcoc(node_child(j)) * abs_charge(node_child(j))) &
             / abs_charge( res_node(i) )
     end do

     ! Shifts and multipole moments
     xs(1:nchild) = xcoc( res_node(i) ) - xshift( node_child(1:nchild) )     ! Shift vector for current node
     ys(1:nchild) = ycoc( res_node(i) ) - yshift( node_child(1:nchild) )
     zs(1:nchild) = zcoc( res_node(i) ) - zshift( node_child(1:nchild) )

     xshift( res_node(i) ) = xcoc( res_node(i) ) ! Shift variable for next level up
     yshift( res_node(i) ) = ycoc( res_node(i) ) 
     zshift( res_node(i) ) = zcoc( res_node(i) ) 

     xdip( res_node(i) ) = 0.
     ydip( res_node(i) ) = 0.
     zdip( res_node(i) ) = 0.
     
     xxquad( res_node(i) ) = 0.
     yyquad( res_node(i) ) = 0.
     zzquad( res_node(i) ) = 0.

     xyquad( res_node(i) ) = 0.  
     yzquad( res_node(i) ) = 0.
     zxquad( res_node(i) ) = 0.

     magmx( res_node(i) ) = 0.
     magmy( res_node(i) ) = 0.
     magmz( res_node(i) ) = 0.

     jx( res_node(i) ) = 0.
     jy( res_node(i) ) = 0.
     jz( res_node(i) ) = 0.

     size_node(res_node(i)) = 0.

     do j = 1,nchild
        ! dipole moment

        xdip( res_node(i) ) = xdip( res_node(i) ) + xdip( node_child(j) ) - charge( node_child(j) )*xs(j) 
        ydip( res_node(i) ) = ydip( res_node(i) ) + ydip( node_child(j) ) - charge( node_child(j) )*ys(j) 
        zdip( res_node(i) ) = zdip( res_node(i) ) + zdip( node_child(j) ) - charge( node_child(j) )*zs(j) 

        ! quadrupole moment
        xxquad( res_node(i) ) = xxquad( res_node(i) ) + xxquad( node_child(j) ) - 2*xdip( node_child(j) )*xs(j) &
             + charge( node_child(j) )*xs(j)**2 
        yyquad( res_node(i) ) = yyquad( res_node(i) ) + yyquad( node_child(j) ) - 2*ydip( node_child(j) )*ys(j) &
             + charge( node_child(j) )*ys(j)**2 
        zzquad( res_node(i) ) = zzquad( res_node(i) ) + zzquad( node_child(j) ) - 2*zdip( node_child(j) )*zs(j) &
             + charge( node_child(j) )*zs(j)**2 

        xyquad( res_node(i) ) = xyquad( res_node(i) ) + xyquad( node_child(j) ) - xdip( node_child(j) )*ys(j) &
             - ydip( node_child(j) )*xs(j) + charge( node_child(j) )*xs(j)*ys(j) 
        yzquad( res_node(i) ) = yzquad( res_node(i) ) + yzquad( node_child(j) ) - ydip( node_child(j) )*zs(j) &
             - zdip( node_child(j) )*ys(j) + charge( node_child(j) )*ys(j)*zs(j) 
        zxquad( res_node(i) ) = zxquad( res_node(i) ) + zxquad( node_child(j) ) - zdip( node_child(j) )*xs(j) &
             - xdip( node_child(j) )*zs(j) + charge( node_child(j) )*zs(j)*xs(j) 

        ! magnetic dipole moment
        magmx( res_node(i) ) = magmx( res_node(i) ) + magmx( node_child(j) ) - 0.5*ys(j)*jz( node_child(j) ) &
		 - 0.5*zs(j)*jy( node_child(j) )
        magmy( res_node(i) ) = magmy( res_node(i) ) + magmy( node_child(j) ) - 0.5*zs(j)*jx( node_child(j) ) &
		 - 0.5*xs(j)*jz( node_child(j) )
        magmz( res_node(i) ) = magmz( res_node(i) ) + magmz( node_child(j) ) - 0.5*xs(j)*jy( node_child(j) ) &
		 - 0.5*ys(j)*jx( node_child(j) )

        jx( res_node(i) ) = jx( res_node(i) ) + jx( node_child(j) )       ! Sum currents of child nodes
        jy( res_node(i) ) = jy( res_node(i) ) + jy( node_child(j) )  
        jz( res_node(i) ) = jz( res_node(i) ) + jz( node_child(j) )  

        ! moments for node size:
     	size_node(res_node(i)) = size_node(res_node(i)) + size_node(node_child(j))

     end do
  end do

  ! Should now have multipole information up to branch list level(s).
  ! By definition, this is complete: each branch node is self-contained.
  ! This information has to be broadcast to the other PEs so that the top levels can be filled in.

  call timer_stop(t_props_twigs)
  call timer_stop(t_local)

end subroutine tree_local







