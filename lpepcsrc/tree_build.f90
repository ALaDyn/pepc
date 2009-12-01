! ===========================================
!
!           TREE_BUILD
!
!   Construct tree using particle keys
!   
!   $Revision$
! ===========================================

subroutine tree_build

  use treevars
  use timings
  implicit none
  include 'mpif.h'

  real*8 :: ts1b=0., ts1e=0., ta1b=0., ta1e=0.

  integer*8, dimension(maxaddress) :: subcell, par_key, res_key        ! All key arrays 64-bit
  integer*8, dimension(nppm+2) :: local_key

 ! integer*8, dimension(2*nppm) :: leaf_key, twig_key
  integer*8, dimension(8) :: sub_key   ! Child partial key

  integer, dimension(maxaddress) :: ix, iy, res_addr, res_node, res_child, res_owner, &
        newentry, treelevel
  integer, dimension(nppm+2)  :: local_plist, local_ind, local_owner
  integer, dimension(0:maxaddress) ::  cell_addr
  logical, dimension(nppm+2) :: part_done            ! finished flag for leaf nodes

  integer*8 ::  parent_key,  key_lo, cell1, cell2

  integer :: i, j, k, nres, ncoll, link, link_addr, & 
       level, ilev, ibit, newleaf, nbound, nlistnew, ipart, ilist, ipoint, &
       childbyte, parent_addr, &
       level_top, level_match, level_diff, iend, &
       i1, i2, ierr

  integer :: key2addr        ! Mapping function to get hash table address from key


  character*1 :: collision(0:maxaddress)
  real :: s
  logical :: resolved

  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  if (tree_debug) write(ipefile,'(/a)') 'TREE BUILD'
  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | LOCAL BUILD'

  ! zero table: need list of 'live' addresses to speed up
  htable%node = 0
  htable%leaves = 0
  htable%childcode = 0
  htable%link = -1
  htable%key = 0_8

  do i=1,npp+2
     local_plist(i) = pelabel(i)       ! Particle (global) label for tracking purposes 
     local_key(i) = pekey(i)           ! Particle key
     local_owner(i) = pepid(i)
  end do
  i=0                     ! List of free (unused #table-addresses)


  collision = " "                               ! Collision flag

  !  Check whether boundary particles are close enough to warrant inclusion
  !  in tree construction.

  nlist = npp
!  level_top = log(1.*npp) / log(8.)   ! Min level for 'same cell' test 
 level_top = 1
! RH neighbour PE

  if ( me /= num_pe-1 ) then
! First find level shared by last particle pair in list
     key_lo = ieor( pekey(npp), pekey(npp-1) )   ! Picks out 1st position where keys differ
     level_diff =  log(1.*key_lo)/log(8.) 
     level_match = max( nlev - level_diff, level_top )    ! Excludes low-level end-cells for discontinuous domains
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
     level_diff =  log(1.0*key_lo)/log(8.) 

     level_match = max( nlev - level_diff, level_top )    ! Excludes low-level end-cells for discontinuous domains
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

  ta1e = MPI_WTIME()  
  t_build_neigh = ta1e-ta1b
  ta1b = MPI_WTIME()

  level = 0                        ! start at root
  local_ind(1:nlist) = (/(k,k=1,nlist)/) ! Sorted local particle index / leaf number: points to particle properties
  part_done(1:nlist) = .false.              ! flag all particles not finished

  !  set up root node

  htable(1)%node = -1                   !  node #
  htable(1)%owner = me                   ! Owner
  htable(1)%key = 1_8                    !  key
  htable(1)%link = -1                   !  collision link
  htable(1)%leaves = npp              ! root contains all leaves, excluding boundary particles
  ntwig = 1
  nleaf = 0
  tablehigh = 0



  do while ( nlist > 0 )              ! While any particle not finished:

     level = level + 1                 ! Next sublevel
     if (level>nlev) then
       write(*,*) 'Problem with tree on PE ',me,' - no more levels '
       write(*,'(a/(i8,o30))') 'Remaining keys: ',(local_plist(i),local_key(i),i=1,nlist)
!       call closefiles
       call MPI_ABORT(MPI_COMM_WORLD,ierr)
       stop
     endif
     ibit = nlev - level               ! bit shift factor (0=highest leaf node, nlev-1 = root)
     newleaf=0                         ! # new leaves at this level
     nbound=0                            ! # boundary particles located (max 2)
     ! Determine subcell # from key
     ! At a given level, these will be unique
     
     do i=1,nlist
        subcell(i) = ishft( local_key(i), -3_8*ibit )    

        ! cell address hash function
        cell_addr(i) = IAND( subcell(i), hashconst)
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
           htable(cell_addr(i))%key = subcell(i)            !   key
           htable(cell_addr(i))%leaves = 1                !   # children
           htable(cell_addr(i))%owner = local_owner(i)       ! Set owner of node equal to owner of particle

        else if ( htable( cell_addr(i) )%node /= 0 .AND. &   ! Entry exists and keys match
             htable(cell_addr(i))%key == subcell(i) ) then       ! => twig cell  

           htable(cell_addr(i))%node = i               ! Overwrite previous (leaf) entry
           htable(cell_addr(i))%leaves =  htable(cell_addr(i))%leaves + 1
           htable(cell_addr(i))%owner = me       ! Set owner of node equal to local PE

        else                                    ! Entry exists and keys don't match
           ncoll = ncoll + 1                    ! Increment collision count
           res_addr(ncoll) = cell_addr(i)        ! Reduced list of addresses 
           res_key(ncoll) = subcell(i)           !     ..    .. of keys
           res_node(ncoll) = i                   !     ..    .. of pointers
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

           link_addr = res_addr(i)                          ! Start of chain
           resolved = .false.                                     ! Resolve flag

           do while (  .not. resolved )       ! Find end of chain  
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
           ipart = pelabel( local_ind( ipoint ) )                ! particle label
           htable( newentry(i) )%leaves = 1                      ! contained leaves
           htable( newentry(i) )%childcode = ipart               ! store label in #-table
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
            htable( newentry(i) )%key = -1_8               ! but retain %link to dummy entry 
	                                                 ! in case it's in the middle of a chain
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
           ! empty - problem with flagging
        endif

     end do


     ! output interim hash table
     if (build_debug) then
        write (ipefile,'(/a,i6/8x,a/8x,a)') 'Hash table at level',level, &
             'entry,    node,    key,     link  children  collision', &
             '----------------------------------------------------- '

        ! flag old & new collisions

        do i=0,tablehigh
           if ( collision(i) == "X") collision(i) = "c" 
           if (htable(i)%node/=0 .and. htable(i)%link/= -1 .and. collision(i)/="c") collision(i)="X"
           if (htable(i)%node /=0 ) write (ipefile,'(2i10,o10,2i10,4x,a1)') &
                i,htable(i)%node,htable(i)%key,htable(i)%link,htable(i)%leaves,collision(i)
        end do
     endif



     ! Make new lists from unfinished particles.
     nlistnew = nlist - newleaf - nbound
     part_done(1:nlist) = (local_plist(1:nlist)==0)
      k=0
     do i=1,nlist
     	if (.not. part_done(i)) then
	  k=k+1
          local_key(k) = local_key(i)
          local_ind(k) = local_ind(i)
	  local_plist(k) = local_plist(i)
	  local_owner(k) = local_owner(i)		  
	endif

     end do

     
     nlist = nlistnew

     if (build_debug) then
        !  Go through list of unresolved entries
        write (ipefile,*) 'So far: ',nleaf,' leaves, ',newleaf, 'new on level,',ntwig,' twigs, ',ncoll,' collisions'
        write (ipefile,*) 'Particles left: ',nlist,', new list:',local_plist(1:nlist)
     endif

  end do

  ta1e = MPI_WTIME()  
  t_build_part = ta1e-ta1b
  ta1b = MPI_WTIME()  

  ! Fill in byte code for children in table (twig entries contain # children at this point)

  nnodes = nleaf + ntwig  ! total number of local tree nodes


  treekey(1:nleaf) = leaf_key(1:nleaf)
  twig_key(ntwig) = 1  ! add root for later
  treekey(nleaf+1:nnodes) = twig_key(1:ntwig)  ! add twigs to list 
!  write (*,*) 'Keys on :',me,nnodes,(treekey(i),i=1,nnodes)

  do i=1,nnodes-1
     childbyte = IAND( treekey(i), hashchild)    ! extract last 3 bits from key
     parent_key = ishft( treekey(i),-3 )      ! parent key
     parent_addr = key2addr(parent_key,'BUILD: childbyte')
     ! Construct children byte-code (8 settable bits indicating which children nodes present)
     htable(parent_addr)%childcode = ibset( htable(parent_addr)%childcode, childbyte )
! This is equivalent to:
!     htable(parent_addr)%childcode =  htable(parent_addr)%childcode + 2**(childbyte)
  end do


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

 ta1e = MPI_WTIME()  
 t_build_byte = ta1e-ta1b
 ts1e = MPI_WTIME()
 call MPI_REDUCE(ts1e-ts1b,t_build,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
end subroutine tree_build







