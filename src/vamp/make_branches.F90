! ===========================================
!
!           MAKE_BRANCHES
!
!   Determine branch lists: 
!     smallest set of twig nodes containing all local particles
!     plus boundary particles if necessary.
!   
!
! ===========================================

subroutine make_branches

  use treevars
  use utils
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none

  integer*8, dimension(size_tree) ::  resolve_key, search_key
  integer*8, dimension(8) :: sub_key   ! Child partial key

  integer, dimension(size_tree) ::newentry, local_node, local_code, local_leaves  ! local branch data
  integer, dimension(size_tree) :: branch_node, branch_code, branch_leaves  ! global htable data for branches
  integer, dimension(size_tree) :: treelevel

  integer, dimension(0:maxaddress) :: cell_addr

  integer*8 ::  keymin, keymax
  integer :: i, j, k, level, nsubset, ncheck, cchild, nchild, newsub, &
       link_addr, ncoll, nres, nbound, newleaf, newtwig, hashaddr
  integer :: nleaf_check, ntwig_check, nleaf_check2
  logical :: resolved

  ! Global arrays used:
  !   treekey

  integer :: key2addr        ! Mapping function to get hash table address from key

  nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
!VAMPINST subroutine_start
       CALL VTENTER(IF_make_branches,VTNOSCL,VTIERR)
!      write(*,*) 'VT: make_branches S>',VTIERR,
!     *    IF_make_branches,ICLASSH
!
  ntwig_me = ntwig
  call check_table('after treebuild     ')

  if (branch_debug) write(ipefile,'(///a)') 'BRANCHES'


  ! Determine minimum set of branch nodes making up local domain

  ncheck = 0
  nbranch = 0
  newsub = 0
  level = 1
  treelevel(1:nnodes)  = log(1.*treekey(1:nnodes))/log(2.**idim)     ! node levels

  nsubset = COUNT( mask = treelevel(1:nnodes) == 1)  ! # nodes at level 1
  search_key(1:nsubset) = pack( treekey(1:nnodes), mask = treelevel(1:nnodes) == 1 )   ! Subset of nodes at same level


  do while ( ncheck < nleaf )

     call sort( search_key(1:nsubset) )  ! Sort keys

     if (nsubset > 0 ) then
        keymin = minval( search_key(1:nsubset),1 )      ! TODO:  This should include cells that have already been added to
        keymax = maxval( search_key(1:nsubset),1 )      ! to the list.  At present, keymin, keymax could lie in middle of 
        ! domain, causing 'gaps' in branch list
        ! =>  Compute from particle keys for given level instead.

        !        keymin = ishft( pekey(1),-idim*(nlev-level) )      ! recover min, max twig keys from particle keys at this level
        !        keymax = ishft( pekey(nlist),-idim*(nlev-level) )      ! recover min, max twig keys from particle keys at this level

        do i=1,nsubset
           if ( (search_key(i) > keymin .and. search_key(i) < keymax ) .or. ( htable( key2addr( search_key(i) ) )%node > 0 )) then
              !  either middle node (complete twig),  or leaf:  so add to domain list
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i) ) )%leaves  ! Augment checksum

           else 
              ! end node: check for complete twigs; otherwise subdivide
              cchild = htable( key2addr( search_key(i) ) )%childcode   !  Children byte-code

              nchild = SUM( (/ (ibits(cchild,j,1),j=0,2**idim-1) /) ) ! # children = sum of bits in byte-code
              sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(cchild,j),j=0,7) /) )  ! Extract sub key from byte code

              resolve_key(newsub+1:newsub+nchild) = IOR( ishft( search_key(i),idim ), sub_key(1:nchild) ) ! Construct keys of children
              newsub = newsub + nchild

           endif
        end do

        if (branch_debug) then
           write (ipefile,'(/a,i7,a,i7/a/(i5,o16,i6,z5,i8))') 'Branches at level:',level, ' Checksum: ',ncheck, &
                '    i      key         node     code     #leaves', &
                (i,search_key(i), &
                htable( key2addr( search_key(i) ) )%node, &                         ! Node #
                htable( key2addr( search_key(i) ) )%childcode, &                         ! Children byte-code 
                htable( key2addr( search_key(i) ) )%leaves, &                           ! # leaves contained in branch 
                i=1,nsubset)
        endif

     endif
     level = level + 1
     search_key(1:newsub) = resolve_key(1:newsub)        ! Put children into search list
     nsubset = newsub
     newsub = 0
  end do

  if (branch_debug) write (ipefile,'(/a/(i6,o16))') 'Domain branch list:',(i,pebranch(i),i=1,nbranch)

  ! Extract info about branches

  do i=1, nbranch
     local_node(i) = htable( key2addr( pebranch(i) ) )%node                       ! Node #
     local_code(i) =  htable( key2addr( pebranch(i) ) )%childcode               ! Children byte-code 
     local_leaves(i) = htable( key2addr( pebranch(i) ) )%leaves             ! # leaves contained in branch 
  end do

  if (ncheck > nleaf) then
     write(*,*) 'Checksum ',ncheck,' /= # leaves on PE ',me
  endif

  call check_table('after local branches')

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize




  ! send copies of branch nodes to all other PEs


  ! first need to find number of branches to be gathered from each PE:
  ! do this by first collecting nbranch values on root.

  call mpi_allgather( nbranch, one, MPI_INTEGER, nbranches, one, MPI_INTEGER, MPI_COMM_WORLD, ierr )


  ! work out stride lengths so that partial arrays placed sequentially in global array

  igap(1) = 0
  igap(2) = nbranches(1)
  do i=3,num_pe
	igap(i) = SUM(nbranches(1:i-1))
  end do
	
! igap(1:num_pe) = (/ 0, nbranches(1), ( SUM( nbranches(1:i-1) ),i=3,num_pe ) /)

  nbranch_sum = SUM( nbranches(1:num_pe) )   ! Total # branches in tree
  igap(num_pe+1) = nbranch_sum


  ! now collect partial arrays together on each PE
  ! using gather-to-all: nbranches and igap already defined locally

  nbuf = nbranch
  recv_counts = nbranches  ! receive buffer lengths and placement
  recv_strides = igap

  call mpi_allgatherv( pebranch, nbuf, MPI_INTEGER8, branch_key, &
                       recv_counts, recv_strides, MPI_INTEGER8, MPI_COMM_WORLD, ierr )

  call mpi_allgatherv( local_node, nbuf, MPI_INTEGER, branch_node, &
                       recv_counts, recv_strides,  MPI_INTEGER, MPI_COMM_WORLD, ierr )

  call mpi_allgatherv( local_code, nbuf, MPI_INTEGER, branch_code, &
                       recv_counts, recv_strides, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  call mpi_allgatherv( local_leaves, nbuf, MPI_INTEGER, branch_leaves, &
                       recv_counts,recv_strides , MPI_INTEGER, MPI_COMM_WORLD, ierr )

  !  need to keep track of owners too
  do i=1,num_pe
     branch_owner(igap(i)+1:igap(i+1) ) = i-1
  end do

  if ( branch_debug .and. me == 0) then
     write (ipefile,'(/a/a/(2i5,o16,i12,z5,i10))') 'Global branch list:', &
          '   #       owner     key         node     code     #leaves', &
          (i,branch_owner(i), branch_key(i), branch_node(i), branch_code(i), branch_leaves(i),i=1,nbranch_sum)

  endif

  if (branch_debug) then
     nleaf_check = count(mask = branch_owner(1:nbranch_sum)/=me .and. branch_node(1:nbranch_sum) >0 )
     nleaf_check2 = count(mask = branch_owner(1:nbranch_sum)/=me .and. branch_leaves(1:nbranch_sum) ==1 )
     ntwig_check = count(mask = branch_owner(1:nbranch_sum)/=me .and. branch_node(1:nbranch_sum) <0 )
     write(ipefile,*) '# branch twigs to be added: ',ntwig_check,' node>0 ',nleaf_check,' branch_leaves:',nleaf_check2 
  endif

  ! Create entries in #table for remote branch nodes

  nres = 0
  do i = 1,nbranch_sum

     if (branch_owner(i) /= me) then
        call make_hashentry( branch_key(i), branch_node(i), branch_leaves(i), branch_code(i), branch_owner(i), hashaddr, ierr )

        nres = nres + 1            ! Count new entries
        newentry(nres) = hashaddr
     endif
  end do


  if (branch_debug ) write (ipefile,'(/a,4(/a20,i6))') 'Check on resolved nodes ','Total # branches: ',nbranch_sum, &
       'local:',nbranch,'new in #table:',nres

  newleaf = 0
  newtwig = 0

  ! Go through new nodes in table and label as leaf/twig
  do i = 1,nres

     if ( htable( newentry(i) )%leaves == 1 ) then
        ! leaf
        newleaf = newleaf + 1 
        htable( newentry(i) )%node = nleaf+newleaf    ! Create local label for remote branch

     else if ( htable( newentry(i) )%leaves > 1 ) then
        ! twig
        newtwig = newtwig + 1
        htable( newentry(i) )%node = -ntwig-newtwig    ! Create local label for remote branch

     else
        ! empty - problem with flagging
     endif

  end do

  if ( branch_debug) write (ipefile,'(2(/a20,i6,a1,i6))') 'New twigs: ',newtwig,'/',ntwig+newtwig, &
       'New leaves:',newleaf,'/',nleaf+newleaf

  nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
  ntwig_me = ntwig
  nleaf = nleaf + newleaf  ! Total # leaves/twigs in local #table
  ntwig = ntwig + newtwig


  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize



!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: make_branches S<',VTIERR,ICLASSH
!
end subroutine make_branches

