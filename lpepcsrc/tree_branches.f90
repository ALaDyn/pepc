! ===========================================
!
!           TREE_BRANCHES
!
!   Determine branch lists: 
!     smallest set of twig nodes containing all local particles
!     plus boundary particles if necessary.
!   
!
! ===========================================

subroutine tree_branches

  use treevars
  use tree_utils

  implicit none
  include 'mpif.h'

!  integer, parameter :: size_t=1000
  integer*8, dimension(nbranch_max) ::  resolve_key, search_key
  integer*8, dimension(8) :: sub_key   ! Child partial key
  integer, dimension(nbranch_local_max) ::  local_node, local_code, local_leaves  ! local branch data
  integer, dimension(nbranch_max) :: newentry, branch_node, branch_code, branch_leaves  ! global htable data for branches
  integer :: treelevel

  integer*8 ::  keymin, keymax
  integer :: i, j, k, level, nsubset, ncheck, cchild, nchild, newsub, &
       link_addr, ncoll, nres, nbound, newleaf, newtwig, hashaddr
  integer :: nleaf_check, ntwig_check, nleaf_check2, ierr
  logical :: resolved

  ! Global arrays used:
  !   treekey

  integer :: key2addr        ! Mapping function to get hash table address from key
  integer :: startlevel = 2  ! Min permitted branch level

!  branch_debug=.true.
  nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
  ntwig_me = ntwig
  if (tree_debug .and. (proc_debug==me .or.proc_debug==-1)) call check_table('after treebuild     ')

  if (tree_debug) write(ipefile,'(a)') 'TREE BRANCHES'
  if (me==0 .and. tree_debug) then
	write(*,'(a)') 'LPEPC | BRANCHES'
!        write(*,'(a,i8)') 'LPEPC | nbranch_local_max = ',nbranch_local_max
  endif


  ! Determine minimum set of branch nodes making up local domain

  ncheck = 0
  nbranch = 0
  newsub = 0
  level = 1
  nsubset=0

  do i=1,nnodes
     treelevel  = log(1.*treekey(i))/log(8.)     ! node levels

     if (treelevel==1) then
        nsubset = nsubset+1  ! # nodes at level 1
        search_key(nsubset) = treekey(i)   ! Subset of nodes at same level
     endif
  end do

  do while ( ncheck < nleaf )

     call sort( search_key(1:nsubset) )  ! Sort keys

     if (nsubset > 0 ) then
        keymin = minval( search_key(1:nsubset),1 )      ! TODO:  This should include cells that have already been added to
        keymax = maxval( search_key(1:nsubset),1 )      ! to the list.  At present, keymin, keymax could lie in middle of 
        ! domain, causing 'gaps' in branch list
        ! =>  Compute from particle keys for given level instead.

        !        keymin = ishft( pekey(1),-3*(nlev-level) )      ! recover min, max twig keys from particle keys at this level
        !        keymax = ishft( pekey(nlist),-3*(nlev-level) )      ! recover min, max twig keys from particle keys at this level

        do i=1,nsubset
           treelevel  = log(1.*search_key(i))/log(8.)     ! node levels
           if ( (search_key(i) > keymin .and. search_key(i) < keymax .and. treelevel>=startlevel) .or. ( htable( key2addr( search_key(i),'BRANCHES: search' ) )%node > 0 )) then
              !  either middle node (complete twig),  or leaf:  so add to domain list
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i),'BRANCHES: ncheck' ) )%leaves  ! Augment checksum

           else 
              ! end node: check for complete twigs; otherwise subdivide
              cchild = htable( key2addr( search_key(i),'BRANCHES: cchild' ) )%childcode   !  Children byte-code

              nchild = SUM( (/ (ibits(cchild,j,1),j=0,7) /) ) ! # children = sum of bits in byte-code
              sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(cchild,j),j=0,7) /) )  ! Extract sub key from byte code

              do j=1,nchild
                 resolve_key(newsub+j) = IOR( ishft( search_key(i),3 ), sub_key(j) ) ! Construct keys of children
              end do
              newsub = newsub + nchild

           endif
        end do

        if (branch_debug) then
           write (ipefile,'(/a,i7,a,i7/a/(i5,o16,i6,z5,i8))') 'Branches at level:',level, ' Checksum: ',ncheck, &
                '    i      key         node     code     #leaves', &
                (i,search_key(i), &
                htable( key2addr( search_key(i),'BRANCHES: debug' ) )%node, &                         ! Node #
                htable( key2addr( search_key(i),'BRANCHES: debug' ) )%childcode, &                         ! Children byte-code 
                htable( key2addr( search_key(i),'BRANCHES:debug' ) )%leaves, &                           ! # leaves contained in branch 
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
     local_node(i) = htable( key2addr( pebranch(i),'BRANCHES: info' ) )%node                       ! Node #
     local_code(i) =  htable( key2addr( pebranch(i),'BRANCHES: info' ) )%childcode               ! Children byte-code 
     local_leaves(i) = htable( key2addr( pebranch(i),'BRANCHES: info' ) )%leaves             ! # leaves contained in branch 
  end do

  if (ncheck > nleaf) then
     write(*,*) 'Checksum ',ncheck,' /= # leaves on PE ',me
  endif

  if (tree_debug .and. (proc_debug==me .or.proc_debug==-1)) call check_table('after local branches     ')

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize




  ! send copies of branch nodes to all other PEs


  ! first need to find number of branches to be gathered from each PE:
  ! do this by first collecting nbranch values on root.

  call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )


  ! work out stride lengths so that partial arrays placed sequentially in global array

  igap(1) = 0
  igap(2) = nbranches(1)
  do i=3,num_pe
	igap(i) = SUM(nbranches(1:i-1))
  end do
	
! igap(1:num_pe) = (/ 0, nbranches(1), ( SUM( nbranches(1:i-1) ),i=3,num_pe ) /)

  nbranch_sum = SUM( nbranches(1:num_pe) )   ! Total # branches in tree
  igap(num_pe+1) = nbranch_sum
  if (me==0 .and. tree_debug) then
    write(*,'(a50,3i8,a3,2i7)') 'LPEPC | BRANCHES: Local, max local, global # branches, (max): ',nbranch,MAXVAL(nbranches),nbranch_sum,'/',nbranch_local_max,nbranch_max
!    write(*,'(a/(3i8))') 'Branch distribution',(i,nbranches(i),igap(i),i=1,num_pe)
  endif

  ! now collect partial arrays together on each PE
  ! using gather-to-all: nbranches and igap already defined locally

!  nbuf = nbranch
!  recv_counts = nbranches  ! receive buffer lengths and placement
!  recv_strides = igap

  call mpi_allgatherv( pebranch, nbranch, MPI_INTEGER8, branch_key, &
                       nbranches, igap, MPI_INTEGER8, MPI_COMM_WORLD, ierr )

  call mpi_allgatherv( local_node, nbranch, MPI_INTEGER, branch_node, &
                       nbranches, igap,  MPI_INTEGER, MPI_COMM_WORLD, ierr )

  call mpi_allgatherv( local_code, nbranch, MPI_INTEGER, branch_code, &
                       nbranches, igap, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  call mpi_allgatherv( local_leaves, nbranch, MPI_INTEGER, branch_leaves, &
                       nbranches, igap, MPI_INTEGER, MPI_COMM_WORLD, ierr )

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

  if ( branch_debug) then
     write (ipefile,'(2(/a20,i6,a1,i6))') 'New twigs: ',newtwig,'/',ntwig+newtwig, 'New leaves:',newleaf,'/',nleaf+newleaf
     write (ipefile,*) 'Local, total # branches = ',nbranch,nbranch_sum
  endif
  nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
  ntwig_me = ntwig
  nleaf = nleaf + newleaf  ! Total # leaves/twigs in local #table
  ntwig = ntwig + newtwig


  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize



end subroutine tree_branches

