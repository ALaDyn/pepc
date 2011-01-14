!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Determime minimum set of local branches and exchange them
!> through all cores.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tree_branches

  use treevars
  use tree_utils
  use timings
  use module_math_tools

  implicit none
  include 'mpif.h'

  ! Timing stuff
  real*8 :: ts1b=0., ts1e=0., ta1b=0., ta1e=0.

  ! Search lists
  integer*8, dimension(nbranch_max) ::  resolve_key, search_key
  integer*8, dimension(8) :: sub_key   ! Child partial key
 
  ! local branch data
  integer, dimension(nbranch_local_max) ::  local_node, local_code, local_leaves

  ! global htable data for branches
  integer, dimension(nbranch_max) :: newentry, branch_node, branch_code, branch_leaves
  integer :: treelevel

  integer :: i, j, k, level, nsubset, ncheck, cchild, nchild, newsub, &
       link_addr, ncoll, nres, nbound, newleaf, newtwig, hashaddr
  integer :: nleaf_check, ntwig_check, nleaf_check2, ierr
  logical :: resolved

  ! stuff for getting a virtual domain
  integer*8 :: right_limit_me, right_limit
  integer*8 :: left_limit_me, left_limit
  integer*8 :: left_virt_limit, right_virt_limit
  integer*8 :: left_cell, right_cell
  integer*8 :: search_key_mod

  ! stuff for estimation
  integer*8 :: branch_level(0:nlev) 
  integer*8 :: branch_level_D1(1:nlev) ! start at level 1
  integer*8 :: branch_level_D2(1:nlev) 
  integer*8 :: D1, D2               ! sub-domains
  integer*8 :: L                    ! inner limit
  integer*8 :: branch_max_local     ! estimation for local branches
  integer*8 :: branch_max           ! global estimation
  integer :: ilevel, pos
  
  ! Mapping function to get hash table address from key
  integer :: key2addr



  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  !  Retain leaves and twigs belonging to local PE
  nleaf_me = nleaf      
  ntwig_me = ntwig
  
  ! Debugging stuff
  if (tree_debug .and. (proc_debug==me .or.proc_debug==-1)) call check_table('after treebuild     ')
  if (tree_debug) write(ipefile,'(a)') 'TREE BRANCHES'
  if (me==0 .and. tree_debug) then
	write(*,'(a)') 'LPEPC | BRANCHES'
  endif

  ! get local key limits
  left_limit_me=pekey(1)
  right_limit_me=pekey(npp)

  ! get key limits for neighbor PE's
  ! and build virtual limits, so that a minimum set a branch nodes comes arround
  ! boundary PE's can access their boundary space fully only need one virtual limit
  if(me.eq.0)then
     right_limit=pekey(npp+1)
     right_virt_limit=bpi_bits(right_limit,right_limit_me,8_8,nlev)
  else if(me.eq.(num_pe-1))then
     left_limit=pekey(npp+1)
     left_virt_limit=bpi_bits(left_limit,left_limit_me,8_8,nlev)
  else
     left_limit=pekey(npp+2)
     right_limit=pekey(npp+1)
     left_virt_limit=bpi_bits(left_limit,left_limit_me,8_8,nlev)
     right_virt_limit=bpi_bits(right_limit,right_limit_me,8_8,nlev)
  end if


  ! Make tough estimation for amount of branches
  branch_level(0)=1;              ! root
  
  ! First find highest power in the Virtual Domain
  L=bpi_bits(left_virt_limit,right_virt_limit,8_8,nlev)
  
  ! divide in two sub-domains
  D1 = L-left_virt_limit
  D2 = right_virt_limit-L+1
  
  ! get estimation: number of branches at all levels
  do ilevel=1,nlev
     pos=3*(nlev-ilevel)
     branch_level_D1(ilevel)=ibits(D1,pos,3)
     branch_level_D2(ilevel)=ibits(D2,pos,3)
     branch_level(ilevel)=branch_level_D1(ilevel)+branch_level_D2(ilevel)
  end do

  ! estimate local number
  branch_max_local = SUM(branch_level(1:nlev))

  ! get global estimation
  Call MPI_ALLREDUCE(branch_max_local, branch_max, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr) 
  

  ! adapt left Virtual Limit 
  if(me.ne.0)then
     left_virt_limit=left_virt_limit-1
  end if

  
  ! Determine minimum set of branch nodes making up local domain
  ncheck = 0
  nbranch = 0
  newsub = 0
  level = 1
  nsubset=0

  ! find all nodes at level 1
  !do i=1,nnodes
  !   treelevel  = log(1.*treekey(i))/log(8.) ! node levels

     !if (treelevel==1) then
      !  nsubset = nsubset+1                  ! # nodes at level 1
       ! search_key(nsubset) = treekey(i)     ! Subset of nodes at same level
     !endif
  !end do

  ! Artificial Startup: mimic first level cells + placeholder bit and check the hashtable
  ! if the key is there. This is done instead of a loop over all nodes where the level were 
  ! checked. Therefore we have O(1) instead of O(N/P). No Sorting is needed because the artificial 
  ! startup is built far-sighted and the keys are already sorted. Note: If the first step is sorted 
  ! the next steps are sorted too. => Sort only first level keys.
  do i=0,7
     if(htable((8+i))%node.ne.0.or.htable((8+i))%key.eq.-1)then ! entry exists
        nsubset = nsubset+1          ! # nodes at level 1
        search_key(nsubset) = (8+i)  ! Subset of nodes at same level
     end if
  end do


  ! while any particle is not in a branch
  do while ( ncheck < nleaf )

     ! found a candidate
     if (nsubset > 0 ) then
        
        ! calculate parent cell for this level from virtual limits 
        left_cell = ibits(left_virt_limit,nlev*3-3*level,3*level)
        right_cell = ibits(right_virt_limit,nlev*3-3*level,3*level)
           
        ! for all nodes at this level
        do i=1,nsubset
           
           ! Important: discount placeholder-bit from search_key-entries
           search_key_mod=search_key(i)-2**(3*level)
           
           ! check if key between limits
           ! Important: discount placeholder-bit from search_key-entries
           if ( (me.eq.0 .and. (search_key_mod < right_cell)) .or. &
		( htable( key2addr( search_key(i),'BRANCHES: search' ) )%node > 0 )) then
              ! twig node between limits or
              ! node is a leaf
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i),'BRANCHES: ncheck' ) )%leaves  ! Augment checksum

           elseif ( (me.eq.num_pe-1 .and. (search_key_mod > left_cell)) .or. &
		( htable( key2addr( search_key(i),'BRANCHES: search' ) )%node > 0 )) then
              ! twig node between limits or
              ! node is a leaf
              nbranch = nbranch + 1
              pebranch(nbranch) = search_key(i)
              ncheck = ncheck +  htable( key2addr( search_key(i),'BRANCHES: ncheck' ) )%leaves  ! Augment checksum

           elseif (  (search_key_mod > left_cell) .and. (search_key_mod < right_cell) .or. &
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

  ! Extract info about branches
  do i=1, nbranch
     local_node(i) = htable( key2addr( pebranch(i),'BRANCHES: info' ) )%node        ! Node #
     local_code(i) =  htable( key2addr( pebranch(i),'BRANCHES: info' ) )%childcode               ! Children byte-code 
     local_leaves(i) = htable( key2addr( pebranch(i),'BRANCHES: info' ) )%leaves             ! # leaves contained in branch 
  end do

  if (ncheck > nleaf) then
     write(*,*) 'Checksum ',ncheck,' /= # leaves on PE ',me
  endif

  if (tree_debug .and. (proc_debug==me .or.proc_debug==-1)) call check_table('after local branches     ')

  ta1e = MPI_WTIME()
  t_branches_find = ta1e-ta1b  
  ta1b = MPI_WTIME()

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
    write(*,'(a50,3i8,a3,2i7)') 'LPEPC | BRANCHES: Local, max local, global # branches, (max): ', &
	nbranch,MAXVAL(nbranches),nbranch_sum,'/',nbranch_local_max,nbranch_max
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

  ta1e = MPI_WTIME()
  t_branches_exchange = ta1e-ta1b  
  ta1b = MPI_WTIME()

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

  ta1e = MPI_WTIME()
  t_branches_integrate = ta1e-ta1b  
  ts1e = MPI_WTIME()
  t_branches = ts1e-ts1b


end subroutine tree_branches

