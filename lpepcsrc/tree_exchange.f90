subroutine tree_exchange

  use treevars
  use timings
  implicit none
  include 'mpif.h'

  real*8 :: ts1b=0., ts1e=0., ta1b=0., ta1e=0.

  integer :: i,ierr, lnode, lcode, lleaves, hashaddr, newleaf, newtwig
  integer*8 :: lnext

  real*8, dimension(nbranch_max) :: pack_size, get_size


  integer,external :: key2addr        ! Mapping function to get hash table address from key
   
  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  if (tree_debug) write(ipefile,'(a)') 'TREE EXCHANGE'
  if (me==0 .and. tree_debug) then
	write(*,'(a)') 'LPEPC | EXCHANGE'
  endif

  ! Pack local branches for shipping
  do i=1,nbranch
     lnode = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%node
     lcode = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%childcode 
     lleaves = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%leaves 
     lnext = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%next
     pack_child(i) = multipole(pebranch(i),lcode,lleaves,me,lnext,&
          charge( lnode ), &
          abs_charge( lnode ), &
          xcoc( lnode), &
          ycoc( lnode), &
          zcoc( lnode), &
          xdip( lnode), &
          ydip( lnode), &
          zdip( lnode), &
          xxquad( lnode), &
          yyquad( lnode), &
          zzquad( lnode), &
          xyquad( lnode), &
          yzquad( lnode), &
          zxquad( lnode), &
          jx( lnode), &
          jy( lnode), &
          jz( lnode), &
          magmx( lnode), &
          magmy( lnode), &
          magmz( lnode) )
     pack_size(i) = size_node(lnode)
  end do
  
  call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  ! work out stride lengths so that partial arrays placed sequentially in global array
  igap(1) = 0
  igap(2) = nbranches(1)
  do i=3,num_pe
	igap(i) = SUM(nbranches(1:i-1))
  end do
  
  nbranch_sum = SUM( nbranches(1:num_pe) )   ! Total # branches in tree
  igap(num_pe+1) = nbranch_sum
  if (me==0 .and. tree_debug) then
    write(*,'(a50,3i8,a3,2i7)') 'LPEPC | BRANCHES: Local, max local, global # branches, (max): ', &
	nbranch,MAXVAL(nbranches),nbranch_sum,'/',nbranch_local_max,nbranch_max
!    write(*,'(a/(3i8))') 'Branch distribution',(i,nbranches(i),igap(i),i=1,num_pe)
  endif

  if (nbranch_sum > size_fetch) then
     write(*,*) 'Too many branches for buffer on ',me,': nbranch_sum=',nbranch_sum,' size_fetch=',size_fetch
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop
  end if
  
  call MPI_ALLGATHERV(pack_child, nbranch, MPI_TYPE_MULTIPOLE, get_child, nbranches, igap, MPI_TYPE_MULTIPOLE, MPI_COMM_WORLD, ierr)
  call MPI_ALLGATHERV(pack_size, nbranch, MPI_REAL8, get_size, nbranches, igap, MPI_REAL8, MPI_COMM_WORLD, ierr)

  newleaf = 0
  newtwig = 0

  ! Integrate remote branches into local tree
  do i = 1,nbranch_sum

     branch_key(i) = get_child(i)%key

     if (get_child(i)%owner /= me) then

        if ( get_child(i)%leaves == 1 ) then
           ! leaf
           newleaf = newleaf + 1 
           lnode = nleaf+newleaf    ! Create local label for remote branch         
        else if ( get_child(i)%leaves > 1 ) then
           ! twig
           newtwig = newtwig + 1
           lnode = -ntwig-newtwig    ! Create local label for remote branch
        else
           write(*,*) 'Problem with flagging on remote branches',me
           call MPI_ABORT(MPI_COMM_WORLD,ierr)
           stop       
        endif
        call make_hashentry( get_child(i)%key, lnode , get_child(i)%leaves, get_child(i)%byte, get_child(i)%owner, hashaddr, ierr )
 
        charge( lnode ) = get_child(i)%q
        abs_charge( lnode ) = get_child(i)%absq
        xcoc( lnode ) = get_child(i)%xcoc
        ycoc( lnode ) = get_child(i)%ycoc
        zcoc( lnode ) = get_child(i)%zcoc
        xdip( lnode ) = get_child(i)%xdip
        ydip( lnode ) = get_child(i)%ydip
        zdip( lnode ) = get_child(i)%zdip
        xxquad( lnode ) = get_child(i)%xxquad
        yyquad( lnode ) = get_child(i)%yyquad
        zzquad( lnode ) = get_child(i)%zzquad
        xyquad( lnode ) = get_child(i)%xyquad
        yzquad( lnode ) = get_child(i)%yzquad
        zxquad( lnode ) = get_child(i)%zxquad
        magmx( lnode ) = get_child(i)%magmx
        magmy( lnode ) = get_child(i)%magmy
        magmz( lnode ) = get_child(i)%magmz
        jx( lnode ) = get_child(i)%jx
        jy( lnode ) = get_child(i)%jy
        jz( lnode ) = get_child(i)%jz

     	size_node( lnode ) =  get_size(i)
       
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
  t_exchange = ts1e - ts1b

end subroutine tree_exchange
