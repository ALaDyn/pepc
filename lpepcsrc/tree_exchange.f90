subroutine tree_exchange

  use treevars
  use timings
  use module_htable
  implicit none
  include 'mpif.h'

  integer :: i,ierr, lnode, lcode, lleaves, hashaddr, newleaf, newtwig
  integer*8 :: lnext
  
  real*8, dimension(nbranch) :: pack_size
  type (multipole), dimension(nbranch) :: pack_mult
  real*8, allocatable :: get_size(:)
  type (multipole),allocatable :: get_mult(:)

  call timer_start(t_exchange_branches)

  if (tree_debug) write(ipefile,'(a)') 'TREE EXCHANGE'
  if (me==0 .and. tree_debug) then
	write(*,'(a)') 'LPEPC | EXCHANGE'
  endif

 
  call timer_start(t_exchange_branches_pack)

  ! Pack local branches for shipping
  do i=1,nbranch
     lnode = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%node
     lcode = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%childcode 
     lleaves = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%leaves 
     lnext = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%next
     pack_mult(i) = multipole(pebranch(i),lcode,lleaves,me,lnext,&
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
          xshift( lnode), &
          yshift( lnode), &
          zshift( lnode) )
     pack_size(i) = size_node(lnode)
  end do

  call timer_stop(t_exchange_branches_pack)
  call timer_start(t_exchange_branches_admininstrative)

  call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  ! work out stride lengths so that partial arrays placed sequentially in global array
  igap(1) = 0
  igap(2) = nbranches(1)
  do i=3,num_pe
	igap(i) = SUM(nbranches(1:i-1))
  end do
  
  nbranch_sum = SUM( nbranches(1:num_pe) )   ! Total # branches in tree
  allocate(get_size(1:nbranch_sum),get_mult(1:nbranch_sum))
  igap(num_pe+1) = nbranch_sum

  if (nbranch_sum > nbranch_max) then
     write(*,*) 'Too many branches for buffer on ',me,': nbranch_sum=',nbranch_sum,' nbranch_max=',nbranch_max
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     stop
  end if
 
  call timer_stop(t_exchange_branches_admininstrative)
  call timer_start(t_exchange_branches_allgatherv)

  call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_MULTIPOLE, get_mult, nbranches, igap, MPI_TYPE_MULTIPOLE, MPI_COMM_WORLD, ierr)
  call MPI_ALLGATHERV(pack_size, nbranch, MPI_REAL8, get_size, nbranches, igap, MPI_REAL8, MPI_COMM_WORLD, ierr)

  call timer_stop(t_exchange_branches_allgatherv)
  call timer_start(t_exchange_branches_integrate)

  newleaf = 0
  newtwig = 0

  ! Integrate remote branches into local tree
  do i = 1,nbranch_sum

     branch_key(i) = get_mult(i)%key

     if (get_mult(i)%owner /= me) then

        if ( get_mult(i)%leaves == 1 ) then
           ! leaf
           newleaf = newleaf + 1 
           lnode = nleaf+newleaf    ! Create local label for remote branch         
        else if ( get_mult(i)%leaves > 1 ) then
           ! twig
           newtwig = newtwig + 1
           lnode = -ntwig-newtwig    ! Create local label for remote branch
        else
           write(*,*) 'Problem with flagging on remote branches',me
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           stop       
        endif
        call make_hashentry( get_mult(i)%key, lnode , get_mult(i)%leaves, get_mult(i)%byte, get_mult(i)%owner, hashaddr, ierr )
 
        htable(hashaddr)%childcode = IBSET( htable(hashaddr)%childcode, CHILDCODE_NODE_TOUCHED ) ! I have touched this node, do not zeor its properties (in tree_global)
        charge( lnode ) = get_mult(i)%q
        abs_charge( lnode ) = get_mult(i)%absq
        xcoc( lnode ) = get_mult(i)%xcoc
        ycoc( lnode ) = get_mult(i)%ycoc
        zcoc( lnode ) = get_mult(i)%zcoc
        xdip( lnode ) = get_mult(i)%xdip
        ydip( lnode ) = get_mult(i)%ydip
        zdip( lnode ) = get_mult(i)%zdip
        xxquad( lnode ) = get_mult(i)%xxquad
        yyquad( lnode ) = get_mult(i)%yyquad
        zzquad( lnode ) = get_mult(i)%zzquad
        xyquad( lnode ) = get_mult(i)%xyquad
        yzquad( lnode ) = get_mult(i)%yzquad
        zxquad( lnode ) = get_mult(i)%zxquad
        xshift( lnode ) = get_mult(i)%xshift
        yshift( lnode ) = get_mult(i)%yshift
        zshift( lnode ) = get_mult(i)%zshift
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
  
  deallocate(get_size,get_mult)

  call timer_stop(t_exchange_branches_integrate)
  call timer_stop(t_exchange_branches)

end subroutine tree_exchange
