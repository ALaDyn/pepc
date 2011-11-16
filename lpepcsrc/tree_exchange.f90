subroutine tree_exchange

  use treevars
  use timings
  use module_htable
  use module_branching
  implicit none
  include 'mpif.h'

  integer :: i,ierr, lnode, lcode, lleaves, hashaddr, newleaf, newtwig
  
  type (t_tree_node), target, dimension(nbranch) :: pack_mult
  type (t_tree_node), pointer :: packm
  type (t_tree_node),allocatable :: get_mult(:)
  integer, allocatable :: igap(:)    !  stride lengths of local branch arrays

  call timer_start(t_exchange_branches)

  if (tree_debug) write(ipefile,'(a)') 'TREE EXCHANGE'
  if (me==0 .and. tree_debug) then
	write(*,'(a)') 'LPEPC | EXCHANGE'
  endif

 
  call timer_start(t_exchange_branches_pack)

  ! Pack local branches for shipping
  do i=1,nbranch
     lnode   = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%node
     lcode   = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%childcode
     lleaves = htable( key2addr( pebranch(i),'EXCHANGE: info' ) )%leaves
     packm=>pack_mult(i)
         packm        = tree_nodes( lnode )
         packm%key    = pebranch(i)   ! TODO: this data is maybe not consistently stored in tree_nodes array
         packm%byte   = lcode  ! therefore, we have to take it directly form the htable --> repair this
         packm%leaves = lleaves
         packm%owner  = me
  end do

  call timer_stop(t_exchange_branches_pack)
  call timer_start(t_exchange_branches_admininstrative)

  call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  ! work out stride lengths so that partial arrays placed sequentially in global array
  allocate (igap(num_pe+3))

  igap(1) = 0
  igap(2) = nbranches(1)
  do i=3,num_pe
	igap(i) = SUM(nbranches(1:i-1))
  end do
  
  nbranch_sum = SUM( nbranches(1:num_pe) )   ! Total # branches in tree
  allocate(get_mult(1:nbranch_sum))
  igap(num_pe+1) = nbranch_sum

  if (nbranch_sum > branch_max_global) then
     write(*,*) 'Too many branches for buffer on ',me,': nbranch_sum=',nbranch_sum,' branch_max_global=',branch_max_global
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     stop
  end if
 
  call timer_stop(t_exchange_branches_admininstrative)
  call timer_start(t_exchange_branches_allgatherv)

  call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_MULTIPOLE, get_mult, nbranches, igap, MPI_TYPE_MULTIPOLE, MPI_COMM_WORLD, ierr)

  deallocate (igap)


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

        !insert received data into local tree
        tree_nodes( lnode )      = get_mult(i)
        tree_nodes( lnode )%byte = htable(hashaddr)%childcode ! TODO: otherwise maybe inconsistent with htable data
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
  
  deallocate(get_mult)

  call timer_stop(t_exchange_branches_integrate)
  call timer_stop(t_exchange_branches)

end subroutine tree_exchange
