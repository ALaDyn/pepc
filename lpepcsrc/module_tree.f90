!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all tree specific helper routines and data fields
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_tree
      implicit none
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public tree_insert_node
      public tree_exchange

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Inserts a given tree node into the next free position in the tree ( -(ntwig+1) or (nleaf+1) )
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine tree_insert_node(tree_node)
        use treevars
        use treetypes
        use module_htable
        implicit none
        include 'mpif.h'

        type(t_tree_node), intent(in) :: tree_node
        integer :: hashaddr, lnode, ierr

        if ( tree_node%leaves == 1 ) then
           nleaf =  nleaf + 1
           lnode =  nleaf
        else if ( tree_node%leaves > 1 ) then
           ! twig
           ntwig =  ntwig + 1
           lnode = -ntwig
        else
           write(*,*) 'Problem with flagging on remote branches', me
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

        call make_hashentry( tree_node%key, lnode , tree_node%leaves, IBSET( tree_node%byte, CHILDCODE_NODE_TOUCHED ), tree_node%owner, hashaddr, ierr )

        if (ierr == 1) then
          write(*,*) 'PE', me, ': Trying to insert already existing node into tree. Updating multipole data instead'
        endif

        !insert received data into local tree
        tree_nodes( lnode ) = tree_node%m

      end subroutine tree_insert_node


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Exchanges tree nodes that are given in local_branch_keys with remote PEs
      !> incoming tree nodes are inserted into tree_nodes array and htable, but the
      !> tree above these nodes is not corrected
      !> outputs keys of new(and own) htable/tree_node entries in branch_keys
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine tree_exchange(local_branch_keys, nbranch, branch_keys, nbranch_sum)

          use treevars, only : me, tree_debug, ipefile, num_pe, tree_nodes, nbranches
          use treetypes
          use timings
          use module_htable
          implicit none
          include 'mpif.h'

          integer*8, intent(in) :: local_branch_keys(1:nbranch)
          integer, intent(in) :: nbranch
          integer*8, intent(inout), allocatable :: branch_keys(:)
          integer, intent(out) :: nbranch_sum

          integer :: i,ierr
          type( t_hash ), pointer :: hbranch
          type (t_tree_node),allocatable :: pack_mult(:), get_mult(:)
          integer, allocatable :: igap(:)    !  stride lengths of local branch arrays

          if (allocated(branch_keys)) deallocate(branch_keys)

          call timer_start(t_exchange_branches)
          call timer_start(t_exchange_branches_pack)

          if (tree_debug) then
              write(ipefile,'(a)') 'TREE EXCHANGE'
              if (me==0) write(*,'(a)') 'LPEPC | EXCHANGE'
          endif

          ! Pack local branches for shipping
          allocate(pack_mult(nbranch))
          do i=1,nbranch
              hbranch      => htable( key2addr( local_branch_keys(i),'EXCHANGE: info' ) )
              pack_mult(i) =  t_tree_node( local_branch_keys(i), hbranch%childcode, hbranch%leaves, me, tree_nodes(hbranch%node) )
          end do

          call timer_stop(t_exchange_branches_pack)
          call timer_start(t_exchange_branches_admininstrative)

          call mpi_allgather( nbranch, 1, MPI_INTEGER, nbranches, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

          ! work out stride lengths so that partial arrays placed sequentially in global array
          allocate (igap(num_pe+3))

          igap(1) = 0
          do i=2,num_pe+1
              igap(i) = igap(i-1) + nbranches(i-1)
          end do

          nbranch_sum = igap(num_pe+1)

          allocate(get_mult(1:nbranch_sum), branch_keys(1:nbranch_sum))

          call timer_stop(t_exchange_branches_admininstrative)
          call timer_start(t_exchange_branches_allgatherv)

          ! actually exchange the branch nodes
          call MPI_ALLGATHERV(pack_mult, nbranch, MPI_TYPE_tree_node, get_mult, nbranches, igap, MPI_TYPE_tree_node, MPI_COMM_WORLD, ierr)

          deallocate(pack_mult)
          deallocate (igap)

          call timer_stop(t_exchange_branches_allgatherv)
          call timer_start(t_exchange_branches_integrate)

          ! Integrate remote branches into local tree
          do i = 1,nbranch_sum
              ! store branch key for later (global tree buildup)
              branch_keys(i) = get_mult(i)%key
              ! insert all remote branches into local data structures (this does *not* prepare the internal tree connections, but only copies multipole properties and creates the htable-entries)
              if (get_mult(i)%owner /= me) call tree_insert_node(get_mult(i))
          end do

          deallocate(get_mult)

          call timer_stop(t_exchange_branches_integrate)
          call timer_stop(t_exchange_branches)

      end subroutine tree_exchange

end module module_tree
