! ===========================================
!
!           GET_NEXT_NODE
!
! $Revision$
!
!   Function to return key of next node in local tree-walk
!   
!
! ===========================================

function next_node(keyin)

  use treevars

  implicit none
  integer*8 :: next_node
  integer*8 :: keyin
  logical :: resolved
  logical, dimension(8) :: keymatch
  integer*8, dimension(8) :: child_key, child_sub
  integer*8 :: node_key, search_key, parent,  child_top
  integer :: node_addr, parent_node, child_byte, nchild
  integer :: jmatch(1), j
  integer :: key2addr        ! Mapping function to get hash table address from key

  search_key = keyin                   
  node_key = keyin                     ! keep key, address of node 
  node_addr = key2addr(keyin)
  resolved = .false.

  !   Search for next sibling, uncle, great-uncle etc

  do while (.not. resolved .and. search_key > 1)
     parent =  ishft(search_key,-3)                 ! parent
     parent_node = htable( key2addr( parent ) )%node   ! parent node pointer

     child_byte = htable( key2addr( parent ) )%childcode                           !  Children byte-code
     nchild = SUM( (/ (ibits(child_byte,j,1),j=0,7) /) )                   ! # children = sum of bits in byte-code
     child_sub(1:nchild) = pack( bitarr, mask=(/ (btest(child_byte,j),j=0,7) /) )  ! Extract child sub-keys from byte code
     child_top = ishft(parent,3)  
     child_key(1:nchild) = IOR( child_top, child_sub(1:nchild) )         ! Construct keys of children

     keymatch=.false.
     keymatch(1:nchild) = (/ (child_key(j) == search_key,j=1,nchild) /)
     jmatch = pack(bitarr, mask = keymatch ) + 1                                  ! Pick out position of current search key in family

     if (jmatch(1) < nchild ) then                                                ! if search_key has 'elder' sibling then
        next_node  = child_key(jmatch(1)+1)                        ! store next_node as sibling of parent/grandparent
        resolved = .true.
     else
        search_key = ishft(search_key, -3)                                     ! Go up one level 
     endif
  end do

  if (.not. resolved .and. search_key == 1) next_node  = 1        ! Top-right corner reached: set pointer=root

end function next_node
