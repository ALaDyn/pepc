!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Function to return key of next node in local tree-walk,
!> i.e. search for next sibling, uncle, great-uncle etc
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function next_node(keyin)

  use treevars

  implicit none
  integer*8 :: next_node
  integer*8, intent(in) :: keyin
  integer :: start_child_idx = 0


  integer*8 :: search_key, parent_key
  integer :: parent_addr
  integer :: parent_child_byte, search_child_idx

  integer :: key2addr        ! Mapping function to get hash table address from key

  search_key = keyin

  ! search for next sibling, uncle, great-uncle etc
  do while (search_key > 1) ! loop through parent nodes up to root
    parent_key        = ishft(search_key,-3)
    parent_addr       = key2addr( parent_key ,"next_node(), get parent_addr" )
    parent_child_byte = ibits( htable( parent_addr ) % childcode, 0, 8)

    search_child_idx  = int(ibits( search_key, 0, 3), kind(search_child_idx) ) ! lower three bits of key

    do ! loop over all siblings
      search_child_idx   = modulo(search_child_idx + 1, 8) ! get next sibling, possibly starting again from first one

      ! if sibling-loop wrapped and reached starting point again --> go up one level
      if ( search_child_idx == start_child_idx ) then
          search_key = parent_key      ! go up one level
          exit
      endif

      ! if sibling exists: next_node has been found
      if ( btest(parent_child_byte, search_child_idx) ) then
          next_node = ior(ishft(parent_key, 3), search_child_idx) ! assemble next_node out of parent-key and new sibling-index
          return
      endif
    end do
  end do

  next_node  = 1 ! nothing has been found, i.e. top-right corner reached: set pointer=root

end function next_node
