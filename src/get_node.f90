! ===========================================
!
!           GET_NODE
!
!   Function to return node number in hash table
!   from key
!   
!
! ===========================================

integer function key2node(keyin)

use treevars
use utils

implicit none
integer*8 :: keyin
integer :: nodeout
integer :: cell_addr, link_addr
logical :: resolved

cell_addr = IAND( keyin, hashconst)     ! cell address hash function

if ( htable( cell_addr )%key == keyin ) then
   key2node = htable( cell_addr )%node      ! Keys match -> found entry
   return
else
   resolved = .false.
   link_addr = cell_addr

   do while (  .not. resolved .or. link_addr >= maxaddress )       ! Repeat until keys match or run out of links
      link_addr = htable(link_addr)%link    ! Next linked entry
      if ( htable( link_addr )%key == keyin ) then
         key2node = htable( link_addr )%node      ! Keys match -> found entry
         return
      endif
   end do
   ! Not resolved - something wrong: invalid key or #-table wrong
   write (*,*) 'Key not resolved in KEY2NODE: check #-table and key list'
   write (*,'(a5,o20,a12,i15)') 'Key #: ',keyin,' Address: ',cell_addr

   ! stop
endif
end function
