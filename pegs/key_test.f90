! ===========================================
!
!           KEY_TEST
!
!   Function to test whether key present in 
!  local hash table 
!   
!
! ===========================================

function key_local(keyin)

  use treevars
  use utils

  implicit none
  integer*8  :: keyin
  integer :: cell_addr, link_addr, ires
  logical :: resolved

  logical :: key_local 

  cell_addr = IAND( keyin, hashconst)     ! cell address hash function

  if ( htable( cell_addr )%key == keyin ) then
     key_local = .true.       ! Keys match -> found entry
     return

  else if (htable(cell_addr)%link == -1) then
     key_local = .false.   ! Keys don't match, and no link here => key not in # table!
     return

  else 
     link_addr = htable(cell_addr)%link  ! Go to link address
     ires = 0

     do while ( link_addr /= -1 )       ! Repeat until keys match or run out of links
        ires = ires + 1
        if ( htable( link_addr )%key == keyin ) then
           key_local = .true.      ! Keys match -> found entry
           return
        else
           link_addr = htable(link_addr)%link    ! Next linked entry
        endif
     end do

     key_local=.false.  ! Not resolved: key not in table

  endif
end function key_local
