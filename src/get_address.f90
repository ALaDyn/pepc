! ===========================================
!
!           GET_ADDRESS
!
!   Function to return hash table entry 
!   address from key
!   
!
! ===========================================

function key2addr(keyin)

  use treevars
  use utils

  implicit none
  integer*8  :: keyin
  integer :: cell_addr, link_addr, ires,i, ierr
  logical :: resolved

  integer :: key2addr 

  cell_addr = IAND( keyin, hashconst)     ! cell address hash function

  if ( htable( cell_addr )%key == keyin ) then
     key2addr = cell_addr       ! Keys match -> found entry
     return
  else
     resolved = .false.
     link_addr = cell_addr
     ires = 0

     do while (  .not. resolved .and. ires <= maxaddress )       ! Repeat until keys match or run out of links
        link_addr = htable(link_addr)%link    ! Next linked entry
        ires = ires + 1
        if ( htable( link_addr )%key == keyin ) then
           key2addr = link_addr      ! Keys match -> found entry
           return
        endif
     end do
     ! Not resolved - something wrong: invalid key or #-table wrong
     write (*,*) 'Key not resolved in KEY2ADDR: check #-table and key list for PE ',me
     write (*,'(a5,o20,a2,i10,a1,a12,i15)') 'Key #: ',keyin,' (',keyin,')',' Address: ',cell_addr
     write (ipefile,*) 'Keys in table:'
     do i=0,maxaddress
        if (htable(i)%key/=0) write(ipefile,'(i8,o21)') i,htable(i)%key
     end do

!     call diagnose_tree
close(75)     
     call closefiles
!     pause
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop
  endif
end function key2addr
