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
!VAMPINST include
!      INCLUDE 'VTcommon.h'
!      INTEGER VTIERR
!

  implicit none
  integer*8  :: keyin
  integer :: cell_addr, link_addr, ires
  logical :: resolved

  integer :: key2addr 

  cell_addr = IAND( keyin, hashconst)     ! cell address hash function

!VAMPINST subroutine_start
!      CALL VTENTER(IF_key2addr,VTNOSCL,VTIERR)
!      write(*,*) 'VT: key2addr S>',VTIERR,
!     *    IF_key2addr,ICLASSH
!
  if ( htable( cell_addr )%key == keyin ) then
     key2addr = cell_addr       ! Keys match -> found entry
!VAMPINST return
!       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: key2addr S<',VTIERR,ICLASSH
!
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
!VAMPINST return
!       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: key2addr S<',VTIERR,ICLASSH
!
           return
        endif
     end do
     ! Not resolved - something wrong: invalid key or #-table wrong
     write (*,*) 'Key not resolved in KEY2ADDR: check #-table and key list for PE ',me
     write (*,'(a5,o20,a2,i10,a1,a12,i15)') 'Key #: ',keyin,' (',keyin,')',' Address: ',cell_addr
!     call diagnose_tree
close(75)     
close(ipefile)
!     pause
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
!VAMPINST stop
!      CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: key2addr S<',VTIERR,ICLASSH
!
     stop
  endif
!VAMPINST subroutine_end
!       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: key2addr S<',VTIERR,ICLASSH
!
end function key2addr
