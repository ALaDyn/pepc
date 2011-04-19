!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Calculates inverse mapping from the hash-key to the hash-address.
!>
!> @param[in] keyin inverse mapping candidate.
!> @param[in] cmark a description.
!> @param[out] key2addr the adress if the key exists
!> @exception if key does not exist, the whole program is aborted
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function key2addr(keyin,cmark)

  use treevars
  implicit none
  include 'mpif.h'

  integer*8, intent(in)  :: keyin
  integer :: cell_addr, link_addr, ires,i, ierr
  logical :: resolved
  character(LEN=*) :: cmark
  integer :: key2addr

  cell_addr = int(IAND( keyin, hashconst))     ! cell address hash function

  if ( htable( cell_addr )%key == keyin ) then
     key2addr = cell_addr       ! Keys match -> found entry
     return
  else
     resolved = .false.
     link_addr = cell_addr
     ires = 0

     do while (  .not. resolved .and. ires <= maxaddress )       ! Repeat until keys match or run out of links
        link_addr = htable(link_addr)%link    ! Next linked entry
        if (link_addr == -1 ) then
           write (*,'(a,a20)') 'Key not resolved in KEY2ADDR at ',cmark
           write (*,*) 'check #-table and key list for PE ',me
           write(*,*) 'Bad address'
           write(*,'(a15,o22)') 'Key = ',keyin
           write(*,*) 'Initial address =',cell_addr
           write(*,*) '# const =',hashconst
           write(*,*) 'ires =',ires
           exit
        endif
        ires = ires + 1
        if ( htable( link_addr )%key == keyin ) then
           key2addr = link_addr      ! Keys match -> found entry
           return
        endif
     end do
     ! Not resolved - something wrong: invalid key or #-table wrong
     write (*,'(a5,o24,a2,i20,a1,a12,i15)') 'Key #: ',keyin,' (',keyin,')',' Address: ',cell_addr
     write (ipefile,*) 'Keys in table:'
     do i=0,maxaddress
        if (htable(i)%key/=0) write(ipefile,'(i8,o25,i10)') i,htable(i)%key,htable(i)%link
     end do

!     call diagnose_tree
close(75)     
!     call closefiles
!     pause
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop
  endif
end function key2addr






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Calculates inverse mapping from the hash-key to the hash-address.
!>
!> @param[in] keyin inverse mapping candidate.
!> @return true if candidate exists
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function testaddr(keyin)

  use treevars
  implicit none
  include 'mpif.h'

  integer*8, intent(in)  :: keyin

  integer :: cell_addr, link_addr, ires
  logical :: resolved
  logical :: testaddr

  ! cell address hash function
  cell_addr = INT(IAND(keyin,hashconst))

  ! Keys match -> found entry
  if ( htable( cell_addr )%key == keyin ) then
     testaddr = .true.
     return
  else
     resolved = .false.
     link_addr = cell_addr
     ires = 0

     do while (.not.resolved .and. ires <= maxaddress )       ! Repeat until keys match or run out of links
       link_addr = htable(link_addr)%link    ! Next linked entry

       ! no more links, thus no entry for keyin
       if (link_addr == -1 ) then
           testaddr = .false.
           return
        endif

        ires = ires + 1

        ! Keys match -> found entry
        if ( htable( link_addr )%key == keyin ) then
           testaddr = .true.
           return
        endif
     end do
  endif
end function testaddr
