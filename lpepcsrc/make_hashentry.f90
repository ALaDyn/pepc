! ===========================================
!
!           MAKE_HASH
!
!  Make entry in hash-table - returns address 'newentry'
!  Resolve collision if necessary
!   
!
! ===========================================

subroutine make_hashentry( keyin, nodein, leavesin, codein, ownerin, newentry, ierror)

  use treevars
  use utils

  implicit none

  integer*8, intent(in) :: keyin
  integer, intent(in) :: nodein, leavesin, codein, ownerin      ! call input parameters
  integer, intent(out) :: newentry  ! address in # table returned to calling routine
  integer :: ierr
  integer :: link_addr, cell_addr, ierror, k, free, addr_count
  logical :: resolved

  ierror = 0
  addr_count = 1 
  free = free_addr(iused)

  ! perform address hashing on key
  cell_addr = IAND( keyin, hashconst)

  tablehigh = max(tablehigh,cell_addr)                 ! Track highest address

  if ( htable( cell_addr )%node == 0 .and. htable( cell_addr )%key /=-1 ) then       ! Is entry empty?
     newentry = cell_addr                 ! Yes, so create new entry:
     htable( cell_addr )%node = nodein          !   local pointer
     htable( cell_addr )%key =  keyin            !   key
     htable( cell_addr )%leaves = leavesin       !   # contained nodes
     htable( cell_addr )%childcode = codein       !  child byte-code or particle #
     htable( cell_addr )%owner = ownerin       ! PE owning branch node

     if (point_free(cell_addr) /= 0) then     ! Check if new address in collision res. list
        free_addr( point_free(cell_addr) ) = free_addr(sum_unused)  ! Replace free address with last on list
        point_free(free_addr(sum_unused)) = point_free(cell_addr)   ! Reset pointer
        point_free(cell_addr) = 0
        sum_unused = sum_unused - 1
     endif

  else if ( htable( cell_addr )%node /= 0 .AND. htable(cell_addr)%key == keyin ) then  
     ! Entry exists and keys match
     ! => local node  so skip
     ierror = 1


  else            ! Entry exists and keys don't match: COLLISION

     if ( htable( cell_addr )%link == -1) then     ! Entry was occupied without existing link

        newentry =  free  ! Pick next free # address from list of unused table positions (computed at end of treebuild)

        if (htable(free)%node /= 0 ) then
           write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
           write (*,*) 'PE ',me,' key ',keyin,' entry',newentry,' used ',iused,'/',sum_unused
           pause
        endif
        htable( free )%node = nodein                     
        htable( free )%key = keyin
        htable( free )%childcode = codein
        htable( free )%leaves = leavesin
        htable( free )%owner = ownerin
        htable( free )%link = -1
        htable( cell_addr )%link = free     ! Create link from 1st duplicate entry to new entry
        iused = iused + 1        ! increment free address counter


     else if ( htable( cell_addr )%link /= -1 ) then     ! Address occupied with link already

        link_addr = cell_addr                          ! Start of chain
        resolved = .false.                                     ! Resolve flag

        do while (  .not. resolved .and. addr_count < maxaddress)       ! Find end of chain  
           link_addr = htable(link_addr)%link

           if ( htable( link_addr )%key == keyin ) then
              ! Occupied with same key -> local or boundary node, so skip
              resolved = .true.
              ierror = 1

           else if ( htable(link_addr)%node == 0 .and. htable (link_addr)%link == -1 ) then
              ! Found end of chain: entry was occupied by boundary node, so reuse entry
              ! link from previous chain entry still valid
              newentry = link_addr
              htable( link_addr )%node = nodein                     
              htable( link_addr )%key = keyin
              htable( link_addr )%childcode = codein
              htable( link_addr )%leaves = leavesin
              htable( link_addr )%owner = ownerin

              if (point_free(link_addr) /= 0) then     ! Check if new address in collision res. list
                 free_addr( point_free(link_addr) ) = free_addr(sum_unused)  ! Replace free address with last on list
                 point_free(free_addr(sum_unused)) = point_free(link_addr)   ! Reset pointer
                 point_free(link_addr) = 0
                 sum_unused = sum_unused - 1
              endif
              resolved = .true.

           else if ( htable(link_addr)%node /= 0 .and. htable (link_addr)%link == -1 ) then
              ! Found end of chain: entry occupied with no link
              newentry = free          

              if (htable(free)%node /= 0 ) then
                 write (*,*) 'Something wrong with address list for collision resolution (free_addr in treebuild)'
                 write (*,*) 'PE ',me,' key ',keyin,' entry',newentry,' used ',iused,'/',sum_unused
                 pause
              endif

              htable( free )%node = nodein                     
              htable( free )%key = keyin
              htable( free )%childcode = codein
              htable( free )%leaves = leavesin
              htable( free )%owner = ownerin
              htable( free )%link = -1
              htable( link_addr )%link = free     ! Create link from 1st duplicate entry to new entry
              iused = iused + 1                  ! Increment free_address counter
              resolved = .true.
           else
	     addr_count = addr_count + 1
              ! not yet resolved - go to next link in chain
           endif
        end do

        if (addr_count >= maxaddress) then
	  write (ipefile,'(a5,o20,a)') 'Key ',keyin,' not resolved in MAKE_HASHENTRY'
	  call MPI_ABORT(ierr)
        endif
     else
        write (ipefile,'(a5,o20,a)') 'Key ',keyin,' not resolved in MAKE_HASHENTRY'
        ierror = 2
     end if
     tablehigh = max(tablehigh,free)                 ! Track highest address

  end if


end subroutine make_hashentry

