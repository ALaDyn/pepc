!  ================================
!
!         TREE_STATS
!
!     Perform  stats on tree comm.
!
!
!  ================================


subroutine tree_stats(timestamp)

  use treevars
  implicit none
  include 'mpif.h'

  integer :: i,ierr, timestamp
  integer, dimension(nppm) :: particles, local_keys, nonlocal_keys, fetches, ships, total_keys
  integer :: nkeys_tot, nkeys_me
  character*6 :: cdump
  character*40 :: cfile

! particle distrib 
  call MPI_GATHER(npp, 1, MPI_INTEGER, particles, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  nkeys_me=nleaf_me+ntwig_me
  nkeys_tot=nleaf+ntwig
  call MPI_GATHER(nkeys_me, 1, MPI_INTEGER, local_keys, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(nkeys_tot, 1, MPI_INTEGER, total_keys, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(sum_fetches, 1, MPI_INTEGER, fetches, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(sum_ships, 1, MPI_INTEGER, ships, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  


if (me.eq.0) then
  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  cfile="log/stats"//"."//cdump(1:6)

  open (60,file=cfile)    
  write (60,'(a60/(7i10))') 'PE  parts local_keys  nl_keys  tot_keys fetches ships', &
      (i-1,particles(i),local_keys(i),total_keys(i)-local_keys(i),total_keys(i),fetches(i),ships(i),i=1,num_pe-1)
  close(60)

endif

end subroutine tree_stats
