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
  integer, dimension(num_pe) :: particles, nonlocal_keys, fetches, ships, total_keys, tot_nleaf, tot_ntwig
  character*6 :: cdump
  character*40 :: cfile

! particle distrib 
  call MPI_GATHER(npp, 1, MPI_INTEGER, particles, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(ntwig_me, 1, MPI_INTEGER, tot_ntwig, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(nleaf_me, 1, MPI_INTEGER, tot_nleaf, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(nkeys_total, 1, MPI_INTEGER, total_keys, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(sum_fetches, 1, MPI_INTEGER, fetches, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(sum_ships, 1, MPI_INTEGER, ships, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  
  call MPI_GATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL, 0,  MPI_COMM_WORLD, ierr )  


if (me.eq.0) then
  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  cfile="log/stats"//"."//cdump(1:6)
  
  open (60,file=cfile)    
  write (60,'(a98/(4i10,F8.4,5i10))') '         PE     parts    nleaf     ntwig   ratio    nl_keys', &
	'   tot_keys   fetches    ships    work', &
      (i-1,particles(i),tot_nleaf(i),tot_ntwig(i),1.0*tot_nleaf(i)/(1.0*tot_ntwig(i)), &
	total_keys(i)-(tot_nleaf(i)+tot_ntwig(i)),total_keys(i),fetches(i),ships(i),int(work_loads(i)),i=1,num_pe)
 close(60)

endif

end subroutine tree_stats
