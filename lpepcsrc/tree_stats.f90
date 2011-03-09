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
  integer, dimension(num_pe) :: particles, fetches, ships, total_keys, tot_nleaf, tot_ntwig
  character*6 :: cdump
  character*40 :: cfile
  integer :: max_fetches, max_nbranch, gmax_leaves, gmax_twigs, total_part
  real :: total_work, average_work

  ! particle distrib
  call MPI_GATHER(npp, 1, MPI_INTEGER, particles, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(ntwig_me, 1, MPI_INTEGER, tot_ntwig, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(nleaf_me, 1, MPI_INTEGER, tot_nleaf, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(nkeys_total, 1, MPI_INTEGER, total_keys, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(sum_fetches, 1, MPI_INTEGER, fetches, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(sum_ships, 1, MPI_INTEGER, ships, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL, 0,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(sum_fetches, max_fetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(nbranch, max_nbranch, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  nnodes=nleaf+ntwig
  call MPI_ALLREDUCE(nleaf, gmax_leaves, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ntwig, gmax_twigs, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )

  call MPI_GATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
  part_imbal_max = MAXVAL(particles)
  part_imbal_min = MINVAL(particles)
  part_imbal = (part_imbal_max-part_imbal_min)/1.0/npart*num_pe
  total_work = SUM(work_loads)
  average_work = total_work/num_pe
  work_imbal_max = MAXVAL(work_loads)/average_work
  work_imbal_min = MINVAL(work_loads)/average_work
  work_imbal = 0.
  do i=1,num_pe
     work_imbal = work_imbal + abs(work_loads(i) - average_work)/average_work/num_pe
  end do

  total_part = sum(particles)

  if (me.eq.0) then
    ! get filename suffix from dump counter
    do i=0,4
       cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)
    end do
    cdump(1:1) = achar(timestamp/10**5 + 48)

    cfile="stats"//"."//cdump(1:6)
  
    open (60,file=cfile)

    write(60,'(a20,i7,a22)') 'Tree stats for CPU ', me, ' and global statistics'
    write(60,'(a50,3i12)') 'new npp, npart, (max): ',npp,npart,nppm
    write(60,'(a50,1i12)') 'total # particles: ', total_part
    write(60,'(a50,2e12.4)') 'total/ave work: ', total_work, average_work
    write(60,'(a50,3i12)') 'local # leaves, twigs, keys: ',nleaf_me,ntwig_me,nleaf_me+ntwig_me
    write(60,'(a50,3i12)') 'non-local # leaves, twigs, keys: ',nleaf-nleaf_me,ntwig-ntwig_me,nleaf+ntwig-nleaf_me-ntwig_me
    write(60,'(a50,3i12,f12.1,a6,i12)') 'final # leaves, twigs, keys, (max): ',nleaf,ntwig,nleaf+ntwig, &
              (nleaf+ntwig)/(.01*maxaddress),' % of ',maxaddress
    write(60,'(a50,3i12)') 'maximum # leaves, twigs, keys: ',maxleaf, maxtwig, maxaddress
    write(60,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # leaves: ',gmax_leaves, gmax_leaves/(.01*maxleaf), ' % of  ', maxleaf
    write(60,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # twigs: ', gmax_twigs, gmax_twigs/(.01*maxtwig), ' % of  ', maxtwig
    write(60,'(a50,3i12,a3,2i12)') 'Local, max local, global # branches, (max): ', &
              nbranch,max_nbranch,nbranch_sum,'/',nbranch_local_max,nbranch_max

    write (60,'(a50,2i12)') 'Max # multipole fetches/ships per cpu: ',maxval(fetches), maxval(ships)
    write (60,'(a50,2i12)') 'Local #  multipole fetches & ships: ',sum_fetches,sum_ships
    write (60,'(a50,3f12.3)') 'Load imbalance percent,min,max: ',work_imbal,work_imbal_min,work_imbal_max
    write (60,'(a50,f12.3,2i12)') 'Particle imbalance ave,min,max: ',part_imbal,part_imbal_min,part_imbal_max
    write (60,'(a50,2i12)') 'cumulative/maximum # of entries in request queue: ', cum_req_list_length, max_req_list_length
    write (60,'(a50,3i12)') '# of comm-loop iterations (tot,send,recv): ', comm_loop_iterations(:)
    write (60,*) '###########################################################################'
    write (60,'(2a/(4i10,F8.4,5i10,F8.4))') '         PE     parts    nleaf     ntwig   ratio    nl_keys', &
              '   tot_keys   fetches    ships    work   rel.work*ncpu', &
              (i-1,particles(i),tot_nleaf(i),tot_ntwig(i),1.0*tot_nleaf(i)/(1.0*tot_ntwig(i)), &
              total_keys(i)-(tot_nleaf(i)+tot_ntwig(i)),total_keys(i),fetches(i),ships(i),int(work_loads(i)),&
              work_loads(i)/average_work,i=1,num_pe)
    close(60)

  endif

end subroutine tree_stats
