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
  use tree_walk_communicator
  use tree_walk_utils
  use module_branching
  implicit none
  include 'mpif.h'

  integer :: i,ierr, timestamp
  integer, dimension(num_pe) :: particles, fetches, ships, total_keys, tot_nleaf, tot_ntwig
  real*8, dimension(num_pe) ::  num_interactions, num_mac_evaluations  ! Load balance arrays
  character*40 :: cfile
  integer :: max_nbranch,min_nbranch, gmax_leaves, gmax_twigs, total_part
  real*8 :: average_interactions, average_mac_evaluations, total_interactions, total_mac_evaluations, max_interactions, max_mac_evaluations
  real :: part_imbal=0.
  real*8 :: work_imbal=0.
  real*8 :: work_imbal_max, work_imbal_min  ! load stats
  integer ::  part_imbal_max, part_imbal_min
  real*8 :: global_thread_workload(-4:4)

  ! particle distrib
  call MPI_GATHER(npp,         1, MPI_INTEGER, particles,  1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(ntwig_me,    1, MPI_INTEGER, tot_ntwig,  1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(nleaf_me,    1, MPI_INTEGER, tot_nleaf,  1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(nkeys_total, 1, MPI_INTEGER, total_keys, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(sum_fetches, 1, MPI_INTEGER, fetches,    1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(sum_ships,   1, MPI_INTEGER, ships,      1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(interactions_local,    1, MPI_REAL8, num_interactions,      1, MPI_REAL8,   0,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(mac_evaluations_local, 1, MPI_REAL8, num_mac_evaluations,   1, MPI_REAL8,   0,  MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(nbranch, max_nbranch,     1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(nbranch, min_nbranch,     1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload( 1), global_thread_workload( 1), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload( 2), global_thread_workload( 2), 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload( 3), global_thread_workload( 3), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload( 4), global_thread_workload( 4), 2, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload(-1), global_thread_workload(-1), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload(-2), global_thread_workload(-2), 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload(-3), global_thread_workload(-3), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(thread_workload(-4), global_thread_workload(-4), 2, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
  nnodes=nleaf+ntwig
  call MPI_REDUCE(nleaf, gmax_leaves, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE(ntwig, gmax_twigs,  1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr )

  part_imbal_max = MAXVAL(particles)
  part_imbal_min = MINVAL(particles)
  part_imbal = (part_imbal_max-part_imbal_min)/1.0/npart*num_pe

  total_interactions       = SUM(num_interactions)
  total_mac_evaluations    = SUM(num_mac_evaluations)
  max_interactions         = MAXVAL(num_interactions)
  max_mac_evaluations      = MAXVAL(num_mac_evaluations)
  average_interactions     = total_interactions    / num_pe
  average_mac_evaluations  = total_mac_evaluations / num_pe
  work_imbal_max = max_interactions/average_interactions
  work_imbal_min = MINVAL(num_interactions)/average_interactions
  work_imbal = 0.
  do i=1,num_pe
     work_imbal = work_imbal + abs(num_interactions(i) - average_interactions)/average_interactions/num_pe
  end do

  total_part = sum(particles)

  if (me.eq.0) then
    call system("mkdir -p " // "stats")
    write(cfile,'("stats/stats.",i6.6)') timestamp
  
    open (60,file=trim(cfile))

    write (60,'(a20,i7,a22)') 'Tree stats for CPU ', me, ' and global statistics'
    write (60,*) '######## GENERAL DATA #####################################################################'
    write (60,'(a50,3i12)') '# procs, walk_threads, max_particles_per_thread: ', num_pe, num_walk_threads, max_particles_per_thread
    write (60,'(a50,i12,f12.2,i12)') 'nintmax, np_mult, size_tree: ',nintmax, np_mult,size_tree
    write (60,'(a50,3i12)') 'npp, npart, nppm(max): ',npp,npart,nppm
    write (60,'(a50,2i12)') 'total # particles, N/P: ',total_part,int(npart/num_pe)
    write (60,*) '######## TREE STRUCTURES ##################################################################'
    write (60,'(a50,3i12)') 'local # leaves, twigs, keys: ',nleaf_me,ntwig_me,nleaf_me+ntwig_me
    write (60,'(a50,3i12)') 'non-local # leaves, twigs, keys: ',nleaf-nleaf_me,ntwig-ntwig_me,nleaf+ntwig-nleaf_me-ntwig_me
    write (60,'(a50,3i12,f12.1,a6,i12)') 'final # leaves, twigs, keys, (max): ',nleaf,ntwig,nleaf+ntwig, &
              (nleaf+ntwig)/(.01*maxaddress),' % of ',maxaddress
    write (60,'(a50,3i12)') 'maximum # leaves, twigs, keys(=maxaddress): ',maxleaf, maxtwig, maxaddress
    write (60,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # leaves: ',gmax_leaves, gmax_leaves/(.01*maxleaf), ' % of  ', maxleaf
    write (60,'(a50,1i12,1f12.1, a6,1i12)') 'Global max # twigs: ', gmax_twigs, gmax_twigs/(.01*maxtwig), ' % of  ', maxtwig
    write (60,*) '######## BRANCHES #########################################################################'
    write (60,'(a50,3i12)') '#branches local, max_global, min_global: ', nbranch,max_nbranch,min_nbranch
    write (60,'(a50,2i12)') '#branches global sum estimated, sum actual: ',branch_max_global,nbranch_sum
    write (60,'(a50,2i12)') 'max res.space for local branches, global br.: ', branch_max_local,branch_max_global
    write (60,*) '######## WALK-COMMUNICATION ###############################################################'
    write (60,'(a50,2i12)') 'Max # multipole fetches/ships per cpu: ',maxval(fetches), maxval(ships)
    write (60,'(a50,2i12)') 'Min # multipole fetches/ships per cpu: ',minval(fetches), minval(ships)
    write (60,'(a50,2i12)') 'Local #  multipole fetches & ships: ',sum_fetches,sum_ships
    write (60,'(a50,2i12)') 'cumulative/maximum # of entries in request queue: ', cum_req_list_length, max_req_list_length
    write (60,'(a50,3i12)') '# of comm-loop iterations (tot,send,recv): ', comm_loop_iterations(:)
    write (60,*) '######## WORKLOAD AND WALK ################################################################'
    write (60,'(a50,3e12.4)')       'total/ave/max_local # interactions(work): ', total_interactions, average_interactions, max_interactions
    write (60,'(a50,3e12.4)')       'total/ave/max_local # mac evaluations: ', total_mac_evaluations, average_mac_evaluations, max_mac_evaluations
    write (60,'(a50,3f12.3)')       'Load imbalance percent,min,max: ',work_imbal,work_imbal_min,work_imbal_max
    write (60,'(a50,f12.3,2i12)')   'Particle imbalance ave,min,max: ',part_imbal,part_imbal_min,part_imbal_max
    write (60,*) '######## WALK-WORKER-THREAD WORKLOAD ######################################################'
    write (60,'(a50)')              'average # processed particles per thread    '
    write (60,'(a50,3f12.3)')       '  threads on exclusive cores, shared cores: ', thread_workload(1), thread_workload(3)
    write (60,'(a50,3f12.3)')       '  maximum relative deviation: ', thread_workload(2), thread_workload(4)
    write (60,'(a50)')              'average wallclocktime per thread    '
    write (60,'(a50,3f12.3)')       '  threads on exclusive cores, shared cores: ', thread_workload(-1) , thread_workload(-3)
    write (60,'(a50,3f12.3)')       '  maximum relative deviation: ', thread_workload(-2), thread_workload(-4)
    write (60,*) '###########################################################################################'
    write (60,*) '######## DETAILED DATA ####################################################################'
    write (60,'(2a/(4i10,F8.4,6i15,F8.4))') '         PE     parts    nleaf     ntwig   ratio    nl_keys', &
              '   tot_keys   fetches    ships    #interactions(work)   #mac_evals   rel.work*ncpu', &
              (i-1,particles(i),tot_nleaf(i),tot_ntwig(i),1.0*tot_nleaf(i)/(1.0*tot_ntwig(i)), &
              total_keys(i)-(tot_nleaf(i)+tot_ntwig(i)),total_keys(i),fetches(i),ships(i),int(num_interactions(i)),int(num_mac_evaluations(i)),&
              num_interactions(i)/average_interactions,i=1,num_pe)
    close(60)

  endif

end subroutine tree_stats
