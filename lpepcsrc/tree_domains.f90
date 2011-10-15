
!  ================================
!
!         TREE_DOMAINS
!
!     Domain decomposition:
!     Share particle keys amoung PEs
!     - weighting according to load incurred on 
!       previous timestep
!
!  ================================


subroutine tree_domains(indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold,weighted,curve_type_)

  use treevars
  use tree_utils
  use timings
  use module_spacefilling
  use module_branching
  implicit none
  include 'mpif.h'

  integer, intent(in) :: weighted
  integer, intent(out) :: indxl(nppm),irnkl(nppm)
  integer, intent(out) :: islen(num_pe),irlen(num_pe)
  integer, intent(out) :: fposts(num_pe+1),gposts(num_pe+1)
  integer :: npnew,npold
  integer, intent(in) :: curve_type_ !< type of space-filling curve

  integer*8, dimension(nppm) :: ixd, iyd, izd
  integer*8, dimension(nppm) :: local_key

  integer :: i, j, ind_recv, inc, prev, next, handle(4)

  real*8 :: s
  real*8 :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
  logical :: boundary_debug=.false. 
  logical :: identical_keys=.false. 

  integer :: status(MPI_STATUS_SIZE), ierr

  ! arrays for parallel sort

  type (particle) :: ship_parts(nppm), get_parts(nppm)

  integer*8 :: w1(nppm)
  integer :: keycheck_pass, ipp

  logical :: sort_debug
  real*8 :: xboxsize, yboxsize, zboxsize

  real*8 imba

  type (particle) :: ship_props, get_props


  interface
     subroutine slsort_keys(nin,nmax,keys,workload,balance_weight,max_imbalance,nout,indxl,irnkl,scounts,rcounts,sdispls,rdispls,keys2,irnkl2,size,rank)
       integer,intent(in) :: nin,nmax,balance_weight,size,rank
       real*8,intent(in) :: max_imbalance
       integer,intent(out) :: nout,indxl(*),irnkl(*),scounts(*),rcounts(*),sdispls(*),rdispls(*),irnkl2(*)
       integer*8,intent(out) :: keys2(*)
       integer*8,intent(inout) :: keys(*)
       real*8,intent(inout) :: workload(*)
     end subroutine slsort_keys
  end interface

  real*8 work2(nppm)
  integer irnkl2(nppm)

  curve_type = curve_type_

  call timer_start(t_domains)
  call timer_start(t_domains_keys)

  sort_debug=domain_debug

  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | DOMAINS..'


  ! Find limits of local simulation region
  xmin_local = minval(particles(1:npp)%x(1))
  xmax_local = maxval(particles(1:npp)%x(1))
  ymin_local = minval(particles(1:npp)%x(2))
  ymax_local = maxval(particles(1:npp)%x(2))
  zmin_local = minval(particles(1:npp)%x(3))
  zmax_local = maxval(particles(1:npp)%x(3))

  ! Find global limits
  call MPI_ALLREDUCE(xmin_local, xmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(xmax_local, xmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ymin_local, ymin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ymax_local, ymax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(zmin_local, zmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(zmax_local, zmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )

  xboxsize = xmax-xmin
  yboxsize = ymax-ymin
  zboxsize = zmax-zmin

  ! Safety margin - put buffer region around particles
  xmax = xmax + xboxsize/10000.
  xmin = xmin - xboxsize/10000.
  ymax = ymax + yboxsize/10000.
  ymin = ymin - yboxsize/10000.
  zmax = zmax + zboxsize/10000.
  zmin = zmin - zboxsize/10000.

  boxsize = max(xmax-xmin, ymax-ymin, zmax-zmin)

  if (domain_debug) write (ipefile,'(4(a15,f12.4/))') &
       'xmin = ',xmin,'xmax = ',xmax, &
       'ymin = ',ymin,'ymax = ',ymax, &
       'zmin = ',zmin,'zmax = ',zmax, &
       'boxsize = ',boxsize

  s=boxsize/2**nlev       ! refinement length

  call compute_particle_keys(local_key)
  particles(1:npp)%key = local_key(1:npp)

  call timer_stop(t_domains_keys)
  call timer_start(t_domains_sort)

  ! Define wraps for ring network  0 -> 1 -> 2 -> ... ... -> num_pe-1 -> 0 ...
  if (me == 0) then
     prev = num_pe - 1
  else
     prev = me-1
  endif

  if (me == num_pe-1 ) then
     next = 0
  else
     next = me+1
  endif

  imba = 0.01

  npold = npp
  npnew = npp

  call timer_start(t_domains_add_sort)

  ! start permutation of local key list
  work2 = particles(:)%work

  call timer_start(t_domains_sort_pure)

  ! perform index sort on keys
  call slsort_keys(npold,nppm-2,local_key,work2,weighted,imba,npnew,indxl,irnkl,islen,irlen,fposts,gposts,w1,irnkl2,num_pe,me)

  ! FIXME: every processor has to have at least one particle
  if (npnew < 2) then
    write(*,'("PE ", I8, " has less than two particles after sorting (had ", I8, " before) - currently this can lead to errors --> aborting")') me, npold
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  endif

  call timer_stop(t_domains_sort_pure)

  call timer_stop(t_domains_sort)
  call timer_start(t_domains_ship)

  npp = npnew

  ! Now permute particle properties
  ! Set up particle structure 
  call timer_start(t_domains_add_pack)

  do i=1,npold
     ship_parts(i) = particles( indxl(i) )
  enddo

  call timer_stop(t_domains_add_pack)

  call timer_start(t_domains_add_alltoallv)

  ! perform permute
  call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_particle, &
          get_parts, irlen, gposts, mpi_type_particle, &
          MPI_COMM_WORLD,ierr)

  call timer_stop(t_domains_add_alltoallv)

  call timer_start(t_domains_add_unpack)

  do i=1,npp
     particles( irnkl(i) ) = get_parts(i)
  enddo

  call timer_stop(t_domains_add_unpack)


  if (npp > nppm) then
   write(*,*) "Something went seriously wrong during sorting: there are more than nppm local particles now."
   write(*,*) "npp = ", npp, ">   nppm = ", nppm
   write(*,*) "All local particle fields are too short. Aborting..."
   call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  endif

  call timer_stop(t_domains_ship)
  call timer_stop(t_domains_add_sort)
  call timer_start(t_domains_bound)

  if (domain_debug) then
     do j=1,npp
        ixd(j) = SUM( (/ (2_8**i*ibits( particles(j)%key-iplace,3*i,1 ), i=0,nlev-1) /) )
        iyd(j) = SUM( (/ (2_8**i*ibits( particles(j)%key-iplace,3*i+1,1 ), i=0,nlev-1) /) )
        izd(j) = SUM( (/ (2_8**i*ibits( particles(j)%key-iplace,3*i+2,1 ), i=0,nlev-1) /) )
     end do
     write (ipefile,'(/a/a/(z21,i6,a2,i8,6f12.4))') 'Particle list after key sort:', &
          '  key,                owner,    from PE  |  label  Fetched coords      derived from key', &
          (particles(i)%key, particles(i)%pid, '|', &
           particles(i)%label, &
           particles(i)%x(1),particles(i)%x(2),particles(i)%x(3),&
           ixd(i)*s+xmin,iyd(i)*s+ymin,izd(i)*s+zmin,&
           i=1,npp)

     write(ipefile,'(/)')
  endif

  particles(1:npp)%pid = me  ! new owner

  ! Each PE now has sorted segment of particles of length npp
  ! Note that now npp /= npart/num_pe, only approx depending on key distribution, or target shape.

  ! Check for identical keys
  keycheck_pass=1

  do i=2,npp-1
     if (particles(i)%key == particles(i-1)%key) then
        particles(  i)%key = particles(  i)%key + 1  ! Augment higher key
        particles(i+1)%key = particles(i+1)%key + 2  ! Augment next higher key to avoid 'triplet'
        !  Tweak momenta and positions to ensure particles drift apart
        particles(i-1)%u(1) = particles(i-1)%u(1) *0.99999
        particles(i+1)%u(1) = particles(i+1)%u(1) *1.00001

        particles(i-1)%x(1) = particles(i-1)%x(1) *0.99999
        particles(i+1)%x(1) = particles(i+1)%x(1) *1.00001

        ! TODO: need 'ripple' here up to next large gap in keys i+1->npp

        write(*,'(a15,i5,a8,i3,a30,2i9,3i10,a25,o25,a12,o25)') 'LPEPC | PE ',me,' pass ',keycheck_pass, &
             ' WARNING: identical keys found for index, npp, label1, label2  ', &
               i,npp,particles(i-1)%label,particles(i)%label,particles(i+1)%label, &
             ' -  to: ',particles(i)%key,' next key: ',particles(i+1)%key
     endif
  end do


  ! Special case for last pair to avoid possibility of identical keys across processor boundary
  ! - work back down from last key

  ipp=npp-1
  identical_keys=.true.
  keycheck_pass=2

  do while (identical_keys .and. ipp.gt.2)
     identical_keys=.false.
     if (particles(ipp+1)%key == particles(ipp)%key) then
        particles(ipp)%key = particles(ipp)%key-1
        write(*,'(a15,i9,a8,i3,a30,2i15/a25,o30)') 'LPEPC | PE ',me,' pass ',keycheck_pass, &
             ' WARNING: identical keys found for particles  ',particles(ipp+1)%label,particles(ipp)%label, &
             'LPEPC | Lower key decreased to:  ',particles(ipp)%key
        identical_keys=.true.
     endif
     ipp=ipp-1
  end do

  ! Copy boundary particles to adjacent PEs to ensure proper tree construction
  !  - if we do not do this, can get two particles on separate PEs 'sharing' a leaf
  ship_props = particles( 1 )

  ! Ship 1st particle data to end of list of LH neighbour PE

  if (me /= 0 ) then
     call MPI_ISEND( ship_props, 1, mpi_type_particle, prev, 1, MPI_COMM_WORLD, handle(1), ierr ) 
     call MPI_REQUEST_FREE(handle(1),ierr) 
  endif

  ! Place incoming data at end of array
  if ( me /= num_pe-1) then
     call MPI_RECV( get_props, 1, mpi_type_particle, next, 1,  MPI_COMM_WORLD, status, ierr )
     particles(npp+1) = get_props
  endif

  ! Ship  end particle data to start of list of RH neighbour PE
  ship_props = particles( npp )

  if (me /= num_pe-1 ) then
     call MPI_ISEND( ship_props, 1, mpi_type_particle, next, 2, MPI_COMM_WORLD, handle(3), ierr )
     call MPI_REQUEST_FREE(handle(3),ierr) 
  endif

  ! Place incoming data at end of array

  if (me == num_pe-1) then
     ind_recv = npp+1   ! PEn array has not yet received boundary value
  else
     ind_recv = npp+2
  endif

  if ( me /= 0) then
     call MPI_RECV( get_props, 1, mpi_type_particle, prev, 2,  MPI_COMM_WORLD, status, ierr )
     particles(ind_recv) = get_props
  endif

  ! Initialize VLD-stuff
  call get_virtual_local_domain()
  ! CSBE - Cross Sum Branch Node Estimator
  call get_local_apriori_est()
  call get_global_apriori_est()


  if (boundary_debug) then
     if (me /= 0 .and. me /= num_pe-1) then
        inc = 1
     else
        inc = 0
     endif

     write (ipefile,'(/a/a/(i5,z21,2i8,3f12.5,1f12.3))') 'Particle list after boundary swap:', &
          ' index   key,     label,   on PE,    x      y     q', &
          (i,particles(i)%key,particles(i)%label,particles(i)%pid,&
           particles(i)%x(1),particles(i)%x(2),particles(i)%x(3),particles(i)%q,i=1,npp+1+inc)

     write(ipefile,'(/)')
  endif

  call timer_stop(t_domains_bound)
  call timer_stop(t_domains)

  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | ..done'

end subroutine tree_domains
