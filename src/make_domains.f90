
!  ================================
!
!         MAKE_DOMAINS
!
!     Domain decomposition:
!     Share particle keys amoung PEs
!     - weighting according to load incurred on 
!       previous timestep
!
!  ================================


subroutine make_domains

  use treevars
  use utils

  implicit none


  integer*8, dimension(nppm) :: ix, iy, iz
  integer*8, dimension(nppm) :: ixd, iyd, izd
  integer*8, dimension(nppm) :: local_key

  integer ::  source_pe(nppm)
  integer :: i, j, ind_recv, inc


  integer :: i1, i2, i3, ibegin, rec_count, iwait, ipack, ipe, source_id, send_count, nplace, insert_loc
  integer :: send_part_count, n_reqs
  integer :: ixc, iyc, izc, nbits


  real :: s
  real :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
  logical :: boundary_debug=.false., key_debug=.false.


  ! arrays for parallel sort
  integer :: npg, mp
  integer*8 :: xarray(nppm),keys(nppm),w1(nppm),wi2(nppm),wi3(nppm)
  integer :: indxl(nppm),irnkl(nppm)

  integer :: islen(num_pe),irlen(num_pe)
  integer :: fposts(num_pe+1),gposts(num_pe+1)

  integer :: npnew,npold
  integer :: iteration, niterations
  integer :: errcount

  integer, dimension(nppm) ::  w2, w3 ! scratch arrays for integer*4 permute
  real, dimension(nppm) :: wr1, wr2, wr3 ! Scratch for real array permute
  integer*8 :: tmp



  ! Find limits of local simulation region
  xmin_local = minval(x(1:npp))
  xmax_local = maxval(x(1:npp))
  ymin_local = minval(y(1:npp))
  ymax_local = maxval(y(1:npp))
  zmin_local = minval(z(1:npp))
  zmax_local = maxval(z(1:npp))

  ! Find global limits

  call MPI_ALLREDUCE(xmin_local, xmin, one, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(xmax_local, xmax, one, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ymin_local, ymin, one, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ymax_local, ymax, one, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(zmin_local, zmin, one, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(zmax_local, zmax, one, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )


  boxsize = max(xmax-xmin, ymax-ymin, zmax-zmin)

  ! Safety margin - put buffer region around particles
  xmax = xmax + boxsize/1000.
  xmin = xmin - boxsize/1000.
  ymax = ymax + boxsize/1000.
  ymin = ymin - boxsize/1000.
  zmax = zmax + boxsize/1000.
  zmin = zmin - boxsize/1000.

  boxsize = max(xmax-xmin, ymax-ymin, zmax-zmin)

  s=boxsize/2**nlev       ! refinement length

  ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)

  ix(1:npp) = ( x(1:npp) - xmin )/s           ! partial keys
  iy(1:npp) = ( y(1:npp) - ymin )/s           !
  iz(1:npp) = ( z(1:npp) - zmin )/s        

  ! construct keys by interleaving coord bits and add placeholder bit
  ! - note use of 64-bit constants to ensure correct arithmetic

  nbits = nlev+1
  do j = 1,npp
     local_key(j) = iplace + &
          SUM( (/ (8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) ),i=0,nbits-1) /) )
  end do


  if (domain_debug) then

     write (ipefile,'(/a/a/(z21,i8,5f12.4))') 'Particle list before key sort:', &
          '  key,             label   coords     q ', &
          (local_key(i),pelabel(i),x(i),y(i),z(i),q(i),work(i),i=1,npp) 

     call blankn(ipefile)
  endif

  ! Use Parallel Sort by Regular Sampling (PSRS) 

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  if (domain_debug) then
     write (ipefile,*) 'MPI psrssort() commencing'
     write (ipefile,*) 'iproc=',me
     write (ipefile,*) 'num_pe=',num_pe
     write (ipefile,*) 'npp=',npp
  endif

  iteration = 0
  niterations = 3  ! Max # iterations
  errcount = 1
  npold = npp


  ! start permutation of local key list
  do i=1,npp
     xarray(i) = local_key(i)
  enddo

  do while (errcount /= 0 .or. iteration == niterations)
     iteration = iteration + 1
     write (ipefile,*) ' '

     if (iteration .lt. niterations) then
        do i=1,npold
           keys(i) = xarray(i)
        enddo
     else
        if (domain_debug .and. errcount .eq. 0) then
           write (ipefile,*) 'This time we are using already sorted data'
        endif
        npold = npnew
        do i=1,npold
           xarray(i) = w1(i)
           keys(i) = w1(i)
        enddo
     endif

  if (domain_debug) then
          write (ipefile,'(a/(i5,z20))') 'input array:  ',(i,keys(i),i=1,npold)
          write (ipefile,*) 'npold=',npold
 endif

     ! perform index sort on keys
     call pll_weightsort(nppm,npold,npnew,num_pe,me,keys,indxl,irnkl, &
          islen,irlen,fposts,gposts,w1,work,load_balance)


     do i=1,npold
        w1(i) = xarray(i)
     enddo

     call pll_permute(nppm,npold,npnew,num_pe,me,w1,wi2,wi3, &
          indxl,irnkl,islen,irlen,fposts,gposts)

  if (domain_debug) then
          write (ipefile,'(a/(i5,z20))') 'output array: ',(i,w1(i),i=1,npnew)
          write (ipefile,*) 'npnew=',npnew
  endif

     ! Check if sort finished
     errcount = 0
     do i=2,npnew
        if (w1(i) .lt. w1(i-1)) then
           errcount = errcount + 1
           if (domain_debug .and. errcount .lt. 10) then
              write (ipefile,'(a,i5,2z20)') 'i,w1(i),w1(i-1)=',i,w1(i),w1(i-1)
           endif
        endif
     enddo


     ! swap end items
     if (me .ne. 0) then
        call MPI_SEND(w1, one, MPI_INTEGER8, me_minus_one, tag1, MPI_COMM_WORLD, ierr)
     endif

     if (me .ne. num_pe-1) then
        call MPI_RECV(tmp, one, MPI_INTEGER8, me_plus_one, tag1, MPI_COMM_WORLD, status, ierr)
     endif

     if (me .ne. num_pe-1) then
        if (domain_debug .and. tmp .lt. w1(npnew)) then          ! still something to sort
           write (ipefile,'(a,i3,a1,2z20)') 'w1(npnew), w1(1) from',me+1, '=',w1(npnew),tmp
           errcount = errcount + 1  
        endif
     endif

     if (errcount .ne. 0 .and. domain_debug) then
        write (ipefile,*) 'errcount=',errcount
     endif

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  enddo



  npp = npnew
  pekey(1:npp) = w1(1:npp)

  ! Check for identical keys

  do i=2,npp
     if (pekey(i) == pekey(i-1)) then
        write(*,*) 'WARNING: identical keys found - increase max # levels'
        pekey(i) = pekey(i) + 1  ! Augment higher key
     endif
  end do

  ! Now permute remaining particle properties : x,y,z; vx,vy,vz; q,m, label, load

  source_pe(1:npold) = pepid(1:npold)   ! where particle came from


  call pll_permute(nppm,npold,npnew,num_pe,me,source_pe,w2,w3, &   ! source PE
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,pelabel,w2,w3, &   ! label
       indxl,irnkl,islen,irlen,fposts,gposts)

  ! Should ship these as one big vector
  call pll_permute(nppm,npold,npnew,num_pe,me,x,wr2,wr3, &       ! coords
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,y,wr2,wr3, &       !  
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,z,wr2,wr3, &       !  
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,ux,wr2,wr3, &       ! velocities
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,uy,wr2,wr3, &       !
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,uz,wr2,wr3, &       ! 
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,q,wr2,wr3, &       ! charge
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,m,wr2,wr3, &       ! mass
       indxl,irnkl,islen,irlen,fposts,gposts)

  call pll_permute(nppm,npold,npnew,num_pe,me,work,wr2,wr3, &       ! workload
       indxl,irnkl,islen,irlen,fposts,gposts)


  pepid(1:npp) = me  ! new owner




  if (domain_debug) then
     do j=1,npp
        ixd(j) = SUM( (/ (2_8**i*ibits( pekey(j)-iplace,3*i,1 ), i=0,nbits-2) /) )
        iyd(j) = SUM( (/ (2_8**i*ibits( pekey(j)-iplace,3*i+1,1 ), i=0,nbits-2) /) )
        izd(j) = SUM( (/ (2_8**i*ibits( pekey(j)-iplace,3*i+2,1 ), i=0,nbits-2) /) )
     end do
     write (ipefile,'(/a/a/(z21,2i6,a2,i8,6f12.4))') 'Particle list after key sort:', &
          '  key,                owner,    from PE  |  label  Fetched coords      derived from key', &
          (pekey(i),pepid(i),source_pe(i),'|', &
          pelabel(i),x(i),y(i),z(i),ixd(i)*s+xmin,iyd(i)*s+ymin,izd(i)*s+zmin,i=1,npp) 

     call blankn(ipefile)
  endif


  ! Each PE now has sorted segment of particles of length npp
  ! Note that now npp /= npart/num_pe, only approx depending on key distribution, or target shape.


  ! Copy boundary particles to adjacent PEs to ensure proper tree construction
  !  - if we don't do this, can get two particles on separate PEs 'sharing' a leaf

  ship_props = particle ( x(1), y(1), z(1), ux(1), uy(1), uz(1), q(1), m(1), work(1), &
   			  pekey(1), pelabel(1), pepid(1) )

!  write (*,'(9f12.3,z20,2i6)') ship_props
  
  ! Ship 1st particle data to end of list of LH neighbour PE

  if (me /= 0 ) then
     call MPI_SEND( ship_props, 1, mpi_type_particle, me_minus_one, tag1, MPI_COMM_WORLD, ierr ) 
  endif

  ! Place incoming data at end of array
  if ( me /= lastpe) then
     call MPI_RECV( get_props, 1, mpi_type_particle, me_plus_one, tag1,  MPI_COMM_WORLD, status, ierr )
     x(npp+1) = get_props%x
     y(npp+1) = get_props%y
     z(npp+1) = get_props%z
     ux(npp+1) = get_props%ux
     uy(npp+1) = get_props%uy
     uz(npp+1) = get_props%uz
     q(npp+1) = get_props%q
     m(npp+1) = get_props%m
     work(npp+1) = get_props%work
     pekey(npp+1) = get_props%key
     pelabel(npp+1) = get_props%label
     pepid(npp+1) = get_props%pid
  endif


  ! Ship  end particle data to start of list of RH neighbour PE
  
  ship_props = particle ( x(npp), y(npp), z(npp), ux(npp), uy(npp), uz(npp), q(npp), m(npp), work(npp), &
   			  pekey(npp), pelabel(npp), pepid(npp) )
  
if (me /= lastpe ) then
     call MPI_SEND( ship_props, 1, mpi_type_particle, me_plus_one, tag1, MPI_COMM_WORLD, ierr ) 
  endif

  ! Place incoming data at end of array
  
  if (me == lastpe) then
     ind_recv = npp+1   ! PEn array hasn't yet received boundary value
  else
     ind_recv = npp+2
  endif

  if ( me /= 0) then
     call MPI_RECV( get_props, 1, mpi_type_particle, me_minus_one, tag1,  MPI_COMM_WORLD, status, ierr )
     x(ind_recv) = get_props%x
     y(ind_recv) = get_props%y
     z(ind_recv) = get_props%z
     ux(ind_recv) = get_props%ux
     uy(ind_recv) = get_props%uy
     uz(ind_recv) = get_props%uz
     q(ind_recv) = get_props%q
     m(ind_recv) = get_props%m
     work(ind_recv) = get_props%work
     pekey(ind_recv) = get_props%key
     pelabel(ind_recv) = get_props%label
     pepid(ind_recv) = get_props%pid
  endif



  if (boundary_debug) then
     if (me /= 0 .and. me /= lastpe) then
        inc = 1
     else
        inc = 0
     endif



     write (ipefile,'(/a/a/(i5,z21,2i8,3f12.5,2f12.3))') 'Particle list after boundary swap:', &
          ' index   key,     label,   on PE,    x      y     q       m', &
          (i,pekey(i),pelabel(i),pepid(i),x(i),y(i),z(i),q(i),m(i),i=1,npp+1+inc) 

     call blankn(ipefile)
  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

end subroutine make_domains
