
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


subroutine tree_domains(xl,yl,zl)


  use treevars
  use tree_utils
  use utils
  implicit none
  include 'mpif.h'

  integer, parameter :: nkmax=100
  real, intent(in) :: xl,yl,zl  ! initial box limits
  integer*8, dimension(nppm) :: ix, iy, iz
  integer*8, dimension(nppm) :: ixd, iyd, izd
  integer*8, dimension(nppm) :: local_key
  integer*8, dimension(2) :: ixbox, iybox, izbox, key_box

  integer ::  source_pe(nppm)
  integer :: i, j, ind_recv, inc, prev, next, handle(4)


  integer :: nbits


  real*8 :: s
  real*8 :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
  logical :: boundary_debug=.false. 

  integer status(MPI_STATUS_SIZE), ierr, tag1

  ! arrays for parallel sort

  type (particle) :: ship_parts(nppm), get_parts(nppm)

  integer*8 :: xarray(nppm),keys(nppm),w1(nppm),wi2(nppm),wi3(nppm)
  integer :: indxl(nppm),irnkl(nppm)

  integer :: islen(num_pe),irlen(num_pe)
  integer :: fposts(num_pe+1),gposts(num_pe+1)

  integer :: npnew,npold
  integer :: iteration, niterations
  integer :: errcount

  integer, dimension(nppm) ::  w2, w3 ! scratch arrays for integer*4 permute
  real, dimension(nppm) :: wr2, wr3 ! Scratch for real array permute
  integer*8 :: tmp
  logical :: sort_debug=.false.
  real*8 xboxsize, yboxsize, zboxsize

  !POMP$ INST BEGIN(keys)

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | DOMAINS'


  ! Find limits of local simulation region
  xmin_local = minval(x(1:npp))
  xmax_local = maxval(x(1:npp))
  ymin_local = minval(y(1:npp))
  ymax_local = maxval(y(1:npp))
  zmin_local = minval(z(1:npp))
  zmax_local = maxval(z(1:npp))

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

  ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)

  ix(1:npp) = ( x(1:npp) - xmin )/s           ! partial keys
  iy(1:npp) = ( y(1:npp) - ymin )/s           !
  iz(1:npp) = ( z(1:npp) - zmin )/s        

  ! construct keys by interleaving coord bits and add placeholder bit
  ! - note use of 64-bit constants to ensure correct arithmetic

  nbits = nlev+1
  do j = 1,npp
     !     local_key(j) = iplace + &
     !          SUM( (/ (8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) ),i=0,nbits-1) /) )
     local_key(j) = iplace
     do i=0,nbits-1
        local_key(j) = local_key(j) &
             + 8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) )
     end do
  end do

  !  Find keys corresponding to corners of graphics box (initial target container)
  if (xmin<0) then
     ixbox(1) = -xmin/s
     ixbox(2) = (xl - xmin)/s
  else
     ixbox(1) = 0
     ixbox(2) = (xmax-xmin)/s
  endif

  if (ymin<0) then
     iybox(1) = -ymin/s
     iybox(2) = (yl - ymin)/s
  else
     iybox(1) = 0
     iybox(2) = (ymax - ymin)/s
  endif

  if (zmin<0) then
     izbox(1) =-zmin/s
     izbox(2) = (zl - zmin)/s
  else
     izbox(1) = 0
     izbox(2) = (zmax - zmin)/s
  endif

  do j=1,2
     key_box(j) = iplace + &
          SUM( (/ (8_8**i*(4_8*ibits( izbox(j),i,1) &
          + 2_8*ibits( iybox(j),i,1 ) &
          + 1_8*ibits( ixbox(j),i,1) ),i=0,nbits-1) /) )
  end do


  if (domain_debug) then
     write (ipefile,'(a,2z20)') 'Box keys:', key_box(1),key_box(2)
     write (ipefile,'(/a/a/(z21,i8,3f12.4,3i8,2f12.4))') 'Particle list before key sort:', &
          '  key,             label   coords     q ', &
          (local_key(i),pelabel(i),x(i),y(i),z(i),ix(i),iy(i),iz(i),q(i),work(i),i=1,npp) 
     !    write (ipefile,'(/a/a/(z21,i8,3f12.4,3i8))') '(last 10):', &
     !         '  key,                  label        coords              q ', &
     !         (local_key(i),pelabel(i),x(i),y(i),z(i),ix(i),iy(i),iz(i),q(i),work(i),i=max(1,npp-10),npp) 

     call blankn(ipefile)
  endif

  !POMP$ INST END(keys)

  !POMP$ INST BEGIN(sort)

  ! Use Parallel Sort by Regular Sampling (PSRS) 

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  if (domain_debug) then
     write (*,*) 'MPI psrssort() commencing'
     write (*,*) 'iproc=',me
     write (*,*) 'num_pe=',num_pe
     write (*,*) 'npp=',npp
  endif

  iteration = 0
  niterations = 3  ! Max # iterations
  errcount = 1
  npold = npp
  npnew = npp


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

     if (domain_debug.and.me==0) then
        write (*,'(a/(i5,z20,f15.3))') 'input array (1st 10):  ',(i,keys(i),work(i),i=1,min(10,npold))
        write (*,'(a/(i5,z20,f15.3))') 'input array (last 10):  ',(i,keys(i),work(i),i=max(1,npold-10),npold)
        write (*,*) 'npold=',npold
     endif

     ! perform index sort on keys

     call pswssort(nppm,npold,npnew,num_pe,me,keys, &
          indxl,irnkl,islen,irlen,fposts,gposts,w1,work,key_box,load_balance,sort_debug)


     do i=1,npold
        w1(i) = xarray(i)
     enddo

     ! permute keys according to sorted indices
     call pll_permute(nppm,npold,npnew,num_pe,me,w1,wi2,wi3, &
          indxl,irnkl,islen,irlen,fposts,gposts)

     if (domain_debug) then
        write (*,*) 'npnew=',npnew,' on ',me
        !         write (ipefile,'(a/(i5,z20))') 'output array (1st 10): ',(i,w1(i),i=1,10)
        !         write (ipefile,'(a/(i5,z20))') 'output array (last 10): ',(i,w1(i),i=npnew-10,npnew)
        if (me.eq.50) then
           write (*,'(a/(i5,z20))') 'output array (1st 10): ',(i,w1(i),i=1,min(10,npnew))
           write (*,'(a/(i5,z20))') 'output array (last 10): ',(i,w1(i),i=max(1,npnew-10),npnew)
        endif
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

     ! swap end items
     if (me .ne. 0) then
        call MPI_ISEND(w1, 1, MPI_INTEGER8, prev, 1, MPI_COMM_WORLD, handle(1), ierr)
        call MPI_REQUEST_FREE(handle(1),ierr)
     endif

     if (me .ne. num_pe-1) then
        call MPI_RECV(tmp, 1, MPI_INTEGER8, next, 1, MPI_COMM_WORLD, status, ierr)
     endif


     if (me .ne. num_pe-1) then
        !    if (me == 50 ) then
        if (domain_debug .and. tmp .lt. w1(npnew)) then          ! still something to sort
           write (*,'(a,i3,a1,2z20)') 'w1(npnew), w1(1) from',me+1, '=',w1(npnew),tmp
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
        write(*,'(a,o21)') 'WARNING: identical keys found:  ',pekey(i)
        pekey(i) = pekey(i) + 1  ! Augment higher key
        write(*,'(a,o21)') 'Upper key increased to:  ',pekey(i)
     endif
  end do

  ! Now permute remaining particle properties : x,y,z; vx,vy,vz; q,m, label, load

  source_pe(1:npold) = pepid(1:npold)   ! where particle came from

  ! Set up particle structure - keys and source_pe are dummies
  ! ( pekey is already sorted)

  do i=1,npold
     ship_parts(i) = particle( x(indxl(i)), y(indxl(i)), z(indxl(i)), &
          ux(indxl(i)), uy(indxl(i)), uz(indxl(i)), &
          q(indxl(i)), m(indxl(i)), work(indxl(i)), &
          ax(indxl(i)), ay(indxl(i)), az(indxl(i)), &
          keys(indxl(i)), pelabel(indxl(i)), source_pe(indxl(i))    )
  enddo


  ! perform permute
  call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_particle, &
       get_parts, irlen, gposts, mpi_type_particle, &
       MPI_COMM_WORLD,ierr)

  do i=1,npp
     x(irnkl(i)) = get_parts(i)%x
     y(irnkl(i)) = get_parts(i)%y
     z(irnkl(i)) = get_parts(i)%z
     ux(irnkl(i)) = get_parts(i)%ux
     uy(irnkl(i)) = get_parts(i)%uy
     uz(irnkl(i)) = get_parts(i)%uz
     q(irnkl(i)) = get_parts(i)%q
     m(irnkl(i)) = get_parts(i)%m
     work(irnkl(i)) = get_parts(i)%work
     ax(irnkl(i)) = get_parts(i)%ax
     ay(irnkl(i)) = get_parts(i)%ay
     az(irnkl(i)) = get_parts(i)%az
     pelabel(irnkl(i)) = get_parts(i)%label
  enddo


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
       ax(1),ay(1),az(1), pekey(1), pelabel(1), pepid(1) )

  !  write (*,'(9f12.3,z20,2i6)') ship_props

  ! Ship 1st particle data to end of list of LH neighbour PE

  if (me /= 0 ) then
     call MPI_ISEND( ship_props, 1, mpi_type_particle, prev, 1, MPI_COMM_WORLD, handle(1), ierr ) 
     call MPI_REQUEST_FREE(handle(1),ierr) 
  endif

  ! Place incoming data at end of array
  if ( me /= num_pe-1) then
     !     call MPI_IRECV( get_props, 1, mpi_type_particle, next, 1,  MPI_COMM_WORLD, handle(2), ierr )
     call MPI_RECV( get_props, 1, mpi_type_particle, next, 1,  MPI_COMM_WORLD, status, ierr )
     x(npp+1) = get_props%x
     y(npp+1) = get_props%y
     z(npp+1) = get_props%z
     ux(npp+1) = get_props%ux
     uy(npp+1) = get_props%uy
     uz(npp+1) = get_props%uz
     q(npp+1) = get_props%q
     m(npp+1) = get_props%m
     work(npp+1) = get_props%work
     ax(npp+1) = get_props%ax
     ay(npp+1) = get_props%ay
     az(npp+1) = get_props%az
     pekey(npp+1) = get_props%key
     pelabel(npp+1) = get_props%label
     pepid(npp+1) = get_props%pid
  endif



  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !  if (me /= num_pe-1) call MPI_WAIT(handle(2),status,ierr)


  ! Ship  end particle data to start of list of RH neighbour PE

  ship_props = particle ( x(npp), y(npp), z(npp), ux(npp), uy(npp), uz(npp), q(npp), m(npp), work(npp), &
       ax(npp),ay(npp),az(npp), pekey(npp), pelabel(npp), pepid(npp) )

  if (me /= num_pe-1 ) then
     call MPI_ISEND( ship_props, 1, mpi_type_particle, next, 2, MPI_COMM_WORLD, handle(3), ierr )
     call MPI_REQUEST_FREE(handle(3),ierr) 
  endif

  ! Place incoming data at end of array

  if (me == num_pe-1) then
     ind_recv = npp+1   ! PEn array hasn't yet received boundary value
  else
     ind_recv = npp+2
  endif

  if ( me /= 0) then
     !     call MPI_IRECV( get_props, 1, mpi_type_particle, prev, 2,  MPI_COMM_WORLD, handle(4), ierr )
     call MPI_RECV( get_props, 1, mpi_type_particle, prev, 2,  MPI_COMM_WORLD, status, ierr )
     x(ind_recv) = get_props%x
     y(ind_recv) = get_props%y
     z(ind_recv) = get_props%z
     ux(ind_recv) = get_props%ux
     uy(ind_recv) = get_props%uy
     uz(ind_recv) = get_props%uz
     q(ind_recv) = get_props%q
     m(ind_recv) = get_props%m
     work(ind_recv) = get_props%work
     ax(ind_recv) = get_props%ax
     ay(ind_recv) = get_props%ay
     az(ind_recv) = get_props%az
     pekey(ind_recv) = get_props%key
     pelabel(ind_recv) = get_props%label
     pepid(ind_recv) = get_props%pid
  endif


  !  if (me /= 0) call MPI_WAIT(handle(4),status,ierr)

  if (boundary_debug) then
     if (me /= 0 .and. me /= num_pe-1) then
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
  !POMP$ INST END(sort)

end subroutine tree_domains
