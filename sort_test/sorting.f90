subroutine sorting

  use physvars
  use tree_utils
  implicit none
  include 'mpif.h'

  logical :: sort_debug

  integer :: indxl(nppm),irnkl(nppm)
  integer :: islen(n_cpu),irlen(n_cpu)
  integer :: fposts(n_cpu+1),gposts(n_cpu+1)
  integer :: npnew,npold,npp,iteration,niterations,errcount

  integer*8 :: xarray(nppm),keys(nppm),w1(nppm),wi2(nppm),wi3(nppm),compare(nppm)
  real*8 :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
  real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: xboxsize, yboxsize, zboxsize, boxsize, s
  integer*8, dimension(nppm) :: ix, iy, iz
  integer :: nbits, nlev
  integer :: j, ierr, prev, next, handle(4), status(MPI_STATUS_SIZE)
  integer*8 :: iplace, tmp
  integer*8, dimension(nppm) :: local_key
  integer*8, dimension(2) :: ixbox, iybox, izbox, key_box
  integer*8, dimension(n_cpu+1)::  pivots
  integer :: load_balance   ! Balances particles in || sort according to work load
  real*8 :: t1, t2
  real :: work_local
  integer*8 :: k,i

!  if ((db_level > 1) .and. (my_rank == 0)) then
!     write(*,*) '-= Sorting =-'
!     write(*,*) 'memory consumption = ',(9*size(indxl)*8+5*size(islen)*8)/1024/1024,'MiB'
!  end if

!  if ((db_level > 1) .and. (my_rank == 0)) then
!     write(*,*) '-= pbalsort =-'
!     write(*,*) 'memory consumption = ',(4*size(local_key)*8+4*size(indxl)*4+14*1000000*4+9*size(pivots)*4)/1024/1024,'MiB'
!  end if

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  if (my_rank==0) write(*,'(a)') ' Generating keys...'
  
  npp = np_local
  sort_debug=.false.
  work_local = 1
  load_balance = 0
  nlev = 20                     ! max refinement level
  iplace = 2_8**(3*nlev)           ! place holder bit

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

  s=boxsize/2**nlev       ! refinement length

  ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)

  ix(1:npp) = ( x(1:npp) - xmin )/s           ! partial keys
  iy(1:npp) = ( y(1:npp) - ymin )/s           !
  iz(1:npp) = ( z(1:npp) - zmin )/s           !

  ! construct keys by interleaving coord bits and add placeholder bit
  ! - note use of 64-bit constants to ensure correct arithmetic

  nbits = nlev+1
  do j = 1,npp
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

  if (my_rank==0) write(*,'(a)') ' Sorting keys...'

  iteration = 0
  niterations = 3  ! Max # iterations
  errcount = 1

  npold = npp
  npnew = npp

  do i=1,npp
     xarray(i) = local_key(i)
!     write(*,'(2i4,o30)') my_rank,i,xarray(i)
  enddo

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first
  t1 = MPI_WTIME()

  do while (errcount /= 0 .or. iteration == niterations)
    
     if (my_rank == 0) write(*,*) 'Starting iteration no.',iteration,'...'

     iteration = iteration + 1

     if (iteration .lt. niterations) then
        do i=1,npold
           keys(i) = xarray(i)
        enddo
     else
        npold = npnew
        do i=1,npold
           xarray(i) = w1(i)
           keys(i) = w1(i)
        enddo
     endif

!     call pswssort(nppm,npold,npnew,n_cpu,my_rank,keys,&
!                   indxl,irnkl,islen,irlen,fposts,gposts,w1,work,key_box,load_balance,sort_debug)
!     call psrssort(nppm,npold,npnew,n_cpu,my_rank,keys,&
!                   indxl,irnkl,islen,irlen,fposts,gposts,w1)
     call pbalsort(nppm,npold,npnew,n_cpu,my_rank,keys,&
                   indxl,irnkl,islen,irlen,fposts,gposts,pivots,w1,work,npp,load_balance,sort_debug,work_local)

     do i=1,npold
        w1(i) = xarray(i)
     enddo

     call pll_permute(nppm,npold,npnew,n_cpu,my_rank,w1,wi2,wi3, &
                      indxl,irnkl,islen,irlen,fposts,gposts)

     ! Check if sort finished
     errcount = 0
     do i=2,npnew
        if (w1(i) .lt. w1(i-1)) then
           errcount = errcount + 1
        endif
     enddo

     ! Define wraps for ring network  0 -> 1 -> 2 -> ... ... -> num_pe-1 -> 0 ...
     if (my_rank == 0) then
        prev = n_cpu - 1
     else
        prev = my_rank-1
     endif

     if (my_rank == n_cpu-1 ) then
        next = 0
     else
        next = my_rank+1
     endif
     
     ! swap end items
     if (my_rank .ne. 0) then
        call MPI_ISEND(w1, 1, MPI_INTEGER8, prev, 1, MPI_COMM_WORLD, handle(1), ierr)
        call MPI_REQUEST_FREE(handle(1),ierr)
     endif


     if (my_rank .ne. n_cpu-1) then
        call MPI_RECV(tmp, 1, MPI_INTEGER8, next, 1, MPI_COMM_WORLD, status, ierr)
     endif


     if (my_rank .ne. n_cpu-1) then
        if (tmp .lt. w1(npnew)) then          ! still something to sort
           errcount = errcount + 1  
        endif
     endif

     if (errcount .ne. 0) then
        if (my_rank == 0) write (*,*) 'errcount=',errcount,'=> Still something to do...'
     else
        if (my_rank == 0) write (*,*) 'errcount=',errcount,'=> Everything is fine! Cleaning up ...\n'
     endif

     call MPI_BARRIER( MPI_COMM_WORLD, ierr) 

  end do

  t2 = MPI_WTIME()

  do i=1,npnew
     keys(i) = w1(i)
  enddo

!  write(*,*) my_rank, k
  if (my_rank == 0) write(*,*) "Time:",t2-t1

end subroutine sorting
