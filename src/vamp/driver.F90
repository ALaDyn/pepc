
program driver

  use utils
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none
  include 'mpif.h'

  integer, parameter :: mprocs=3, &
       !       nbits = 6, &
!VAMPINST program_start
      CALL VTdefs()
      CALL VTENTER(IF_driver,VTNOSCL,VTIERR)
      write(*,*) 'VT: driver P>',VTIERR, &
         IF_driver,ICLASSH
!
  npg = 36, &
  mp=npg+64

  integer :: np=1+(npg-1)/mprocs

  integer :: xarray(mp),keys(mp),w1(mp)
  integer :: indxl(mp),irnkl(mp)

  integer :: islen(mprocs),irlen(mprocs)
  integer :: fposts(mprocs+1),gposts(mprocs+1)

  integer :: iseed
  real :: x

  integer :: i,npnew,npold,nprocs,iproc
  integer :: iteration, niterations, fd, ierr
  integer :: errcount

  integer :: stat(MPI_STATUS_SIZE), w2(mp), w3(mp)
  integer :: tmp
  logical :: example=.true.
  save


  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  fd = iproc+10

  if (nprocs .ne. mprocs) then
     write (fd,*) 'error: nprocs != mprocs'
     write (fd,*) 'nprocs,mprocs=', nprocs, mprocs
     call MPI_ABORT(MPI_COMM_WORLD, 136, ierr)
  endif
  np = np - 1 + iproc

  write (fd,*) 'MPI psrssort() commencing'
  write (fd,*) 'iproc=',iproc
  write (fd,*) 'nprocs=',nprocs
  write (fd,*) 'np=',np

  if (nprocs .lt. mprocs) then
     write (fd,*) ' Whoops! nprocs .lt. mprocs'
     write (fd,*) ' mprocs, nprocs =', mprocs, nprocs
     call MPI_ABORT(MPI_COMM_WORLD, 138, ierr)
  endif


  if (mp*nprocs .lt. npg) then
     write (fd,*) 'Whoops! Storage space per processor', &
          ' is too small'
     write (fd,*) 'npg, mp, nprocs =',npg,mp,nprocs
     call MPI_ABORT(MPI_COMM_WORLD, 139, ierr)
  endif

  niterations = 2
  npold = np

  iseed = 2019
  iseed = -iseed + iproc

  if (example) then  ! Data from example in psrs paper
     if (iproc==0) then
        xarray(1:np) = (/ 16,2,17,24,33,28,30,1,0,27,9 /)
     else if (iproc==1) then
        xarray(1:np) = (/ 7,8,11,12,13,18,19,21,23,29,34,35 /)
     else if (iproc ==2) then
        xarray(1:np) = (/ 3,4,5,6,10,14,15,20,22,26,31,32,25 /)        
     endif

  else
! random start array
     do i=1,np
        x = rano(iseed)
        write(fd,*) x
        xarray(i) = np*nprocs*x
     enddo
  endif

  do iteration=1,niterations
     write (fd,*) ' '

     if (iteration .lt. niterations) then
        do i=1,npold
           keys(i) = xarray(i)
        enddo
     else
        if (errcount .eq. 0) then
           write (fd,*) 'This time we are using already sorted data'
        endif
        npold = npnew
        do i=1,npold
           xarray(i) = w1(i)
           keys(i) = w1(i)
        enddo
     endif

     write (fd,*) 'input array:  ',(keys(i),i=1,npold)
     write (fd,*) 'npold=',npold


     call pll_indexsort(npold,npnew,nprocs,iproc,keys,indxl,irnkl, &
          islen,irlen,fposts,gposts,w1)


     write (fd,*) 'npold= ',npold, 'npnew=',npnew
     write (fd,*) 'keys:',keys(1:npold)
     write (fd,*) 'local index: ',indxl(1:npold)
     write (fd,*) '3rd part',irnkl(1:npold)
     write (fd,*) 'fposts: ',fposts(1:nprocs+1)
     write (fd,*) 'gposts: ',gposts(1:nprocs+1)
     write (fd,*) 'slen:',islen(1:nprocs)
     write (fd,*) 'rlen:',irlen(1:nprocs)
 
     if (npnew .gt. mp) then
        write (fd,*) 'error: npnew > mp; npnew,mp=',npnew,mp
        call MPI_ABORT(MPI_COMM_WORLD, 238, ierr)
     endif


     do i=1,npold
        w1(i) = xarray(i)
     enddo

     call pll_permute(npold,npnew,nprocs,iproc,w1,w2,w3, &
          indxl,irnkl,islen,irlen,fposts,gposts)

     write (fd,*) 'output array: ',(w1(i),i=1,npnew)


! Check if sort finished
     errcount = 0
     do i=1,npnew
        if (w1(i) .lt. w1(i-1)) then
           errcount = errcount + 1
           if (errcount .lt. 10) then
              write (fd,*) 'i,w1(i),w1(i-1)=',i,w1(i),w1(i-1)
           endif
        endif
     enddo


     if (iproc .ne. 0) then
        call MPI_SEND(w1, 1, MPI_INTEGER, iproc-1, 0, MPI_COMM_WORLD, ierr)
     endif

     if (iproc .ne. nprocs-1) then
        call MPI_RECV(tmp, 1, MPI_INTEGER, iproc+1, 0, MPI_COMM_WORLD, stat, ierr)
     endif

     if (iproc .ne. nprocs-1) then
        if (tmp .lt. w1(npnew)) then          ! still something to sort
           write (fd,*) 'w1(npnew), w1(1) from',iproc+1, '=',w1(npnew),tmp
           errcount = errcount + 1  
        endif
     endif

     if (errcount .ne. 0) then
        write (fd,*) 'errcount=',errcount
     endif

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  enddo

  call MPI_FINALIZE(MPI_COMM_WORLD, ierr)

!VAMPINST program_end
      CALL VTLEAVE(ICLASSH,VTIERR)
      write(*,*) 'VT: driver P<',VTIERR,ICLASSH
!
end program driver
