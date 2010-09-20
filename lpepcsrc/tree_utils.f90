!
!                          Tree sort utility module
!                       *   Parallel Sorting
!

module tree_utils

  ! Interfaces
  interface pll_regsort
     module procedure psrssort
  end interface

  interface pll_weightsort
     module procedure pswssort
  end interface

  interface pll_balsort
     module procedure pbalsort
  end interface

  interface pll_permute
     module procedure psrsperm_i4
     module procedure psrsperm_i8
     module procedure psrsperm_r8
  end interface

  interface nwaymerge
     module procedure nwaymrg
  end interface

  interface sort
     module procedure sort_i
  end interface

  interface indexsort
     module procedure indsort_i8
     module procedure indsort_i4
  end interface

  interface swap
     module procedure swap_ab
  end interface

  interface unique
     module procedure uniq8
  end interface

  interface locaddress
     module procedure locaddress_int4_8
     module procedure locaddress_int8_8
     module procedure locaddress_real8_8
  end interface

  interface bpi
     module procedure bpi_int8_8
  end interface

  interface bpi_bits
     module procedure bpi_bits_int8_8
  end interface

contains

! ================
! permute integer*4
! ================

  subroutine psrsperm_i4(nppm,np,npnew,nprocs,iproc,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)

    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs,iproc

    integer, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) :: indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts

    integer :: i,ierr

    do i=1,np
       w1(i) = array(indxl(i))
    enddo

    call MPI_ALLTOALLV(  w1, islen, fposts, MPI_INTEGER, &
                         w2, irlen, gposts, MPI_INTEGER, &
                         MPI_COMM_WORLD,ierr)

    do i=1,npnew
       array(irnkl(i)) = w2(i)
    enddo

  end subroutine psrsperm_i4


! ================
! permute integer*8
! ================

  subroutine psrsperm_i8(nppm,np,npnew,nprocs,iproc,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)


    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs,iproc

    integer*8, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) :: indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts

    integer :: i,ierr


    do i=1,np
       w1(i) = array(indxl(i))
    enddo

    call MPI_ALLTOALLV(  w1, islen, fposts, MPI_INTEGER8, &
                         w2, irlen, gposts, MPI_INTEGER8, &
                         MPI_COMM_WORLD,ierr)

    do i=1,npnew
       array(irnkl(i)) = w2(i)
    enddo

  end subroutine psrsperm_i8


! ================
  !  Permute REAL*8 array()
! ================

  subroutine psrsperm_r8(nppm,np,npnew,nprocs,iproc,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)


    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs,iproc

    real*8, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) ::  indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts


    integer :: i,ierr


    do i=1,np
       w1(i) = array(indxl(i))
    enddo


    call MPI_alltoallv(  w1,islen, fposts, MPI_REAL8, &
         w2, irlen, gposts, MPI_REAL8, &
         MPI_COMM_WORLD,ierr)

    do i=1,npnew
       array(irnkl(i)) = w2(i)
    enddo

  end subroutine psrsperm_r8





! unweighted || sort

  subroutine psrssort(nppm,np,npnew,nprocs,iproc,keys,indxl,irnkl,islen,irlen,fposts,gposts,w1)


    ! Xiaobo Li, Paul Lu, Jonathan Schaeffer,
    ! John Shillington, Pok Sze Wong, Hanmao Shi,
    ! "On the versatility of parallel sorting by regular sampling",
    ! Parallel Computing, volume 19, 1079--1103, Oct 1993.
    ! http://www.cs.utoronto.ca/~paullu/Papers/psrs.ps.Z
    !
    ! sorting articles:
    ! http://www.lpac.ac.uk/SEL-HPC/Articles/GeneratedHtml/hpc.sort.html
    !
    !
    !     (C. Mobarry/GSFC, J. Crawford/VSEP)
    !
    !     Convex SPP-1000 Exemplar MPI implementation
    !     Written by C. Mobarry and J. Crawford
    !     NASA/GSFC Code 934, Greenbelt MD, 20771
    !     Clark.Mobarry@gsfc.nasa.gov, {mobarry,crawford}@maxwell.gsfc.nasa.gov

    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs,iproc
    integer*8, dimension(nppm) ::  keys, &      ! array of keys to be sorted.
                                   w1       ! work array
    integer, dimension(nppm) ::    indxl, irnkl ! origin locations of the keys 


    integer, dimension(nprocs) :: islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle

    integer, parameter :: maxprocs = 1024
    integer :: itabr(maxprocs), itabl(maxprocs+1)

    integer*8 :: work(maxprocs*(maxprocs+1))
    integer*8 :: fpval(maxprocs+1)
    integer*8 :: lmax, ak

    integer*8 :: buf(maxprocs)

    integer :: i,j,k,step, tag1, ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: fd, nsamp
    character(13) :: cfmt

    fd = 29

    tag1=0


    !     Independent s  !     Note that indx() is a local index on the process.


    call indexsort( keys,indxl, np, nppm )   ! Index sort from Num. Rec.

    do i=1,np
       w1(i) = keys(indxl(i))
    enddo

    !     w1() is now the sorted keys.
    !     indxl() is now the local indexes for the sort.

    lmax = w1(np)

    !     Choose nproc evenly spaced values out of every bin.
    !     Store them all in work().

    k = 1
    step=np/nprocs

    do i=1,nprocs
       work(i) = w1(k)
       k = k + step
    enddo

    !     work(1:nprocs) are the sampled key values for the fenceposts.

    tag1 = tag1 + 1

    if (iproc .ne. 0) then
       call MPI_SEND(work, nprocs, MPI_INTEGER8, 0, tag1, MPI_COMM_WORLD, ierr)
    else
       do i=1,nprocs-1
          call MPI_RECV(buf, nprocs, MPI_INTEGER8, MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, status, ierr)
          do j=1,nprocs
             work(status(MPI_SOURCE)*nprocs+j) = buf(j)
          enddo
       enddo
    endif


    if (iproc .eq. 0) then
       !        Insertion sort the fencepost values of keys and indexes.
       do i=2,nprocs*nprocs
          ak = work(i)
          do j=i,2,-1
             if (work(j-1) .le. ak) goto 300
             work(j) = work(j-1)
          enddo
          j = 1
300       continue
          work(j) = ak
       enddo
       !        After the insertion sort.
       !        work() are the sorted sampled key values for the fenceposts.
       nsamp=nprocs*nprocs
       cfmt = "(/a15,"//achar(mod(nsamp/10,10)+48) // achar(mod(nsamp,10)+48) // "(i4))"
!       write(fd,cfmt) 'Sorted sample: ',(work(i),i=1,nsamp)

       k = 1
       do i=1,nprocs*nprocs,nprocs
          fpval(k) = work(i)
          k = k + 1
       enddo
    endif

    call MPI_BCAST(fpval, nprocs+1, MPI_INTEGER8, 0,  MPI_COMM_WORLD, ierr)


    fpval(nprocs+1) = lmax+1
!    write (fd,*) 'Pivots: ',fpval(1:nprocs+1)

    !     Determine segment boundaries. Within each bin, fposts(i) is the
    !     start of the ith shuffle segment.
    fposts(1) = 1
    k = 2
    do i=1,np
       !        The first element may be greater than several fencepost values,
       !        so we must use a do-while loop.
       do while (w1(i) .ge. fpval(k))
          fposts(k) = i
          k = k + 1
       enddo
    enddo

    !     The last element may not be greater than the last fencepost value, so
    !     we must assign an appropriate value to every fencepost past the last.
    do i=k,nprocs+1
       fposts(i) = np+1
    enddo

    !     Every process needs fposts() values from every other process, so we will
    !     give each process a copy of all the fposts()s in work().

    do i=1,nprocs
       islen(i) = fposts(i+1) - fposts(i)
    enddo

    call MPI_ALLTOALL( islen,1,MPI_INTEGER, &
                       irlen,1,MPI_INTEGER, &
                       MPI_COMM_WORLD,ierr)


    !     Make sure that "fposts" and "gposts" are zero based for MPI_ALLTOALLV.
    !     fposts and gposts are the addresses of the segment boundaries.
    fposts(1) = 0
    gposts(1) = 0
    do i=1,nprocs
       fposts(i+1) = fposts(i) + islen(i)
       gposts(i+1) = gposts(i) + irlen(i)
    enddo

    npnew = gposts(nprocs+1)

    call MPI_ALLTOALLV(  w1  ,islen,fposts,MPI_INTEGER8, &
                         keys,irlen,gposts,MPI_INTEGER8, &
                         MPI_COMM_WORLD,ierr)


    !     Set up the information for the merge:
    do i=1,nprocs+1
       itabl(i) = gposts(i)
    enddo

    !     Merge the segments within each bin.
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl,itabr,iproc)

  end subroutine psrssort


!  ===================================================================================================



  subroutine pbalsort(nppm,np,npnew,nprocs,iproc,keys, &
       indxl,irnkl, islen,irlen,fposts,gposts,pivot,kw1,wload,nkeys_total,balance,debug,work_local)


    ! P. Gibbon
    !
    ! Adapted from: Xiaobo Li, et al.,
    ! "On the versatility of parallel sorting by regular sampling",
    ! Parallel Computing, volume 19, 1079--1103, Oct 1993.
    ! http://www.cs.utoronto.ca/~paullu/Papers/psrs.ps.Z
    !
    !  Selects pivots according to key distribution across whole range
    !  Keys subdivided according to oct-tree structure, avoiding jumps when crossing
    !  octant boundaries.
    !
    !  Use integer arithmetic for particle weights & local, global key distribution
    !

    implicit  none
    include 'mpif.h'

    integer, intent(in) :: nppm,np,nprocs,iproc,nkeys_total
    real*8, intent(in) :: wload(nppm)  ! particle work loads
    integer, intent(out) :: npnew
    real, intent(in) :: work_local  ! total local work load
    integer*8, dimension(nppm) ::  keys, &      ! array of keys to be sorted.
                                   kw1,kw2,kred      ! work arrays
    integer :: pbin(nppm)
    integer, dimension(nppm) ::  indxl, irnkl ! origin locations of the keys 
    logical, intent(in) :: debug
    integer, intent(in) :: balance
    integer*8 :: pivot(nprocs+1)
    integer*8 :: work1(nppm), work2(nppm)  ! integer arrays for particle loads
    integer, dimension(nprocs) :: islen, irlen  ! send/receive segment lengths
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle
    integer :: itabr(nprocs), itabl(nprocs+1)  ! send/receive strides for A2A sort
    integer*8 :: work_loads(nprocs) ! Aggregate work loads
    integer :: total_keys(nprocs)  ! Total # keys in local tree(s) 
    real*8 :: ave_nkeys,finc  ! Average # keys

    integer, parameter :: maxbin=1000000  ! Max # bins for key distrib
    integer*8, dimension(maxbin)  :: f_local, f_global, f_final  ! Key distribution functions
    integer*8, dimension(maxbin) ::  search_list, retain_list, bin_list
    integer, dimension(maxbin) :: index_bin
    logical :: finished(maxbin)

    integer*8 :: lmax, lmin, key_min, key_max, gkey_min, gkey_max, step ! Key mins and maxes and step size
    integer*8 :: step_reduced, key_reduce, pshift
    integer*8, parameter :: iplace=8_8**20
    integer ::  ibin, itag
    integer :: status(MPI_STATUS_SIZE),ierr
    integer :: alpha  ! alpha now in range 1-100
    real :: ave_work
    integer*8 :: f_integral, work_average, work_threshold
    integer*8 :: total_work, checksum
    integer ::  i,j,k,ipe, iset,fd, nfill, proc_debug, nbin
    integer :: level_dist=8, nlev=20, level_strip
    integer, parameter :: lev_map=15 ! max level for key distrib
    integer*8 :: maxbox

    integer :: ilev, p, np_left, new_bins, nfbins
    character(13) :: cfmt
    integer, save :: icall
    data icall/0/

    double PRECISION :: ts, ttotal, twhile, twhile_allreduce, talltoall

    ttotal = 0.0
    twhile = 0.0
    twhile_allreduce = 0.0
    talltoall = 0.0

    ttotal = MPI_Wtime()

!    fd = iproc+10
    fd=6
    itag=0
    proc_debug = 0
!  Make key map for binning - need to put in setup routine 
search_list = 8_8**lev_map  ! place holder

! Sort local keys
    !     Independent s  !     Note that indx() is a local index on the process.
    call indexsort( keys,indxl, np, nppm )   ! Index sort from Num. Rec.

    do i=1,np
       kw1(i) = keys(indxl(i))
       work1(i) = wload(indxl(i))  ! apply sort to work loads too & convert to integer
       kw2(i) = kw1(i)   ! duplicate work arrays
       work2(i) = work1(i)
    enddo

    !     kw1 now contains the sorted keys; work1 the sorted loads
    !     kw2 the reduced keys
    !     indxl is now the local indexes for the sort.

    lmax = kw1(np)   ! local max
    lmin = kw1(1)    ! local min



    ! Determine global min,max

    call MPI_ALLREDUCE(lmax, gkey_max, 1, MPI_INTEGER8, MPI_MAX,  MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE(lmin, gkey_min, 1, MPI_INTEGER8, MPI_MIN,  MPI_COMM_WORLD, ierr )
! Total # keys
    call MPI_GATHER(nkeys_total, 1, MPI_INTEGER, total_keys, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
    ave_nkeys = SUM(total_keys)/1./nprocs

   if (debug .and. iproc==0) then
	write (*,'(a,i12/(2i10))') "Ave # keys:", int(ave_nkeys),(i,total_keys(i),i=1,nprocs)
   endif

!  Get actual work loads from previous step
!    call MPI_ALLGATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
!    ave_work=SUM(work_loads(1:nprocs))


! Set min/max limits
    key_min = gkey_min
    key_max = gkey_max
    
! Shortcut parallel sort if only 1 CPU by raising particles/bin threshold
    if (nprocs==1) then
       alpha=1
    else
       alpha=25
    endif

    if (debug.and.iproc==proc_debug) then
	write (*,'(2(a12,o30,a12,o30/))') &
         'local min: ',lmin,' local max: ',lmax, &
         'global min: ',gkey_min,' global max: ',gkey_max

    endif

!  place 8 level 1 boxes on search list
   do k = 1,8
      search_list(k) = search_list(k) + 8_8**(lev_map-1)*(k-1)
      finished(k)=.false.
   end do

   np_left = np
   nbin = 8
   nfbins = 0  ! Final # bins
   ilev = 1

   twhile = MPI_Wtime()

   do while (ilev <= lev_map)


! Find local key distribution within current list of boxes (to level lev_map)
! Keys reduced to same level as boxes, so just need to count # with same key.

    level_strip = nlev-ilev

    if (debug .and. iproc==proc_debug ) then
!    if (debug ) then
	write(*,'(a4,i4,a40,4i10)') 'PE',iproc,'ilev, level_strip, nbin, np_left: ',ilev, level_strip,nbin,np_left
    endif
     f_local(1:nbin) = 0  ! reuse work array
     ibin = 1    
     p=1
     step = 1  ! bin size - keys have been reduced to same level as boxes 

     do while (p .le. np_left .and. ibin .le. nbin)

  ! Strip off lower order bits of particle key to match box level 
        key_reduce = 8_8**(lev_map-ilev)*ishft( kw2(p), -3_8*level_strip ) 


	if (key_reduce == search_list(ibin)) then

!    if (debug .and. iproc==proc_debug ) then
!       write(fd,'(a30,2o30)') 'key, reduced key ',kw2(p), key_reduce
!    endif
! keys match - update bin count and move on to next particle

	  if (nprocs.gt.1) then
            loadbal: select case(balance) 
   	    case(1)   ! Equal aggregate interaction lists
	      f_local(ibin) = f_local(ibin) + work2(p) 
 
            case(2)  ! Equal int lists corrected by # local keys 
	      finc = 1.*work2(p)*ave_nkeys/nkeys_total 
	      f_local(ibin) = f_local(ibin) + finc

	    case default   ! Equal # particles
	      f_local(ibin) = f_local(ibin) + 1
 	    end select loadbal

  	   else
	    f_local(ibin) = f_local(ibin) + 1
	   endif


        if (debug .and. iproc==proc_debug) then
	   write(fd,'(a30,o30,2i16,1pe12.1,i16,1pe12.1,i16)') 'reduced key, bin, f_local ', &
		key_reduce,ibin,work2(p),ave_nkeys,nkeys_total,finc,f_local(ibin)
        endif
	  pbin(p) = ibin ! remember bin number 
	  p=p+1
	else
! no match - try next bin
	  ibin=ibin+1
	endif
      end do


    ! Global distrib - must make sure all CPUs participate, even if locally finished

    ts = MPI_Wtime();

    call MPI_ALLREDUCE(f_local, f_global, nbin, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr )

    twhile_allreduce = twhile_allreduce + (MPI_Wtime() - ts);
    
    if (ilev==1) then 
       ave_work=1.0*SUM(f_global(1:nbin))/nprocs        ! need more particles than procs here
       work_threshold = ave_work/alpha
    endif

!    if (debug .and. iproc==proc_debug ) then
!       write(fd,*) 'Particles left:',np_left
!       write(fd,*) 'Bins left:',nbin
!       write(fd,'(a30/(i8,o30,f12.1))') 'search list distrib: ',(i,search_list(i),f_global(i),i=1,nbin)
!       write(fd,*) 'Average bin weight',ave_work
!       write(fd,*) 'Threshold bin weight',ave_work*alpha
!    endif


!  Analyse global distribution on search list
    new_bins=0
    ilev=ilev+1  ! New refinement level

    do ibin=1,nbin

	if (f_global(ibin) > work_threshold .and. ilev<lev_map) then
	   ! bin too heavy: put 8 sub-boxes on new search list
	    do k=0,7
	      new_bins=new_bins+1
              retain_list(new_bins) = search_list(ibin) + 8_8**(lev_map-ilev)*k
	    end do

	else if (f_global(ibin) /= 0) then
	   nfbins=nfbins+1   ! Copy bin onto 'final' list here 
	   bin_list(nfbins) = search_list(ibin)
	   f_final(nfbins) = f_global(ibin)
	   finished(ibin)=.true.
        else
!  discard bin if f_global = 0
	endif
    end do

! Make new key list from incompleted bins: pack operation
    k=0
    do p=1,np_left
	if ( .not. finished(pbin(p)) ) then
	   k=k+1
	   kw2(k)=kw2(p)
	   work2(k)=work2(p)
        endif
    end do
    np_left = k

!  Make new search list
    search_list(1:new_bins) = retain_list(1:new_bins)
    finished(1:new_bins)=.false.
    nbin = new_bins
  end do 

  twhile = MPI_Wtime() - twhile

  nbin = nfbins

    call indexsort( bin_list,index_bin, nbin, maxbin )   ! Index sort from Num. Rec.

    do i=1,nbin
       retain_list(i) = bin_list(index_bin(i))
       f_global(i) = f_final(index_bin(i))
    enddo

    !     kw1 now contains the sorted keys; work1 the sorted loads
    checksum = SUM(f_final(1:nbin))
    total_work=ave_work*nprocs

    if (debug .and. iproc==proc_debug ) then
!    if (iproc==proc_debug ) then
       write(fd,*) 'Final # bins/max/level:',nbin,maxbin,lev_map
       write(fd,*) 'Average bin weight',ave_work
       write(fd,*) 'Threshold bin weight',ave_work*alpha
       write(fd,*) 'Total work',total_work
       write(fd,*) 'Final checksum:',checksum
    endif

       if (balance==0 .and. checksum/=total_work) then
	 if (iproc==proc_debug) then
	   write(fd,*) 'Problem with binning - writing out list to bin.out'
	   open(90,file='bin.out')
	   write(90,'(a30/(i8,o30,f12.1))') 'final bin list ',(i,retain_list(i),f_global(i),i=1,nbin)
	   close(90)
         endif
!	 call closefiles
	 call MPI_FINALIZE(ierr)
         stop
       endif

! Do cumulative integral of f(ibin) and set pivots where int(f) = iproc/nprocs*N

    f_integral = 0
!    ave_work = SUM(f_global(1:nbin))/nprocs  ! Average # particles/PE

    ipe=1
    pivot(1) = ishft(gkey_min,-3_8*(nlev-lev_map))  ! absolute lowest pivot

    do ibin = 1,nbin
       f_integral = f_integral + f_global(ibin)
       if (f_integral >= ipe*ave_work ) then  ! If int(f) exceeds multiple of work average, set pivot
          ipe=ipe+1
          pivot(ipe) = retain_list(ibin)  ! Set next pivot
	endif
   end do

    pivot(nprocs+1) = ishft(gkey_max,-3_8*(nlev-lev_map))+1  ! absolute highest pivot


 !   if (debug .and. iproc==proc_debug ) then
 !      write (fd,'(a20/(10x,i5,o20))') 'Pivots: ',(i,pivot(i),i=1,nprocs+1)
 !   endif

!  if (debug .and. icall==0 .and. iproc==proc_debug) then
  if (debug .and. icall==0 ) then
!	write(*,*) 'PE ',iproc,'Writing key dist'
        open(90,file='fglobal.data')
       write(90,'(a30/(i6,2f12.3))') '! Local & global key distributions: ',(i,f_local(i),f_global(i),i=1,nbin)
	close(90)
  endif


    !     Determine segment boundaries. Within each bin, fposts(i) is the
    !     start of the ith shuffle segment.
    ! need reduced keys here

    fposts(1) = 1
    k = 2

    do i=1,np
       key_reduce = ishft(kw1(i),-3_8*(nlev-lev_map))  ! reduced keys
       !        The first element may be greater than several fencepost values,
       !        so we must use a do-while loop.
       do while (key_reduce .ge. pivot(k) .and. key_reduce .le. pivot(nprocs+1))
          fposts(k) = i
          k = k + 1
	  if (k>nprocs+1) then
	    write(*,*) 'post not found'
	    write(*,'(a5,o20)') 'key = ',i
	    write(*,'(a5,o20)') 'fp = ',pivot(k)
	  endif
       enddo
    enddo

    !     The last element may not be greater than the last fencepost value, so
    !     we must assign an appropriate value to every fencepost past the last.

    do i=k,nprocs+1
       fposts(i) = np+1
    enddo

    !     Every process needs fposts() values from every other process, so we will
    !     give each process a copy of all the fposts()s

    do i=1,nprocs
       islen(i) = fposts(i+1) - fposts(i)
    enddo


    talltoall = MPI_Wtime()

    call MPI_ALLTOALL( islen,1,MPI_INTEGER, &
                       irlen,1,MPI_INTEGER, &
                       MPI_COMM_WORLD,ierr)


    !     Make sure that "fposts" and "gposts" are zero based for MPI_ALLTOALLV.
    !     fposts and gposts are the addresses of the segment boundaries.
    fposts(1) = 0
    gposts(1) = 0
    do i=1,nprocs
       fposts(i+1) = fposts(i) + islen(i)
       gposts(i+1) = gposts(i) + irlen(i)
    enddo

    npnew = gposts(nprocs+1)
    if (npnew.gt.nppm) then
      write(*,*) 'Problem with particle balance on proc',iproc,' np,npnew,npmm:',np,npnew,nppm
    endif
!  Use full keys for swap

    call MPI_ALLTOALLV(  kw1  ,islen,fposts,MPI_INTEGER8, &
                         keys,irlen,gposts,MPI_INTEGER8, &
                         MPI_COMM_WORLD,ierr)
!    if (debug .and. iproc==proc_debug) then
 !      write (*,'(a10,i7/a10/(4i12))') 'npnew: ',npnew, &
 !      'posts: ',(fposts(i),gposts(i),islen(i),irlen(i),i=1,nprocs)
 !   endif

    talltoall = MPI_Wtime() - talltoall


    !     Set up the information for the merge:
    do i=1,nprocs+1
       itabl(i) = gposts(i)
    enddo


    !     Merge the segments within each bin.
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl(1:nprocs+1),itabr(1:nprocs),iproc)

!    if (debug .and. iproc==proc_debug ) then
!       write (fd,'(a20/(10x,5i8))') 'fp, is, gp, ir ',(i,fposts(i),islen(i),gposts(i),irlen(i),i=1,nprocs+1)
!    endif
    icall = icall + 1          ! update call count
    
    ttotal = MPI_Wtime() - ttotal

    if (iproc .eq. -1) then
      write(*,*) 'pbalsort: total:', ttotal
      write(*,*) 'pbalsort: while:', twhile
      write(*,*) 'pbalsort:   allreduce:', twhile_allreduce
      write(*,*) 'pbalsort: alltoall:', talltoall
    endif
    
  end subroutine pbalsort


!  ===================================================================================================



  subroutine pbalsortr(nppm,np,npnew,nprocs,iproc,keys, &
       indxl,irnkl, islen,irlen,fposts,gposts,pivot,kw1,wload,key_box,balance,debug,work_local)


    ! P. Gibbon
    !
    ! Adapted from: Xiaobo Li, et al.,
    ! "On the versatility of parallel sorting by regular sampling",
    ! Parallel Computing, volume 19, 1079--1103, Oct 1993.
    ! http://www.cs.utoronto.ca/~paullu/Papers/psrs.ps.Z
    !
    !  Selects pivots according to key distribution across whole range
    !  Keys subdivided according to oct-tree structure, avoiding jumps when crossing
    !  octant boundaries.
    !
    !  Older version using real distrib. function
    !

    implicit  none
    include 'mpif.h'

    integer, intent(in) :: nppm,np,nprocs,iproc
    real*8, intent(in) :: wload(nppm)  ! particle work loads
    integer, intent(out) :: npnew
    real, intent(in) :: work_local  ! total local work load
    integer*8, dimension(nppm) ::  keys, &      ! array of keys to be sorted.
                                   kw1,kw2,kred      ! work arrays
    integer :: pbin(nppm)
    integer, dimension(nppm) ::  indxl, irnkl ! origin locations of the keys 
    logical, intent(in) :: debug
    integer, intent(in) :: balance
    integer*8 :: pivot(nprocs+1)
    real :: work1(nppm), work2(nppm)
    integer*8, dimension(2) :: key_box
    integer, dimension(nprocs) :: islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle
    integer :: itabr(nprocs), itabl(nprocs+1)
    real :: work_loads(nprocs)
    real :: load_correct(nprocs)
    integer*8 :: lmax, lmin, key_min, key_max, gkey_min, gkey_max, step ! Key mins and maxes and step size
    integer*8 :: step_reduced, key_reduce, pshift
    integer*8, parameter :: iplace=8_8**20
!    integer*8, parameter :: iplace=0
    integer ::  ibin, itag
    integer :: status(MPI_STATUS_SIZE),ierr
    real :: alpha, load_sum, ave_work, f_integral, work_average, load_adjust
    integer :: total_work, checksum
    integer ::  i,j,k,ipe, iset,fd, nfill, proc_debug, nbin
    integer :: level_dist=8, nlev=20, level_strip
    integer, parameter :: lev_map=15 ! max level for key distrib
    integer*8 :: maxbox
    integer, parameter :: maxbin=1500000  ! Max # bins for key distrib
    real, dimension(maxbin)  :: f_local, f_global, f_final
    integer*8, dimension(maxbin) ::  search_list, retain_list, bin_list
    integer, dimension(maxbin) :: index_bin
    logical :: finished(maxbin)
    integer :: ilev, p, np_left, new_bins, nfbins
    character(13) :: cfmt
    integer, save :: icall
    data icall/0/


!    fd = iproc+10
    fd=6
    itag=0
    proc_debug = 0

!  Make key map for binning - need to put in setup routine 
search_list = 8_8**lev_map  ! place holder

! Sort local keys
    !     Independent s  !     Note that indx() is a local index on the process.
    call indexsort( keys,indxl, np, nppm )   ! Index sort from Num. Rec.

    do i=1,np
       kw1(i) = keys(indxl(i))
       work1(i) = wload(indxl(i))  ! apply sort to work loads too
	kw2(i) = kw1(i)   ! duplicate work arrays
	work2(i) = work1(i)
	
    enddo

    !     kw1 now contains the sorted keys; work1 the sorted loads
    !     kw2 the reduced keys
    !     indxl is now the local indexes for the sort.

    lmax = kw1(np)   ! local max
    lmin = kw1(1)    ! local min



    ! Determine global min,max

    call MPI_ALLREDUCE(lmax, gkey_max, 1, MPI_INTEGER8, MPI_MAX,  MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE(lmin, gkey_min, 1, MPI_INTEGER8, MPI_MIN,  MPI_COMM_WORLD, ierr )

!  Get actual work loads from previous step
!    call MPI_ALLGATHER(work_local, 1, MPI_REAL, work_loads, 1, MPI_REAL,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
!    ave_work=SUM(work_loads(1:nprocs))


! Set min/max limits
    key_min = gkey_min
    key_max = gkey_max
    
! Shortcut parallel sort if only 1 CPU by raising particles/bin threshold
    if (nprocs==1) then
       alpha=1
    else
       alpha=0.05
    endif

    if (debug.and.iproc==proc_debug) then
	write (*,'(2(a12,o30,a12,o30/))') &
         'local min: ',lmin,' local max: ',lmax, &
         'global min: ',gkey_min,' global max: ',gkey_max

    endif

!  place 8 level 1 boxes on search list
   do k = 1,8
      search_list(k) = search_list(k) + 8_8**(lev_map-1)*(k-1)
      finished(k)=.false.
   end do

   np_left = np
   nbin = 8
   nfbins = 0  ! Final # bins
   ilev = 1


   do while (ilev <= lev_map)


! Find local key distribution within current list of boxes (to level lev_map)
! Keys reduced to same level as boxes, so just need to count # with same key.

    level_strip = nlev-ilev

    if (debug .and. iproc==proc_debug ) then
!    if (debug ) then
	write(*,'(a4,i4,a40,4i10)') 'PE',iproc,'ilev, level_strip, nbin, np_left: ',ilev, level_strip,nbin,np_left
    endif
     f_local(1:nbin) = 0.  ! reuse work array
     ibin = 1    
     p=1
     step = 1  ! bin size - keys have been reduced to same level as boxes 

     do while (p .le. np_left .and. ibin .le. nbin)

  ! Strip off lower order bits of particle key to match box level 
        key_reduce = 8_8**(lev_map-ilev)*ishft( kw2(p), -3_8*level_strip ) 


	if (key_reduce == search_list(ibin)) then

!    if (debug .and. iproc==proc_debug ) then
!       write(fd,'(a30,2o30)') 'key, reduced key ',kw2(p), key_reduce
!    endif
! keys match - update bin count and move on to next particle
 	if (balance==1 .and. nprocs.gt.1) then
	  f_local(ibin) = f_local(ibin) + work2(p) 
	else
	  f_local(ibin) = f_local(ibin) + 1
	endif
	  pbin(p) = ibin ! remember bin number 
	  p=p+1
	else
! no match - try next bin
	  ibin=ibin+1
	endif
      end do


    ! Global distrib - must make sure all CPUs participate, even if locally finished

    call MPI_ALLREDUCE(f_local, f_global, nbin, MPI_REAL, MPI_SUM,  MPI_COMM_WORLD, ierr )

    if (ilev==1) ave_work=SUM(f_global(1:nbin))/nprocs

!    if (debug .and. iproc==proc_debug ) then
!       write(fd,*) 'Particles left:',np_left
!       write(fd,*) 'Bins left:',nbin
!       write(fd,'(a30/(i8,o30,f12.1))') 'search list distrib: ',(i,search_list(i),f_global(i),i=1,nbin)
!       write(fd,*) 'Average bin weight',ave_work
!       write(fd,*) 'Threshold bin weight',ave_work*alpha
!    endif


!  Analyse global distribution on search list
    new_bins=0
    ilev=ilev+1  ! New refinement level

    do ibin=1,nbin

	if (f_global(ibin) > ave_work*alpha .and. ilev<lev_map) then
	   ! bin too heavy: put 8 sub-boxes on new search list
	    do k=0,7
	      new_bins=new_bins+1
              retain_list(new_bins) = search_list(ibin) + 8_8**(lev_map-ilev)*k
	    end do

	else if (f_global(ibin) /= 0) then
	   nfbins=nfbins+1   ! Copy bin onto 'final' list here 
	   bin_list(nfbins) = search_list(ibin)
	   f_final(nfbins) = f_global(ibin)
	   finished(ibin)=.true.
        else
!  discard bin if f_global = 0
	endif
    end do

! Make new key list from incompleted bins: pack operation
    k=0
    do p=1,np_left
	if ( .not. finished(pbin(p)) ) then
	   k=k+1
	   kw2(k)=kw2(p)
	   work2(k)=work2(p)
        endif
    end do
    np_left = k

!  Make new search list
    search_list(1:new_bins) = retain_list(1:new_bins)
   finished(1:new_bins)=.false.
    nbin = new_bins
  end do 

  nbin = nfbins

    call indexsort( bin_list,index_bin, nbin, maxbin )   ! Index sort from Num. Rec.

    do i=1,nbin
       retain_list(i) = bin_list(index_bin(i))
       f_global(i) = f_final(index_bin(i))
    enddo

    !     kw1 now contains the sorted keys; work1 the sorted loads
    checksum = SUM(f_final(1:nbin))
    total_work=ave_work*nprocs

    if (debug .and. iproc==proc_debug ) then
!    if (iproc==proc_debug ) then
       write(fd,*) 'Final # bins/max/level:',nbin,maxbin,lev_map
       write(fd,*) 'Average bin weight',ave_work
       write(fd,*) 'Threshold bin weight',ave_work*alpha
       write(fd,*) 'Total work',total_work
       write(fd,*) 'Final checksum:',checksum
    endif

       if (balance==0 .and. checksum/=total_work) then
	 if (iproc==proc_debug) then
	   write(fd,*) 'Problem with binning - writing out list to bin.out'
	   open(90,file='bin.out')
	   write(90,'(a30/(i8,o30,f12.1))') 'final bin list ',(i,retain_list(i),f_global(i),i=1,nbin)
	   close(90)
         endif
!	 call closefiles
	 call MPI_FINALIZE(ierr)
         stop
       endif

! Do cumulative integral of f(ibin) and set pivots where int(f) = iproc/nprocs*N

    f_integral = 0.
!    ave_work = SUM(f_global(1:nbin))/nprocs  ! Average # particles/PE

    ipe=1
    pivot(1) = ishft(gkey_min,-3_8*(nlev-lev_map))  ! absolute lowest pivot

    do ibin = 1,nbin
       f_integral = f_integral + f_global(ibin)
       if (f_integral >= ipe*ave_work ) then  ! If int(f) exceeds multiple of work average, set pivot
          ipe=ipe+1
          pivot(ipe) = retain_list(ibin)  ! Set next pivot
	endif
   end do

    pivot(nprocs+1) = ishft(gkey_max,-3_8*(nlev-lev_map))+1  ! absolute highest pivot


 !   if (debug .and. iproc==proc_debug ) then
 !      write (fd,'(a20/(10x,i5,o20))') 'Pivots: ',(i,pivot(i),i=1,nprocs+1)
 !   endif

!  if (debug .and. icall==0 .and. iproc==proc_debug) then
  if (debug .and. icall==0 ) then
!	write(*,*) 'PE ',iproc,'Writing key dist'
        open(90,file='fglobal.data')
       write(90,'(a30/(i6,2f12.3))') '! Local & global key distributions: ',(i,f_local(i),f_global(i),i=1,nbin)
	close(90)
  endif


    !     Determine segment boundaries. Within each bin, fposts(i) is the
    !     start of the ith shuffle segment.
    ! need reduced keys here

    fposts(1) = 1
    k = 2

    do i=1,np
       key_reduce = ishft(kw1(i),-3_8*(nlev-lev_map))  ! reduced keys
       !        The first element may be greater than several fencepost values,
       !        so we must use a do-while loop.
       do while (key_reduce .ge. pivot(k) .and. key_reduce .le. pivot(nprocs+1))
          fposts(k) = i
          k = k + 1
	  if (k>nprocs+1) then
	    write(*,*) 'post not found'
	    write(*,'(a5,o20)') 'key = ',i
	    write(*,'(a5,o20)') 'fp = ',pivot(k)
	  endif
       enddo
    enddo

    !     The last element may not be greater than the last fencepost value, so
    !     we must assign an appropriate value to every fencepost past the last.

    do i=k,nprocs+1
       fposts(i) = np+1
    enddo

    !     Every process needs fposts() values from every other process, so we will
    !     give each process a copy of all the fposts()s

    do i=1,nprocs
       islen(i) = fposts(i+1) - fposts(i)
    enddo



    call MPI_ALLTOALL( islen,1,MPI_INTEGER, &
                       irlen,1,MPI_INTEGER, &
                       MPI_COMM_WORLD,ierr)


    !     Make sure that "fposts" and "gposts" are zero based for MPI_ALLTOALLV.
    !     fposts and gposts are the addresses of the segment boundaries.
    fposts(1) = 0
    gposts(1) = 0
    do i=1,nprocs
       fposts(i+1) = fposts(i) + islen(i)
       gposts(i+1) = gposts(i) + irlen(i)
    enddo

    npnew = gposts(nprocs+1)


!  Use full keys for swap

    call MPI_ALLTOALLV(  kw1  ,islen,fposts,MPI_INTEGER8, &
                         keys,irlen,gposts,MPI_INTEGER8, &
                         MPI_COMM_WORLD,ierr)
!    if (debug .and. iproc==proc_debug) then
 !      write (*,'(a10,i7/a10/(4i12))') 'npnew: ',npnew, &
 !      'posts: ',(fposts(i),gposts(i),islen(i),irlen(i),i=1,nprocs)
 !   endif


    !     Set up the information for the merge:
    do i=1,nprocs+1
       itabl(i) = gposts(i)
    enddo


    !     Merge the segments within each bin.
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl(1:nprocs+1),itabr(1:nprocs),iproc)

!    if (debug .and. iproc==proc_debug ) then
!       write (fd,'(a20/(10x,5i8))') 'fp, is, gp, ir ',(i,fposts(i),islen(i),gposts(i),irlen(i),i=1,nprocs+1)
!    endif
    icall = icall + 1          ! update call count
  end subroutine pbalsortr



! ========================================================================


  subroutine pswssort(nppm,np,npnew,nprocs,iproc,keys, &
       indxl,irnkl, islen,irlen,fposts,gposts,kw1,wload,key_box,balance,debug)


    ! P. Gibbon
    !
    ! Adapted from: Xiaobo Li, et al.,
    ! "On the versatility of parallel sorting by regular sampling",
    ! Parallel Computing, volume 19, 1079--1103, Oct 1993.
    ! http://www.cs.utoronto.ca/~paullu/Papers/psrs.ps.Z
    !
    !  Selects pivots according to key distribution across whole range
    !
    ! sorting articles:
    ! http://www.lpac.ac.uk/SEL-HPC/Articles/GeneratedHtml/hpc.sort.html
    !
    !

    implicit  none
    include 'mpif.h'

    integer, intent(in) :: nppm,np,nprocs,iproc
    real*8, intent(in) :: wload(nppm)  ! particle work loads
    integer, intent(out) :: npnew
    integer, parameter :: binmult=1000000   !TODO: need to reduce size of f() arrays
!    integer, parameter :: binmult=100000   !TODO: need to reduce size of f() arrays
    integer*8, dimension(nppm) ::  keys, &      ! array of keys to be sorted.
                                   kw1       ! work array
    integer, dimension(nppm) ::  indxl, irnkl ! origin locations of the keys 
    logical :: debug
    integer, intent(in) :: balance
    real :: w2(nppm)
    integer*8, dimension(2) :: key_box
    integer, dimension(nprocs) :: islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle
    integer :: itabr(nprocs), itabl(nprocs+1)
    real, save, dimension(binmult)  :: f_local, f_global
    integer*8 :: fpval(nprocs+1)
    integer*8 :: lmax, lmin, key_min, key_max, gkey_min, gkey_max, step ! Key mins and maxes and step size
    integer*8 :: step_reduced
    integer :: nbin, ibin, itag
    integer :: status(MPI_STATUS_SIZE),ierr
    real :: ave_work, f_integral
    integer ::  i,j,k, fd, nfill, proc_debug
    character(13) :: cfmt

    nbin = binmult  ! must correspond to array size

!    fd = iproc+10
    fd=6
    itag=0
    proc_debug = 0

    !     Independent s  !     Note that indx() is a local index on the process.
    call indexsort( keys,indxl, np, nppm )   ! Index sort from Num. Rec.

    do i=1,np
       kw1(i) = keys(indxl(i))
       w2(i) = wload(indxl(i))  ! apply sort to work loads too
    enddo

    !     kw1 now contains the sorted keys; w2 the sorted loads
    !     indxl is now the local indexes for the sort.

    lmax = kw1(np)   ! local max
    lmin = kw1(1)    ! local min

    ! Determine global min,max

    call MPI_ALLREDUCE(lmax, gkey_max, 1, MPI_INTEGER8, MPI_MAX,  MPI_COMM_WORLD, ierr )

    call MPI_ALLREDUCE(lmin, gkey_min, 1, MPI_INTEGER8, MPI_MIN,  MPI_COMM_WORLD, ierr )


    step=(gkey_max - gkey_min)/nbin + 1
    step_reduced = (key_box(2) - key_box(1))/nbin + 1

! Set min/max limits
!    key_min = key_box(1)
!    key_max = key_box(2)
    key_min = gkey_min
    key_max = gkey_max
    
    if (debug.and.iproc==proc_debug) write (*,'(3(a12,z20,a12,z20/),a12,i8,a12,z20)') &
         'local min: ',lmin,' local max: ',lmax, &
         'global min: ',gkey_min,' global max: ',gkey_max, &
         ' box_min: ',key_min,' box_max: ',key_max, &
         ' nbin ',nbin,' step', step

    f_local(1:nbin) = 0.0

    ! Find local key distribution
    do k=1,np
          ! bin inside container limits
          ibin = (kw1(k)-key_min)/step + 1
          ibin = max(min(ibin,nbin),1)
       if (balance==1) then
          f_local(ibin) = f_local(ibin) + w2(k)  ! Can include actual force load on particle here
       else
          f_local(ibin) = f_local(ibin) + 1  ! No load balancing - try to get equal # particles  
       endif
    enddo

    if (debug .and. iproc==proc_debug) then
       cfmt = "(/a15,"//achar(mod(nbin/100,10)+48)//achar(mod(nbin/10,10)+48) // achar(mod(nbin,10)+48) // "(f12.4))"
       write(fd,'(a15/(f12.3))') 'Local key distrib: ',(f_local(i),i=1,nbin,nbin/10)
    endif

    ! Global distrib
    call MPI_ALLREDUCE(f_local, f_global, nbin, MPI_REAL, MPI_SUM,  MPI_COMM_WORLD, ierr )

    if (debug .and. iproc==proc_debug) then
       write(fd,'(a15/(f12.3))') 'Global key distrib: ',(f_global(i),i=1,nbin,nbin/10)
       write(fd,*) 'Checksum, N=SUM(f)=',SUM(f_global(1:nbin))
    endif

    ! Do cumulative integral of f(ibin) and set pivots where int(f) = iproc/nprocs*N

    f_integral = 0.
    ave_work = SUM(f_global(1:nbin))/nprocs  ! Average # particles/PE
    fpval(1) = gkey_min  ! absolute lowest pivot
    i=1

    do ibin = 1,nbin
       f_integral = f_integral + f_global(ibin)
       if (f_integral >= i*ave_work ) then  ! If int(f) exceeds multiple of work average, set pivot
	  nfill = (f_integral-i*ave_work)/ave_work
	  if (nfill>0) then 
	        write (*,'(a20,i8/a40,i8/a10,f12.3,a10,f12.3)') 'Problem on Processor ',iproc &
		, 'integral jumps by more than work average: nfill= ',nfill &
                , 'increment = ', f_global(ibin), ' ave= ',ave_work
!		open(90,file='fglobal.data')
!		write (90,'((i8,f12.2))') (j,f_global(j),j=1,nbin)
!	       close(90)
!	       call closefiles
	       call MPI_ABORT(MPI_COMM_WORLD,ierr)
         	stop
          endif
          i=i+1
          fpval(i) = key_min + ibin*step  ! Set next highest pivot
       endif
    end do

    fpval(nprocs+1) = gkey_max+1  ! Set absolute highest pivot

    if (debug .and. iproc==proc_debug) then
       write (fd,*) 'Ave work: ',ave_work
       write (fd,*) 'f_Integral: ',f_integral
       write (fd,'(a10,z20/)') 'key_min:',gkey_min
       write (fd,'(a10/(10x,i5,z20))') 'Pivots: ',(i,fpval(i),i=1,nprocs+1)
       write (fd,'(/a10,z20)') 'key_max:',gkey_max
    endif

    !     Determine segment boundaries. Within each bin, fposts(i) is the
    !     start of the ith shuffle segment.
    fposts(1) = 1
    k = 2
    do i=1,np
       !        The first element may be greater than several fencepost values,
       !        so we must use a do-while loop.
       do while (kw1(i) .ge. fpval(k) .and. kw1(i).le.key_max)
          fposts(k) = i
          k = k + 1
	  if (k>nprocs+1) then
	    write(*,*) 'post not found'
	    write(*,'(a5,o20)') 'key = ',i
	    write(*,'(a5,o20)') 'fp = ',fpval(k)
	  endif
       enddo
    enddo

    !     The last element may not be greater than the last fencepost value, so
    !     we must assign an appropriate value to every fencepost past the last.
    do i=k,nprocs+1
       fposts(i) = np+1
    enddo

    !     Every process needs fposts() values from every other process, so we will
    !     give each process a copy of all the fposts()s

    do i=1,nprocs
       islen(i) = fposts(i+1) - fposts(i)
    enddo

    call MPI_ALLTOALL( islen,1,MPI_INTEGER, &
                       irlen,1,MPI_INTEGER, &
                       MPI_COMM_WORLD,ierr)


    !     Make sure that "fposts" and "gposts" are zero based for MPI_ALLTOALLV.
    !     fposts and gposts are the addresses of the segment boundaries.
    fposts(1) = 0
    gposts(1) = 0
    do i=1,nprocs
       fposts(i+1) = fposts(i) + islen(i)
       gposts(i+1) = gposts(i) + irlen(i)
    enddo

    npnew = gposts(nprocs+1)


    call MPI_ALLTOALLV(  kw1  ,islen,fposts,MPI_INTEGER8, &
                         keys,irlen,gposts,MPI_INTEGER8, &
                         MPI_COMM_WORLD,ierr)
    if (debug .and. iproc==57) then

       write (*,'(a10,i7/a10/(4i12))') 'npnew: ',npnew, &
       'posts: ',(fposts(i),gposts(i),islen(i),irlen(i),i=1,nprocs)
    endif
    !     Set up the information for the merge:
    do i=1,nprocs+1
       itabl(i) = gposts(i)
    enddo

    !     Merge the segments within each bin.
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl(1:nprocs+1),itabr(1:nprocs),iproc)

  end subroutine pswssort


  subroutine nwaymrg(nppm,np,nprocs,keys,irnkl,itabl,itabr,me)

    implicit none
    integer :: nppm,np,nprocs,me
    integer :: irnkl(nppm)
    integer*8 :: keys(nppm)
    integer :: itabl(nprocs+1),itabr(nprocs)

    integer, parameter :: maxprocs=1024
    integer*8 :: ik,itabk(nprocs+1)
    integer :: i,j,k,il,ir,ij,irnkj(nprocs+1),jprocs

    !     Make sure that the indicies are "one" baised.
    if (itabl(1) .ne. 1) then
       do i=nprocs+1,1,-1
          itabl(i) = itabl(i) - itabl(1) + 1
       enddo
    endif

    do i=1,nprocs
       itabr(i) = itabl(i+1) - 1
    enddo

!    if (me.eq.57) then
!       write(*,'(2i10)') (i,itabl(i),i=1,nprocs+1)
!    endif

    do j=1,nprocs
       itabk(j) = keys(itabl(j))
       irnkj(j) = j
    enddo

    !     Check for empty segments and remove them.
    jprocs = nprocs
    i = 1
    do j=1,nprocs
       if (itabl(irnkj(i)) .gt. itabr(irnkj(i))) then
          jprocs = jprocs - 1
          do k=i,jprocs
             irnkj(k) = irnkj(k+1)
          enddo
       else
          i = i + 1
       endif
    enddo

    !     Sort the N-way merging table:
    do j=2,jprocs
       !        Consider each of the original elements in turn.
       ij = irnkj(j)
       ik = itabk(ij)
       !           and look for a place to insert it
       !           The slot "j" is now empty.
       do i=j-1,1,-1
          if(itabk(irnkj(i)) .le. ik) goto 202
          irnkj(i+1) = irnkj(i)
       enddo
       i=0
202    continue
       irnkj(i+1) = ij
    enddo
    !     The merging table is now in sorted order.


    !     proceed with the merge
    do i=1,np

       !        Remove the smallest element from the merging list.
       ij = irnkj(1)

       !        refresh the merge table
       il = itabl(ij) + 1
       ir = itabr(ij)
       ik = keys(il)
       itabk(ij) = ik
       itabl(ij) = il
       irnkl(il-1) = i

       !        pick out each element in turn

       !        The first slot is now empty.
       if( ir .ge. il ) then
          !           look for the slot to insert the new data.
          do j=1,jprocs-1
             if(itabk(irnkj(j+1)) .ge. ik) goto 203
             irnkj(j) = irnkj(j+1)
          enddo
          j=jprocs
203       continue
          irnkj(j) = ij
       else
          !           retire a slot
          jprocs = jprocs-1
          do j=1,jprocs
             irnkj(j) = irnkj(j+1)
          enddo
       endif
    enddo

  end subroutine nwaymrg


  !  ================================
  !
  !         SORT_HEAP
  !
  !     Sort 64-bit integer array into ascending order
  !     using heap-sort algorithm
  !    (Numerical Recipes f90, p1171)
  !
  !  ================================

  subroutine sort_i(iarr)
    implicit none
    integer*8, intent(inout) :: iarr(:)
    integer :: i,n

    n = size(iarr)


    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call swap( iarr(1),iarr(i) )      ! Clear space at end of array and retire top of heap into it
       call sift_down( 1,i-1)
    end do

  contains
    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain 
      integer, intent(in) :: l,r     ! the heap structure
      integer :: j,jold    ! index
      integer*8 :: a

      a = iarr(l)

      jold = l
      j = l + l
      do                   ! do while j <= r
         if (j > r) exit
         if (j < r) then
            if (iarr(j) < iarr(j+1)) j = j+1
         endif
         if (a >= iarr(j)) exit       ! Found a's level, so terminate sift-down
         iarr(jold) = iarr(j)
         jold = j                    ! Demote a and continue
         j = j+j
      end do
      iarr(jold) = a                  ! Put a into its slot

    end subroutine sift_down

  end subroutine sort_i


!  Index sort for 8-byte integer list

  subroutine indsort_i8(iarr,list,n,nppm)
    implicit none

    integer, intent(in) :: n,nppm

    integer*8, dimension(nppm), intent(in) :: iarr
    integer, dimension(nppm), intent(inout) :: list
    integer :: i, indxt, ir, l, j
    integer*8 :: q

    list(1:n) =  (/ (i,i=1,n) /)  

    if(n==1) return
    l= n/2 + 1
    ir = n

    do
       if (l>1) then
          l=l-1
          indxt = list(l)
          q = iarr(indxt)
       else
          indxt = list(ir)
          q = iarr(indxt)
          list(ir) = list(1)
          ir = ir - 1
          if (ir == 1) then
             list(1) = indxt
             return
          endif
       endif

       i = l
       j = l+l
       do while (j <= ir)
          if (j < ir) then
             if (iarr(list(j)) < iarr(list(j+1)) ) j=j+1
          endif
          if (q < iarr(list(j)) ) then
             list(i) = list(j)
             i=j
             j=j+j
          else
             j = ir+1
          endif
       end do
       list(i) = indxt
    end do

  end subroutine indsort_i8


!  Index sort for 4-byte integer list

  subroutine indsort_i4(iarr,list,n,nppm)
    implicit none

    integer, intent(in) :: n,nppm

    integer, dimension(nppm), intent(in) :: iarr
    integer, dimension(nppm), intent(inout) :: list
    integer :: i, indxt, ir, l, j
    integer :: q

    list(1:n) =  (/ (i,i=1,n) /)  

    if(n==1) return
    l= n/2 + 1
    ir = n

    do
       if (l>1) then
          l=l-1
          indxt = list(l)
          q = iarr(indxt)
       else
          indxt = list(ir)
          q = iarr(indxt)
          list(ir) = list(1)
          ir = ir - 1
          if (ir == 1) then
             list(1) = indxt
             return
          endif
       endif

       i = l
       j = l+l
       do while (j <= ir)
          if (j < ir) then
             if (iarr(list(j)) < iarr(list(j+1)) ) j=j+1
          endif
          if (q < iarr(list(j)) ) then
             list(i) = list(j)
             i=j
             j=j+j
          else
             j = ir+1
          endif
       end do
       list(i) = indxt
    end do

  end subroutine indsort_i4



  subroutine swap_ab(p,q)
    integer*8 :: p,q, dum
    dum = p
    p=q
    q = dum
  end subroutine swap_ab




! Removes duplicate entries from key list
 
  subroutine uniq8(key, ntot)
    integer :: ntot, nu
    integer*8, dimension(ntot+1) :: key
    nu=0
    do i=1,ntot
       if (key(i) /= key(i+1) ) then
          nu=nu+1
          key(nu) = key(i) ! pack
       endif
    end do
    ntot=nu
  end subroutine uniq8


  ! Get address of variable
  subroutine locaddress_int4_8(location, addr, ierr)
    
    implicit none
    include 'mpif.h'
    
    integer,intent(in) :: location
    integer(kind=MPI_ADDRESS_KIND),intent(out) :: addr 
    integer :: ierr
    
    ierr=0
    call MPI_GET_ADDRESS( location, addr, ierr )
    
  end subroutine locaddress_int4_8


  subroutine locaddress_int8_8(location, addr, ierr)
    
    implicit none
    include 'mpif.h'
    
    integer*8,intent(in):: location
    integer(kind=MPI_ADDRESS_KIND),intent(out) :: addr 
    integer :: ierr
    
    ierr=0
    call MPI_GET_ADDRESS( location, addr, ierr )
    
  end subroutine locaddress_int8_8


  subroutine locaddress_real8_8(location, addr, ierr)
    
    implicit none
    include 'mpif.h'
    
    real*8,intent(in) :: location
    integer(kind=MPI_ADDRESS_KIND),intent(out) :: addr 
    integer :: ierr
    
    ierr=0
    call MPI_GET_ADDRESS( location, addr, ierr )
    
  end subroutine locaddress_real8_8

  
  subroutine bpi_int8_8(a, b, base, res)
    
    implicit none    
    
    integer*8,intent(in) :: a, b, base
    integer*8,intent(out) :: res
    integer*8 :: k, ka, kb, temp
    integer*8 :: i 
    integer*8 :: pot
    
    if (a.eq.b)then
       res=a
       return
    end if
    
    ! swap for correct interval 
    if (b .le. a) then
       call swap(a,b)
    end if
    
    pot = floor(log(REAL(b))/log(REAL(base)))
    do i = pot, 0, -1
       k = b/(base**i)
       res = k * base**i
       if (a.lt.res) then
          return
       end if
    end do
    
  end subroutine bpi_int8_8
  
  subroutine bpi_bits_int8_8(a, b, base, res)
    
    implicit none   

    integer*8,intent(in) :: a, b, base
    integer*8,intent(out) :: res
    integer*8 :: l_cell,r_cell
    integer*8 :: k, pot
    integer*8 :: i 
    integer*8 :: nbits
    

    nbits = 3

    ! swap for correct interval 
    if (b .le. a) then
       call swap(a,b)
    end if

    
    do i=1,63/nbits
     
       l_cell = ibits(a,63-i*nbits+1,nbits)
       r_cell = ibits(b,63-i*nbits+1,nbits)
     
       if(l_cell.ne.r_cell)then
          pot=(63/nbits)-i+1
          exit
       end if
       
    end do

    k=b/base**pot
    res=k*base**pot
    
    
  end subroutine bpi_bits_int8_8

end module tree_utils
