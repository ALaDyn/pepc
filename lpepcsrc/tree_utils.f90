!
!                          Tree sort utility module
!                       *   Parallel Sorting
!

module tree_utils

  ! Interfaces
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

  subroutine psrsperm_i4(nppm,np,npnew,nprocs,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)

    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs

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

  subroutine psrsperm_i8(nppm,np,npnew,nprocs,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)


    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs

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

  subroutine psrsperm_r8(nppm,np,npnew,nprocs,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)


    implicit none
    include 'mpif.h'

    integer :: nppm,np,npnew,nprocs

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
                                   kw1,kw2      ! work arrays
    integer :: pbin(nppm)
    integer, dimension(nppm) ::  indxl, irnkl ! origin locations of the keys 
    logical, intent(in) :: debug
    integer, intent(in) :: balance
    integer*8 :: pivot(nprocs+1)
    integer*8 :: work1(nppm), work2(nppm)  ! integer arrays for particle loads
    integer, dimension(nprocs) :: islen, irlen  ! send/receive segment lengths
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle
    integer :: itabr(nprocs), itabl(nprocs+1)  ! send/receive strides for A2A sort
    integer :: total_keys(nprocs)  ! Total # keys in local tree(s) 
    real*8 :: ave_nkeys,finc  ! Average # keys

    integer, parameter :: maxbin=1000000  ! Max # bins for key distrib
    integer*8, dimension(maxbin)  :: f_local, f_global, f_final  ! Key distribution functions
    integer*8, dimension(maxbin) ::  search_list, retain_list, bin_list
    integer, dimension(maxbin) :: index_bin
    logical :: finished(maxbin)

    integer*8 :: lmax, lmin, key_min, key_max, gkey_min, gkey_max, step ! Key mins and maxes and step size
    integer*8 :: key_reduce
    integer*8, parameter :: iplace=8_8**20
    integer ::  ibin, itag
    integer :: ierr
    integer :: alpha  ! alpha now in range 1-100
    real*8 :: ave_work
    integer*8 :: f_integral, work_threshold
    integer*8 :: total_work, checksum
    integer ::  i,k,ipe, fd, proc_debug, nbin
    integer :: nlev=20, level_strip
    integer, parameter :: lev_map=15 ! max level for key distrib

    integer :: ilev, p, np_left, new_bins, nfbins
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
       work1(i) = int(wload(indxl(i)),kind(work1(i)))  ! apply sort to work loads too & convert to integer
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
	      f_local(ibin) = f_local(ibin) + int(finc,kind(f_local))

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
       ave_work=1.0_8*SUM(f_global(1:nbin))/nprocs        ! need more particles than procs here
       work_threshold = int(ave_work/alpha,kind(work_threshold))
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
    total_work=int(ave_work*nprocs,kind(total_work))

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
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl(1:nprocs+1),itabr(1:nprocs))

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

! ========================================================================

  subroutine nwaymrg(nppm,np,nprocs,keys,irnkl,itabl,itabr)

    implicit none
    integer :: nppm,np,nprocs
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
         if (a >= iarr(j)) exit       ! Found a`s level, so terminate sift-down
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
    integer   :: location
    integer*8 :: addr
    integer   :: ierr
! IBM machines powerX, bg/l use internal function
    addr = LOC(location)
    ierr=0
! Linux: use mpi1 address function
!   call MPI_ADDRESS( location, addr4, ierr )
!   addr=addr4
  end subroutine locaddress_int4_8

  subroutine locaddress_int8_8(location, addr, ierr)
    integer*8 :: location
    integer*8 :: addr
    integer   :: ierr
    addr = LOC(location)
    ierr=0
!   call MPI_ADDRESS( location, addr4, ierr )
!   addr=addr4
  end subroutine locaddress_int8_8


  subroutine locaddress_real8_8(location, addr, ierr)
    real*8    :: location
    integer*8 :: addr
    integer   :: ierr
    addr = LOC(location)
    ierr=0
!   call MPI_ADDRESS( location, addr4, ierr )
!   addr=addr4
  end subroutine locaddress_real8_8

  subroutine bpi_int8_8(a, b, base, res)
    
    implicit none    
    include 'mpif.h'

    integer :: ierr
    
    integer*8,intent(in) :: a, b, base
    integer*8,intent(out) :: res
    integer*8 :: k
    integer*8 :: i 
    integer*8 :: pot
    
  
    pot = floor(log(REAL(b))/log(REAL(base)))
    do i = pot, 0, -1
       
       k = b/base**i
       res = k * base**i
       if (a.lt.res) then
          return
       end if
    end do

    res=-1
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)


  end subroutine bpi_int8_8

  subroutine bpi_bits_int8_8(a, b, base, res, levels)
 
    implicit none
    include 'mpif.h'

    integer :: ierr
 
    integer*8,intent(in) :: a, b, base
    integer,intent(in) :: levels
    integer*8,intent(out) :: res
    integer*8 :: i
    integer*8 :: bn, pos
  
    do i=1,levels
       pos=3*(levels-i)
       if(ibits(a,pos,3).ne.ibits(b,pos,3))then
          bn=8**(levels-i)
          res=b/bn*bn ! integer division
          return
       end if
    end do

    res=-1
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

  end subroutine bpi_bits_int8_8

end module tree_utils
