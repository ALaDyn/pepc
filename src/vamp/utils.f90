!
!                          General utility module
!                       *   I/O routines
!                       *   Sorting
!                       *   Random number gen.
!

module utils

  ! Interfaces
  interface pll_regsort
     module procedure psrssort
  end interface

  interface pll_weightsort
     module procedure pswssort
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
     module procedure indsort_i
  end interface

  interface swap
     module procedure swap_ab
  end interface

  interface blank
     module procedure blankn, blank6
  end interface

  interface cputime
     module procedure cput
  end interface

contains

! ================
! permute integer*4
! ================

  subroutine psrsperm_i4(nppm,np,npnew,nprocs,iproc,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)

    use my_mpidefs
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

    integer :: nppm,np,npnew,nprocs,iproc

    integer, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) :: indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts

    integer :: i

    save

!VAMPINST subroutine_start
       CALL VTENTER(IF_psrsperm_i4,VTNOSCL,VTIERR)
!      write(*,*) 'VT: psrsperm_i4 S>',VTIERR,
!     *    IF_psrsperm_i4,ICLASSH
!
    do i=1,np
       w1(i) = array(indxl(i))
    enddo

    call MPI_ALLTOALLV(  w1, islen, fposts, MPI_INTEGER, &
                         w2, irlen, gposts, MPI_INTEGER, &
                         MPI_COMM_WORLD,ierr)

    do i=1,npnew
       array(irnkl(i)) = w2(i)
    enddo

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: psrsperm_i4 S<',VTIERR,ICLASSH
!
  end subroutine psrsperm_i4


! ================
! permute integer*8
! ================

  subroutine psrsperm_i8(nppm,np,npnew,nprocs,iproc,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)

    use my_mpidefs
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

    integer :: nppm,np,npnew,nprocs,iproc

    integer*8, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) :: indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts

    integer :: i

    save

!VAMPINST subroutine_start
       CALL VTENTER(IF_psrsperm_i8,VTNOSCL,VTIERR)
!      write(*,*) 'VT: psrsperm_i8 S>',VTIERR,
!     *    IF_psrsperm_i8,ICLASSH
!
    do i=1,np
       w1(i) = array(indxl(i))
    enddo

    call MPI_ALLTOALLV(  w1, islen, fposts, MPI_INTEGER8, &
                         w2, irlen, gposts, MPI_INTEGER8, &
                         MPI_COMM_WORLD,ierr)

    do i=1,npnew
       array(irnkl(i)) = w2(i)
    enddo

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: psrsperm_i8 S<',VTIERR,ICLASSH
!
  end subroutine psrsperm_i8


! ================
  !  Permute REAL*8 array()
! ================

  subroutine psrsperm_r8(nppm,np,npnew,nprocs,iproc,array,w1,w2, indxl,irnkl,islen,irlen,fposts,gposts)

    use my_mpidefs
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

    integer :: nppm,np,npnew,nprocs,iproc

    real, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) ::  indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts


    integer :: i

    save

!VAMPINST subroutine_start
       CALL VTENTER(IF_psrsperm_r8,VTNOSCL,VTIERR)
!      write(*,*) 'VT: psrsperm_r8 S>',VTIERR,
!     *    IF_psrsperm_r8,ICLASSH
!
    do i=1,np
       w1(i) = array(indxl(i))
    enddo


    call MPI_alltoallv(  w1,islen, fposts, MPI_REAL8, &
         w2, irlen, gposts, MPI_REAL8, &
         MPI_COMM_WORLD,ierr)

    do i=1,npnew
       array(irnkl(i)) = w2(i)
    enddo

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: psrsperm_r8 S<',VTIERR,ICLASSH
!
  end subroutine psrsperm_r8


  subroutine psrssort(nppm,np,npnew,nprocs,iproc,keys,indxl,irnkl,islen,irlen,fposts,gposts,w1)

    use my_mpidefs
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

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
    integer :: lmax, ak

    integer*8 :: buf(maxprocs)

    integer :: i,j,k,step

    integer :: fd, nsamp
    character(13) :: cfmt

!VAMPINST subroutine_start
       CALL VTENTER(IF_psrssort,VTNOSCL,VTIERR)
!      write(*,*) 'VT: psrssort S>',VTIERR,
!     *    IF_psrssort,ICLASSH
!
    fd = 20

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

    nbuf = nprocs

    if (iproc .ne. 0) then
       call MPI_SEND(work, nbuf, MPI_integer8, root, tag1, MPI_COMM_WORLD, ierr)
    else
       do i=1,nprocs-1
          nbuf = nprocs
          call MPI_RECV(buf, nbuf, MPI_integer8, MPI_ANY_SOURCE, tag1, MPI_COMM_WORLD, status, ierr)
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
       write(fd,cfmt) 'Sorted sample: ',(work(i),i=1,nsamp)

       k = 1
       do i=1,nprocs*nprocs,nprocs
          fpval(k) = work(i)
          k = k + 1
       enddo
    endif

    nbuf = nprocs+1
    call MPI_BCAST(fpval, nbuf, MPI_integer8, root,  MPI_COMM_WORLD, ierr)


    fpval(nprocs+1) = lmax+1
    write (fd,*) 'Pivots: ',fpval(1:nprocs+1)

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

    call MPI_ALLTOALL( islen,one,MPI_INTEGER, &
                       irlen,one,MPI_INTEGER, &
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
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl,itabr)

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: psrssort S<',VTIERR,ICLASSH
!
  end subroutine psrssort

  subroutine pswssort(nppm,np,npnew,nprocs,iproc,keys,indxl,irnkl, islen,irlen,fposts,gposts,w1,wload,balance)

   use my_mpidefs
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

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
    !     (C. Mobarry/GSFC, J. Crawford/VSEP)
    !
    !     Convex SPP-1000 Exemplar MPI implementation
    !     Written by C. Mobarry and J. Crawford
    !     NASA/GSFC Code 934, Greenbelt MD, 20771
    !     Clark.Mobarry@gsfc.nasa.gov, {mobarry,crawford}@maxwell.gsfc.nasa.gov

    integer :: nppm,np,npnew,nprocs,iproc
    integer*8, dimension(nppm) ::  keys, &      ! array of keys to be sorted.
                                   w1       ! work array
    integer, dimension(nppm) ::  indxl, irnkl ! origin locations of the keys 
    real :: wload(nppm), w2(nppm)

    integer, dimension(nprocs) :: islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle
    integer, parameter :: maxprocs = 1024
    integer :: itabr(maxprocs), itabl(maxprocs+1)
    real, dimension(8*nppm)  :: f_local, f_global
    integer*8 :: fpval(maxprocs+1)
    integer*8 :: lmax, lmin, key_min, key_max, step  ! Key mins and maxes and step size
    integer :: ak, nbin, ibin, npost
    real :: ave_work, f_integral
    integer :: buf(maxprocs)
    integer ::  i,j,k, fd, nsamp
    character(13) :: cfmt
    logical :: debug=.false., balance

!    fd = iproc+10
!VAMPINST subroutine_start
       CALL VTENTER(IF_pswssort,VTNOSCL,VTIERR)
!      write(*,*) 'VT: pswssort S>',VTIERR,
!     *    IF_pswssort,ICLASSH
!
    fd=20
    itag=0

    !     Independent s  !     Note that indx() is a local index on the process.
    call indexsort( keys,indxl, np, nppm )   ! Index sort from Num. Rec.

    do i=1,np
       w1(i) = keys(indxl(i))
       w2(i) = wload(indxl(i))  ! apply sort to work loads too
    enddo

    !     w1() is now the sorted keys, w2() the sorted loads
    !     indxl() is now the local indexes for the sort.

    lmax = w1(np)   ! local max
    lmin = w1(1)    ! local min

    ! Determine global min,max

    call MPI_ALLREDUCE(lmax, key_max, one, MPI_INTEGER8, MPI_MAX,  MPI_COMM_WORLD, ierr )
    call MPI_ALLREDUCE(lmin, key_min, one, MPI_INTEGER8, MPI_MIN,  MPI_COMM_WORLD, ierr )

    !     Choose bin size for key distribution.
    nbin = 8*nppm  ! # bins - may want to refine this for very large N
    step=(key_max - key_min)/nbin + 1

    if (debug) write (fd,*) 'keymax: ',key_max,' keymin: ',key_min,' nbin ',nbin+1,' step', step

    f_local(1:nbin) = 0.0

    ! Find local key distribution
    do k=1,np
       ibin = (w1(k)-key_min)/step + 1
       if (balance) then
          f_local(ibin) = f_local(ibin) + w2(k)  ! Can include actual force load on particle here
       else
          f_local(ibin) = f_local(ibin) + 1  ! No load balancing - try to get equal # particles  
       endif
    enddo

    if (debug) then
       cfmt = "(/a15,"//achar(mod(nbin/10,10)+48) // achar(mod(nbin,10)+48) // "(f12.4))"
       write(fd,cfmt) 'Local key distrib: ',(f_local(i),i=1,nbin)
    endif

    ! Global distrib
    nbuf = nbin
    call MPI_ALLREDUCE(f_local, f_global, nbuf, MPI_REAL8, MPI_SUM,  MPI_COMM_WORLD, ierr )

    if (debug) then
       write(fd,cfmt) 'Global key distrib: ',(f_global(i),i=1,nbin)
       write(fd,*) 'Checksum, N=SUM(f)=',SUM(f_global(1:nbin))
    endif

    ! Do cumulative integral of f(ibin) and set pivots where int(f) = iproc/nprocs*N

    f_integral = 0.
    ave_work = SUM(f_global(1:nbin))/nprocs  ! Average # particles/PE
    fpval(1) = key_min  ! lowest pivot
    i=1

    do ibin = 1,nbin
       f_integral = f_integral + f_global(ibin)
       if (f_integral >= i*ave_work ) then  ! If int(f) exceeds multiple of work average, set pivot
          i=i+1
          fpval(i) = fpval(1)+ibin*step  ! Set next highest pivot
       endif
    end do

    if (debug) write (fd,*) 'Pivots: ',fpval(1:nprocs+1)


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
    !     give each process a copy of all the fposts()s

    do i=1,nprocs
       islen(i) = fposts(i+1) - fposts(i)
    enddo

    call MPI_ALLTOALL( islen,one,MPI_INTEGER, &
                       irlen,one,MPI_INTEGER, &
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
    call nwaymerge(nppm,npnew,nprocs,keys,irnkl,itabl,itabr)

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: pswssort S<',VTIERR,ICLASSH
!
  end subroutine pswssort


  subroutine nwaymrg(nppm,np,nprocs,keys,irnkl,itabl,itabr)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

    implicit none
    integer :: nppm,np,nprocs
    integer :: irnkl(nppm)
    integer*8 :: keys(nppm)
    integer :: itabl(nprocs+1),itabr(nprocs)

    integer, parameter :: maxprocs=1024
    integer*8 :: ik,itabk(maxprocs)
    integer :: i,j,k,il,ir,ij,irnkj(maxprocs),jprocs

    !     Make sure that the indicies are "one" baised.
!VAMPINST subroutine_start
       CALL VTENTER(IF_nwaymrg,VTNOSCL,VTIERR)
!      write(*,*) 'VT: nwaymrg S>',VTIERR,
!     *    IF_nwaymrg,ICLASSH
!
    if (itabl(1) .ne. 1) then
       do i=nprocs+1,1,-1
          itabl(i) = itabl(i) - itabl(1) + 1
       enddo
    endif

    do i=1,nprocs
       itabr(i) = itabl(i+1) - 1
    enddo

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

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: nwaymrg S<',VTIERR,ICLASSH
!
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
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    implicit none
    integer*8, intent(inout) :: iarr(:)
    integer :: i,n

!VAMPINST subroutine_start
       CALL VTENTER(IF_sort_i,VTNOSCL,VTIERR)
!      write(*,*) 'VT: sort_i S>',VTIERR,
!     *    IF_sort_i,ICLASSH
!
    n = size(iarr)


    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call swap( iarr(1),iarr(i) )      ! Clear space at end of array and retire top of heap into it
       call sift_down( 1,i-1)
    end do

!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: sort_i S<',VTIERR,ICLASSH
!
  contains
    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain 
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
      integer, intent(in) :: l,r     ! the heap structure
      integer :: j,jold    ! index
      integer*8 :: a

!VAMPINST subroutine_start
       CALL VTENTER(IF_sift_down,VTNOSCL,VTIERR)
!      write(*,*) 'VT: sift_down S>',VTIERR,
!     *    IF_sift_down,ICLASSH
!
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

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: sift_down S<',VTIERR,ICLASSH
!
    end subroutine sift_down


  end subroutine sort_i


  subroutine indsort_i(iarr,list,n,nppm)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    implicit none

    integer, intent(in) :: n,nppm

    integer*8, dimension(nppm), intent(in) :: iarr
    integer, dimension(nppm), intent(inout) :: list
    integer :: i, indxt, ir, l, j
    integer*8 :: q

!VAMPINST subroutine_start
       CALL VTENTER(IF_indsort_i,VTNOSCL,VTIERR)
!      write(*,*) 'VT: indsort_i S>',VTIERR,
!     *    IF_indsort_i,ICLASSH
!
    list(1:n) =  (/ (i,i=1,n) /)  

    if(n==1)  Then
!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: indsort_i S<',VTIERR,ICLASSH
!
      RETURN
      END IF
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
!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: indsort_i S<',VTIERR,ICLASSH
!
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

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: indsort_i S<',VTIERR,ICLASSH
!
  end subroutine indsort_i

  subroutine swap_ab(p,q)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    integer*8 :: p,q, dum
!VAMPINST subroutine_start
       CALL VTENTER(IF_swap_ab,VTNOSCL,VTIERR)
!      write(*,*) 'VT: swap_ab S>',VTIERR,
!     *    IF_swap_ab,ICLASSH
!
    dum = p
    p=q
    q = dum
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: swap_ab S<',VTIERR,ICLASSH
!
  end subroutine swap_ab


  subroutine blankn(ichan)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
!VAMPINST subroutine_start
       CALL VTENTER(IF_blankn,VTNOSCL,VTIERR)
!      write(*,*) 'VT: blankn S>',VTIERR,
!     *    IF_blankn,ICLASSH
!
    write(ichan,'(/)')
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: blankn S<',VTIERR,ICLASSH
!
  end subroutine blankn

  subroutine blank6
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
!VAMPINST subroutine_start
       CALL VTENTER(IF_blank6,VTNOSCL,VTIERR)
!      write(*,*) 'VT: blank6 S>',VTIERR,
!     *    IF_blank6,ICLASSH
!
    write(6,'(/)')
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: blank6 S<',VTIERR,ICLASSH
!
  end subroutine blank6

  subroutine cput(sec)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    real :: sec
    integer :: ic1, ir1, im1
    !    CALL SYSTEM_CLOCK(COUNT=IC1, COUNT_RATE=IR1, COUNT_MAX=IM1)
    !    sec = 1.*ic1/ir1
!VAMPINST subroutine_start
       CALL VTENTER(IF_cput,VTNOSCL,VTIERR)
!      write(*,*) 'VT: cput S>',VTIERR,
!     *    IF_cput,ICLASSH
!
    call cpu_time(sec)
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: cput S<',VTIERR,ICLASSH
!
  end subroutine cput

  ! Random number scrambler
  ! =======================
  !
  !  called with: 
  !               x=rano(iseed)
  !
  !  returns real number in the interval (0.0 - 1.0)
  !  set iseed = -ve integer before using call
  !
  !  Example of calling routine:
  !
  !      subroutine coords
  !      include 'common.h'
  !
  !      iseed1 = -11
  !      iseed2 = -7
  !
  !
  !      do i=1,n
  !	x(i)=xlen*rano(iseed1)
  !       y(i)=ylen*rano(iseed2)
  !      end do
  !
  !
  !      end


  ! - routine taken from Numerical Recipies, p195
  !
  real function rano(idum)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    implicit none
    integer :: idum
    real, save :: dseed, dum
    real, save :: v(97), y
    integer, save :: iff, icall, i, j
    data iff,icall/0,0/
!VAMPINST subroutine_start
       CALL VTENTER(IF_rano,VTNOSCL,VTIERR)
!      write(*,*) 'VT: rano S>',VTIERR,
!     *    IF_rano,ICLASSH
!
    if (idum.lt.0.or.iff.eq.0) then
       !  first call generates new sequence
       iff = 1
       dseed=abs(idum)*1.0
       idum=1
       do  j=1,97
          dum=genran(dseed)
       end do
       do j=1,97
          v(j)=genran(dseed)
       end do
       y=genran(dseed)
    endif

    !  next index - make sure we don't overstep array bounds if
    !  generator returns a 0.0 or 1.0

    j=max(mod(1+int(97.*y),98),1)
    if(j.gt.97.or.j.lt.1) then
       write (6,*) 'Call: ',icall
       write (6,*) 'idum = ',idum,'j = ',j,' y= ',y
       write (6,*) 'Random No. generator not initialised properly'
       write (6,*) 'dummy =',dum,' dseed=',dseed
       write (6,100) (i,v(i),i=1,97)
100    format (i4,f10.6)
!VAMPINST stop
!      CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: rano S<',VTIERR,ICLASSH
!
       stop
    endif
    !  get next variate and generate new one to fill gap

    y=v(j)
    rano=y
    v(j)=genran(dseed)
    icall = icall + 1

!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: rano S<',VTIERR,ICLASSH
!
    return
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: rano S<',VTIERR,ICLASSH
!
  end function rano


  real function genran (dseed)                                      
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    !                                  specifications for arguments         
    real ::  dseed                                          
    !                                  specifications for local variables   
    !    real ::  d2p31m,d2p31                                   
    real ::  d2p31m,d2p31    
    !                                  d2p31m=(2**31) - 1                   
    !                                  d2p31 =(2**31)(or an adjusted value) 
    data               d2p31m/2147483647.0/                          
    data               d2p31 /2147483648.0/                          
    !                                  first executable statement           
!VAMPINST subroutine_start
       CALL VTENTER(IF_genran,VTNOSCL,VTIERR)
!      write(*,*) 'VT: genran S>',VTIERR,
!     *    IF_genran,ICLASSH
!
    dseed = mod(16807.0*dseed,d2p31m)                               
    genran = dseed / d2p31                                            
!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: genran S<',VTIERR,ICLASSH
!
    return                                                            
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: genran S<',VTIERR,ICLASSH
!
  end function genran

  !  ========================================
  !
  !      function   PHASE
  !
  !  Returns -pi -> pi (4-quadrant) phase of complex pair (x,y)
  !
  !  ========================================

  real function phase(x,y)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
    implicit none
    real, parameter :: pi=3.1415926536
    real, intent(in) :: x,y
    integer :: sx,sy,itx

!VAMPINST subroutine_start
       CALL VTENTER(IF_phase,VTNOSCL,VTIERR)
!      write(*,*) 'VT: phase S>',VTIERR,
!     *    IF_phase,ICLASSH
!
    sx=sign(1.,x)
    sy=sign(1.,y)
    itx = (1-sx)/2  ! 0 or 1

    ! special cases first
    if (x.eq.0 .and. y.eq.0) then
       phase = 0.
    else if (x.eq.0.) then
       phase = sy*pi/2.
    else if (y.eq.0) then
       phase = itx*pi
    else
       phase = atan(y/x)+itx*sy*pi
    endif

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: phase S<',VTIERR,ICLASSH
!
  end function phase


end module utils
