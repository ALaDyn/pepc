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

  interface unique
     module procedure uniq8
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

    save

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

    save

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

    real, dimension(nppm) ::  array, w1, w2
    integer, dimension(nppm) ::  indxl, irnkl
    integer, dimension(nprocs) ::  islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts


    integer :: i,ierr

    save

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
    integer :: lmax, ak

    integer*8 :: buf(maxprocs)

    integer :: i,j,k,step, tag1, ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: fd, nsamp
    character(13) :: cfmt

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
       write(fd,cfmt) 'Sorted sample: ',(work(i),i=1,nsamp)

       k = 1
       do i=1,nprocs*nprocs,nprocs
          fpval(k) = work(i)
          k = k + 1
       enddo
    endif

    call MPI_BCAST(fpval, nprocs+1, MPI_INTEGER8, 0,  MPI_COMM_WORLD, ierr)


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



  subroutine pswssort(nppm,np,npnew,nprocs,iproc,keys, &
       indxl,irnkl, islen,irlen,fposts,gposts,w1,wload,key_box,balance,debug)


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

    implicit  none
    include 'mpif.h'

    integer, intent(in) :: nppm,np,nprocs,iproc
    integer, intent(out) :: npnew
    integer, parameter :: binmult=800000   !TODO: need to reduce size of f() arrays
    integer*8, dimension(nppm) ::  keys, &      ! array of keys to be sorted.
                                   w1       ! work array
    integer, dimension(nppm) ::  indxl, irnkl ! origin locations of the keys 
    logical :: debug, balance
    real :: wload(nppm), w2(nppm)
    integer*8, dimension(2) :: key_box
    integer, dimension(nprocs) :: islen, irlen
    integer, dimension(nprocs+1) :: fposts, gposts !  fencepost index and key values for shuffle
!    integer, parameter :: maxprocs = 1024
    integer :: itabr(nprocs), itabl(nprocs+1)
    real, dimension(binmult)  :: f_local, f_global
    integer*8 :: fpval(nprocs+1)
    integer*8 :: lmax, lmin, key_min, key_max, gkey_min, gkey_max, step ! Key mins and maxes and step size
    integer*8 :: step_reduced
    integer :: nbin, ibin, itag
    integer :: status(MPI_STATUS_SIZE),ierr
    real :: ave_work, f_integral
    integer ::  i,k, fd, nfill
    character(13) :: cfmt

    nbin = binmult  ! must correspond to array size

!    fd = iproc+10
    fd=20
    itag=0

    !     Independent s  !     Note that indx() is a local index on the process.
    call indexsort( keys,indxl, np, nppm )   ! Index sort from Num. Rec.

    do i=1,np
       w1(i) = keys(indxl(i))
       w2(i) = wload(indxl(i))  ! apply sort to work loads too
    enddo

    !     w1 now contains the sorted keys; w2 the sorted loads
    !     indxl is now the local indexes for the sort.

    lmax = w1(np)   ! local max
    lmin = w1(1)    ! local min

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
    
    if (debug.and.iproc==0) write (*,'(3(a12,z20,a12,z20/),a12,i8,a12,z20)') &
         'local min: ',lmin,' local max: ',lmax, &
         'global min: ',gkey_min,' global max: ',gkey_max, &
         ' box_min: ',key_min,' box_max: ',key_max, &
         ' nbin ',nbin,' step', step

    f_local(1:nbin) = 0.0

    ! Find local key distribution
    do k=1,np
          ! bin inside container limits
          ibin = (w1(k)-key_min)/step + 1
          ibin = max(min(ibin,nbin),1)
       if (balance) then
          f_local(ibin) = f_local(ibin) + w2(k)  ! Can include actual force load on particle here
       else
          f_local(ibin) = f_local(ibin) + 1  ! No load balancing - try to get equal # particles  
       endif
    enddo

    if (debug) then
       cfmt = "(/a15,"//achar(mod(nbin/100,10)+48)//achar(mod(nbin/10,10)+48) // achar(mod(nbin,10)+48) // "(f12.4))"
       write(fd,'(a15/(f12.3))') 'Local key distrib: ',(f_local(i),i=1,nbin,nbin/10)
    endif

    ! Global distrib
    call MPI_ALLREDUCE(f_local, f_global, nbin, MPI_REAL8, MPI_SUM,  MPI_COMM_WORLD, ierr )

    if (debug) then
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
	        write (*,*) 'Problem on Processor ',iproc
		write (*,*) 'integral jumps by more than work average: nfill= ',nfill
                write (*,*) 'increment = ', f_global(ibin), ' ave= ',ave_work
          endif
          i=i+1
          fpval(i) = key_min + ibin*step  ! Set next highest pivot
       endif
    end do

    if (debug) then
       write (fd,*) 'Ave work: ',ave_work
       write (fd,*) 'f_Integral: ',f_integral
       write (fd,'(a10/(z20))') 'Pivots: ',fpval(1:nprocs+1)
    endif

    !     Determine segment boundaries. Within each bin, fposts(i) is the
    !     start of the ith shuffle segment.
    fposts(1) = 1
    k = 2
    do i=1,np
       !        The first element may be greater than several fencepost values,
       !        so we must use a do-while loop.
       do while (w1(i) .ge. fpval(k) .and. w1(i).le.key_max)
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


  subroutine indsort_i(iarr,list,n,nppm)
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

  end subroutine indsort_i

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



end module tree_utils
