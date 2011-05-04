!
!                          Tree sort utility module
!                       *   Parallel Sorting
!

module tree_utils

  interface pll_permute
     module procedure psrsperm_i4
     module procedure psrsperm_i8
     module procedure psrsperm_r8
  end interface

  interface sort
     module procedure sort_i
  end interface

  interface swap
     module procedure swap_ab
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


  subroutine swap_ab(p,q)
    integer*8 :: p,q, dum
    dum = p
    p=q
    q = dum
  end subroutine swap_ab

end module tree_utils
