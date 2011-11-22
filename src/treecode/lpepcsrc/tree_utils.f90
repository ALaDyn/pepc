!
!                          Tree sort utility module
!                       *   Parallel Sorting
!

module tree_utils

  interface sort
     module procedure sort_i
  end interface

  interface swap
     module procedure swap_ab
  end interface

contains


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
