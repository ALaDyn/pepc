!
!                          Tree sort utility module
!                       *   Parallel Sorting
!

module module_utils

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
  !    Modified 2011.12.02 by Andreas Breslau
  !    to optionally return an integer array containing the original positions of the sorted array (iarr).
  !    This can be used by the calling code to sort several other arrays according to iarr.
  !
  !  ================================

  subroutine sort_i(iarr, map)
    use module_debug
    implicit none
    integer*8, intent(inout) :: iarr(:)
    integer, optional, intent(out) :: map(:)           !< the optional integer array filled with positions for sorting other arrays according to iarr
    integer :: i,n

    if (present(map)) then
      if (size(map) .ne. size(iarr)) then
         DEBUG_ERROR(*, "Error in tree_utils sort_i: optional second parameter has not the same size as first argument. Ignoring second argument.")
      endif
    end if

    n = size(iarr)

    ! only sort map, if map present
    if(present(map)) map = [ (i, i=1,n) ]

    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call sort_swap_ab( 1,i )      ! Clear space at end of array and retire top of heap into it
       call sift_down( 1,i-1)
    end do

  contains
    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain the heap structure
      ! Modified 2011.12.02 by Andreas Breslau to sort the map array according to iarr (see above description)
      integer, intent(in) :: l,r
      integer :: j,jold    ! index
      integer*8 :: a
      integer :: b
      
      a = iarr(l)

      ! only sort map, if map present
      if(present(map)) b= map(l)

      jold = l
      j = l + l
      do                   ! do while j <= r
         if (j > r) exit
         if (j < r) then
            if (iarr(j) < iarr(j+1)) j = j+1
         endif
         if (a >= iarr(j)) exit       ! Found a`s level, so terminate sift-down
         iarr(jold) = iarr(j)

         ! only sort map, if map present
         if(present(map)) map(jold) = map(j)

         jold = j                    ! Demote a and continue
         j = j+j
      end do
      iarr(jold) = a                  ! Put a into its slot

      ! only sort map, if map present
      if(present(map)) map(jold) = b

    end subroutine sift_down

    subroutine sort_swap_ab(p,q)
      integer :: p,q, dum
      integer*8 dum8

      dum8 = iarr(p)
      iarr(p)=iarr(q)
      iarr(q) = dum8

      ! only sort map, if map present
      if(present(map)) then
         dum = map(p)
         map(p)=map(q)
         map(q) = dum
      end if

    end subroutine sort_swap_ab

  end subroutine sort_i


  subroutine swap_ab(p,q)
    integer*8 :: p,q, dum
    dum = p
    p=q
    q = dum
  end subroutine swap_ab

end module module_utils
