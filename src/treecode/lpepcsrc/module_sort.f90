! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!>
!> Provides `sort()`, a local heap sort for int64 arrays.
!>
module module_sort
  implicit none
  private

  interface sort
     module procedure sort_i
  end interface

  public sort

  contains

  !>
  !> Sort 64-bit integer array into ascending order using heap-sort algorithm. (Numerical Recipes f90, p1171)
  !>
  !> @note to optionally return an integer array containing the original positions of the sorted array (iarr).
  !> This can be used by the calling code to sort several other arrays according to iarr.
  !>
  subroutine sort_i(iarr, map)
    use module_pepc_types
    use module_debug
    implicit none
    integer(kind_key), intent(inout) :: iarr(:)
    integer(kind_particle), optional, intent(out) :: map(:) !< the optional integer array filled with positions for sorting other arrays according to iarr
    integer(kind_particle) :: i,n

    if (present(map)) then
      DEBUG_ASSERT_MSG(size(map)==size(iarr), *, "Error in tree_utils sort_i: optional second parameter has not the same size as first argument. Ignoring second argument.")
    end if

    n = size(iarr, kind=kind(n))

    ! only sort map, if map present
    if(present(map)) map = [ (i, i=1,n) ]

    do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
       call sift_down(i,n)
    end do

    do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
       ! during retirement/promotion (heap selection) phase.
       call sort_swap_ab( 1_kind_particle, i )      ! Clear space at end of array and retire top of heap into it
       call sift_down(    1_kind_particle, i-1_kind_particle)
    end do

    contains

    subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain the heap structure
      ! Modified 2011.12.02 by Andreas Breslau to sort the map array according to iarr (see above description)
      integer(kind_particle), intent(in) :: l,r
      integer(kind_particle) :: j,jold    ! index
      integer(kind_key) :: a
      integer(kind_particle) :: b
      
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
      integer(kind_particle) :: p,q,dum
      integer(kind_key) ::  dum8

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

end module module_sort
