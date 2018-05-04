!
! Example illustrating the use of the frand123 Fortran interface for the
! generation of 32-bit signed integer random numbers.
! The signed integer random numbers ar used as experiment for coin flip.
!

program integer32
   use, intrinsic :: iso_c_binding, only: c_int64_t, c_int32_t
   use frand123, only: frand123_state_kind, frand123_state_size, &
                       frand123Init, frand123Integer32
   implicit none

   ! state for frand123
   integer(kind=frand123_state_kind), dimension(frand123_state_size) :: state

   ! init state
   call frand123Init( state, 0, 0 )

   ! use coinFlipScalar
   write(*,*) 'Generate one random number at a time'
   call coinFlipScalar( state )

   ! use coinFlip10
   write(*,*) 'Generate all 10 random numbers afront'
   call coinFlip10( state )

contains

   ! simulate coin flip by generating one random number at a time
   subroutine coinFlipScalar( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! result of coin flip
      integer( kind = c_int32_t ) :: coin

      ! iteration variable
      integer :: i

      ! simulate 10 coin flips
      ! > 0 heads
      ! < 0 tails
      ! = 0 what a luck
      do i = 1, 10
         ! create 32-bit signed random integer
         call frand123Integer32( state, coin )
         ! result
         if( coin .eq. 0 ) then
            write(*,*) 'What a luck'
         elseif( coin .gt. 0 ) then
            write(*,*) 'heads'
         else
            write(*,*) 'tails'
         endif
      enddo
   end subroutine coinFlipScalar

   ! simulate coin flip by generating all random numbers afront
   subroutine coinFlip10( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! result of coin flip
      integer( kind = c_int32_t ), dimension( 10 ) :: coin

      ! iteration variable
      integer :: i

      ! create 32-bit signed random integers
      call frand123Integer32( state, coin )

      ! simulate 10 coin flips
      ! > 0 heads
      ! < 0 tails
      ! = 0 what a luck
      do i = 1, 10
         ! result
         if( coin( i ) .eq. 0 ) then
            write(*,*) 'What a luck'
         elseif( coin( i ) .gt. 0 ) then
            write(*,*) 'heads'
         else
            write(*,*) 'tails'
         endif
      enddo
   end subroutine coinFlip10
end program integer32
