program testOdd
   use, intrinsic :: iso_c_binding
   use frand123
   implicit none

   ! state
   integer(kind=frand123_state_kind), dimension(frand123_state_size) :: state

   ! seed for frand123
   integer( kind = frand123_state_kind ), dimension( 2 ) :: seed

   ! success variable
   logical :: passed

   ! initialize state
   seed = (/ 12, 94 /)
   call frand123Init( state, 9, 5, seed )

   ! test frand123Double
   passed = oddDoubles( state )
   if( passed ) then
      write(*,*) 'oddDoubles passed'
   else
      write(*,*) 'oddDoubles failed'
      stop( 1 )
   endif

   ! test frand123Single
   passed = oddSingles( state )
   if( passed ) then
      write(*,*) 'oddSingle passed'
   else
      write(*,*) 'oddSingle failed'
      stop( 1 )
   endif

   ! test frand123NormDouble
   passed = oddNormDoubles( state )
   if( passed ) then
      write(*,*) 'oddNormDoubles passed'
   else
      write(*,*) 'oddNormDoubles failed'
      stop( 1 )
   endif

   ! test frand123NormSingle
   passed = oddNormSingles( state )
   if( passed ) then
      write(*,*) 'oddNormSingle passed'
   else
      write(*,*) 'oddNormSingle failed'
      stop( 1 )
   endif

contains

   ! test odd number of doubles
   logical function oddDoubles( state )
      implicit none
      integer(kind=frand123_state_kind), dimension(frand123_state_size), &
                                         intent( inout ) :: state
      
      real( kind = c_double ), dimension( 5 ) :: r
      integer :: i

      ! assume positive
      oddDoubles = .true.

      ! init to zero to find errors
      do i = 1, 5
         r( i ) = -2000.d0
      enddo
      ! fill with random numbers
      call frand123Rand( state, r )
      ! check whether all are set
      do i = 1, 5
         if( r( i ) .lt. -1.d0 ) then
            write(*, '( "oddDoubles: index ", I1, " not set by frand123" )' ) i
            oddDoubles = .false.
         endif
      enddo
   end function oddDoubles

   ! test odd number of singles
   logical function oddSingles( state )
      implicit none
      integer(kind=frand123_state_kind), dimension(frand123_state_size), &
                                         intent( inout ) :: state
      
      real( kind = c_float ), dimension( 5 ) :: r
      integer :: i

      ! assume positive
      oddSingles = .true.

      ! init to zero to find errors
      do i = 1, 5
         r( i ) = -2000.
      enddo
      ! fill with random numbers
      call frand123Rand( state, r )
      ! check whether all are set
      do i = 1, 5
         if( r( i ) .lt. -1. ) then
            write(*, '( "oddDoubles: index ", I1, " not set by frand123" )' ) i
            oddSingles = .false.
         endif
      enddo
   end function oddSingles

   ! test odd number of normal doubles
   logical function oddNormDoubles( state )
      implicit none
      integer(kind=frand123_state_kind), dimension(frand123_state_size), &
                                         intent( inout ) :: state
      
      real( kind = c_double ), dimension( 5 ) :: r
      integer :: i

      ! assume positive
      oddNormDoubles = .true.

      ! init to zero to find errors
      do i = 1, 5
         r( i ) = -2000.d0
      enddo
      ! fill with random numbers
      call frand123Rand( state, 10.d0, 1.d0, r )
      ! check whether all are set
      do i = 1, 5
         if( r( i ) .lt. -1.d0 ) then
            write(*, '( "oddDoubles: index ", I1, " not set by frand123" )' ) i
            oddNormDoubles = .false.
         endif
      enddo
   end function oddNormDoubles

   ! test odd number of normal singles
   logical function oddNormSingles( state )
      implicit none
      integer(kind=frand123_state_kind), dimension(frand123_state_size), &
                                         intent( inout ) :: state
      
      real( kind = c_float ), dimension( 5 ) :: r
      integer :: i

      ! assume positive
      oddNormSingles = .true.

      ! init to zero to find errors
      do i = 1, 5
         r( i ) = -2000.
      enddo
      ! fill with random numbers
      call frand123Rand( state, 10., 1., r )
      ! check whether all are set
      do i = 1, 5
         if( r( i ) .lt. -1. ) then
            write(*, '( "oddDoubles: index ", I1, " not set by frand123" )' ) i
            oddNormSingles = .false.
         endif
      enddo
   end function oddNormSingles

end program testOdd
