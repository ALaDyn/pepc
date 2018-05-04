!
! Example illustrating the use of the frand123 C interface for the generation
! of single precision uniformly distributed random numbers.
! The example uses a naive Monte-Carlo approach to approximate pi.
!
! It generates points within the unit square and computes the ratio of points
! within the upper right quadrant of the unit circle.
! This ratio approximates pi / 4
!

program piSingle
   use, intrinsic :: iso_c_binding, only: c_float, c_double
   use omp_lib, only: omp_get_wtime
   use frand123, only: frand123_state_kind, frand123_state_size, &
                       frand123Init, frand123Single
   implicit none

   ! integer kind for large integers
   integer, parameter :: ik = selected_int_kind( 10 )

   ! number of points used for the approximation
   integer( kind = ik ), parameter :: numberOfPoints = 1000 * 1000 * 100   

   ! state for frand123
   integer(kind=frand123_state_kind), dimension(frand123_state_size) :: state

   ! timing variables
   real( kind = c_double ) :: startTime, endTime

   ! result variables
   real( kind = c_float ) :: varPiScalar, varPiTwo, varPiFour, varPi2000

   ! initialize state
   call frand123Init( state, 0, 0 )

   ! use piScalar
   startTime = omp_get_wtime()
   varPiScalar = piScalar( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating  one random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) varPiScalar, endTime - startTime

   ! use piTwo
   startTime = omp_get_wtime()
   varPiTwo = piTwo( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating  two random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) varPiTwo, endTime - startTime

   ! use piFour
   startTime = omp_get_wtime()
   varPiFour = piFour( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating four random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) varPiFour, endTime - startTime

   ! use pi2000
   startTime = omp_get_wtime()
   varPi2000 = pi2000( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating 2000 random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) varPiTwo, endTime - startTime

contains
   
   ! comput pi drawing a single single precision uniformly distributed random
   ! number at a time
   real( kind = c_float ) function piScalar( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      integer( kind = ik ) :: inside ! counter for points inside unit circle
      integer( kind = ik ) :: i ! loop variable
      real( kind = c_float ) :: x, y ! coordinates of the point
      real( kind = c_float ) :: r2 ! square of the distance to (0,0)

      ! approximate pi
      inside = 0
      do i = 1, numberOfPoints
         ! generate point in the unit square
         call frand123Single( state, x )
         call frand123Single( state, y )
         ! compute square of radius and check whether inside unit circle
         r2 = x * x + y * y
         if( r2 .lt. 1.d0 ) then
            inside = inside + 1
         endif
      enddo
      ! pi = 4 * ratio between points inside unit circle and points generated
      piScalar = 4.d0 * real( inside, c_float ) / &
                        real( numberOfPoints, c_float )
   end function piScalar
   
   ! comput pi drawing two single precision uniformly distributed random
   ! number at a time
   real( kind = c_float ) function piTwo( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      integer( kind = ik ) :: inside ! counter for points inside unit circle
      integer( kind = ik ) :: i ! loop variable
      real( kind = c_float ), dimension( 2 ) :: pos ! coordinates of the point
      real( kind = c_float ) :: r2 ! square of the distance to (0,0)

      ! approximate pi
      inside = 0
      do i = 1, numberOfPoints
         ! generate point in the unit square
         call frand123Single( state, pos )
         ! compute square of radius and check whether inside unit circle
         r2 = sum( pos ** 2 )
         if( r2 .lt. 1.d0 ) then
            inside = inside + 1
         endif
      enddo
      ! pi = 4 * ratio between points inside unit circle and points generated
      piTwo = 4.d0 * real( inside, c_float ) / &
                     real( numberOfPoints, c_float )
   end function piTwo
   
   ! comput pi drawing four single precision uniformly distributed random
   ! number at a time
   real( kind = c_float ) function piFour( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      integer( kind = ik ) :: inside ! counter for points inside unit circle
      integer( kind = ik ) :: i ! loop variable
      real( kind = c_float ), dimension( 4 ) :: pos ! coordinates of the point
      real( kind = c_float ) :: r2 ! square of the distance to (0,0)

      ! approximate pi
      inside = 0
      do i = 1, numberOfPoints, 2
         ! generate point in the unit square
         call frand123Single( state, pos )
         ! compute square of radius and check whether inside unit circle
         ! first point
         r2 = pos( 1 )**2 + pos(2)**2
         if( r2 .lt. 1.d0 ) then
            inside = inside + 1
         endif
         ! second point
         r2 = pos(3)**2 + pos(4)**2
         if( r2 .lt. 1.d0 ) then
            inside = inside + 1
         endif
      enddo
      ! pi = 4 * ratio between points inside unit circle and points generated
      piFour = 4.d0 * real( inside, c_float ) / &
                      real( numberOfPoints, c_float )
   end function piFour
   
   ! comput pi drawing 2000 single precision uniformly distributed random
   ! number at a time
   real( kind = c_float ) function pi2000( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      integer( kind = ik ) :: inside ! counter for points inside unit circle
      integer( kind = ik ) :: i, j ! loop variables
      real( kind = c_float ), dimension( 2000 ) :: buffer ! buffer for RNG
      real( kind = c_float ) :: r2 ! square of the distance to (0,0)

      ! approximate pi
      inside = 0
      do i = 1, numberOfPoints, 1000
         ! generate point in the unit square
         call frand123Single( state, buffer )
         do j = 1, 1000
            ! compute square of radius and check whether inside unit circle
            r2 = sum( buffer( 2 * j - 1:2 * j ) ** 2 )
            if( r2 .lt. 1.d0 ) then
               inside = inside + 1
            endif
         enddo
      enddo
      ! pi = 4 * ratio between points inside unit circle and points generated
      pi2000 = 4.d0 * real( inside, c_float ) / &
                      real( numberOfPoints, c_float )
   end function pi2000
end program piSingle
