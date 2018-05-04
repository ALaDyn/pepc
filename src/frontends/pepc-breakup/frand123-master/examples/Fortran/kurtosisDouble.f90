!
! Example illustrating the use of the frand123 Fortran interface for the
! generation of double precision normally distributed random numbers.
! The example estimates the excess kurtosis of the underlying normal
! distribution with expectation mu = 0.7 and variance sigma = 1.6
!

program kurtosisDouble
   use, intrinsic :: iso_c_binding, only: c_double
   use omp_lib, only: omp_get_wtime
   use frand123, only: frand123_state_kind, frand123_state_size, &
                       frand123Init, frand123NormDouble
   implicit none

   ! integer kind for large integers
   integer, parameter :: ik = selected_int_kind( 10 )

   ! number of draws used to estimate the excess kurtosis
   integer( kind = ik ), parameter :: numberOfSamples = 1000 * 1000 * 100

   ! expectation and variance of the normal distribution to draw from
   real( kind = c_double ), parameter :: mu = real( 0.7, c_double )
   real( kind = c_double ), parameter :: sigma = real( 1.6, c_double )

   ! state for frand123
   integer(kind=frand123_state_kind), dimension(frand123_state_size) :: state

   ! timing variables
   real( kind = c_double ) :: startTime, endTime

   ! result variables
   real( kind = c_double ) :: varScalar, varTwo, var1000

   ! initialize state
   call frand123Init( state, 0, 0 )

   ! use excessKurtosisScalar
   startTime = omp_get_wtime()
   varScalar = excessKurtosisScalar( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating  one random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) varScalar, endTime - startTime

   ! use excessKurtosisTwo
   startTime = omp_get_wtime()
   varTwo = excessKurtosisTwo( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating  two random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) varTwo, endTime - startTime

   ! use excessKurtosis1000
   startTime = omp_get_wtime()
   var1000 = excessKurtosis1000( state )
   endTime = omp_get_wtime()
   write(*, '( "Result generating 1000 random number at a time: ", ES11.4, &
          " took ", ES11.4, " seconds" )' ) var1000, endTime - startTime

contains

   ! estimate excess kurtosis by drawing a single double precision sample at a
   ! time
   real( kind = c_double ) function excessKurtosisScalar( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      real( kind = c_double ) :: excessKurtosis
      real( kind = c_double ) :: draw
      integer( kind = ik ) :: i

      ! iterate over samples
      excessKurtosis = real( 0, c_double )
      do i = 1, numberOfSamples
         ! draw random number
         call frand123NormDouble( state, mu, sigma, draw )
         ! centralize
         draw = draw - mu
         ! sum up contributions
         excessKurtosis = excessKurtosis + draw ** 4
      enddo
      ! apply remaining operations
      excessKurtosisScalar = excessKurtosis / &
          ( real( numberOfSamples, c_double) * sigma ** 4 ) - real( 3., c_double )
   end function excessKurtosisScalar

   ! estimate excess kurtosis by drawing two double precision sample at a time
   real( kind = c_double ) function excessKurtosisTwo( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      real( kind = c_double ) :: excessKurtosis
      real( kind = c_double ), dimension( 2 ) :: draw
      integer( kind = ik ) :: i

      ! iterate over samples
      excessKurtosis = real( 0, c_double )
      do i = 1, numberOfSamples, 2
         ! draw random number
         call frand123NormDouble( state, mu, sigma, draw )
         ! centralize
         draw = draw - mu
         ! sum up contributions
         excessKurtosis = excessKurtosis + sum( draw ** 4 )
      enddo
      ! apply remaining operations
      excessKurtosisTwo = excessKurtosis / &
          ( real( numberOfSamples, c_double) * sigma ** 4 ) - real( 3., c_double )
   end function excessKurtosisTwo

   ! estimate excess kurtosis by drawing 1000 double precision sample at a time
   real( kind = c_double ) function excessKurtosis1000( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      ! local variables
      real( kind = c_double ) :: excessKurtosis
      real( kind = c_double ), dimension( 1000 ) :: draw
      integer( kind = ik ) :: i

      ! iterate over samples
      excessKurtosis = real( 0, c_double )
      do i = 1, numberOfSamples, 1000
         ! draw random number
         call frand123NormDouble( state, mu, sigma, draw )
         ! centralize
         draw = draw - mu
         ! sum up contributions
         excessKurtosis = excessKurtosis + sum( draw ** 4 )
      enddo
      ! apply remaining operations
      excessKurtosis1000 = excessKurtosis / &
          ( real( numberOfSamples, c_double) * sigma ** 4 ) - real( 3., c_double )
   end function excessKurtosis1000
end program kurtosisDouble
