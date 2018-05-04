program testNormSinglePerformance
   use frand123
   use omp_lib
   use, intrinsic :: iso_c_binding, only: c_float
   implicit none

   integer( kind = frand123_state_kind ), dimension( frand123_state_size ) :: state
   integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
   integer, parameter :: ctr_kind = selected_int_kind( 14 )
   integer( kind = ctr_kind ), parameter :: rounds = 1000 * 1000 * 100
   integer( kind = ctr_kind ) :: i
   real( kind = c_float ), parameter :: mu = 0.d0
   real( kind = c_float ), parameter :: sigma = 1.d0
   real( kind = c_float ), dimension( 2 ) :: buffer
   real( kind = c_float ), dimension( 2 ) :: resArr
   real( kind = c_float ) :: res
   double precision :: startTime, stopTime

   ! serial
   startTime = omp_get_wtime()
   ! initialize state
   seed = (/ 0, 0 /)
   call frand123Init( state, 0, 0, seed)
   ! run for rounds
   resArr = 0.d0
   do i = 1, rounds
      ! generate pair of random variates
      call frand123NormSingle( state, mu, sigma, buffer )
      ! summation
      resArr = resArr + buffer
   enddo
   ! sum up two entries of resArr
   res = sum( resArr ) / rounds
   stopTime = omp_get_wtime()
   write(*,'( "Serial version: result: ", E14.7, ", runtime: ", E11.4 )' ) &
      res, stopTime - startTime

   ! parallel
   startTime = omp_get_wtime();
   !$OMP parallel default( none ) private( buffer, state ) shared( seed, resArr )
   call frand123Init( state, 0, omp_get_thread_num(), seed )
   resArr = 0.d0
   !$OMP do reduction( +: resArr )
   do i = 1, rounds
      ! generate pair of random variates
      call frand123NormSingle( state, mu, sigma, buffer )
      ! summation
      resArr = resArr + buffer
   enddo
   !$OMP end do
   !$OMP end parallel
   ! sum up two entries of resArr
   res = sum( resArr ) / rounds
   stopTime = omp_get_wtime()
   write(*,'( "Parallel version: result: ", E14.7, ", runtime: ", E11.4 )' ) &
      res, stopTime - startTime
end program
