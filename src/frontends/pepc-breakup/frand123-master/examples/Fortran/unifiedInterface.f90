!
! Example illustrating the use of the unified interface to all random number
! generator subroutines offered by frand123.
!

program unifiedInterface
   use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, &
                                          c_double, c_float
   use frand123, only: frand123_state_kind, frand123_state_size, &
                       frand123Init, frand123Rand

   ! state
   integer(kind=frand123_state_kind), dimension(frand123_state_size) :: state

   ! random numbers of different types
   integer( kind = c_int64_t ) :: scalarInt64
   integer( kind = c_int64_t ), dimension( 3 ) :: vecInt64
   integer( kind = c_int32_t ) :: scalarInt32
   integer( kind = c_int32_t ), dimension( 3 ) :: vecInt32
   real( kind = c_float ) :: scalarSingle
   real( kind = c_float ), dimension( 3 ) :: vecSingle
   real( kind = c_double ) :: scalarDouble
   real( kind = c_double ), dimension( 3 ) :: vecDouble
   real( kind = c_float ) :: scalarNormSingle
   real( kind = c_float ), dimension( 3 ) :: vecNormSingle
   real( kind = c_double ) :: scalarNormDouble
   real( kind = c_double ), dimension( 3 ) :: vecNormDouble

   ! initialize state
   call frand123Init( state, 0, 0 )
   
   ! generate 64-bit integers
   call frand123Rand( state, scalarInt64 )
   write(*,*) 'Scalar 64-bit integer', scalarInt64
   call frand123Rand( state, vecInt64 )
   write(*,*) 'Vector of 64-bit integers', vecInt64
   
   ! generate 32-bit integers
   call frand123Rand( state, scalarInt32 )
   write(*,*) 'Scalar 32-bit integer', scalarInt32
   call frand123Rand( state, vecInt32 )
   write(*,*) 'Vector of 32-bit integers', vecInt32

   ! generate double precision reals
   call frand123Rand( state, scalarDouble )
   write(*,*) 'Scalar double precision real', scalarDouble
   call frand123Rand( state, vecDouble )
   write(*,*) 'Vector of double precision reals', vecDouble

   ! generate single precision reals
   call frand123Rand( state, scalarSingle )
   write(*,*) 'Scalar single precision real', scalarSingle
   call frand123Rand( state, vecSingle )
   write(*,*) 'Vector of single precision reals', vecSingle

   ! generate double precision normal reals
   call frand123Rand( state, 0.1d0, 1.6d0, scalarNormDouble )
   write(*,*) 'Scalar double precision normal real', scalarNormDouble
   call frand123Rand( state, 0.1d0, 1.6d0, vecNormDouble )
   write(*,*) 'Vector of double precision normal reals', vecNormDouble

   ! generate single precision normal reals
   call frand123Rand( state, scalarNormalSingle )
   write(*,*) 'Scalar single precision normal real', scalarNormalSingle
   call frand123Rand( state, vecNormalSingle )
   write(*,*) 'Vector of single precision normal reals', vecNormalSingle
end program
