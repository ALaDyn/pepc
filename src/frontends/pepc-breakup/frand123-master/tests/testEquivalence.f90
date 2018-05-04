program testEquivalence
   use, intrinsic :: iso_c_binding
   use frand123
   implicit none

   ! state
   integer(kind=frand123_state_kind), dimension(frand123_state_size) :: state

   ! didTestPass
   logical :: passed

   ! test frand123Double
   passed = equivalenceDouble( state )
   if( passed ) then
      write(*,*) 'Test frand123Double passed'
   else
      write(*,*) 'Test frand123Double failed'
      exit( 1 )
   endif

   ! test frand123Single
   passed = equivalenceSingle( state )
   if( passed ) then
      write(*,*) 'Test frand123Single passed'
   else
      write(*,*) 'Test frand123Single failed'
      exit( 1 )
   endif

   ! test frand123NormDouble
   passed = equivalenceNormDouble( state )
   if( passed ) then
      write(*,*) 'Test frand123NormDouble passed'
   else
      write(*,*) 'Test frand123NormDouble failed'
      exit( 1 )
   endif

   ! test frand123NormSingle
   passed = equivalenceNormSingle( state )
   if( passed ) then
      write(*,*) 'Test frand123NormSingle passed'
   else
      write(*,*) 'Test frand123NormSingle failed'
      exit( 1 )
   endif

   ! test frand123Integer64
   passed = equivalenceInteger64( state )
   if( passed ) then
      write(*,*) 'Test frand123Integer64 passed'
   else
      write(*,*) 'Test frand123Integer64 failed'
      exit( 1 )
   endif

   ! test frand123Integer32
   passed = equivalenceInteger32( state )
   if( passed ) then
      write(*,*) 'Test frand123Integer32 passed'
   else
      write(*,*) 'Test frand123Integer32 failed'
      exit( 1 )
   endif

contains

   ! test equivalence of generating 10 x 2 doubles and 1 x 20 doubles
   logical function equivalenceDouble( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
      real( kind = c_double ), dimension( 20 ) :: r10x2, r1x20
      integer :: i

      ! overwrite in case of failure
      equivalenceDouble = .true.

      ! generate random numbers starting from same state
      seed = (/ 1, 2 /)
      call frand123Init( state, 4, 5, seed )
      do i = 1, 10
         call frand123Double( state, r10x2( 2*i-1:2*i ) )
      enddo
      call frand123Init( state, 4, 5, seed )
      call frand123Double( state, r1x20 )

      ! compare random numbers
      do i = 1, 20
         if( abs( r10x2( i ) - r1x20( i ) ) .gt. 1e-15 ) then
            write(*, '( "Entries ", I2, ": absolute difference: ", ES11.4 )' ) &
                        i, abs( r10x2( i ) - r1x20( i ) )
            equivalenceDouble = .false.
         endif
      enddo
   end function equivalenceDouble

   ! test equivalence of generating 5 x 4 singles and 1 x 20 singles
   logical function equivalenceSingle( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
      real( kind = c_float ), dimension( 20 ) :: r5x4, r1x20
      integer :: i

      ! overwrite in case of failure
      equivalenceSingle = .true.

      ! generate random numbers starting from same state
      seed = (/ 1, 2 /)
      call frand123Init( state, 4, 5, seed )
      do i = 1, 5
         call frand123Single( state, r5x4( 4*i-3:4*i ) )
      enddo
      call frand123Init( state, 4, 5, seed )
      call frand123Single( state, r1x20 )

      ! compare random numbers
      do i = 1, 20
         if( abs( r5x4( i ) - r1x20( i ) ) .gt. 1e-7 ) then
            write(*, '( "Entries ", I2, ": absolute difference: ", ES11.4 )' ) &
                        i, abs( r5x4( i ) - r1x20( i ) )
            equivalenceSingle = .false.
         endif
      enddo
   end function equivalenceSingle

   ! test equivalence of generating 10 x 2 normal doubles and 1 x 20 normal doubles
   logical function equivalenceNormDouble( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
      real( kind = c_double ), dimension( 20 ) :: r10x2, r1x20
      integer :: i

      ! overwrite in case of failure
      equivalenceNormDouble = .true.

      ! generate random numbers starting from same state
      seed = (/ 1, 2 /)
      call frand123Init( state, 4, 5, seed )
      do i = 1, 10
         call frand123NormDouble( state, 2.d0, 3.d0, r10x2( 2*i-1:2*i ) )
      enddo
      call frand123Init( state, 4, 5, seed )
      call frand123NormDouble( state, 2.d0, 3.d0, r1x20 )

      ! compare random numbers
      do i = 1, 20
         if( abs( r10x2( i ) - r1x20( i ) ) .gt. 1e-15 ) then
            write(*, '( "Entries ", I2, ": absolute difference: ", ES11.4 )' ) &
                        i, abs( r10x2( i ) - r1x20( i ) )
            equivalenceNormDouble = .false.
         endif
      enddo
   end function equivalenceNormDouble

   ! test equivalence of generating 5 x 4 normal singles and 1 x 20 normal singles
   logical function equivalenceNormSingle( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
      real( kind = c_float ), dimension( 20 ) :: r5x4, r1x20
      integer :: i

      ! overwrite in case of failure
      equivalenceNormSingle = .true.

      ! generate random numbers starting from same state
      seed = (/ 1, 2 /)
      call frand123Init( state, 4, 5, seed )
      do i = 1, 5
         call frand123NormSingle( state, 2., 3., r5x4( 4*i-3:4*i ) )
      enddo
      call frand123Init( state, 4, 5, seed )
      call frand123NormSingle( state, 2., 3., r1x20 )

      ! compare random numbers
      do i = 1, 20
         if( abs( r5x4( i ) - r1x20( i ) ) .gt. 1e-7 ) then
            write(*, '( "Entries ", I2, ": absolute difference: ", ES11.4 )' ) &
                        i, abs( r5x4( i ) - r1x20( i ) )
            equivalenceNormSingle = .false.
         endif
      enddo
   end function equivalenceNormSingle

   ! test equivalence of generating 10 x 2 64-bit integers and 1 x 20 64-bit integers
   logical function equivalenceInteger64( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
      integer( kind = c_int64_t ), dimension( 20 ) :: r10x2, r1x20
      integer :: i

      ! overwrite in case of failure
      equivalenceInteger64 = .true.

      ! generate random numbers starting from same state
      seed = (/ 1, 2 /)
      call frand123Init( state, 4, 5, seed )
      do i = 1, 10
         call frand123Integer64( state, r10x2( 2*i-1:2*i ) )
      enddo
      call frand123Init( state, 4, 5, seed )
      call frand123Integer64( state, r1x20 )

      ! compare random numbers
      do i = 1, 20
         if( r10x2( i ) .ne. r1x20( i ) ) then
            write(*, '( "Entries ", I2, ": difference: ", I20 )' ) &
                        i, r10x2( i ) - r1x20( i )
            equivalenceInteger64 = .false.
         endif
      enddo
   end function equivalenceInteger64

   ! test equivalence of generating 5 x 4 32-bit integers and 1 x 20 32-bit integers
   logical function equivalenceInteger32( state )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), &
                                             intent( inout ) :: state

      integer( kind = frand123_state_kind ), dimension( 2 ) :: seed
      integer( kind = c_int32_t ), dimension( 20 ) :: r5x4, r1x20
      integer :: i

      ! overwrite in case of failure
      equivalenceInteger32 = .true.

      ! generate random numbers starting from same state
      seed = (/ 1, 2 /)
      call frand123Init( state, 4, 5, seed )
      do i = 1, 5
         call frand123Integer32( state, r5x4( 4*i-3:4*i ) )
      enddo
      call frand123Init( state, 4, 5, seed )
      call frand123Integer32( state, r1x20 )

      ! compare random numbers
      do i = 1, 20
         if( r5x4( i ) .ne. r1x20( i ) ) then
            write(*, '( "Entries ", I2, ": difference: ", I20 )' ) &
                        i, r5x4( i ) - r1x20( i )
            equivalenceInteger32 = .false.
         endif
      enddo
   end function equivalenceInteger32

end program testEquivalence
