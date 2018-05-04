module frand123CInterfaces
   use, intrinsic :: iso_c_binding, only: c_float, c_double, c_int64_t, c_int32_t
   implicit none

   ! kind parameter for state variables
   integer, parameter, public :: frand123_state_kind = c_int64_t
   ! size of the state (might be different if additional RNGs be supported
   integer, parameter, public :: frand123_state_size = 4
      
   ! interfaces to C functions frand123Double_scalar and frand123Double
   ! generating random double precision numbers uniformly distributed in (0,1)
   interface
      function frand123Double_scalar_C( state ) &
               bind( C, name='frand123Double_scalar' )
         use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         real( kind = c_double ) :: frand123Double_scalar_C
      end function
      subroutine frand123Double_C( state, lenRes, res ) &
                 bind( C, name='frand123Double' )
         use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         integer( kind = c_int64_t ), value, intent( in ) :: lenRes
         real( kind = c_double ), dimension( lenRes ), intent( inout ) :: res
      end subroutine
   end interface

   ! interfaces to C functions frand123Single_scalar and frand123Single
   ! generating random single precision numbers uniformly distributed in (0,1)
   interface
      function frand123Single_scalar_C( state ) &
               bind( C, name='frand123Single_scalar' )
         use, intrinsic :: iso_c_binding, only: c_float, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         real( kind = c_float ) :: frand123Single_scalar_C
      end function
      subroutine frand123Single_C( state, lenRes, res ) &
                 bind( C, name='frand123Single' )
         use, intrinsic :: iso_c_binding, only: c_float, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         integer( kind = c_int64_t ), value, intent( in ) :: lenRes
         real( kind = c_float ), dimension( lenRes ), intent( inout ) :: res
      end subroutine
   end interface

   ! interfaces to C functions frand123NormDouble_scalar and frand123NormDouble
   ! generating rand double precision numbers normally distributed with
   ! expectation mu and variance sigma
   interface
      function frand123NormDouble_scalar_C( state, mu, sigma ) &
               bind( C, name='frand123NormDouble_scalar' )
         use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         real( kind = c_double ), value, intent( in ) :: mu
         real( kind = c_double ), value, intent( in ) :: sigma
         real( kind = c_double ) :: frand123NormDouble_scalar_C
      end function
      subroutine frand123NormDouble_C( state, mu, sigma, lenRes, res ) &
                 bind( C, name='frand123NormDouble' )
         use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         real( kind = c_double ), value, intent( in ) :: mu
         real( kind = c_double ), value, intent( in ) :: sigma
         integer( kind = c_int64_t ), value, intent( in ) :: lenRes
         real( kind = c_double ), dimension( lenRes ), intent( inout ) :: res
      end subroutine
   end interface

   ! interfaces to C functions frand123NormSingle_scalar and frand123NormSingle
   ! generating random single precision numbers normally distributed with
   ! expectation mu and variance sigma
   interface
      function frand123NormSingle_scalar_C( state, mu, sigma ) &
               bind( C, name='frand123NormSingle_scalar' )
         use, intrinsic :: iso_c_binding, only: c_float, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         real( kind = c_float ), value, intent( in ) :: mu
         real( kind = c_float ), value, intent( in ) :: sigma
         real( kind = c_float ) :: frand123NormSingle_scalar_C
      end function
      subroutine frand123NormSingle_C( state, mu, sigma, lenRes, res ) &
                 bind( C, name='frand123NormSingle' )
         use, intrinsic :: iso_c_binding, only: c_float, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         real( kind = c_float ), value, intent( in ) :: mu
         real( kind = c_float ), value, intent( in ) :: sigma
         integer( kind = c_int64_t ), value, intent( in ) :: lenRes
         real( kind = c_float ), dimension( lenRes ), intent( inout ) :: res
      end subroutine
   end interface

   ! interfaces to C functions frand123Integer64_scalar and frand123Integer64
   ! generating random 64-bit numbers discretely uniformly distributed in
   ! INT64_MIN,..,INT64_MAX
   interface
      function frand123Integer64_scalar_C( state ) &
               bind( C, name='frand123Integer64_scalar' )
         use, intrinsic :: iso_c_binding, only: c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         integer( kind = c_int64_t ) :: frand123Integer64_scalar_C
      end function
      subroutine frand123Integer64_C( state, lenRes, res ) bind( C, name='frand123Integer64' )
         use, intrinsic :: iso_c_binding, only: c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         integer( kind = c_int64_t ), value, intent( in ) :: lenRes
         integer( kind = c_int64_t ), dimension( lenRes ), intent( inout ) :: res
      end subroutine
   end interface

   ! interfaces to C functions frand123Integer32_scalar and frand123Integer32
   ! generating random 32-bit numbers discretely uniformly distributed in
   ! INT32_MIN,..,INT32_MAX
   interface
      function frand123Integer32_scalar_C( state ) &
               bind( C, name='frand123Integer32_scalar' )
         use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         integer( kind = c_int32_t ) :: frand123Integer32_scalar_C
      end function
      subroutine frand123Integer32_C( state, lenRes, res ) bind( C, name='frand123Integer32' )
         use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t, c_int64_t
         implicit none
         integer( kind = c_int64_t ), dimension( 4 ), intent( inout ) :: state
         integer( kind = c_int64_t ), value, intent( in ) :: lenRes
         integer( kind = c_int32_t ), dimension( lenRes ), intent( inout ) :: res
      end subroutine
   end interface

end module frand123CInterfaces
