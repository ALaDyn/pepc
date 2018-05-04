module frand123
   use frand123CInterfaces
   implicit none

   ! interfaces allowing use with scalar and vector arguments
   interface frand123Double
      module procedure frand123Double_scalar
      module procedure frand123Double_vector
   end interface frand123Double

   interface frand123Single
      module procedure frand123Single_scalar
      module procedure frand123Single_vector
   end interface frand123Single

   interface frand123NormDouble
      module procedure frand123NormDouble_scalar
      module procedure frand123NormDouble_vector
   end interface frand123NormDouble

   interface frand123NormSingle
      module procedure frand123NormSingle_scalar
      module procedure frand123NormSingle_vector
   end interface frand123NormSingle

   interface frand123Integer64
      module procedure frand123Integer64_scalar
      module procedure frand123Integer64_vector
   end interface frand123Integer64

   interface frand123Integer32
      module procedure frand123Integer32_scalar
      module procedure frand123Integer32_vector
   end interface frand123Integer32

   ! allowing one unified interface
   interface frand123Rand
      module procedure frand123Double_scalar
      module procedure frand123Double_vector
      module procedure frand123Single_scalar
      module procedure frand123Single_vector
      module procedure frand123NormDouble_scalar
      module procedure frand123NormDouble_vector
      module procedure frand123NormSingle_scalar
      module procedure frand123NormSingle_vector
      module procedure frand123Integer64_scalar
      module procedure frand123Integer64_vector
      module procedure frand123Integer32_scalar
      module procedure frand123Integer32_vector
   end interface frand123Rand

   public :: frand123_state_kind
   public :: frand123_state_size
   public :: frand123Double
   public :: frand123Single
   public :: frand123NormDouble
   public :: frand123NormSingle
   public :: frand123Integer32
   public :: frand123Integer64
   public :: frand123Init
   public :: frand123Rand

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                           !!!!! 
   !!!!!  uniform double precision !!!!!
   !!!!!                           !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! generate single random double precision number
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   scalar for random double precision reals in (0,1)
   subroutine frand123Double_scalar( state, res )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_double ), intent( out ) :: res

      res = frand123Double_scalar_C( state )
   end subroutine frand123Double_scalar

   ! generate size(res) random double precision numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   array to be filled with random double precision reals in (0,1)
   subroutine frand123Double_vector( state, res )
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_double ), dimension(:), intent( inout ) :: res

      integer( kind = c_int64_t ) :: len_res

      ! get length of res
      len_res = size( res )

      ! hand over to C implementation
      call frand123Double_C( state, len_res, res )
   end subroutine frand123Double_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                           !!!!! 
   !!!!!  uniform single precision !!!!!
   !!!!!                           !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! generate a single random single precision numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:  scalar for random single precision reals in (0,1)
   subroutine frand123Single_scalar( state, res )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_float ), intent( out ) :: res

      res = frand123Single_scalar_C( state )
   end subroutine frand123Single_scalar

   ! generate size(res) random single precision numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   array to be filled with random single precision reals in (0,1)
   subroutine frand123Single_vector( state, res )
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_float ), dimension(:), intent( inout ) :: res

      integer( kind = c_int64_t ) :: len_res

      ! get length of res
      len_res = size( res )

      ! hand over to C implementation
      call frand123Single_C( state, len_res, res )
   end subroutine frand123Single_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                          !!!!! 
   !!!!!  normal double precision !!!!!
   !!!!!                          !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! generate single normally distributedrandom double precision numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            mu:    expected value
   !            sigma: variance
   !            res:   scalar random double precision reals
   subroutine frand123NormDouble_scalar( state, mu, sigma, res )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_double ), intent( in ) :: mu
      real( kind = c_double ), intent( in ) :: sigma
      real( kind = c_double ), intent( out ) :: res

      res = frand123NormDouble_scalar_C( state, mu, sigma )
   end subroutine frand123NormDouble_scalar

   ! generate size(res) normally distributedrandom double precision numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            mu:    expected value
   !            sigma: variance
   !            res:   array to be filled with random double precision reals
   subroutine frand123NormDouble_vector( state, mu, sigma, res )
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_double ), intent( in ) :: mu
      real( kind = c_double ), intent( in ) :: sigma
      real( kind = c_double ), dimension(:), intent( inout ) :: res

      integer( kind = c_int64_t ) :: len_res

      ! get length of res
      len_res = size( res )

      ! call c implementation
      call frand123NormDouble_C( state, mu, sigma, len_Res, res )
   end subroutine frand123NormDouble_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                          !!!!! 
   !!!!!  normal single precision !!!!!
   !!!!!                          !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! generate single normally distributedrandom single precision numbers
   ! 
   ! Note: The single precision version always uses the polar form of Box-Muller
   !       as this turned out to be on par with Wichura's PPND7 w.r.t.
   !       performance and excels w.r.t. quality of the random numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            mu:    expected value
   !            sigma: variance
   !            res:   scalar for random single precision reals
   subroutine frand123NormSingle_scalar( state, mu, sigma, res )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_float ), intent( in ) :: mu
      real( kind = c_float ), intent( in ) :: sigma
      real( kind = c_float ), intent( out ) :: res

      res = frand123NormSingle_scalar_C( state, mu, sigma )
   end subroutine frand123NormSingle_scalar

   ! generate size(res) normally distributedrandom single precision numbers
   ! 
   ! Note: The single precision version always uses the polar form of Box-Muller
   !       as this turned out to be on par with Wichura's PPND7 w.r.t.
   !       performance and excels w.r.t. quality of the random numbers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            mu:    expected value
   !            sigma: variance
   !            res:   array to be filled with random single precision reals
   subroutine frand123NormSingle_vector( state, mu, sigma, res )
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      real( kind = c_float ), intent( in ) :: mu
      real( kind = c_float ), intent( in ) :: sigma
      real( kind = c_float ), dimension(:), intent( inout ) :: res

      integer( kind = c_int64_t ) :: len_res

      ! get length of res
      len_res = size( res )

      ! hand over to C implementation
      call frand123NormSingle_C( state, mu, sigma, len_res, res )
   end subroutine frand123NormSingle_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                                !!!!! 
   !!!!!  random 64 bit signed integers !!!!!
   !!!!!                                !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! generate single random 64 bit signed integers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   scalar random 64 bit signed integers
   subroutine frand123Integer64_scalar( state, res )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout) :: state
      integer( kind = c_int64_t ), intent( out ) :: res

      res = frand123Integer64_scalar_C( state )
   end subroutine frand123Integer64_scalar

   ! generate size(res) random 64 bit signed integers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   array to be filled with random 64 bit signed integers
   subroutine frand123Integer64_vector( state, res )
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout) :: state
      integer( kind = c_int64_t ), dimension(:), intent( inout ) :: res

      integer( kind = c_int64_t ) :: len_res

      ! get length of res
      len_res = size( res )

      ! leave rest to C implementation
      call frand123Integer64_C( state, len_res, res )
   end subroutine frand123Integer64_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                                !!!!! 
   !!!!!  random 32 bit signed integers !!!!!
   !!!!!                                !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! generate single random 32 bit signed integers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   scalar for random 32 bit signed integers
   subroutine frand123Integer32_scalar( state, res )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      integer( kind = c_int32_t ), intent( out ) :: res

      res = frand123Integer32_scalar_C( state )
   end subroutine frand123Integer32_scalar

   ! generate size(res) random 32 bit signed integers
   !
   ! Arguments: state: state of the random number generator
   !                   the counter in the state is incremented in every call
   !            res:   array to be filled with random 32 bit signed integers
   subroutine frand123Integer32_vector( state, res )
      use, intrinsic :: iso_c_binding, only: c_int64_t
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      integer( kind = c_int32_t ), dimension(:), intent( inout ) :: res

      integer( kind = c_int64_t ) :: len_res

      ! get length of res
      len_res = size( res )

      ! leave rest to C implementation
      call frand123Integer32_C( state, len_res, res )
   end subroutine frand123Integer32_vector

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!                  !!!!! 
   !!!!!  Initialization  !!!!!
   !!!!!                  !!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! initialize the state as follows:
   ! counter: first  64 bits: first entry of seed
   !          second 64 bits: second entry of seed
   ! key:     first  64 bits: rank
   !          second 64 bits: threadID
   !
   ! Arguments: state:    storage to hold the state of the random number generator
   !                      the state is handed over in each call to the random
   !                      number generators
   !            rank:     rank of the process using the random number generator
   !                      allows MPI-parallel use of the random number generator
   !                      If not in MPI situation, choose freely
   !            threadID: ID of the thread using the random number generator
   !                      allows thread-parallel use of the random number generator
   !                      If not in threaded situation, choose freely
   !            seed:     Seed for the random number generator to allow for
   !                      different or same random numbers with each run
   subroutine frand123Init( state, rank, threadID, seed )
      implicit none
      integer( kind = frand123_state_kind ), dimension( frand123_state_size ), intent( inout ) :: state
      integer, intent( in ) :: rank
      integer, intent( in ) :: threadID
      integer( kind = frand123_state_kind ), dimension( 2 ), intent( in ), optional :: seed

      ! set counter if present
      if(present(seed)) then
         state( 1 ) = seed( 1 )
         state( 2 ) = seed( 2 )
      end if

      ! set key
      state( 3 ) = rank
      state( 4 ) = threadID
   end subroutine frand123Init
end module
