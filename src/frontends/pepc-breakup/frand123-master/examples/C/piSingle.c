/*
 * Example illustrating the use of the frand123 C interface for the generation
 * of single precision uniformly distributed random numbers.
 * The example uses a naive Monte-Carlo approach to approximate pi.
 * 
 * It generates points within the unit square and computes the ratio
 * of points within the upper right quadrant of the unit circle.
 * This ratio approximates pi / 4
 */

#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <omp.h>
#include "../../wrapper/frand123.h"

// number of points for the approximation
const int64_t numberOfPoints = UINT64_C( 1000 ) * UINT64_C( 1000 ) * UINT64_C( 100 );

// compute pi using calls to frand123Single_scalar drawing a single single precision
// uniformly distributed random number at a time
float piScalar( int64_t *state )
{
   // counter for points inside unit circle
   int64_t inside = 0;

   for( int64_t i = 0; i < numberOfPoints; i++ )
   {
      // generate point in unit square
      float x = frand123Single_scalar( state );
      float y = frand123Single_scalar( state );
      // compute square of radius and see whether inside unit circle
      float r2 = x * x + y * y;
      if( r2 <= 1. )
      {
         // inside unit circle
         inside++;
      }
   }
   // pi = 4 * ratio between points inside unit circle and points generated in total
   return 4.f * (float)inside / (float)numberOfPoints;
}

// compute pi using calls to frand123Single drawing both single precision coordinates
// a time from the uniform distribution
float piTwo( int64_t *state )
{
   // counter for points inside unit circle
   int64_t inside = 0;

   for( int64_t i = 0; i < numberOfPoints; i++ )
   {
      // generating point in unit square
      float pos[ 2 ];
      frand123Single( state, 2, pos );
      // compute square of radius and see whether inside unit circle
      float r2 = pos[ 0 ] * pos[ 0 ] + pos[ 1 ] * pos[ 1 ];
      if( r2 <= 1.f )
      {
         // inside unit circle
         inside++;
      }
   }
   // pi = 4 * ratio between points inside unit circle and points generated in total
   return 4.f * (float)inside / (float)numberOfPoints;
}

// compute pi using calls to frand123Single drawing four single precision coordinates
// a time from the uniform distribution
float piFour( int64_t *state )
{
   // counter for points inside unit circle
   int64_t inside = 0;

   for( int64_t i = 0; i < numberOfPoints; i += 2 )
   {
      // generating point in unit square
      float pos[ 4 ];
      frand123Single( state, 4, pos );
      // compute square of radius and see whether inside unit circle
      float r2 = pos[ 0 ] * pos[ 0 ] + pos[ 1 ] * pos[ 1 ];
      if( r2 <= 1.f )
      {
         // inside unit circle
         inside++;
      }
      r2 = pos[ 2 ] * pos[ 2 ] + pos[ 3 ] * pos[ 3 ];
      if( r2 <= 1.f )
      {
         // inside unit circle
         inside++;
      }
   }
   // pi = 4 * ratio between points inside unit circle and points generated in total
   return 4.f * (float)inside / (float)numberOfPoints;
}

// compute pi using calls to frand123Single drawing 2000 single precision coordinates
// a time from the uniform distribution
float pi2000( int64_t *state )
{
   // counter for points inside unit circle
   int64_t inside = 0;
   // buffer for random numbers
   float *buffer;
   buffer = (float*)malloc( 2000 * sizeof( float ) );
   if( buffer == NULL )
   {
      perror( "Malloc of buffer failed" );
   }

   for( int64_t i = 0; i < numberOfPoints; i += 1000 )
   {
      // generating  buffer of random numbers
      frand123Single( state, 2000, buffer );

      // iterate over the associated 1000 points
      for( int j = 0; j < 1000; j++ )
      {
         // ease access to coordinates
         float x = buffer[ UINT64_C( 2 ) * j ];
         float y = buffer[ UINT64_C( 2 ) * j + UINT64_C( 1 ) ];
         // compute square of radius and see whether inside unit circle
         float r2 = x * x + y * y;
         if( r2 <= 1.f )
         {
            // inside unit circle
            inside++;
         }
      }
   }
   // pi = 4 * ratio between points inside unit circle and points generated in total
   return 4.f * (float)inside / (float)numberOfPoints;
}

int main()
{
   double startTime, endTime;

   printf( "Approximating pi using %" PRId64 " points in the unit square.\n", numberOfPoints );

   // initialize state (do not provide seed, instead use whats in the memory location already)
   int64_t state[ 4 ];
   frand123Init( state, 0, 0, NULL );

   // use piScalar
   startTime = omp_get_wtime();
   float varPiScalar = piScalar( state );
   endTime = omp_get_wtime();
   printf( "Result generating one random number at a time:  %e took %e seconds\n", varPiScalar, endTime - startTime );

   // use piTwo
   startTime = omp_get_wtime();
   float varPiTwo = piTwo( state );
   endTime = omp_get_wtime();
   printf( "Result generating two random number at a time:  %e took %e seconds\n", varPiTwo, endTime - startTime );

   // use piFour
   startTime = omp_get_wtime();
   float varPiFour = piFour( state );
   endTime = omp_get_wtime();
   printf( "Result generating four random number at a time: %e took %e seconds\n", varPiFour, endTime - startTime );

   // use pi2000
   startTime = omp_get_wtime();
   float varPi2000 = pi2000( state );
   endTime = omp_get_wtime();
   printf( "Result generating 2000 random number at a time: %e took %e seconds\n", varPi2000, endTime - startTime );

   return 0;
}

