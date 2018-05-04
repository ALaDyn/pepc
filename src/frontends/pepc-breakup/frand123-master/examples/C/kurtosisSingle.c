/* 
 * Example illustrating the use of the frand123 C interface for the generation
 * of single precision normally distributed random numbers.
 * The example estimates the excess kurtosis of the underlying normal distribution
 * with expectation mu = 0.7 and variance sigma = 1.6
 */

#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <omp.h>
#include "../../wrapper/frand123.h"

// number of draws used to estimate the excess kurtosis
const int64_t numberOfSamples = INT64_C( 1000 ) * INT64_C( 1000 ) * INT64_C( 100 );

// expectation and variance of the normal distribution to draw from
const float mu = 0.7;
const float sigma = 1.6;

// estimate excess kurtosis using calls to frand123NormSingle_scalar drawing a single single
// precision sample at a time
float excessKurtosisScalar( int64_t *state )
{
   float excessKurtosis = 0.f;

   for( int64_t i = 0; i < numberOfSamples; i++ )
   {
      // draw random number
      float draw = frand123NormSingle_scalar( state, mu, sigma );
      // centralize
      draw = draw - mu;
      // sum up contributions
      excessKurtosis += draw * draw * draw * draw;
   }
   // apply remaining operations
   return excessKurtosis / ( (float)numberOfSamples * sigma * sigma * sigma * sigma ) - 3.f;
}

// estimate excess kurtosis using calls to frand123NormSingle drawing a two single
// precision samples at a time
float excessKurtosisTwo( int64_t *state )
{
   float excessKurtosis = 0.;

   for( int64_t i = 0; i < numberOfSamples; i += 2 )
   {
      // draw random numbers
      float draw[ 2 ];
      frand123NormSingle( state, mu, sigma, 2, draw );
      for( int64_t j = 0; j < 2; j++ )
      {
         // centralize
         draw[ j ] = draw[ j ] - mu;
         // sum up contributions
         excessKurtosis += draw[ j ] * draw[ j ] * draw[ j ] * draw[ j ];
      }
   }
   // apply remaining operations
   return excessKurtosis / ( (float)numberOfSamples * sigma * sigma * sigma * sigma ) - 3.f;
}

// estimate excess kurtosis using calls to frand123NormSingle drawing dour single
// precision samples at a time
float excessKurtosisFour( int64_t *state )
{
   float excessKurtosis = 0.;

   for( int64_t i = 0; i < numberOfSamples; i += 4 )
   {
      // draw random numbers
      float draw[ 4 ];
      frand123NormSingle( state, mu, sigma, 4, draw );
      for( int64_t j = 0; j < 4; j++ )
      {
         // centralize
         draw[ j ] = draw[ j ] - mu;
         // sum up contributions
         excessKurtosis += draw[ j ] * draw[ j ] * draw[ j ] * draw[ j ];
      }
   }
   // apply remaining operations
   return excessKurtosis / ( (float)numberOfSamples * sigma * sigma * sigma * sigma ) - 3.f;
}

// estimate excess kurtosis using calls to frand123NormSingle drawing a 1000 single
// precision samples at a time
float excessKurtosis1000( int64_t *state )
{
   float excessKurtosis = 0.f;
   float *draw = malloc( 1000 * sizeof( float ) );
   if( draw == NULL )
   {
      perror( "Malloc of draw failed" );
   }

   for( int64_t i = 0; i < numberOfSamples; i += 1000 )
   {
      // draw random numbers
      frand123NormSingle( state, mu, sigma, 1000, draw );
      for( int64_t j = 0; j < 1000; j++ )
      {
         // centralize
         float centralized = draw[ j ] - mu;
         // sum up contributions
         excessKurtosis += centralized * centralized * centralized * centralized;
      }
   }
   // apply remaining operations
   return excessKurtosis / ( (float)numberOfSamples * sigma * sigma * sigma * sigma ) - 3.f;
}

int main()
{
   double startTime, endTime;

   printf( "Approximating excess kurtosis using %" PRId64 " samples.\n", numberOfSamples );

   // initialize state (do not provide seed, instead use whats in the memory location already)
   int64_t state[ 4 ];
   frand123Init( state, 0, 0, NULL );

   // use piScalar
   startTime = omp_get_wtime();
   float varExcessKurtosisScalar = excessKurtosisScalar( state );
   endTime = omp_get_wtime();
   printf( "Result generating one random number at a time:  %e took %e seconds\n", varExcessKurtosisScalar, endTime - startTime );

//   // use piTwo
//   startTime = omp_get_wtime();
//   float varExcessKurtosisTwo = excessKurtosisTwo( state );
//   endTime = omp_get_wtime();
//   printf( "Result generating two random number at a time:  %e took %e seconds\n", varExcessKurtosisTwo, endTime - startTime );
//
//   // use piFour
//   startTime = omp_get_wtime();
//   float varExcessKurtosisFour = excessKurtosisFour( state );
//   endTime = omp_get_wtime();
//   printf( "Result generating four random number at a time: %e took %e seconds\n", varExcessKurtosisFour, endTime - startTime );
//
//   // use pi2000
//   startTime = omp_get_wtime();
//   float varExcessKurtosis1000 = excessKurtosis1000( state );
//   endTime = omp_get_wtime();
//   printf( "Result generating 1000 random number at a time: %e took %e seconds\n", varExcessKurtosis1000, endTime - startTime );

   return 0;
}
