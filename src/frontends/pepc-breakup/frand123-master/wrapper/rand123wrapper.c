#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include "rand123wrapper.h"

#if USE_ARS
   #include <ars.h>

   /*
    * Constants used to first map onto interval [0,1) and then onto (0,1)
    */
   static const double factor_double  = 1.  / ( (double)UINT64_MAX + 1. );
   static const double summand_double = 0.5 / ( (double)UINT64_MAX + 1. );

   /*
    * Function ars2x64_u01 calculates two double precision random numbers
    * uniformly distributed in (0,1).
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 2 double precision reals
    */
   void ars2x64_u01( int64_t *state, double *res )
   {
      // extract counter and key from state
      ars4x32_ctr_t *ctr_ars = (ars4x32_ctr_t*)&state[0];
      ars4x32_key_t *key_ars = (ars4x32_key_t*)&state[2];
      // calc uniformly distributed integers
      ars4x32_ctr_t resArs = ars4x32_R( 6, *ctr_ars, *key_ars );
      // convert to uint64_t
      uint64_t resInts[2];
      resInts[0] = *((uint64_t*)&resArs.v[0]);
      resInts[1] = *((uint64_t*)&resArs.v[2]);
      // convert to uniformly distributed doubles in (0,1)
      res[0] = (double)resInts[0] * factor_double + summand_double;
      res[1] = (double)resInts[1] * factor_double + summand_double;
      // advance counter
      if( ctr_ars->v[0] < UINT32_MAX )
         ctr_ars->v[0]++;
      else
      {
         ctr_ars->v[0] = 0;
         if( ctr_ars->v[1] < UINT32_MAX )
            ctr_ars->v[1]++;
         else
         {
            ctr_ars->v[1] = 0;
            if( ctr_ars->v[2] < UINT32_MAX )
               ctr_ars->v[2]++;
            else
            {
               ctr_ars->v[2] = 0;
               ctr_ars->v[3]++;
            }
         }
      }
      return;
   }

   /*
    * Constants used to first map onto (-0.5,0.5) and then (0,1)
    *
    * *******************
    * **** IMPORTANT ****
    * *******************
    * The value of enlarger was determined experimentally to ensure mapping onto (0,-1).
    * It is included from frand123enlarger.h to allow validation
    * of the choice in testAccuracyFloats
    */
#include "frand123enlarger.h"
   static const float  factor_float   = 1.f / ( (float)UINT32_MAX + 257.f );
   static const float  summand_float  = 0.5f;

   /*
    * Function ars4x32_u01 calculates four single precision random numbers
    * uniformly distributed in [0,1].
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 4 single precision reals
    */
   void ars4x32_u01( int64_t *state, float *res )
   {
      // extract counter and key from state
      ars4x32_ctr_t *ctr_ars = (ars4x32_ctr_t*)&state[0];
      ars4x32_key_t *key_ars = (ars4x32_key_t*)&state[2];
      // move to vector registers
      ars1xm128i_ctr_t c128;
      ars1xm128i_key_t k128;
      c128.v[0].m = _mm_set_epi32( ctr_ars->v[3], ctr_ars->v[2], ctr_ars->v[1], ctr_ars->v[0] );
      k128.v[0].m = _mm_set_epi32( key_ars->v[3], key_ars->v[2], key_ars->v[1], key_ars->v[0] );
      // calc uniformly distributed integers in SEE vector registers
      c128 = ars1xm128i_R( 6, c128, k128 );
      // convert signed integers to signed floats
      __m128 asSignedFloats = _mm_cvtepi32_ps( c128.v[0].m );
      // normalize to [-0.5,0.5] and add shift to end in [0,1]
      __m128 normFactor = _mm_load_ps1( &factor_float );
      __m128 summand    = _mm_load_ps1( &summand_float );
#ifdef USE_FMA
      __m128 restrictedToUnitInterval = _mm_fmadd_ps( asSignedFloats, normFactor, summand );
#else
      __m128 restrictedToUnitCircle = _mm_mul_ps( asSignedFloats, normFactor );
      __m128 restrictedToUnitInterval = _mm_add_ps( restrictedToUnitCircle, summand );
#endif
      // store result in memory
      _mm_storeu_ps( res, restrictedToUnitInterval );
      // advance counter
      if( ctr_ars->v[0] < UINT32_MAX )
         ctr_ars->v[0]++;
      else
      {
         ctr_ars->v[0] = 0;
         if( ctr_ars->v[1] < UINT32_MAX )
            ctr_ars->v[1]++;
         else
         {
            ctr_ars->v[1] = 0;
            if( ctr_ars->v[2] < UINT32_MAX )
               ctr_ars->v[2]++;
            else
            {
               ctr_ars->v[2] = 0;
               ctr_ars->v[3]++;
            }
         }
      }
      return;
   }
   
   /*
    * Function ars2x64_int calculates two 64 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   adress to storage for 2 64 bit signed integers
    */
   void ars2x64_int( int64_t *state, int64_t *res )
   {
      // extract counter and key from state
      ars4x32_ctr_t *ctr_ars = (ars4x32_ctr_t*)&state[0];
      ars4x32_key_t *key_ars = (ars4x32_key_t*)&state[2];
      // calc uniformly distributed integers
      ars4x32_ctr_t resArs = ars4x32_R( 6, *ctr_ars, *key_ars );
      // store in res and reinterpret as int64_t
      res[0] = *((int64_t*)&resArs.v[0]);
      res[1] = *((int64_t*)&resArs.v[2]);
      // advance counter
      if( ctr_ars->v[0] < UINT32_MAX )
         ctr_ars->v[0]++;
      else
      {
         ctr_ars->v[0] = 0;
         if( ctr_ars->v[1] < UINT32_MAX )
            ctr_ars->v[1]++;
         else
         {
            ctr_ars->v[1] = 0;
            if( ctr_ars->v[2] < UINT32_MAX )
               ctr_ars->v[2]++;
            else
            {
               ctr_ars->v[2] = 0;
               ctr_ars->v[3]++;
            }
         }
      }
      return;
   }
   
   /*
    * Function ars4x32_int calculates four 32 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   adress to storage for 4 32 bit signed integers
    */
   void ars4x32_int( int64_t *state, int32_t *res )
   {
      // extract counter and key from state
      ars4x32_ctr_t *ctr_ars = (ars4x32_ctr_t*)&state[0];
      ars4x32_key_t *key_ars = (ars4x32_key_t*)&state[2];
      // calc uniformly distributed integers
      ars4x32_ctr_t resArs = ars4x32_R( 6, *ctr_ars, *key_ars );
      // store in res and reinterpret as int64_t
      res[0] = *((int32_t*)&resArs.v[0]);
      res[1] = *((int32_t*)&resArs.v[1]);
      res[2] = *((int32_t*)&resArs.v[2]);
      res[3] = *((int32_t*)&resArs.v[3]);
      // advance counter
      if( ctr_ars->v[0] < UINT32_MAX )
         ctr_ars->v[0]++;
      else
      {
         ctr_ars->v[0] = 0;
         if( ctr_ars->v[1] < UINT32_MAX )
            ctr_ars->v[1]++;
         else
         {
            ctr_ars->v[1] = 0;
            if( ctr_ars->v[2] < UINT32_MAX )
               ctr_ars->v[2]++;
            else
            {
               ctr_ars->v[2] = 0;
               ctr_ars->v[3]++;
            }
         }
      }
      return;
   }
#else
   #include <threefry.h>

   /*
    * Union conv_union_t simplifies access to bits generated using threefry2x64
    */
   typedef union
   {
      threefry2x64_ctr_t ctr;
      uint64_t ints[2];
   } conv_union_t;

   /*
    * Constants used to first map onto interval [0,1) and then onto (0,1)
    */
   static const double factor_double  = 1.  / ( (double)UINT64_MAX + 1. );
   static const double summand_double = 0.5 / ( (double)UINT64_MAX + 1. );

   /*
    * Function threefry2x64_u01 returns 2 double precision reals
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 2 double precision reals
    */
   void threefry2x64_u01( int64_t *state, double *res )
   {
      // extract counter and key from state
      threefry2x64_ctr_t *ctr_threefry = (threefry2x64_ctr_t*)&state[0];
      threefry2x64_key_t *key_threefry = (threefry2x64_key_t*)&state[2];
      // calc uniformly distributed integers
      conv_union_t resInt;
      resInt.ctr = threefry2x64_R( 13, *ctr_threefry, *key_threefry );
      // convert to uniformly distributed doubles in (0,1)
      res[0] = (double)resInt.ints[0] * factor_double + summand_double;
      res[1] = (double)resInt.ints[1] * factor_double + summand_double;
      // advance counter
      if( ctr_threefry->v[0] < UINT64_MAX )
         ctr_threefry->v[0]++;
      else
      {
         ctr_threefry->v[0] = 0;
         ctr_threefry->v[1]++;
      }
      return;
   }

   /*
    * Constants used to first map onto interval [0,1) and then onto (0,1)
    */
   static const float factor_float   = 1.  / ( (float)UINT32_MAX + 400.f );
   static const float summand_float  = 0.5 / ( (float)UINT32_MAX + 400.f );

   /*
    * Function threefry4x32_u01 returns 4 single precision reals
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *                   res:     address to storage for 4 single precision reals
    */
   void threefry4x32_u01( int64_t *state, float *res )
   {
      // extract counter and key from state
      threefry4x32_ctr_t *ctr_threefry = (threefry4x32_ctr_t*)&state[0];
      threefry4x32_key_t *key_threefry = (threefry4x32_key_t*)&state[2];
      // calc uniformly distributed integers
      threefry4x32_ctr_t resInt = threefry4x32_R( 12, *ctr_threefry, *key_threefry );
      // convert to uniformly distributed floats in (0,1)
      res[0] = (float)resInt.v[0] * factor_float + summand_float;
      res[1] = (float)resInt.v[1] * factor_float + summand_float;
      res[2] = (float)resInt.v[2] * factor_float + summand_float;
      res[3] = (float)resInt.v[3] * factor_float + summand_float;
      // advance counter
      if( ctr_threefry->v[0] < UINT32_MAX )
         ctr_threefry->v[0]++;
      else
      {
         ctr_threefry->v[0] = 0;
         if( ctr_threefry->v[1] < UINT32_MAX )
            ctr_threefry->v[1]++;
         else
         {
            ctr_threefry->v[1] = 0;
            if( ctr_threefry->v[2] < UINT32_MAX )
               ctr_threefry->v[2]++;
            else
            {
               ctr_threefry->v[2] = 0;
               ctr_threefry->v[3]++;
            }
         }
      }
      return;
   }

   /*
    * Function threefry2x64_int returns 2 64 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 2 64 bit signed integers
    */
   void threefry2x64_int( int64_t *state, int64_t *res )
   {
      // extract counter and key from state
      threefry2x64_ctr_t *ctr_threefry = (threefry2x64_ctr_t*)&state[0];
      threefry2x64_key_t *key_threefry = (threefry2x64_key_t*)&state[2];
      // calc uniformly distributed integers
      threefry2x64_ctr_t resThreefry = threefry2x64_R( 13, *ctr_threefry, *key_threefry );
      // reinterprete as signed 64 bit integer
      res[0] = *((int64_t*)&resThreefry.v[0]);
      res[1] = *((int64_t*)&resThreefry.v[1]);
      // advance counter
      if( ctr_threefry->v[0] < UINT64_MAX )
         ctr_threefry->v[0]++;
      else
      {
         ctr_threefry->v[0] = 0;
         ctr_threefry->v[1]++;
      }
      return;
   }

   /*
    * Function threefry4x32 returns 4 32 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:     address to storage for 4 32 bit signed integers
    */
   void threefry4x32_int( int64_t *state, int32_t *res )
   {
      // extract counter and key from state
      threefry4x32_ctr_t *ctr_threefry = (threefry4x32_ctr_t*)&state[0];
      threefry4x32_key_t *key_threefry = (threefry4x32_key_t*)&state[2];
      // calc uniformly distributed integers
      threefry4x32_ctr_t resInt = threefry4x32_R( 12, *ctr_threefry, *key_threefry );
      // reinterprete as 32 bit signed integer
      res[0] = *((int32_t*)&resInt.v[0]);
      res[1] = *((int32_t*)&resInt.v[1]);
      res[2] = *((int32_t*)&resInt.v[2]);
      res[3] = *((int32_t*)&resInt.v[3]);
      // advance counter
      if( ctr_threefry->v[0] < UINT32_MAX )
         ctr_threefry->v[0]++;
      else
      {
         ctr_threefry->v[0] = 0;
         if( ctr_threefry->v[1] < UINT32_MAX )
            ctr_threefry->v[1]++;
         else
         {
            ctr_threefry->v[1] = 0;
            if( ctr_threefry->v[2] < UINT32_MAX )
               ctr_threefry->v[2]++;
            else
            {
               ctr_threefry->v[2] = 0;
               ctr_threefry->v[3]++;
            }
         }
      }
   }

#endif

#ifdef USE_POLAR
   /*
    * Function polar2x64 calculates two double precision random numbers
    * normally distributed with expectation mu and variance sigma using the
    * polar rejection method by Box and Muller
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            mu:    expectation
    *            sigma: variance
    *            res:   address to storage for 2 double precision reals
    */
   void polar2x64( int64_t *state, const double mu, const double sigma, double *res )
   {
      double u[ 2 ];
      double x[ 2 ];
      double r2;
      double f;
      // generate coordinates until within unit circle
      // at least first try successful: probability ~79%
      // at least second try successful: ~95%
      // at least third try successful: ~99%
      do
      {
#if USE_ARS
         ars2x64_u01( state, u );
#else
         threefry2x64_u01( state, u );
#endif
         x[ 0 ] = 2. * u[ 0 ] - 1.;
         x[ 1 ] = 2. * u[ 1 ] - 1.;
         r2 = x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ];
      } while( ( r2 >= 1. ) || ( r2 == 0. ) );
      // compute random numbers
      f = sqrt( -2. * log( r2 ) / r2 );
      res[ 0 ] = mu + sigma * f * x[ 0 ];
      res[ 1 ] = mu + sigma * f * x[ 1 ];
      return;
   }
#else
   /*
    * Function wichura2x64kernel represents the kernel of the AS 241 algorithm extracted for testing
    */
   void wichura2x64kernel( const double *p, const double mu, const double sigma, double *res )
   {
      // constants for the polynomials
      const double A[ 8 ] = 
      { 
         3.3871328727963666080e0,
         1.3314166789178437745e2,
         1.9715909503065514427e3,
         1.3731693765509461125e4,
         4.5921953931549871457e4,
         6.7265770927008700853e4,
         3.3430575583588128105e4,
         2.5090809287301226727e3
      };
      const double B[ 8 ] =
      {
         1.0E0,
         4.2313330701600911252e1,
         6.8718700749205790830e2,
         5.3941960214247511077e3,
         2.1213794301586595867e4,
         3.9307895800092710610e4,
         2.8729085735721942674e4,
         5.2264952788528545610e3
      };
      const double C[ 8 ] =
      {
         1.42343711074968357734e0,
         4.63033784615654529590e0,
         5.76949722146069140550e0,
         3.64784832476320460504e0,
         1.27045825245236838258e0,
         2.41780725177450611770e-1,
         2.27238449892691845833e-2,
         7.74545014278341407640e-4
      };
      const double D[ 8 ] = 
      {
         1.0E0,
         2.05319162663775882187e0,
         1.67638483018380384940e0,
         6.89767334985100004550e-1,
         1.48103976427480074590e-1,
         1.51986665636164571966e-2,
         5.47593808499534494600e-4,
         1.05075007164441684324e-9
      };
      const double E[ 8 ] =
      {
         6.65790464350110377720e0,
         5.46378491116411436990e0,
         1.78482653991729133580e0,
         2.96560571828504891230e-1,
         2.65321895265761230930e-2,
         1.24266094738807843860e-3,
         2.71155556874348757815e-5,
         2.01033439929228813265e-7
      };
      const double F[ 8 ] =
      {
         1.0E0,
         5.99832206555887937690e-1,
         1.36929880922735805310e-1,
         1.48753612908506148525e-2,
         7.86869131145613259100e-4,
         1.84631831751005468180e-5,
         1.42151175831644588870e-7,
         2.04426310338993978564e-15
      };
      // algorithmic constants
      const double split1 = 0.425;
      const double split2 = 5.;
      const double const1 = 0.180625;
      const double const2 = 1.6;
      // variables
      double q[ 2 ], r, prodA[ 2 ], prodB[ 2];
      int i, j;
      bool cond[ 2 ];
      // compute both q at once
      q[ 0 ] = p[ 0 ] - 0.5;
      q[ 1 ] = p[ 1 ] - 0.5;
      // compute abs eval
      cond[ 0 ] = ( fabs( q[ 0 ] ) <= split1 );
      cond[ 1 ] = ( fabs( q[ 1 ] ) <= split1 );
      // run vectorized version of the main route
      // in the following cases
      // both values of |q| <= split1: probability ~72%
      // one |q| <= split1, one |q| > split1: probability ~26%
      if( cond[ 0 ] || cond[ 1 ] )
      {
         #pragma omp simd private( r )
         for( i = 0; i <= 1; i++ )
         {
            r = const1 - q[ i ] * q[ i ];
            prodA[ i ] = A[7];
            prodB[ i ] = B[7];
            for( j = 6; j >= 0; j-- )
            {
               prodA[ i ] = prodA[ i ] * r + A[ j ];
               prodB[ i ] = prodB[ i ] * r + B[ j ];
            }
            // transform to required mean and variance
            res[ i ] = mu + sigma * q[ i ] * prodA[ i ] / prodB[ i ];
         }
      }
      // return in 72.25% of calls
      if( cond[ 0 ] && cond[ 1 ] )
      {
         return;
      }
      // handle first non-standard case
      // probability for need to handle single case: ~13%
      // probability for need to handle both cases: ~2%
      for( i = 0; i <= 1; i++ )
      {
         if( ! cond[ i ] )
         {
            if( q[ i ] < 0 )
            {
               r = p[ i ];
            }
            else
            {
               r = 1. - p[ i ];
            }
            r = sqrt( -log( r ) );
            // not too far within the tail
            if( r <= split2 )
            {
               r = r - const2;
               double prodC = C[ 7 ];
               double prodD = D[ 7 ];
               for( j = 6; j >= 0; j-- )
               {
                  prodC = prodC * r + C[ j ];
                  prodD = prodD * r + D[ j ];
               }
               res[ i ] = prodC / prodD;
            }
            // far within the tail
            else
            {
               r = r - split2;
               double prodE = E[ 7 ];
               double prodF = F[ 7 ];
               for( j = 6; j >= 0; j-- )
               {
                  prodE = prodE * r + E[ j ];
                  prodF = prodF * r + F[ j ];
               }
               res[ i ] = prodE / prodF;
            }
            if( q[ i ] < 0. )
            {
               res[ i ] = -res[ i ];
            }
            // transform to required mean and variance
            res[ i ] = mu + sigma * res[ i ];
         }
      }
   }

   /*
    * Function wichura2x64 calculates two double precision random numbers
    * normally distributed with expectation mu and variance sigma using the
    * inverse transform sampling introduced by Wichura 1988 in algorithm AS241 (PPND16)
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            mu:    expectation
    *            sigma: variance
    *            res:   address to storage for 2 double precision reals
    */
   void wichura2x64( int64_t *state, const double mu, const double sigma, double *res )
   {
      // get uniform random numbers
      double p[ 2 ];
#ifdef USE_ARS
      ars2x64_u01( state, p );
#else
      threefry2x64_u01( state, p );
#endif
      // rest is computed within the kernel
      wichura2x64kernel( p, mu, sigma, res );
      return;
   }
#endif

/*
 * Function polar4x32 calculates four single precision random numbers
 * normally distributed with expectation mu and variance sigma using the
 * polar rejection method by Box and Muller
 *
 * Arguments: state: four elements holding
 *                   counter: first  128 bit
 *                   key:     second 128 bit
 *            mu:    expectation
 *            sigma: variance
 *            res:   address to storage for 4 single precision reals
 */
void polar4x32( int64_t *state, const float mu, const float sigma, float *res )
{
   float u[ 4 ];
   float x[ 4 ];
   float r2[ 2 ];
   float f[ 2 ];
   // generate coordinates until within unit circle
   // at least first try successful: probability ~62%
   // at least second try successful: probability ~85%
   // at least third try successful: probability ~94%
   // this implementation should require ~0.8 calls of ars4x32_u01 or
   // threefry4x32_u01 to generate a single random number but fewer
   // conditionals than one that requires only ~0.7 calls per random number
   do
   {
#if USE_ARS
      ars4x32_u01( state, u );
#else
      threefry4x32_u01( state, u );
#endif
      x[ 0 ] = 2.f * u[ 0 ] - 1.f;
      x[ 1 ] = 2.f * u[ 1 ] - 1.f;
      x[ 2 ] = 2.f * u[ 2 ] - 1.f;
      x[ 3 ] = 2.f * u[ 3 ] - 1.f;
      r2[ 0 ] = x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ];
      r2[ 1 ] = x[ 2 ] * x[ 2 ] + x[ 3 ] * x[ 3 ];
   } while( ( r2[ 0 ] >= 1.f ) || ( r2[ 0 ] == 0.f ) || ( r2[ 1 ] >= 1.f ) || ( r2[ 1 ] == 0.f ) );
   // compute random numbers
   f[ 0 ] = sqrt( -2.f * log( r2[ 0 ] ) / r2[ 0 ] );
   f[ 1 ] = sqrt( -2.f * log( r2[ 1 ] ) / r2[ 1 ] );
   res[ 0 ] = mu + sigma * f[ 0 ] * x[ 0 ];
   res[ 1 ] = mu + sigma * f[ 0 ] * x[ 1 ];
   res[ 2 ] = mu + sigma * f[ 1 ] * x[ 2 ];
   res[ 3 ] = mu + sigma * f[ 1 ] * x[ 3 ];
   return;
}

/*
 * Function polar4x32_two calculates two single precision random numbers
 * normally distributed with expectation mu and variance sigma using the
 * polar rejection method by Box and Muller
 *
 * Note: restriction to 2 random numbers allows for fewer calls to RNG due
 * the specifics of the rejection method
 *
 * Arguments: state: four elements holding
 *                   counter: first  128 bit
 *                   key:     second 128 bit
 *            mu:    expectation
 *            sigma: variance
 *            res:   address to storage for 2 single precision reals
 */
void polar4x32_two( int64_t *state, const float mu, const float sigma, float *res )
{
   float u[ 4 ];
   float x[ 4 ];
   float r2[ 2 ];
   float f;
   bool in[ 2 ];
   // generate coordinates until within unit circle
   // at least first try successful: probability ~95%
   // at least second try successful: probability ~100%
   // this implementation requires ~ 0.52 calls of the
   // RNG per random number
   do
   {
#if USE_ARS
      ars4x32_u01( state, u );
#else
      threefry4x32_u01( state, u );
#endif
      // compute distance from center
      x[ 0 ] = 2.f * u[ 0 ] - 1.f;
      x[ 1 ] = 2.f * u[ 1 ] - 1.f;
      x[ 2 ] = 2.f * u[ 2 ] - 1.f;
      x[ 3 ] = 2.f * u[ 3 ] - 1.f;
      r2[ 0 ] = x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ];
      r2[ 1 ] = x[ 2 ] * x[ 2 ] + x[ 3 ] * x[ 3 ];
      // check which lies within the unit circle
      in[ 0 ] = ( r2[ 0 ] > 0.f ) && ( r2[ 0 ] < 1.f );
      in[ 1 ] = ( r2[ 1 ] > 0.f ) && ( r2[ 1 ] < 1.f );
   } while( ! ( in[ 0 ] || in[ 1 ] ) );
   // compute two random numbers
   if( in[ 0 ] )
   {
      f = sqrt( -2.f * log( r2[ 0 ] ) / r2[ 0 ] );
      res[ 0 ] = mu + sigma * f * x[ 0 ];
      res[ 1 ] = mu + sigma * f * x[ 1 ];
   }
   else
   {
      f = sqrt( -2.f * log( r2[ 1 ] ) / r2[ 1 ] );
      res[ 0 ] = mu + sigma * f * x[ 2 ];
      res[ 1 ] = mu + sigma * f * x[ 3 ];
   }
   return;
}
