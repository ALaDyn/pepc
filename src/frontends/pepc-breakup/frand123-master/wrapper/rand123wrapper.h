#include <stdint.h>

#if USE_ARS
   /*
    * Function ars2x64_u01 calculates two double precision random numbers
    * uniformly distributed in (0,1).
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 2 double precision reals
    */
   void ars2x64_u01( int64_t *state, double *res );

   /*
    * Function ars4x32_u01 calculates four single precision random numbers
    * uniformly distributed in [0,1].
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 4 single precision reals
    */
   void ars4x32_u01( int64_t *state, float *res );
   
   /*
    * Function ars2x64_int calculates two 64 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   adress to storage for 2 64 bit signed integers
    */
   void ars2x64_int( int64_t *state, int64_t *res );
   
   /*
    * Function ars4x32_int calculates four 32 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   adress to storage for 4 32 bit signed integers
    */
   void ars4x32_int( int64_t *state, int32_t *res );
#else
   /*
    * Function threefry2x64_u01 returns 2 double precision reals
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 2 double precision reals
    */
   void threefry2x64_u01( int64_t *state, double *res );

   /*
    * Function threefry4x32_u01 returns 4 single precision reals
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *                   res:     address to storage for 4 single precision reals
    */
   void threefry4x32_u01( int64_t *state, float *res );

   /*
    * Function threefry2x64_int returns 2 64 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:   address to storage for 2 64 bit signed integers
    */
   void threefry2x64_int( int64_t *state, int64_t *res );

   /*
    * Function threefry4x32 returns 4 32 bit signed integers
    *
    * Arguments: state: four elements holding
    *                   counter: first  128 bit
    *                   key:     second 128 bit
    *            res:     address to storage for 4 32 bit signed integers
    */
   void threefry4x32_int( int64_t *state, int32_t *res );
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
   void polar2x64( int64_t *state, const double mu, const double sigma, double *res );
#else
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
   void wichura2x64( int64_t *state, const double mu, const double sigma, double *res );
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
void polar4x32( int64_t *state, const float mu, const float sigma, float *res );

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
void polar4x32_two( int64_t *state, const float mu, const float sigma, float *res );
