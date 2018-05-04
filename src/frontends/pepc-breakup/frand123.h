#ifndef __FRAND123_H__
#define __FRAND123_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/*
 * generate random double precision numbers uniformly distributed in (0,1)
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
double frand123Double_scalar( int64_t *state );
// vectorial version
void frand123Double( int64_t *state, const int64_t lenRes, double *res );

/*
 * generate random single precision numbers uniformly distributed in (0,1)
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            lenRes: number of random variables to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
float frand123Single_scalar( int64_t *state );
// vectorial version
void frand123Single( int64_t *state, const int64_t lenRes, float *res );

/*
 * generate random double precision numbers normally distributed with expectation mu and variance sigma
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            mu:     expectation of the normal distribution
 *            sigma:  variance of the normal distribution
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
double frand123NormDouble_scalar( int64_t *state, const double mu, const double sigma );
// vectorial version
void frand123NormDouble( int64_t *state, const double mu, const double sigma, const int64_t lenRes, double *res );

/*
 * generate random single precision numbers normally distributed with expectation mu and variance sigma
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            mu:     expectation of the normal distribution
 *            sigma:  variance of the normal distribution
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
float frand123NormSingle_scalar( int64_t *state, const float mu, const float sigma );
//vectorial version
void frand123NormSingle( int64_t *state, const float mu, const float sigma, const int64_t lenRes, float *res );

/*
 * generate random 64-bit signed integers uniformly distributed over INT64_MIN,..,INT64_MAX
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
int64_t frand123Integer64_scalar( int64_t *state );
//vectorial version
void frand123Integer64( int64_t *state, const int64_t lenRes, int64_t *res );

/*
 * generate random 32-bit signed integers uniformly distributed over INT32_MIN,..,INT32_MAX
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
int32_t frand123Integer32_scalar( int64_t *state );
// vectorial version
void frand123Integer32( int64_t *state, const int64_t lenRes, int32_t *res );

/*
 * generate random 64-bit unsigned integers uniformly distributed over UINT64_MIN,..,UINT64_MAX
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
uint64_t frand123UnsignedInteger64_scalar( int64_t *state );
//vectorial version
void frand123UnsignedInteger64( int64_t *state, const int64_t lenRes, uint64_t *res );

/*
 * generate random 32-bit unsigned integers uniformly distributed over UINT32_MIN,..,UINT32_MAX
 *
 * Arguments: state:  state of the random number generator
 *                    the counter in the state is incremented appropriately
 *            lenRes: number of random variates to generate
 *            res:    array to be filled with random numbers
 *                    length of array: lenRes
 */
// scalar version returns one random number
uint32_t frand123UnsignedInteger32_scalar( int64_t *state );
// vectorial version
void frand123UnsignedInteger32( int64_t *state, const int64_t lenRes, uint32_t *res );

/*
 * initialize the state for the random number generators used (Threefry or ARS)
 * rank and threadID determine the stream of random numbers used (by determining the key used in Threefry or ARS)
 * seed determines the initial value of the counter and by this the position within the stream of random numbers
 *
 * Arguments: state:    memory for state of the random number generator to initialize
 *            rank:     rank of the program (MPI)/number of the image (PGAS) using this state for a random number generator
 *            threadID: id of the thread using this state with a random number generator (pthreads/OpenMP)
 *            seed:     
 */
void frand123Init( int64_t *state, const int64_t rank, const int64_t threadID, const int64_t *seed );

#ifdef __cplusplus
} // extern "C"
#endif

#endif
