# The frand123 Fortran and C wrapper for the low-level Random123 C library

The frand123 Fortran wrapper for the Random123 C library provides random number generators (RNG) to Fortran using the random bits generated using functions of the Random123 C library.

The Random123 C library (see [publication](http://dx.doi.org/10.1145/2063384.2063405) and [code](https://www.deshawresearch.com/resources_random123.html)) offer the following features:
* Perfect parallelization due to counter-based design with minimal state
* At least 2^64 unique streams with a period of at least 2^128
* Resistent to the BigCrush randomness test in TestU01 (see [publication](http://dx.doi.org/10.1145/1268776.1268777) and [code](http://simul.iro.umontreal.ca/testu01/tu01.html))
* High-performance ARS RNG implemented using Intel AES-NI instructions
* High-performance Threefry RNG implemented in plain C
* Permissive 3-clause BSD-type license (see file "License")

The frand123 Fortran wrapper provides the following capabilities:
* Perfect parallelization due to counter-based design with minimal state
* Generation of uniformly distributed:
    * signed 32-bit integers
    * signed 64-bit integers
    * single precision reals in (0,1)
    * double precision reals in (0,1)
* Comfortable interface for generating arbitrarily sized vectors of random numbers
* Interface for initializing the state of the RNG

## Reference Fortran

### frand123Single( state, res ) / frand123Rand( state, res )
#### Description
* Fill vector res with single precision real random numbers uniformly distributed in (0,1)
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __res__: memory to which random numbers are stored to
    * dimension: scalar or arbitrary length array
    * kind: real of kind _c_float_
    * intent: _inout_

### frand123Double( state, res ) / frand123Rand( state, res )
#### Description
* Fill vector res with double precision real random numbers uniformly distributed in (0,1)
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __res__: memory to which random numbers are stored to
    * dimension: scalar or arbitrary length array
    * kind: real of kind _c_double_
    * intent: _inout_

### frand123NormDouble( state, mu, sigma, res ) / frand123Rand( state, mu, sigma, res )
#### Description
* Fill vector res with double precision real random numbers normally distributed with mean mu and variance sigma
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __mu__: mean of the normal distribution
    * dimension: scalar
    * kind: real of kind _c_double_
    * intent: in
* __sigma__: variance of the normal distribution
    * dimension: scalar
    * kind: real of kind _c_double_
    * intent: in
* __res__: memory to which random numbers are stored to
    * dimension: scalar or arbitrary length array
    * kind: real of kind _c_double_
    * intent: _inout_

### frand123NormSingle( state, mu, sigma, res ) / frand123Rand( state, mu, sigma, res )
#### Description
* Fill vector res with single precision real random numbers normally distributed with mean mu and variance sigma
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __mu__: mean of the normal distribution
    * dimension: scalar
    * kind: real of kind _c_float_
    * intent: in
* __sigma__: variance of the normal distribution
    * dimension: scalar
    * kind: real of kind _c_float_
    * intent: in
* __res__: memory to which random numbers are stored to
    * dimension: scalar or arbitrary length array
    * kind: real of kind _c_float_
    * intent: _inout_

### frand123Integer32( state, res ) / frand123Rand( state, res )
#### Description
* Fill vector res with signed 32-bit integer random numbers uniformly distributed between (iclusive) INT32_MIN and INT32_MAX
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __res__: memory to which random numbers are stored to
    * dimension: scalar or arbitrary length array
    * kind: real of kind _c_int32_t
    * intent: _inout_

### frand123Integer64( state, res ) / frand123Rand( state, res )
#### Description
* Fill vector res with signed 64-bit integer random numbers uniformly distributed between (iclusive) INT64_MIN and INT64_MAX
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __res__: memory to which random numbers are stored to
    * dimension: scalar or arbitrary length array
    * kind: real of kind _c_int64_t
    * intent: _inout_

### frand123Init( state, rank, threadID, seed )
#### Description
* Initialize state for use in serial, MPI-parallel, thread-parallel, and MPI- and thread-parallel settings
* The key is initialized as follows:
    * first  64 bits: rank
    * second 64 bits: threadID
* The counter is initialized as follows:
    * first  64 bits: first element of seed
    * second 64 bits: second element of seed

#### Arguments
* __state__: memory in which to initialize state in
    * dimension: _frand123_state_size_
    * kind: integer of kind _frand123_state_kind_
    * intent: _inout_
* __rank__: MPI rank of the caller
    * dimension: scalar
    * kind: integer
    * intent: _in_
* __threadID__: thread ID of the thread using the RNG
    * dimension: scalar
    * kind: integer
    * intent: _in_
* __seed__: seed to be used to initialize the counter
    * dimension: 2
    * kind: integer of kind _frand123_state_kind_
    * intent: _in_

## Reference C

### void frand123Single( int64_t *state, const long long lenRes, float *res )
#### Description
* Fill vector res with single precision real random numbers uniformly distributed in (0,1)
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### float frand123Single_scalar( int64_t *state )
#### Description
* Return single random single precision number uniformly distributed in (0,1)
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG

### void frand123Double( int64_t *state, const long long lenRes, double *res )
#### Description
* Fill vector res with double precision real random numbers uniformly distributed in (0,1)
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### double frand123Double_scalar( int64_t *state )
#### Description
* Return single random double precision number uniformly distributed in (0,1)
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG

### void frand123NormSingle( int64_t *state, const float mu, const float sigma, const long long lenRes, float *res )
#### Description
* Fill vector res with single precision real random numbers normally distributed with expectation mu and variance sigma
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
* __mu__: expectation of the normal distribution to draw from
* __sigma__: variance of the normal distribution to draw from
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### float frand123NormSingle_scalar( int64_t *state, const float mu, const float sigma )
#### Description
* Return single random single precision number normally distributed with expectation mu and variance sigma
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG
* __mu__: expectation of the normal distribution to draw from
* __sigma__: variance of the normal distribution to draw from

### void frand123NormDouble( int64_t *state, const double mu, const double sigma, const long long lenRes, double *res )
#### Description
* Fill vector res with double precision real random numbers normally distributed with expectation mu and variance sigma
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
* __mu__: expectation of the normal distribution to draw from
* __sigma__: variance of the normal distribution to draw from
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### double frand123NormDouble_scalar( int64_t *state, const double mu, const double sigma )
#### Description
* Return single random double precision number normally distributed with expectation mu and variance sigma
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG
* __mu__: expectation of the normal distribution to draw from
* __sigma__: variance of the normal distribution to draw from

### void frand123Integer32( int64_t *state, const long long lenRes, int32_t *res )
#### Description
* Fill vector res with 32-bit integer random numbers discretely uniformly distributed on INT32_MIN,..,INT32_MAX
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### int32_t frand123Integer32_scalar( int64_t *state )
#### Description
* Return single random 32-bit integer number discretely uniformly distributed on INT32_MIN,..,INT32_MAX
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG

### void frand123Integer64( int64_t *state, const long long lenRes, int64_t *res )
#### Description
* Fill vector res with 64-bit integer random numbers discretely uniformly distributed on INT64_MIN,..,INT64_MAX
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### int64_t frand123Integer64_scalar( int64_t *state )
#### Description
* Return single random 64-bit integer number discretely uniformly distributed on INT64_MIN,..,INT64_MAX
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

### void frand123UnsignedInteger32( int64_t *state, const long long lenRes, uint32_t *res )
#### Description
* Fill vector res with 32-bit unsigned integer random numbers discretely uniformly distributed on 0,..,2^32-1
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 4

#### Arguments
* __state__: state used by the RNG
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### uint32_t frand123UnsignedInteger32_scalar( int64_t *state )
#### Description
* Return single random 32-bit unsigned integer number discretely uniformly distributed on 0,..,2^32-1
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG

### void frand123UnsignedInteger64( int64_t *state, const long long lenRes, uint64_t *res )
#### Description
* Fill vector res with 64-bit unsigned integer random numbers discretely uniformly distributed on 0,..,2^64-1
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits
* Random numbers are always generated in chunks of 2

#### Arguments
* __state__: state used by the RNG
* __lenRes__: length of array res
* __res__: memory to which random numbers are stored to

### uint64_t frand123Integer64_scalar( int64_t *state )
#### Description
* Return single random 64-bit unsigned integer number discretely uniformly distributed on 0,..,2^64-1
* The counter within the state is incremented appropriately
* Uses either Threefry or ARS for generation of random bits

#### Arguments
* __state__: state used by the RNG

### void frand123Init( int64_t *state, int64_t rank, int64_t threadID, int64_t *seed )
#### Description
* Initialize state for use in serial, MPI-parallel, thread-parallel, and MPI- and thread-parallel settings
* The key is initialized as follows:
    * first  64 bits: rank
    * second 64 bits: threadID
* The counter is initialized as follows:
    * first  64 bits: first element of seed
    * second 64 bits: second element of seed

#### Arguments
* __state__: memory in which to initialize state in
    * dimension: 4
* __rank__: MPI rank of the caller
* __threadID__: thread ID of the thread using the RNG
* __seed__: seed to be used to initialize the counter
    * dimension: 2

## Installation
The Makefile was tested with:
* gcc + gfortran 7.3.0
* icc + ifortran 17.0.2

### Build options
* static library only: _make lib64/libfrand123.a_
* dynamic library only: _make lib64/libfrand123.so_
* static and dynamic libraries: _make all_
* run tests: _make tests_

### Enabling features
* use ARS: add _ars=y_ to _make_ command
* use FMA3 in ARS: add _ars=y fma=y_ to _make_ command
* use gcc: add _gcc=y_ to _make_ command
* use Polar Box-Muller instead of Wichura''s AS 241 PPND16 algorithm for double precision normally distributed random numbers: add _use_polar=y_ to _make_ command

## Examples
For examples, please consult the tests subdirectory

## Tests
### testRandSingle
#### Description
This test consists of three stages:
1. Write out an array of 10^8 single precision uniformly distributed real random numbers to disk
2. Read random numbers into octave
3. Test the null hypothesis on randomness of the numbers against the alternative hypothesis that they are not random with significance level 10^{-5} (runstest, cf. MATLAB documentation)

### testRandDouble
#### Description
This test consists of three stages:
1. Write out an array of 10^8 double precision uniformly distributed real random numbers to disk
2. Read random numbers into octave
3. Test the null hypothesis on randomness of the numbers against the alternative hypothesis that they are not random with significance level 10^{-5} (runstest, cf. MATLAB documentation)

### testAccuracyFloats
#### Description
This test assesses that the value added onto UINT32_MAX in the mapping onto the interval (0,1) is sufficiently large.

### testMomentsSingle
#### Description
This test compares the first 75 moments of the generated single precision real random numbers uniformly distributed in (0,1) to the respective moments of the uniform distribution.
The test is passed if the relative error is below 10^-3

### testMomentsDouble
#### Description
This test compares the first 75 moments of the generated double precision real random numbers uniformly distributed in (0,1) to the respective moments of the uniform distribution.
The test is passed if the relative error is below 10^-3

### testCentralMomentsSingle
#### Description
This test compares the first 75 central moments of the generated single precision real random numbers uniformly distributed in (0,1) to the respective central moments of the uniform distibution in (0,1).
The test is passed if:
* for even moments, the relative error is below 10^-3
* for odd moments, the error of the central moment is below 10^-3

### testCentralMomentsDouble
#### Description
This test compares the first 75 central moments of the generated double precision real random numbers uniformly distributed in (0,1) to the respective central moments of the uniform distibution in (0,1).
The test is passed if:
* for even moments, the relative error is below 10^-3
* for odd moments, the error of the central moment is below 10^-3

### testWichura2x64Kernel
#### Description
This test generates 10^8 uniformly distributed douple precision real random numbers and applies the PPND16 implementation in as241.c (taken from GRASS GIS) to generate reference normally distributed random numbers.
Then, these reference values are compared to those computed from the exact same uniformly distributed numbers using the function wichura2x64Kernel.
Deviations up to 1e-14 are tolerated.

### testRandNormDoublePython
#### Description
This test generates 10^8 normally distributed random numbers and applies the hypothesis-tests skewtest, kurtosistest and normtest from the scipy stats package.
Bounds on the p-values are defined on a per-algorithm basis.

### testNormDoublePerformance
#### Description
Prints out serial and OpenMP parallel timings for the generation and summation of 10^9 random numbers

### testCentralMomentsNormDouble
#### Description
This test generates 10^9 normally distributed random numbers and compute the 2nd, 4th, 6th, 8th, 10th, 12th, 14th, 16th, 18th, 20th central moments and compares these to the central moments of the standard normal distribution.

### testRandNormSinglePython
#### Description
This test generates 10^8 normally distributed random numbers and applies the hypothesis-tests skewtest, kurtosistest and normtest from the scipy stats package.
Bounds on the p-values are defined on a per-RNG basis.

### testNormSinglePerformance
#### Description
Prints out serial and OpenMP parallel timings for the generation and summation of 10^9 random numbers.

### testEquivalence
#### Description
Generates random numbers in pieces of 2/4 and 20 at a time and checks that both versions return the same random numbers when starting with the same state.

### testOdd
#### Description
Tests that generation of odd number of random variables works correctly.

## License
still to do
note that as241.c is taken from GRASS GIS and the original license (most likely GPL >= 2) applies
