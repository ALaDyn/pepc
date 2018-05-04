/*
 * Example illustrating the use of the frand123 C interface for the generation
 * of 32-bit signed and unsigned integer random numbers.
 * The unsigned integer random numbers are used as containers for ranom bits.
 * The signed integer random numbers are used as experiment for coin flip.
 */

#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include "../../wrapper/frand123.h"

// check if the position-th bit is set in bitset
bool isBitSet( uint64_t bitset, short position )
{
   // create a mask by shifting 1 left by position, i.e. calculating 2^position
   uint32_t mask = UINT64_C( 1 ) << position;
   // compare to bitset
   uint32_t masked = bitset & mask;
   // return true if it is set
   return masked > 0;
}

// check whether certain bits are set in random unsigned 32-bit integers
void checkBitSets( int64_t *state )
{
   // generate five uint64_t as containers for random bits
   uint32_t randomBits[ 5 ];
   frand123UnsignedInteger32( state, 5, randomBits );
   
   // check each case
   for( int i = 1; i < 5; i++ )
   {
      printf( "isBitSet( %" PRIu32 ", 31): %s\n", randomBits[ i ], isBitSet( randomBits[ i ], 31 ) ? "true" : "false" );
      printf( "isBitSet( %" PRIu32 ",  0): %s\n", randomBits[ i ], isBitSet( randomBits[ i ], 0 ) ? "true" : "false" );
   }

   return;
}

// simulates a coin flip:
// > 0 heads
// < 0 tails
// = 0 what a luck
void naiveCoinFlip( int64_t *state )
{
   char heads[ 6 ] = "heads";
   char tails[ 6 ] = "tails";
   char zero[ 12 ] = "what a luck";
   for( int i = 0; i < 10; i++ )
   {
      int32_t flip = frand123Integer32_scalar( state );
      char *text;
      if( flip == 0 )
      {
         text = zero;
      }
      else if( flip > 0 )
      {
         text = heads;
      }
      else
      {
         text = tails;
      }
      printf( "Flip %d: %s\n", i, text );
   }
   return;
}

int main()
{
   // build state
   int64_t state[ 4 ];
   frand123Init( state, 0, 0, NULL );

   checkBitSets( state );

   printf( "\n\n" );

   naiveCoinFlip( state );

   return 0;
}
