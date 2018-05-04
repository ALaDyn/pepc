#include <stdio.h>
#include <ars.h>
#include <limits.h>
#include <gdef.h>
#include <bbattery.h>
#include <sstring.h>

const double fac1 = 1. / ( (double)UINT32_MAX + 1. );
const double fac2 = 0.5 / ( (double)UINT32_MAX + 1. );

typedef struct state
{
   ars4x32_key_t key;
   ars4x32_ctr_t ctr;
   short inBuffer;
   uint32_t buffer[4];
} state_t;

void generateRandomBits( state_t *state )
{
   ars4x32_ctr_t int_res = ars4x32_R( 5, state->ctr, state->key );
   state->buffer[0] = int_res.v[0];
   state->buffer[1] = int_res.v[1];
   state->buffer[2] = int_res.v[2];
   state->buffer[3] = int_res.v[3];
   state->inBuffer = 4;
   if( state->ctr.v[0] < UINT_MAX )
      state->ctr.v[0]++;
   else
   {
      state->ctr.v[0] = 0;
      if( state->ctr.v[1] < UINT_MAX )
         state->ctr.v[1]++;
      else
      {
         state->ctr.v[1] = 0;
         if( state->ctr.v[2] < UINT_MAX )
            state->ctr.v[2]++;
         else
         {
            state->ctr.v[2] = 0;
            if( state->ctr.v[3] < UINT_MAX )
               state->ctr.v[3]++;
            else
               state->ctr.v[3] = 0;
         }
      }
   }
   return;
}

uint32_t getRandomBitsFromState( state_t *state )
{
   state->inBuffer = state->inBuffer - 1;
   return state->buffer[state->inBuffer];
}

double getFloat( void *param, void *state_void )
{
   state_t *state = (state_t*)state_void;
   uint32_t bits;
   if( state->inBuffer == 0 )
      generateRandomBits( state );
   bits = getRandomBitsFromState( state );
   return (double)bits * fac1 + fac2;
}

unsigned long getBits( void *param, void *state_void )
{
   state_t *state = (state_t*)state_void;
   if( state->inBuffer == 0 )
      generateRandomBits( state );
   return getRandomBitsFromState( state );
//   state_t *state = (state_t*)state_void;
//   ars4x32_ctr_t int_res = ars4x32_R( 5, state->ctr, state->key );
//   if( state->ctr.v[0] < UINT_MAX )
//      state->ctr.v[0]++;
//   else
//   {
//      state->ctr.v[0] = 0;
//      if( state->ctr.v[1] < UINT_MAX )
//         state->ctr.v[1]++;
//      else
//      {
//         state->ctr.v[1] = 0;
//         if( state->ctr.v[2] < UINT_MAX )
//            state->ctr.v[2]++;
//         else
//         {
//            state->ctr.v[2] = 0;
//            if( state->ctr.v[3] < UINT_MAX )
//               state->ctr.v[3]++;
//            else
//               state->ctr.v[3] = 0;
//         }
//      }
//   }
//   return (unsigned long)int_res.v[0];
}

void write( void *state )
{
   return;
}

int main()
{
   unif01_Gen gen;
   state_t state;
   char name[4] = "ars";
   char a[100];
   gen.state = (void*)&state;
   gen.GetBits = getBits;
   gen.GetU01 = getFloat;
   gen.Write = write;
   gen.name = name;
   gen.param = (void*)&a;
   bbattery_BigCrush( &gen );
   return 0;
}

