/*
 * Test case assessing that the static const float enlarger 257.f is sufficiently large
 * to restrict floats in the range [INT32_MIN,INT32_MAX] onto the open set (-0.5,0.5)
 */

#include <stdio.h>
#include <stdint.h>
#include "../wrapper/frand123enlarger.h"

int main()
{
   const float rec = 1.f / ( (float)UINT32_MAX + enlarger );
   float res_max = (float)INT32_MAX * rec;
   float res_min = (float)INT32_MIN * rec;
   if( ( res_max < 0.5f ) && ( res_min > -0.5f ) )
   {
      printf( "passed accuracy test for floats\n" );
      return 0;
   }
   else
   {
      printf( "failed accuracy test for floats\n" );
      return 1;
   }
}
