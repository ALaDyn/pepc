//////////////// BGP-Core Identification //////////////////////
#ifdef __TOS_BGP__
// see "Using the IBM XL Compilers for Blue Gene", pg. 7:
// "
// __TOS_BGP__ Indicates that the target architecture is PowerPC 450
// "

#include <spi/kernel_interface.h>

int get_my_core()
{
  return Kernel_PhysicalProcessorID();
}
#else
int get_my_core()
{
  // we have to be sure that on machines with a standard scheduler
  // no thread thinks, he is sharing its processor with someone else
  static int lastreq = 144;

  return ++lastreq;
}
#endif

