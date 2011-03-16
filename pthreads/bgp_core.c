#include "fortran2c_types.h"

//////////////// BGP-Core Identification //////////////////////
FINT_TYPE_C get_my_core()
{
  return get_my_core_();
}

#ifdef _BGP

#include <spi/kernel_interface.h>

FINT_TYPE_C get_my_core_()
{
  return Kernel_PhysicalProcessorID();
}
#else
FINT_TYPE_C get_my_core_()
{
  // we have to be sure that on machines with a standard scheduler
  // no thread thinks, he is sharing its processor with someone else
  static FINT_TYPE_C lastreq = 144;

  return ++lastreq;
}
#endif

