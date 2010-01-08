/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_basic_auto.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


slint pepcparts_merge2_basic_auto_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* pepcparts_sl_proto, sl_func pepcparts_merge2_basic_auto_01_x */
{
  if (xmin(s0->size, s1->size) <= sqrt(s0->size + s1->size)) return pepcparts_merge2_basic_sbin_01_x(s0, s1, sx);

  return pepcparts_merge2_basic_straight_01_x(s0, s1, sx);
}
