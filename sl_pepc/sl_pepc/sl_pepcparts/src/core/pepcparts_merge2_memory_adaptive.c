/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepcparts_merge2_memory_adaptive.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


slint pepcparts_merge2_memory_adaptive(elements_t *s0, elements_t *s1, elements_t *sx) /* pepcparts_sl_proto, sl_func pepcparts_merge2_memory_adaptive */
{
  slint v;

  v = pepcparts_merge2_basic_auto_01_x(s0, s1, sx);

  if (v >= 0) return v;

  v = pepcparts_merge2_compo_tridgell(s0, s1, sx);

  if (v >= 0) return v;

  return pepcparts_merge2_compo_hula(s0, s1, sx);
  
  if (v >= 0) return v;

  return pepcparts_merge2_compo_hula(s0, s1, NULL);
}
