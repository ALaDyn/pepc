/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepckeys_merge2_memory_adaptive.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


slint pepckeys_merge2_memory_adaptive(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_memory_adaptive */
{
  slint v;

  v = pepckeys_merge2_basic_auto_01_x(s0, s1, sx);

  if (v >= 0) return v;

  v = pepckeys_merge2_compo_tridgell(s0, s1, sx);

  if (v >= 0) return v;

  return pepckeys_merge2_compo_hula(s0, s1, sx);
  
  if (v >= 0) return v;

  return pepckeys_merge2_compo_hula(s0, s1, NULL);
}
