/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepckeys_sort_insert.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


slint pepckeys_sort_insert_(elements_t *s, elements_t *sx) /* sl_func pepckeys_sort_insert_ */
{
  slint i, j;

  for (i = 1; i < s->size; i++)
  {
    if (key_cmp_lt(s->keys[i], s->keys[i - 1]))
    {
      j = i - 1;
      elem_copy_at(s, i, sx, 0);

      do
      {
        elem_copy_at(s, j, s, j + 1);
        if (--j < 0) break;

      } while (key_cmp_lt(*sx->keys, s->keys[j]));

      elem_copy_at(sx, 0, s, j + 1);
    }
  }

  return 0;
}


slint pepckeys_sort_insert(elements_t *s, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_sort_insert */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_insert);

  if (sx == NULL)
  {
    sx = &_sx;
    pepckeys_elements_alloc(sx, 1, 1, 1);
  }

  pepckeys_sort_insert_(s, sx);

  if (sx == &_sx) pepckeys_elements_free(sx);

  rti_tstop(rti_tid_sort_insert);

  return 0;
}
