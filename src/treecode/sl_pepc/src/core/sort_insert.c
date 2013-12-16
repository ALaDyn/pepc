/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_insert.c
 *  
 */


#include "sl_common.h"


#ifdef key_integer

slint_t sort_insert_bmask_kernel(elements_t *s, elements_t *sx, slkey_pure_t bmask) /* sl_proto, sl_func sort_insert_bmask_kernel */
{
  slint_t i, j;

  for (i = 1; i < s->size; i++)
  {
    if (key_pure_cmp_lt(key_purify(s->keys[i]) & bmask, key_purify(s->keys[i - 1]) & bmask))
    {
      j = i - 1;
      elem_copy_at(s, i, sx, 0);

      do
      {
        elem_copy_at(s, j, s, j + 1);
        if (--j < 0) break;

      } while (key_pure_cmp_lt(key_purify(*sx->keys) & bmask, key_purify(s->keys[j]) & bmask));

      elem_copy_at(sx, 0, s, j + 1);
    }
  }

  return 0;
}

#endif


slint_t sort_insert_kernel(elements_t *s, elements_t *sx) /* sl_func sort_insert_kernel */
{
  slint_t i, j;

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


slint_t sort_insert(elements_t *s, elements_t *sx) /* sl_proto, sl_func sort_insert */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_insert);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  sort_insert_kernel(s, sx);

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_insert);

  return 0;
}
