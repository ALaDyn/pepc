/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_radix_1bit.c
 *  
 */

/* simple MSDF 1-bit radix-sort */

#include "sl_common.h"


#ifdef key_integer


slint_t split2_b_1brs(elements_t *s, elements_t *sx, slkey_pure_t bmask) /* sl_func split2_b_1brs */
{
  elements_t xl, xh;

  elem_assign(s, &xl);
  elem_assign_at(s, s->size - 1, &xh);

  while (1)
  {
    while (xl.keys < xh.keys)
    if (key_purify(*xl.keys) & bmask) break; else elem_inc(&xl);

    while (xl.keys < xh.keys)
    if (key_purify(*xh.keys) & bmask) elem_dec(&xh); else break;

    if (xl.keys >= xh.keys) break;

    elem_xchange(&xl, &xh, sx);
    elem_inc(&xl);
    elem_dec(&xh);
  }

  return xl.keys - s->keys + ((key_purify(*xl.keys) & bmask) == 0);
}


slint_t sort_radix_1bit_kernel(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow) /* sl_proto, sl_func sort_radix_1bit_kernel */
{
  slkey_pure_t bmask;

  elements_t xl, xh;

  slint_t n0, n1;

  elem_assign(s, &xl);

  while (xl.size > 1)
  {
    bmask = z_powof2_typed(rhigh, slkey_pure_t);

    n0 = split2_b_1brs(&xl, sx, bmask);
    n1 = xl.size - n0;

    if (rhigh <= rlow) break;

    rhigh--;

    xl.size = n0;
    
    if (n0 <= n1)
    {
      sort_radix_1bit_kernel(&xl, sx, rhigh, rlow);

      elem_add(&xl, n0);
      xl.size = n1;

    } else
    {
      elem_assign_at(&xl, n0, &xh);
      xh.size = n1;

      sort_radix_1bit_kernel(&xh, sx, rhigh, rlow);
    }
  }

  return 0;
}


slint sort_radix_1bit(elements_t *s, elements_t *sx, slint_t rhigh, slint_t rlow) /* sl_proto, sl_func sort_radix_1bit */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;

  sort_radix_1bit_kernel(s, sx, rhigh, rlow);

  if (sx == &_sx) elements_free(sx);

  return 0;
}


#endif /* key_integer */
