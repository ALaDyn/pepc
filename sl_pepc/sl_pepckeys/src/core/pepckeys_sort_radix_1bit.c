/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepckeys_sort_radix_1bit.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* simple MSDF 1-bit radix-sort */

#include "sl_common.h"


#ifdef key_integer


#define radix_low    0
#define radix_high   (sizeof(slkey_pure_t) * 8 - 1)


inline slint pepckeys_split2_b_1brs(elements_t *s, elements_t *sx, slkey_pure_t bmask) /* sl_func pepckeys_split2_b_1brs */
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


slint pepckeys_rs_1bit_rec(elements_t *s, elements_t *sx, slint rhigh, slint rlow) /* sl_func pepckeys_rs_1bit_rec */
{
  slkey_pure_t bmask = powof2_typed(rhigh, slkey_pure_t);

  elements_t xl, xh;

  slint n0, n1;

  elem_assign(s, &xl);

  while (xl.size > 1)
  {
    n0 = pepckeys_split2_b_1brs(&xl, sx, bmask);
    n1 = xl.size - n0;
    
    if (rhigh <= rlow) break;

    rhigh--;
    bmask >>= 1;

    xl.size = n0;
    
    if (n0 <= n1)
    {
      pepckeys_rs_1bit_rec(&xl, sx, rhigh, rlow);

      elem_add(&xl, n0);
      xl.size = n1;

    } else
    {
      elem_assign_at(&xl, n0, &xh);
      xh.size = n1;

      pepckeys_rs_1bit_rec(&xh, sx, rhigh, rlow);
    }
  }

  return 0;
}


slint pepckeys_sort_radix_1bit(elements_t *s, elements_t *sx, slint rhigh, slint rlow) /* pepckeys_sl_proto, sl_func pepckeys_sort_radix_1bit */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  if (sx == NULL)
  {
    sx = &_sx;
    pepckeys_elements_alloc(sx, 1, 1, 1);
  }

  if (rhigh < 0) rhigh = radix_high;
  if (rlow < 0) rlow = radix_low;

  pepckeys_rs_1bit_rec(s, sx, rhigh, rlow);

  if (sx == &_sx) pepckeys_elements_free(sx);

  return 0;
}


#endif /* key_integer */
