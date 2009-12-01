/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_basic.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


slint pepcparts_merge2_basic_01_x(elements_t *s0, elements_t *s1, elements_t *sx, m2x_func _x0_1, m2x_func _0x_1) /* pepcparts_sl_proto, sl_func pepcparts_merge2_basic_01_x */
{
  elements_t s, e, d, x;

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx == NULL)
  {
    pepcparts_elements_alloc(&x, xmin(s0->size, s1->size), 1, 1);
    sx = &x;

  } else if (sx->size < xmin(s0->size, s1->size)) return -1;

  if (s0->size <= s1->size)
  {
    /* initialize */
    elem_assign(s0, &s); elem_assign_at(s0, s0->size, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s0 */
    while (s.keys != e.keys)
    if (key_cmp_le(*s.keys, *e.keys)) elem_inc(&s); else break;

    /* evacuate wrong-placed elements of s0 */
    d.size = s.size = e.keys - s.keys;
    elem_ncopy(&s, &d, d.size);

    /* call merge2 subroutine */
    (_x0_1)(s1, &d, &s);

  } else
  {
    /* initialize */
    elem_assign_at(s1, s1->size - 1, &s); elem_assign_at(s1, -1, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s1 */
    while (s.keys != e.keys)
    if (key_cmp_gt(*s.keys, *e.keys)) elem_dec(&s); else break;

    /* evacuate wrong-placed elements of s1 */
    d.size = e.size = s.keys - e.keys;
    elem_inc(&e);
    elem_ncopy(s1, &d, d.size);

    /* call merge2 subroutine */
    (_0x_1)(s0, &d, &e);
  }

  if (sx == &x) pepcparts_elements_free(&x);

  return 0;
}


slint pepcparts_merge2_basic_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t, m2X_func _X0_1, m2X_func _0X_1) /* pepcparts_sl_proto, sl_func pepcparts_merge2_basic_01_X */
{
  elements_t s, e, d, x;

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx->size < xmin(s0->size, s1->size)) return -1;

  if (t == NULL)
  {
    pepcparts_elements_alloc(&x, 1, 1, 1);
    t = &x;
  }

  if (s0->size <= s1->size)
  {
    /* initialize */
    elem_assign(s0, &s); elem_assign_at(s0, s0->size, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s0 */
    while (s.keys != e.keys)
    if (key_cmp_le(*s.keys, *e.keys)) elem_inc(&s); else break;

    /* evacuate wrong-placed elements of s0 */
    d.size = s.size = e.keys - s.keys;
    pepcparts_elem_nxchange_ro0(&s, &d, d.size, t);

    /* call merge2 subroutine */
    (_X0_1)(s1, &d, &s, t);

  } else
  {
    /* initialize */
    elem_assign_at(s1, s1->size - 1, &s); elem_assign_at(s1, -1, &e);
    elem_assign(sx, &d);

    /* skip already well-placed elements of s1 */
    while (s.keys != e.keys)
    if (key_cmp_gt(*s.keys, *e.keys)) elem_dec(&s); else break;

    /* evacuate wrong-placed elements of s1 */
    d.size = e.size = s.keys - e.keys;
    elem_inc(&e);
    pepcparts_elem_nxchange_ro0(s1, &d, d.size, t);

    /* call merge2 subroutine */
    (_0X_1)(s0, &d, &e, t);
  }

  if (t == &x) pepcparts_elements_free(&x);

  return 0;
}
