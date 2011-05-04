/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_common.c
 *  timestamp: 2011-02-08 21:32:53 +0100
 *  
 */


#include "sl_common.h"


slint merge2_basic_auto_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_auto_01_x */
{
  if (z_min(s0->size, s1->size) <= sqrt(s0->size + s1->size)) return merge2_basic_sbin_01_x(s0, s1, sx);

  return merge2_basic_straight_01_x(s0, s1, sx);
}


slint merge2_basic_01_x(elements_t *s0, elements_t *s1, elements_t *sx, m2x_func _x0_1, m2x_func _0x_1) /* sl_proto, sl_func merge2_basic_01_x */
{
  elements_t s, e, d, x;

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx == NULL)
  {
    elements_alloc(&x, z_min(s0->size, s1->size), SLCM_ALL);
    sx = &x;

  } else if (sx->size < z_min(s0->size, s1->size)) return -1;

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

  if (sx == &x) elements_free(&x);

  return 0;
}


slint merge2_basic_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t, m2X_func _X0_1, m2X_func _0X_1) /* sl_proto, sl_func merge2_basic_01_X */
{
  elements_t s, e, d, x;

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx->size < z_min(s0->size, s1->size)) return -1;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
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
    elem_nxchange_ro0(&s, &d, d.size, t);

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
    elem_nxchange_ro0(s1, &d, d.size, t);

    /* call merge2 subroutine */
    (_0X_1)(s0, &d, &e, t);
  }

  if (t == &x) elements_free(&x);

  return 0;
}


#define the_merge2_01(s0, s1, sx)  merge2_basic_sseq_01(s0, s1, sx)


/*slint merge2_simplify_s0(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements)
{
  return 0;
}*/


slint merge2_simplify_s1(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements) /* sl_proto, sl_func merge2_simplify_s1 */
{
  slint m, m0, m1;

  elements_t _s0, _s1, x;


  s1elements = z_min(s1elements, s1->size);

/*  printf("simplifying %d elements from s1\n", s1elements);*/

  if (s1elements == 0) return 0;

  /** find the s1elements highest elements of s0 and s1 **/

  m = s1elements;
  m0 = m1 = 0;

  elem_assign_at(s0, s0->size - 1, &_s0);
  elem_assign_at(s1, s1->size - 1, &_s1);

  while (m-- > 0 && _s0.keys >= s0->keys)
  if (key_cmp_ge(*_s0.keys, *_s1.keys))
  {
    elem_dec(&_s0);
    ++m0;
  } else
  {
    elem_dec(&_s1);
  }

  elem_inc(&_s0);

  m1 = s1elements - m0;
  elem_assign_at(s1, s1->size - m1, &_s1);

/*  printf("highest %d elements, %d from s0, %d from s1\n", s1elements, m0, m1);*/

  /** bring the s1elements highest elements to the end of s1 **/

  elem_assign_at(s1, s1->size - s1elements, &x);

  elem_nxchange(&_s0, &x, m0, sx);

  /* merge the highest elements */
  x.size = m0;
  _s1.size = m1;

/*  printf("merging x(%d) @ %p & _s1(%d) @ %p\n", x.size, x.keys, _s1.size, _s1.keys);*/

  the_merge2_01(&x, &_s1, sx);

/*  elements_print_keys(s1);*/

  /* merge the highest elements of s0 */
  elem_assign(s0, &x); x.size -= m0;

  _s0.size = m0;

/*  printf("merging x(%d) @ %p & _s0(%d) @ %p\n", x.size, x.keys, _s0.size, _s0.keys);*/

  the_merge2_01(&x, &_s0, sx);

/*  elements_print_keys(s0);*/

  s1->size -= s1elements;

  return 0;
}


#undef the_merge2_01


slint merge2_memory_adaptive(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_memory_adaptive */
{
  slint v;

  v = merge2_basic_auto_01_x(s0, s1, sx);

  if (v >= 0) return v;

  v = merge2_compo_tridgell(s0, s1, sx);

  if (v >= 0) return v;

  return merge2_compo_hula(s0, s1, sx);
  
  if (v >= 0) return v;

  return merge2_compo_hula(s0, s1, NULL);
}
