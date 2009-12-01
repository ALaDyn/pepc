/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_basic_sseq.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


#define the_search_lt  pepckeys_sl_search_sequential_lt
#define the_search_le  pepckeys_sl_search_sequential_le
#define the_search_gt  pepckeys_sl_search_sequential_gt
#define the_search_ge  pepckeys_sl_search_sequential_ge

#define the_nmove      elem_nmove
#define the_rotate     pepckeys_elem_rotate


slint pepckeys_merge2_basic_sseq_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_sseq_x0_1 */
{
  slint n;
  elements_t src0, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  /* process until one of the srcs is empty */
  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_le(&src0, src1.keys);

    if (n > 0)
    {
      the_nmove(&src0, &dst, n);
      elem_add(&src0, n);
      elem_add(&dst, n);
      src0.size -= n;
    }

    elem_copy(&src1, &dst);
    elem_inc(&src1);
    elem_inc(&dst);
  }

  /* copy the remaining elements of s1 to dst */
  src1.size = src1e.keys - src1.keys;
  if (src1.size > 0) elem_ncopy(&src1, &dst, src1.size);

  return 0;
}


slint pepckeys_merge2_basic_sseq_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_sseq_0x_1 */
{
  slint n;
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, sx->size, &dst);

  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_gt(&src0, src1.keys);

    if (n > 0)
    {
      elem_sub(&dst, n);
      elem_sub(&src0e, n);
      the_nmove(&src0e, &dst, n);
      src0.size -= n;
    }

    elem_dec(&dst);
    elem_copy(&src1, &dst);
    elem_dec(&src1);
  }

  /* copy the remaining elements of s1 to the front */
  src1.size = src1.keys - src1e.keys;
  if (src1.size > 0) elem_ncopy(&src1, s0, src1.size);

  return 0;
}


slint pepckeys_merge2_basic_sseq_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_sseq_01_x */
{
  return pepckeys_merge2_basic_01_x(s0, s1, sx, pepckeys_merge2_basic_sseq_x0_1, pepckeys_merge2_basic_sseq_0x_1);
}


slint pepckeys_merge2_basic_sseq_01(elements_t *s0, elements_t *s1, elements_t *t) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_sseq_01 */
{
  slint k;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    pepckeys_elements_alloc(&x, 1, 1, 1);
    t = &x;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  elem_assign_at(s1, s1->size - 1, &last);

  while (_s0.size > 0 && _s1.size > 0)
  if (_s0.size <= _s1.size)
  {
    k = the_search_lt(&_s1, _s0.keys);

    the_rotate(&_s0, _s0.size, k, t);

    elem_add(&_s0, k + 1);
    _s0.size -= 1;
    elem_add(&_s1, k);
    _s1.size -= k;

  } else
  {
    k = the_search_gt(&_s0, last.keys);

    elem_sub(&_s1, k);

    the_rotate(&_s1, k, _s1.size, t);

    elem_sub(&last, k + 1);
    _s0.size -= k;
    _s1.size -= 1;
  }

  if (t == &x) pepckeys_elements_free(&x);

  return 0;
}
