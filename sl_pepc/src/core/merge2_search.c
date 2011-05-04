/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_search.c
 *  timestamp: 2011-02-08 21:32:53 +0100
 *  
 */


#include "sl_common.h"


#define the_ncopy      elem_ncopy
#define the_nmove      elem_nmove
#define the_rotate     elem_rotate


#define the_search_lt  sl_search_sequential_lt
#define the_search_le  sl_search_sequential_le
#define the_search_gt  sl_search_sequential_gt
#define the_search_ge  sl_search_sequential_ge


slint merge2_basic_sseq_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sseq_x0_1 */
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


slint merge2_basic_sseq_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sseq_0x_1 */
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


slint merge2_basic_sseq_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sseq_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_sseq_x0_1, merge2_basic_sseq_0x_1);
}


slint merge2_basic_sseq_01(elements_t *s0, elements_t *s1, elements_t *t) /* sl_proto, sl_func merge2_basic_sseq_01 */
{
  slint k;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
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

  if (t == &x) elements_free(&x);

  return 0;
}


#undef the_search_lt
#undef the_search_le
#undef the_search_gt
#undef the_search_ge


#define the_search_lt  sl_search_binary_lt
#define the_search_le  sl_search_binary_le
#define the_search_gt  sl_search_binary_gt
#define the_search_ge  sl_search_binary_ge


slint merge2_basic_sbin_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sbin_x0_1 */
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
  if (src1.size > 0) the_ncopy(&src1, &dst, src1.size);

  return 0;
}


slint merge2_basic_sbin_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sbin_0x_1 */
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
  elem_sub(&src1, src1.size - 1);
  if (src1.size > 0) the_ncopy(&src1, s0, src1.size);

  return 0;
}


slint merge2_basic_sbin_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_sbin_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_sbin_x0_1, merge2_basic_sbin_0x_1);
}


slint merge2_basic_sbin_01(elements_t *s0, elements_t *s1, elements_t *t) /* sl_proto, sl_func merge2_basic_sbin_01 */
{
  slint k;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
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

  if (t == &x) elements_free(&x);

  return 0;
}


#undef the_search_lt
#undef the_search_le
#undef the_search_gt
#undef the_search_ge


/* using "hybrid-search" for binary-merge from [Hwang, Lin] (Knuth v3) */

#define the_search_lt  sl_search_hybrid_lt
#define the_search_le  sl_search_hybrid_le
#define the_search_gt  sl_search_hybrid_gt
#define the_search_ge  sl_search_hybrid_ge


slint merge2_basic_shyb_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_shyb_x0_1 */
{
  slint n, l;
  elements_t src0, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  l = z_powof2((slint) (log((double) s0->size / (double) s1->size) / log(2.0)));

  /* process until one of the srcs is empty */
  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_le(&src0, src1.keys, l);

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


slint merge2_basic_shyb_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_shyb_0x_1 */
{
  slint n, l;
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, sx->size, &dst);

  l = z_powof2((slint) (log((double) s0->size / (double) s1->size) / log(2.0)));

  while (src0.size > 0 && src1.keys != src1e.keys)
  {
    n = the_search_gt(&src0, src1.keys, l);

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


slint merge2_basic_shyb_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_basic_shyb_01_x */
{
  return merge2_basic_01_x(s0, s1, sx, merge2_basic_shyb_x0_1, merge2_basic_shyb_0x_1);
}


slint merge2_basic_shyb_01(elements_t *s0, elements_t *s1, elements_t *t) /* sl_proto, sl_func merge2_basic_shyb_01 */
{
  slint k, l;
  elements_t x, _s0, _s1, last;

  if (t == NULL)
  {
    elements_alloc(&x, 1, SLCM_ALL);
    t = &x;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  elem_assign_at(s1, s1->size - 1, &last);

  while (_s0.size > 0 && _s1.size > 0)
  if (_s0.size <= _s1.size)
  {
    l = z_powof2((slint) (log((double) _s1.size / (double) _s0.size) / log(2.0)));

    k = the_search_lt(&_s1, _s0.keys, l);

    the_rotate(&_s0, _s0.size, k, t);

    elem_add(&_s0, k + 1);
    _s0.size -= 1;
    elem_add(&_s1, k);
    _s1.size -= k;

  } else
  {
    l = z_powof2((slint) (log((double) _s0.size / (double) _s1.size) / log(2.0)));

    k = the_search_gt(&_s0, last.keys, l);

    elem_sub(&_s1, k);

    the_rotate(&_s1, k, _s1.size, t);

    elem_sub(&last, k + 1);
    _s0.size -= k;
    _s1.size -= 1;
  }

  if (t == &x) elements_free(&x);

  return 0;
}


#undef the_search_lt
#undef the_search_le
#undef the_search_gt
#undef the_search_ge
