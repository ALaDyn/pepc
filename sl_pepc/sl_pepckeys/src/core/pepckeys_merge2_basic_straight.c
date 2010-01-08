/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_basic_straight.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


/* front2back, overriding the content of the extraspace */
slint pepckeys_merge2_basic_straight_x0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_x0_1 */
{
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst);

  /* process until one of the srcs is empty */
  while (src0.keys != src0e.keys && src1.keys != src1e.keys)
  {
    if (key_cmp_le(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_inc(&src0);
    } else
    {
      elem_copy(&src1, &dst);
      elem_inc(&src1);
    }
    elem_inc(&dst);
  }

  /* copy the remaining elements of the 1st src */
  src1.size = src1e.keys - src1.keys;
  if (src1.size > 0) elem_ncopy(&src1, &dst, src1.size);

  return 0;
}


/* back2front, overriding the content of the extraspace */
slint pepckeys_merge2_basic_straight_0x_1(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_0x_1 */
{
  elements_t src0, src0e, src1, src1e, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign_at(s0, s0->size - 1, &src0); elem_assign_at(s0, -1, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign_at(s1, -1, &src1e);
  elem_assign_at(sx, s1->size - 1, &dst);

  /* process until one of the srcs is empty */
  while (src0.keys != src0e.keys && src1.keys != src1e.keys)
  {
    if (key_cmp_gt(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_dec(&src0);
    } else
    {
      elem_copy(&src1, &dst);
      elem_dec(&src1);
    }
    elem_dec(&dst);
  }

  /* copy the remaining elements of the 1st src */
  src1.size = src1.keys - src1e.keys;
  if (src1.size > 0) elem_ncopy(s1, s0, src1.size);

  return 0;
}


slint pepckeys_merge2_basic_straight_01_x(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_01_x */
{
  return pepckeys_merge2_basic_01_x(s0, s1, sx, pepckeys_merge2_basic_straight_x0_1, pepckeys_merge2_basic_straight_0x_1);
}


slint pepckeys_merge2_basic_straight_x_0_1(elements_t *s0, elements_t *s1, elements_t *sx) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_x_0_1 */
{
  elements_t dst, dste, src0e, src1e;

  elem_assign_at(s0, s0->size, &src0e);
  elem_assign_at(s1, s1->size, &src1e);
  elem_assign(sx, &dst); elem_assign_at(sx, sx->size, &dste);

  while (s0->keys != src0e.keys && s1->keys != src1e.keys && dst.keys != dste.keys)
  {
    if (key_cmp_le(*s0->keys, *s1->keys))
    {
      elem_copy(s0, &dst);
      elem_inc(s0);

    } else
    {
      elem_copy(s1, &dst);
      elem_inc(s1);
    }
    elem_inc(&dst);
  }

  s0->size = src0e.keys - s0->keys;
  s1->size = src1e.keys - s1->keys;

  return dst.keys - sx->keys;
}


slint pepckeys_merge2_basic_straight_X0_1(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_X0_1 */
{
  elements_t src0, src0e, src1, src1l, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size - 1, &src1l);
  elem_assign(sx, &dst);

  /* create the hole */
  elem_copy(&dst, t);

  while (src0.keys != src0e.keys && src1.keys != src1l.keys)
  {
    if (key_cmp_le(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src0);
      elem_inc(&src0);

    } else
    {
      elem_copy(&src1, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src1);
      elem_inc(&src1);
    }
  }

  /* now: either src0 is empty or src1 contains just one element */

  while (src0.keys != src0e.keys)
  if (key_cmp_le(*src0.keys, *src1.keys))
  {
    /* remove one element from src0 */
    elem_copy(&src0, &dst);
    elem_inc(&dst);
    elem_inc(&src0);

  } else break; /* the last of src1 has to move */

  /* close the hole at dst */
  elem_copy(t, &dst);

  /* exchange the remaining elements at src1 with dst */
  src1.size = src1l.keys - src1.keys + 1;
  pepckeys_elem_nxchange_ro0(&src1, &dst, src1.size, t);

  return 0;
}


slint pepckeys_merge2_basic_straight_0X_1(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_0X_1 */
{
  elements_t src0, src1, src0e, src1l, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign_at(s0, s0->size - 1, &src0); elem_assign_at(s0, -1, &src0e);
  elem_assign_at(s1, s1->size - 1, &src1); elem_assign(s1, &src1l);
  elem_assign_at(sx, s1->size - 1, &dst);

  /* create the hole */
  elem_copy(&dst, t);

  while (src0.keys != src0e.keys && src1.keys != src1l.keys)
  {
    if (key_cmp_gt(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_dec(&dst);
      elem_copy(&dst, &src0);
      elem_dec(&src0);

    } else
    {
      elem_copy(&src1, &dst);
      elem_dec(&dst);
      elem_copy(&dst, &src1);
      elem_dec(&src1);
    }
  }

  /* now: either src0 is empty or src1 contains just one element */

  while (src0.keys != src0e.keys)
  if (key_cmp_gt(*src0.keys, *src1.keys))
  {
    /* remove one (maybe the last) from src0 */
    elem_copy(&src0, &dst);
    elem_dec(&dst);
    elem_dec(&src0);

  } else break; /* the last of src1 has to move */

  /* close the hole at dst */
  elem_copy(t, &dst);

  /* xchange the remaining elements at s1 with s0 */
  src1.size = src1.keys - src1l.keys + 1;
  elem_inc(&src0);
  pepckeys_elem_nxchange_ro0(s1, &src0, src1.size, t);

  return 0;
}


slint pepckeys_merge2_basic_straight_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_01_X */
{
  return pepckeys_merge2_basic_01_X(s0, s1, sx, t, pepckeys_merge2_basic_straight_X0_1, pepckeys_merge2_basic_straight_0X_1);
}


slint pepckeys_merge2_basic_straight_X0_1u(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t) /* pepckeys_sl_proto, sl_func pepckeys_merge2_basic_straight_X0_1u */
{
  elements_t src0, src0e, src1, src1l, dst;

  /* if the second list is empty, there is nothing to do */
  if (s1->size == 0) return 0;

  /* initialize */
  elem_assign(s0, &src0); elem_assign_at(s0, s0->size, &src0e);
  elem_assign(s1, &src1); elem_assign_at(s1, s1->size - 1, &src1l);
  elem_assign(sx, &dst);

  /* create the hole */
  elem_copy(&dst, t);

  while (src0.keys != src0e.keys && src1.keys != src1l.keys)
  {
    if (key_cmp_le(*src0.keys, *src1.keys))
    {
      elem_copy(&src0, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src0);
      elem_inc(&src0);

    } else
    {
      elem_copy(&src1, &dst);
      elem_inc(&dst);
      elem_copy(&dst, &src1);
      elem_inc(&src1);
    }
  }

  /* now: either src0 is empty or src1 contains just one element */

  while (src0.keys != src0e.keys)
  if (key_cmp_le(*src0.keys, *src1.keys))
  {
    /* remove one element from src0 */
    elem_copy(&src0, &dst);
    elem_inc(&dst);
    elem_inc(&src0);

  } else
  {
    /* the last of src1 has to move */
    elem_copy(&src1, &dst);

    /* let the hole (dst) point to the last element of s1 */
    elem_assign(&src1, &dst);

    /* let src1 point behind src1l for computing the right number of unfinished elements of s1 (0!) */
    elem_inc(&src1);

    break;
  }

  /* close the hole at dst */
  elem_copy(t, &dst);

  /* return the number of unfinished elements of s1 */
  return src1l.keys - src1.keys + 1;
}
