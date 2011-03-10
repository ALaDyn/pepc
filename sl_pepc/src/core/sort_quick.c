/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_quick.c
 *  timestamp: 2011-02-08 21:32:53 +0100
 *  
 */

/* quicksort from "Introduction to Algorithms"
   - iterative version with own stack requires pre-allocated stack with max. depth -> bad
   - half-recursive version (qs_halfrec) using "tail recursion" (pg. 162)
   - reducing max. stack depth to \Theta(lg n) proposed by problem 7-4c (pg. 162)
   - hybrid version using own stack of fixed size and recursive calls if the stack exceeds
*/


#include "sl_common.h"


slint qs_halfrec_i(elements_t *s, slint l, slint h, elements_t *xs) /* sl_func qs_halfrec_i */
{
  slint q, i, j;

  elements_t x;

  while (l < h)
  {
    /* partitioning the array */
    elem_assign_at(s, h, &x);

    for (i = j = l; j < h; j++)
    if (key_cmp_le(s->keys[j], *x.keys))
    {
      /* exchange i and j */
      elem_xchange_at(s, i, s, j, xs);
      i++;
    }
    /* exchange i and h */
    elem_xchange_at(s, i, s, h, xs);

    q = i;

    /* call recursive with the smaller part to reduce max. stack depth */
    if (q - l < h - q)
    {
      q--;
      if (l < q) qs_halfrec_i(s, l, q, xs);
      l = q + 2;

    } else
    {
      q++;
      if (q < h) qs_halfrec_i(s, q, h, xs);
      h = q - 2;
    }
  }

  return 0;
}


slint sort_quick_i(elements_t *s, elements_t *xs) /* sl_func sort_quick_i */
{
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  qs_halfrec_i(s, 0, s->size - 1, xs);

  if (xs == &txs) elements_free(xs);

  return 0;
}


slint qs_halfrec_p(elements_t xl, elements_t xh, elements_t *xs) /* sl_func qs_halfrec_p */
{
  elements_t xi, xj;

  while (xl.keys < xh.keys)
  {
    /* partitioning the array */
    elem_assign(&xl, &xi);
    elem_assign(&xl, &xj);

    for (; xj.keys < xh.keys;)
    {
      if (key_cmp_le(*xj.keys, *xh.keys))
      {
        /* exchange i and j */
        elem_xchange(&xi, &xj, xs);
        elem_inc(&xi);
      }
      elem_inc(&xj);
    }

    /* exchange i and h */
    elem_xchange(&xi, &xh, xs);

    /* call recursive with the smaller part to reduce max. stack depth */
    if (xi.keys - xl.keys < xh.keys - xi.keys)
    {
      elem_dec(&xi);
      if (xl.keys < xi.keys) qs_halfrec_p(xl, xi, xs);
      elem_assign_at(&xi, 2, &xl);

    } else
    {
      elem_inc(&xi);
      if (xi.keys < xh.keys) qs_halfrec_p(xi, xh, xs);
      elem_assign_at(&xi, -2, &xh);
    }
  }

  return 0;
}


slint sort_quick_p(elements_t *s, elements_t *xs) /* sl_func sort_quick_p */
{
  elements_t xh, txs;

  if (s == NULL) return -1;

  rti_tstart(rti_tid_sort_quick);

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size > 1)
  {
    elem_assign_at(s, s->size - 1, &xh);
    qs_halfrec_p(*s, xh, xs);
  }

  if (xs == &txs) elements_free(xs);

  rti_tstop(rti_tid_sort_quick);

  return 0;
}


#define stack_size 32

slint qs_hybrid(elements_t xl, elements_t xh, elements_t *xs) /* sl_func qs_hybrid */
{
  elements_t xi, xj;
  struct { elements_t xl, xh; } stack[stack_size], *sp = stack, *lsp = stack + stack_size;

#define push(l, h)  ((sp < lsp)?(sp->xl = l, sp->xh = h, ++sp, 1):0)
#define pop(l, h)   ((sp > stack)?(--sp, l = sp->xl, h = sp->xh, 1):0)

  while (1)
  {
    /* if the current part is to small */
    if (xl.keys >= xh.keys)
    if (!pop(xl, xh)) break; /* pop a new part from the stack, or break if empty */

    /* partitioning the array */
    elem_assign(&xl, &xi);
    elem_assign(&xl, &xj);

    for (; xj.keys < xh.keys;)
    {
      if (key_cmp_le(*xj.keys, *xh.keys))
      {
        /* exchange i and j */
        elem_xchange(&xi, &xj, xs);
        elem_inc(&xi);
      }
      elem_inc(&xj);
    }

    /* put the pivot in the middle */
    elem_xchange(&xi, &xh, xs);

    /* push smaller part on the stack (call recursive if full) to reduce max. stack depth */
    if (xi.keys - xl.keys < xh.keys - xi.keys)
    {
      elem_dec(&xi);

      if (xl.keys < xi.keys)
      if (!push(xl, xi)) qs_hybrid(xl, xi, xs); /* call if push fails */

      elem_assign_at(&xi, 2, &xl);

    } else
    {
      elem_inc(&xi);

      if (xi.keys < xh.keys)
      if (!push(xi, xh)) qs_hybrid(xi, xh, xs); /* call if push fails */

      elem_assign_at(&xi, -2, &xh);
    }
  }

#undef push
#undef pop

  return 0;
}


slint sort_quick_h(elements_t *s, elements_t *xs) /* sl_func sort_quick_h */
{
  elements_t xh, txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  elem_assign_at(s, s->size - 1, &xh);
  qs_hybrid(*s, xh, xs);

  if (xs == &txs) elements_free(xs);

  return 0;
}

slint sort_quick(elements_t *s, elements_t *xs) /* sl_proto, sl_func sort_quick */
{
  return sort_quick_p(s, xs);
}
