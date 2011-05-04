/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_heap.c
 *  timestamp: 2010-12-19 13:24:42 +0100
 *  
 */


#include "sl_common.h"



inline slint hs_heapify_i(elements_t *s, slint i, slint size, elements_t *xs) /* sl_func hs_heapify_i */
{
  slint l, r, largest;

  elem_copy_at(s, i, xs, 0);

  while ((l = 2 * i + 1) < size)
  {
    largest = l;
    r = 2 * i + 2;

    if (r < size)
    if (key_cmp_gt(s->keys[r], s->keys[l])) largest = r;

    if (key_cmp_gt(*xs->keys, s->keys[largest])) break;

    elem_copy_at(s, largest, s, i);

    i = largest;
  }

  elem_copy_at(xs, 0, s, i);

  return 0;
}


slint sort_heap_i(elements_t *s, elements_t *xs) /* sl_func sort_heap_i */
{
  slint i;
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  /* build the heap */
  for (i = s->size / 2 - 1; i >= 0; --i) hs_heapify_i(s, i, s->size, xs);

  /* extract the maxima */
  for (i = s->size - 1; i > 0; --i)
  {
    elem_xchange_at(s, 0, s, i, xs);
    hs_heapify_i(s, 0, i, xs);
  }

  if (xs == &txs) elements_free(xs);

  return 0;
}


slint sort_heap_p(elements_t *s, elements_t *xs) /* sl_func sort_heap_p */
{
  elements_t txs, last, si, sj, l;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  if (s->size < 2) return 0;

  elem_assign_at(s, s->size, &last);

  /* build the heap */
  for (elem_assign_at(s, s->size / 2 - 1, &si); si.keys >= s->keys; elem_dec(&si))
  {
    elem_assign(&si, &sj);

    elem_copy(&sj, xs);

    while (1)
    {
      elem_assign_at(s, (sj.keys - s->keys) * 2 + 1, &l); /* l = i * 2 + 1; */

      if (l.keys > last.keys) break;

      if (l.keys < last.keys)
      if (key_cmp_gt(l.keys[1], *l.keys)) elem_inc(&l);

      if (key_cmp_gt(*xs->keys, *l.keys)) break;

      elem_copy(&l, &sj);

      elem_assign(&l, &sj);
    }

    elem_copy(xs, &sj);
  }

  /* extract the maxima */
  for (elem_assign_at(s, s->size - 1, &si); si.keys > s->keys;)
  {
    elem_copy(&si, xs);
    elem_copy(s, &si);

    /* si becomes the last element in the heap */
    elem_dec(&si);

    elem_assign(s, &sj);

    while (1)
    {
      elem_assign_at(s, (sj.keys - s->keys) * 2 + 1, &l); /* l = i * 2 + 1; */

      if (l.keys > si.keys) break;

      if (l.keys < si.keys)
      if (key_cmp_gt(l.keys[1], *l.keys)) elem_inc(&l);

      if (key_cmp_gt(*xs->keys, *l.keys)) break;

      elem_copy(&l, &sj);

      elem_assign(&l, &sj);
    }

    elem_copy(xs, &sj);
  }

  if (xs == &txs) elements_free(xs);

  return 0;
}


slint sort_heap(elements_t *s, elements_t *xs) /* sl_proto, sl_func sort_heap */
{
  return sort_heap_p(s, xs);
}
