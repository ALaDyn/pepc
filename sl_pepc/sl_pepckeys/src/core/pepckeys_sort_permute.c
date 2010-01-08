/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_permute.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


inline slint pepckeys_sort_permute_forward_(elements_t *s, elements_t *sx, slint *perm, slint offset) /* sl_func pepckeys_sort_permute_forward_ */
{
  elements_t src, end;

  slint j, i, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  i = offset;
  ia = perm;
  while (src.keys != end.keys)
  {
    while (*ia != i)
    {
      j = *ia - offset;

      elem_xchange_at(&src, 0, s, j, sx);

      ja = perm + j;
      *ia = *ja;
      *ja = j + offset;
    }

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


inline slint pepckeys_sort_permute_forward_masked(elements_t *s, elements_t *sx, slint *perm, slint offset, slint bit_mask) /* sl_func pepckeys_sort_permute_forward_masked */
{
  elements_t src, end;

  slint inv_bit_mask = bit_mask ^ -1L;
  slint j, i, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  i = 0;
  ia = perm;
  while (src.keys != end.keys)
  {
    if ((*ia & bit_mask) == 0)
    {
      j = *ia - offset;

      while (i != j)
      {
        elem_xchange_at(&src, 0, s, j, sx);

        ja = perm + j;
        j = *ja - offset;
        *ja |= bit_mask;
      }
    }

    *ia &= inv_bit_mask;

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


slint pepckeys_sort_permute_forward(elements_t *s, elements_t *sx, slint *perm, slint offset, slint mask_bit) /* pepckeys_sl_proto, sl_func pepckeys_sort_permute_forward */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_permute_forward);

  if (sx == NULL)
  {
    sx = &_sx;
    pepckeys_elements_alloc(sx, 1, 1, 1);
  }

  if (mask_bit < 0) pepckeys_sort_permute_forward_(s, sx, perm, offset);
  else pepckeys_sort_permute_forward_masked(s, sx, perm, offset, 1L << xmin(mask_bit, (sizeof(slint) * 8) - 1));

  if (sx == &_sx) pepckeys_elements_free(sx);

  rti_tstop(rti_tid_sort_permute_forward);

  return 0;
}


inline slint pepckeys_sort_permute_backward_(elements_t *s, elements_t *sx, slint *perm, slint offset) /* sl_func pepckeys_sort_permute_backward_ */
{
  elements_t src, end, e0, e1, *from, *to, *t;

  slint i, j, k, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  from = &e0;
  to = &e1;

  i = offset;
  ia = perm;
  while (src.keys != end.keys)
  {
    if (*ia != i)
    {
      elem_copy(&src, sx);
      elem_assign(&src, to);

      ja = ia;
      j = i;

      while (i != (k = *ja))
      {
        elem_assign_at(s, k - offset, from);
        elem_copy(from, to);

        t = to;
        to = from;
        from = t;

        *ja = j;
        ja = perm + (j = k) - offset;
      }

      elem_copy(sx, to);
      *ja = j;
    }

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


inline slint pepckeys_sort_permute_backward_masked(elements_t *s, elements_t *sx, slint *perm, slint offset, slint bit_mask) /* sl_func pepckeys_sort_permute_backward_masked */
{
  elements_t src, end, e0, e1, *from, *to, *t;

  slint inv_bit_mask = bit_mask ^ -1L;

  slint i, k, *ia, *ja;

  elem_assign(s, &src);
  elem_assign_at(s, s->size, &end);

  from = &e0;
  to = &e1;

  i = offset;
  ia = perm;
  while (src.keys != end.keys)
  {
    if ((*ia & bit_mask) == 0)
    {
      elem_copy(&src, sx);
      elem_assign(&src, to);

      ja = ia;

      while (i != (k = *ja))
      {
        elem_assign_at(s, k - offset, from);
        elem_copy(from, to);

        t = to;
        to = from;
        from = t;

        *ja |= bit_mask;
        ja = perm + k - offset;
      }

      elem_copy(sx, to);

      *ja |= bit_mask;
    }

    *ia &= inv_bit_mask;

    ia++;
    i++;
    elem_inc(&src);
  }

  return 0;
}


slint pepckeys_sort_permute_backward(elements_t *s, elements_t *sx, slint *perm, slint offset, slint mask_bit) /* pepckeys_sl_proto, sl_func pepckeys_sort_permute_backward */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_permute_backward);

  if (sx == NULL)
  {
    sx = &_sx;
    pepckeys_elements_alloc(sx, 1, 1, 1);
  }

  if (mask_bit < 0) pepckeys_sort_permute_backward_(s, sx, perm, offset);
  else pepckeys_sort_permute_backward_masked(s, sx, perm, offset, 1L << xmin(mask_bit, (sizeof(slint) * 8) - 1));

  if (sx == &_sx) pepckeys_elements_free(sx);

  rti_tstop(rti_tid_sort_permute_backward);

  return 0;
}
