/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sl_elem.c
 *  timestamp: 2010-12-19 13:24:41 +0100
 *  
 */


#include "sl_common.h"

#include <stdarg.h>


slint_t elem_set_data(elements_t *e, ...) /* sl_proto, sl_func elem_set_data */
{
  va_list ap;

  va_start(ap, e);

#define xelem_call_data  xelem_buf(e) = va_arg(ap, xelem_sltype_t *);
#include "sl_xelem_call.h"

  va_end(ap);

  return 0;
}


inline slint_t elem_reverse_aio(elements_t *e, elements_t *t) /* sl_func elem_reverse_aio */
{
  elements_t front, back, end;

  elem_assign(e, &front);
  elem_assign_at(e, e->size - 1, &back);
  elem_assign_at(e, e->size / 2, &end);

  while (front.keys < end.keys)
  {
    elem_xchange(&front, &back, t);
    elem_inc(&front);
    elem_dec(&back);
  }

  return 0;
}


inline slint_t elem_reverse_obo(elements_t *e, elements_t *t) /* sl_func elem_reverse_obo */
{
  elements_t front, back, end;

  elem_assign(e, &front);
  elem_assign_at(e, e->size - 1, &back);
  elem_assign_at(e, e->size / 2, &end);

#define xelem_call \
  while (xelem_buf(&front) > xelem_buf(&end)) \
  { \
    xelem_xchange(&front, &back, t); \
    xelem_inc(&front); \
    xelem_dec(&back); \
  }
#include "sl_xelem_call.h"

  return 0;
}


slint_t elem_reverse(elements_t *e, elements_t *t) /* sl_proto, sl_func elem_reverse */
{
  return elem_reverse_aio(e, t);
}


slint_t elem_nxchange_at(elements_t *e0, slint_t at0, elements_t *e1, slint_t at1, slint_t n, elements_t *t) /* sl_proto, sl_func elem_nxchange_at */
{
  elements_t _e0, _e1, end;

  elem_assign_at(e0, at0, &_e0);
  elem_assign_at(e1, at1, &_e1);
  elem_assign_at(&_e0, n, &end);

  while (_e0.keys < end.keys)
  {
    elem_xchange(&_e0, &_e1, t);
    elem_inc(&_e0);
    elem_inc(&_e1);
  }

  return 0;
}


slint_t elem_nxchange(elements_t *e0, elements_t *e1, slint_t n, elements_t *t) /* sl_proto, sl_func elem_nxchange */
{
  return elem_nxchange_at(e0, 0, e1, 0, n, t);
}


slint_t elem_nxchange_ro0(elements_t *e0, elements_t *e1, slint_t n, elements_t *t) /* sl_proto, sl_func elem_nxchange_ro0 */
{
  elements_t _e0, _e1, end;

  elem_assign(e0, &_e0);
  elem_assign(e1, &_e1);
  elem_assign_at(e1, n, &end);

  elem_copy(&_e1, t); /* create the hole */
  elem_copy(&_e0, &_e1);
  elem_inc(&_e1);

  while (_e1.keys < end.keys)
  {
    elem_copy(&_e1, &_e0);
    elem_inc(&_e0);
    elem_copy(&_e0, &_e1);
    elem_inc(&_e1);
  }

  elem_copy(t, &_e0); /* close the hole */

  return 0;
}


/* the easy way, with 3x reverse */
inline slint_t elem_rotate_3rev_aio(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_3rev_aio */
{
  elements_t _e;

  if (m == 0 || n == 0) return 0;

  /* reverse 2nd part */
  elem_assign_at(e, m, &_e);
  _e.size = n;
  elem_reverse_aio(&_e, t);

  /* reverse 1st part */
  elem_assign(e, &_e);
  _e.size = m;
  elem_reverse_aio(&_e, t);

  /* reverse all */
  _e.size = m + n;
  elem_reverse_aio(&_e, t);

  return 0;
}


inline slint_t elem_rotate_3rev_obo(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_3rev_obo */
{
  elements_t _e;

  if (m == 0 || n == 0) return 0;

  /* reverse 2nd part */
  elem_assign_at(e, m, &_e);
  _e.size = n;
  elem_reverse_obo(&_e, t);

  /* reverse 1st part */
  elem_assign(e, &_e);
  _e.size = m;
  elem_reverse_obo(&_e, t);

  /* reverse all */
  _e.size = m + n;
  elem_reverse_obo(&_e, t);

  return 0;
}


inline slint_t elem_rotate_cycles_aio(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_cycles_aio */
{
  slint_t k;
  elements_t start, half, current, next;

  if (m == 0 || n == 0) return 0;

  elem_assign(e, &start);
  elem_assign_at(e, n, &half);

  k = m + n;
  while (k > 0)
  {
    elem_copy(&start, t);
    elem_assign(&start, &current);

    while (1)
    {
      if (current.keys < half.keys) elem_assign_at(&current, m, &next);
      else elem_assign_at(&current, -n, &next);

      if (next.keys == start.keys) break;

      elem_copy(&next, &current);
      elem_assign(&next, &current);
      --k;
    }

    elem_copy(t, &current);
    --k;

    elem_inc(&start);
  }

  return 0;
}


/* FIXME: hier fehlt noch etwas vom xcall! */
inline slint_t elem_rotate_cycles_obo(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_func elem_rotate_cycles_obo */
{
/*  slint_t k;*/
  elements_t start, half/*, current, next*/;

  if (m == 0 || n == 0) return 0;

  elem_assign(e, &start);
  elem_assign_at(e, n, &half);

#define xelem_call \
  k = m + n; \
  while (k > 0) \
  { \
    xelem_copy(&start, t); \
    xelem_assign(&start, &current); \
\
    while (1) \
    { \
      if (xelem_buf(&current) < xelem_buf(&half)) xelem_assign_at(&current, m, &next); \
      else xelem_assign_at(&current, -n, &next); \
\
      if (xelem_buf(&next) == xelem_buf(&start)) break; \
\
      xelem_copy(&next, &current); \
      xelem_assign(&next, &current); \
      --k; \
    } \
\
    xelem_copy(t, &current); \
    --k; \
\
    xelem_inc(&start); \
  }

  return 0;
}


slint_t elem_rotate(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_proto, sl_func elem_rotate */
{
  return elem_rotate_3rev_aio(e, m, n, t);
}


/* retain order of the 1st part, do it back2front (for the sake of simplicity) */
slint_t elem_rotate_ro0(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_proto, sl_func elem_rotate_ro0 */
{
  elements_t e0, e1;

  elem_assign_at(e, m, &e0);
  elem_assign_at(e, m + n, &e1);

  elem_copy(&e1, t);
  elem_copy(&e0, &e1);

  while (e0.keys > e->keys)
  {
    elem_dec(&e1);
    elem_copy(&e1, &e0);
    elem_dec(&e0);
    elem_copy(&e0, &e1);
  }

  elem_copy(t, &e0);

  return 0;
}


/* retain order of the 2nd part, do it front2back (for the sake of simplicity) */
slint_t elem_rotate_ro1(elements_t *e, slint_t m, slint_t n, elements_t *t) /* sl_proto, sl_func elem_rotate_ro1 */
{
  elements_t e0, e1, end;

  elem_assign(e, &e0);
  elem_assign_at(e, m, &e1);
  elem_assign_at(e, m + n, &end);

  elem_copy(&e0, t);
  elem_copy(&e1, &e0);

  while (e1.keys < end.keys)
  {
    elem_inc(&e0);
    elem_copy(&e0, &e1);
    elem_inc(&e1);
    elem_copy(&e1, &e0);
  }

  elem_copy(t, &e1);

  return 0;
}
