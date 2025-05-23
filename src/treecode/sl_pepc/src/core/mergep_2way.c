/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/mergep_2way.c
 *  
 */


/* sl_macro MP2W_TRACE_IF */


#include "sl_common.h"


#ifndef MP2W_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MP2W_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MP2W_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mergep_2way_ip_int(elements_t *s, elements_t *sx, slint_t p, int *displs, merge2x_f m2x) /* sl_proto, sl_func mergep_2way_ip_int */
{
  slint_t i, step, counts[p];
  elements_t s0, s1;


  Z_TRACE_IF(MP2W_TRACE_IF, "merging %" slint_fmt " sub-sequences", p);

  for (i = 0; i < p - 1; ++i) counts[i] = displs[i + 1] - displs[i];
  counts[p - 1] = s->size - displs[p - 1];

  Z_TRACE_ARRAY_IF(MP2W_TRACE_IF, i, p, " %d", displs[i], "displs =");
  Z_TRACE_ARRAY_IF(MP2W_TRACE_IF, i, p, " %" slint_fmt, counts[i], "counts =");
  
  for (step = 1; step < p; step *= 2)
  {
    for (i = 0; i < p; i += 2 * step)
    {
      if (i + step < p)
      {
        elem_assign_at(s, displs[i], &s0);
        s0.size = counts[i];
        elem_assign_at(s, displs[i + step], &s1);
        s1.size = counts[i + step];
        
        Z_TRACE_IF(MP2W_TRACE_IF, "merging %" slint_fmt " and %" slint_fmt, s0.size, s1.size);

        if (s0.size > 0 && s1.size > 0) m2x(&s0, &s1, sx);
        
        counts[i] = s0.size + s1.size;
      }
    }
  }
  
  return 0;
}


static inline void rec(elements_t *s, elements_t *sx, slint_t p, int offset, int *displs, merge2x_f m2x)
{
  slint_t p0, p1;
  elements_t s0, s1;

  p0 = p / 2;
  p1 = p - p0;

  elem_assign(s, &s0); s0.size = displs[p0] - offset;
  if (p0 > 1) rec(&s0, sx, p0, offset, displs, m2x);

  elem_assign_at(s, displs[p0] - offset, &s1); s1.size = s->size - (displs[p0] - offset);
  if (p1 > 1) rec(&s1, sx, p1, displs[p0], displs + p0, m2x);

  m2x(&s0, &s1, sx);
}


slint_t mergep_2way_ip_int_rec(elements_t *s, elements_t *sx, slint_t p, int *displs, merge2x_f m2x) /* sl_proto, sl_func mergep_2way_ip_int_rec */
{
  if (p > 1) rec(s, sx, p, 0, displs, m2x);
  
  return 0;
}
