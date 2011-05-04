/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/mergep_2way.c
 *  timestamp: 2011-02-10 21:20:49 +0100
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

  Z_TRACE_ARRAY_IF(MP2W_TRACE_IF, "displs =", " %d", i, p, displs);
  Z_TRACE_ARRAY_IF(MP2W_TRACE_IF, "counts =", " %" slint_fmt, i, p, counts);
  
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
