/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_sort_merge.c
 *  timestamp: 2009-11-20 11:30:44 +0100
 *  
 */


#include "sl_common.h"


slint_t pepckeys_mpi_sort_merge(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_sort_merge */
{
  merge2x_f m2x = pepckeys_merge2_memory_adaptive;
  sortnet_f sn = pepckeys_sn_batcher;

#ifdef key_integer
  pepckeys_sort_radix(s0, xs, -1, -1, -1);
#else
  pepckeys_sort_quick(s0, xs);
#endif

  pepckeys_mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

  return 0;
}
