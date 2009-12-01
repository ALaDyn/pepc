/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_sort_merge.c
 *  timestamp: 2009-11-20 11:30:44 +0100
 *  
 */


#include "sl_common.h"


slint_t pepcparts_mpi_sort_merge(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_sort_merge */
{
  merge2x_f m2x = pepcparts_merge2_memory_adaptive;
  sortnet_f sn = pepcparts_sn_batcher;

#ifdef key_integer
  pepcparts_sort_radix(s0, xs, -1, -1, -1);
#else
  pepcparts_sort_quick(s0, xs);
#endif

  pepcparts_mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

  return 0;
}
