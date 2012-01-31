/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_partition_exact_radix.c
 *  
 */


/* sl_macro MPER_TRACE_IF */


#include "sl_common.h"


#ifndef MPER_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPER_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPER_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_partition_exact_radix(elements_t *s, partcond_t *pcond, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_radix */
{
  binning_t bm;


  binning_radix_create(&bm, rhigh, rlow, rwidth, sorted|SL_SORTED_OUT);

  mpi_partition_exact_generic(s, pcond, &bm, scounts, rcounts, size, rank, comm);

  binning_radix_destroy(&bm);

  return 0;
}
