/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_sample_select_qs.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* */

#include "sl_common.h"


slint pepckeys_mpi_sample_select_qs(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_sample_select_qs */
{
  slint i;
  slint all_local_sizes[size], iths[size - 1];

  rti_tstart(rti_tid_mpi_sample_select_qs);
  rti_tstart(rti_tid_mpi_sample_select_qs_pre);

  /* prepare the classification_info */
  ci->nclasses = size;
  ci->keys = sl_alloc(size - 1, sizeof(slkey_pure_t));
  ci->counts = NULL;
  ci->masks = NULL;

  if (size > 1)
  {
    /* distribute the local sizes */
    MPI_Allgather(&s->size, 1, sl_int_type_mpi, all_local_sizes, 1, sl_int_type_mpi, comm);

    /* computing the desired ith values (prefix-sums) */
    iths[0] = all_local_sizes[0];
    for (i = 1; i < size - 1; ++i) iths[i] = iths[i - 1] + all_local_sizes[i];
  }

  rti_tstop(rti_tid_mpi_sample_select_qs_pre);
  rti_tstart(rti_tid_mpi_sample_select_qs_select);
  
  pepckeys_mpi_select_qs(s, size - 1, iths, pepckeys_pivot_random, threshold, ci->keys, size, rank, comm);

  rti_tstop(rti_tid_mpi_sample_select_qs_select);
  rti_tstop(rti_tid_mpi_sample_select_qs);

  return 0;
}
