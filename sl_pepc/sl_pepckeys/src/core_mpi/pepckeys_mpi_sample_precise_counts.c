/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_sample_precise_counts.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* - präzisiert/ermittelt die ci->counts für gegebene ci->keys
   - nutzt eventuell schon vorhandene Informationen aus ci
*/


#include "sl_common.h"


slint pepckeys_mpi_sample_precise_counts(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_sample_precise_counts */
{
  slint i, j, k, l;
  elements_t _s, _e;
  
  slkey_t *ks, *ke;
  
  slint n = ci->nclasses - 1;

  slint all_local_sizes[size];
  slint *als = all_local_sizes;

  slint local_lt_eq_counts[2 * n], all_local_lt_eq_counts[size * 2 * n];
  slint *llec = local_lt_eq_counts, *allec = all_local_lt_eq_counts;
  
  slint prefix_sums[size];

  /* prepare ci */
  ci->counts = sl_alloc(n, sizeof(slint));

  if (ci->nclasses < 2) return 0;

  rti_tstart(rti_tid_mpi_sample_precise);
  rti_tstart(rti_tid_mpi_sample_precise_llec);

  /* prepare to compute the local_lt_eq counts */
  for (i = 0; i < n; ++i) llec[i * 2] = llec[i * 2 + 1] = 0;
  
  /* every process computes the necessary numbers of keys less than and equal to the split-keys */
  key_assign(s->keys, ks);
  key_assign_at(s->keys, s->size, ke);
  while (ks != ke)
  {
    for (i = 0; i < n; ++i)
    if (key_pure_cmp_le(key_purify(*ks), ci->keys[i]))
    {
      if (key_pure_cmp_eq(key_purify(*ks), ci->keys[i])) ++llec[2 * i + 1];
      else ++llec[2 * i];

      break;
    }
    
    key_inc(ks);
  }

  rti_tstop(rti_tid_mpi_sample_precise_llec);
  rti_tstart(rti_tid_mpi_sample_precise_gather);

  /* distribute or reuse the local sizes */
  if (ci->all_local_sizes) als = ci->all_local_sizes;
  else MPI_Allgather(&s->size, 1, sl_int_type_mpi, all_local_sizes, 1, sl_int_type_mpi, comm);

  /* distribute the local_lt_eq numbers */
  MPI_Allgather(llec, 2 * n, sl_int_type_mpi, allec, 2 * n, sl_int_type_mpi, comm);

  rti_tstop(rti_tid_mpi_sample_precise_gather);
  rti_tstart(rti_tid_mpi_sample_precise_detect);

  /* compute the prefix_sums from all_local_sizes */
  prefix_sums[0] = 0;
  for (i = 1; i < size; ++i) prefix_sums[i] = prefix_sums[i - 1] + als[i - 1];

  /* compute the precise counts */
  k = 0;
  for (i = 0; i < n; ++i)
  {
    for (j = 0; j < size; ++j) k += allec[j * 2 * n + (2 * i)];

    for (j = 0; j < size; ++j)
    {
      if (i > 0 && key_pure_cmp_eq(ci->keys[i - 1], ci->keys[i])) allec[j * 2 * n + (2 * i) + 1] = allec[j * 2 * n + (2 * (i - 1)) + 1];
      
      l = xmax(0, xmin(allec[j * 2 * n + (2 * i) + 1], prefix_sums[i + 1] - k));

      if (rank == j) ci->counts[i] = l;

      k += allec[j * 2 * n + (2 * i) + 1];
    }
  }

  rti_tstop(rti_tid_mpi_sample_precise_detect);
  rti_tstop(rti_tid_mpi_sample_precise);

  return 0;
}
