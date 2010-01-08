/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_splitk_dummy.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* same as a normal mpi_split, but send_buffer-counts never get full, therefore elements are only copied locally, no communication */


#include "sl_common.h"


#define K2C_ONLY


slint pepcparts_mpi_splitk_dummy(elements_t *s, k2c_func k2c, void *ci, elements_t *sx, slint *send_stats, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_splitk_dummy */
{
  slint i, j, k, t;

  slint local_sb_counts[size];
  
  slint _send_stats[size];

  elements_t sb[size], sb_current[size];
  elements_t src, dst, end;

  
  if (s == NULL || ci == NULL || sx == NULL) return -1;

  /* need send_buffers with at least one element per foreign process */
  if (sx->size < size - 1) return -2;

  rti_tstart(rti_tid_mpi_splitk_dummy);
  rti_tstart(rti_tid_mpi_splitk_dummy_init);

  if (send_stats == NULL) send_stats = _send_stats;

  /* initials */
  j = sx->size;
  k = size - 1;
  for (i = 0; i < size; ++i)
  {
    /* init the local send_buffer counters */
    local_sb_counts[i] = 0;

    /* prepare the send_buffers */
    if (i != rank)
    {
      elem_assign_at(sx, sx->size - j, &sb[i]);
      sb[i].size = (j / k) + (j % k != 0);
      j -= sb[i].size;
      --k;

    } else elem_null(&sb[i]);
    elem_assign(&sb[i], &sb_current[i]);
    
    send_stats[i] = 0;
  }

  elem_assign(s, &src);
  elem_assign(s, &dst);
  elem_assign_at(s, s->size, &end);

  rti_tstop(rti_tid_mpi_splitk_dummy_init);
  rti_tstart(rti_tid_mpi_splitk_dummy_loop);

  while (1)
  {
    /* distribute the elements to the send_buffer, as long as possible (elements left and target send_buffer not full) */
    while (src.keys != end.keys)
    {
      /* compute the target-process of the current element */
      t = (k2c)(src.keys, src.keys - s->keys, ci);

      ++send_stats[t];

#ifndef K2C_ONLY

      /* is the local process the target? */
      if (t == rank)
      {
        /* if necessary, move the element on the local process */
        if (src.keys != dst.keys) elem_copy(&src, &dst);

        /* update the dst-position */
        elem_inc(&dst);

      } else /* the target is another process (need to send the element) */
      {
        /* break, if the according send_buffer is full */
        if (local_sb_counts[t] >= sb[t].size) break;

        /* copy the element to the according send_buffer */
        elem_copy(&src, &sb_current[t]);
        elem_inc(&sb_current[t]);

        ++local_sb_counts[t];
        
        if (local_sb_counts[t] >= sb[t].size)
        {
          elem_sub(&sb_current[t], local_sb_counts[t]);
          local_sb_counts[t] = 0;
        }
      }

#endif
      
      /* update the src-position */
      elem_inc(&src);
    }

    break;
  }

  rti_tstop(rti_tid_mpi_splitk_dummy_loop);
  rti_tstop(rti_tid_mpi_splitk_dummy);
  
  return 0;
}
