/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_sample_complete.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* - very intensiv way
   - mpi_gather key-exchange + local sorting + split-key-detection + broadcast
   - just for the sake of simplicity and for producing a valid classification
   - assuming a communicator-wide constant threshold (when used)
*/

#include "sl_common.h"


#define the_sort_di(s, sx)  sort_radix_di(s, sx, -1, -1, -1)


slint pepckeys_mpi_sample_complete(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_sample_complete */
{
  const slint root = 0;

  slint i, j, global_size;

  int local_size, recv_counts[size], recv_displs[size];
  elements_t all_keys, e, t;

  rti_tstart(rti_tid_mpi_sample_complete);

  /* prepare the ci structure */
  ci->nclasses = size;
  ci->keys = sl_alloc(size - 1, sizeof(slkey_pure_t));
  ci->counts = NULL;
  ci->masks = NULL;

  rti_tstart(rti_tid_mpi_sample_complete_gather);

  /* collect the local sizes */
  local_size = s->size;
  MPI_Gather(&local_size, 1, MPI_INT, recv_counts, 1, MPI_INT, root, comm);

  /* root only */
  if (rank == root)
  {
    /* compute recv_counts/displs */
    recv_displs[0] = 0;
    for (i = 1; i < size; ++i)
    {
      global_size += recv_counts[i];
      recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
    }
    global_size = recv_displs[size - 1] + recv_counts[size - 1];

    /* alloc space for keys */
    elements_alloc_di(&all_keys, global_size, 1, 0);
  }

  /* gatherv the keys at root */
  MPI_Gatherv(s->keys, s->size, pepckeys_key_mpi_datatype, all_keys.keys, recv_counts, recv_displs, pepckeys_key_mpi_datatype, root, comm);

  rti_tstop(rti_tid_mpi_sample_complete_gather);
  rti_tstart(rti_tid_mpi_sample_complete_detect);

  if (rank == root)
  {
    /* local sort the keys */
    the_sort_di(&all_keys, NULL);

    /* detect split-keys (using threshold?) */
    for (i = 0; i < size - 1; ++i) ci->keys[i] = key_purify(all_keys.keys[recv_displs[i + 1] / key_size_mpi]);

/*    pepckeys_elements_print_keys(&all_keys);*/
  }

  rti_tstop(rti_tid_mpi_sample_complete_detect);
  rti_tstart(rti_tid_mpi_sample_complete_bcast);

  /* distribute the split-keys */
  MPI_Bcast(ci->keys, size - 1, pepckeys_key_pure_mpi_datatype, root, comm);

  rti_tstop(rti_tid_mpi_sample_complete_bcast);
  rti_tstop(rti_tid_mpi_sample_complete);

  return 0;
}
