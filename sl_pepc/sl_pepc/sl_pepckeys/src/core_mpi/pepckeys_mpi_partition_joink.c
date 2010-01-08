/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_partition_joink.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* cycle nicht immer zu erst bei Prozess 0 suchen?
   - 1. nach und nach bei anderem Startknoten beginnen (mit offset+modulo oder Permutation)
   - 2. zufällige Reihenfolge der Knoten (Permutation)
   
   Vorteile:
   - kein sendbuffer
   - kein alltoall für sendbuffer-sizes
   - gegenüber pepckeys_mpi_splitk.c klare Aufteilung in local und Komm.
   - schneller?
   - Verwendung für PSES? mit abschliessendem lokalen sortieren anstatt p-way merge, Schlüsselverteilung (min, max) für radixsort verwenden
   - space (1) und time O(n) optimal p-way merge bekannt? nö -> mit 2-way merge und sorting networks
*/

#include "sl_common.h"


#define white -1
#define black -2

slint pepckeys_dfs_visit_dcrbe(slint v, slint *max_flow, slint *pred, slint *m, slint n, slint *complete, slint *sendrecv_size, slint *send_rank, slint *recv_rank, int rank) /* sl_func pepckeys_dfs_visit_dcrbe */
{
  slint i, j, k, f;
  
  for (i = 0; i < n; i++)
  if ((k = m[v * n + i]) > 0)
  {
    if (pred[i] == white) /* not visited, yet */
    {
      pred[i] = v;  /* mark as grey */

      /* visit */
      f = pepckeys_dfs_visit_dcrbe(i, max_flow, pred, m, n, complete, sendrecv_size, send_rank, recv_rank, rank);

      pred[i] = black;  /* mark as black */

      /* if the visit found a back-edge? */
      if (f >= 0)
      {
        *complete -= *max_flow;
        m[v * n + i] -= *max_flow;

        if (rank == v)
        {
          *sendrecv_size = *max_flow;
          *send_rank = i;
        }

        if (rank == i) *recv_rank = v;

        /* if the back-edge was ending at me? */
        if (f == v) return -1;

        return f;
      }

    } else if (pred[i] != black) /* not white, not black -> it's grey, it's a back edge */
    {
      /* compute the max flow (checking all edges of the located cycle) */
      *max_flow = k;
      j = v;

      while (j != i)
      {
        *max_flow = xmin(*max_flow, m[pred[j] * n + j]);
        j = pred[j];
      }

      *complete -= *max_flow;
      m[v * n + i] -= *max_flow;

      if (rank == v)
      {
        *sendrecv_size = *max_flow;
        *send_rank = i;
      }

      /* set the source for the first process of this cycle */
      if (rank == i) *recv_rank = v;

      return i; /* signal which grey vertex was found */
    }
  }

  return -1;
}


slint pepckeys_disjoint_cycle_reduce_back_edges(slint *m, slint n, slint *complete, slint *sendrecv_size, slint *send_rank, slint *recv_rank, int rank) /* sl_func pepckeys_disjoint_cycle_reduce_back_edges */
{
  slint i, pred[n], flow;

  /* initialize predecessors */
  for (i = 0; i < n; ++i) pred[i] = white;

  for (i = 0; i < n; i++)
  if (pred[i] == white)
  {
    pred[i] = i;
    
    pepckeys_dfs_visit_dcrbe(i, &flow, pred, m, n, complete, sendrecv_size, send_rank, recv_rank, rank);
    
    pred[i] = black;
  }

  return 0;
}


slint pepckeys_mpi_partition_joink(elements_t *s, slint *sizes, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_partition_joink */
{
  slint i, j, k, l, m, n;

  slint *local_sizes = sizes, all_local_sizes[size * size];

  slint sendrecv_size, send_rank, recv_rank;

  slint complete, former_complete;

  elements_t parts[size];

  MPI_Status status;

/*  printf("%d here: size = %d\n", rank, size);*/

  if (size < 2) return 0;

/*  MPI_Barrier(comm);*/

  rti_treset(rti_tid_mpi_partition_joink_loop_flow);
  rti_treset(rti_tid_mpi_partition_joink_loop_dist);

  rti_tstart(rti_tid_mpi_partition_joink);
  rti_tstart(rti_tid_mpi_partition_joink_init);

  /* initials */
  elem_assign(s, &parts[0]);
  for (i = 1; i < size; i++) elem_assign_at(&parts[i - 1], sizes[i - 1], &parts[i]);

  /* distribute the local_sizes */
  MPI_Allgather(local_sizes, size, sl_int_type_mpi, all_local_sizes, size, sl_int_type_mpi, comm);

  complete = 0;
  for (i = 0; i < size; ++i)
  {
    all_local_sizes[i * size + i] = 0;
    for (j = i + 1; j < size; ++j) complete += all_local_sizes[i * size + j] + all_local_sizes[j * size + i];
  }

  rti_tstop(rti_tid_mpi_partition_joink_init);
  rti_tstart(rti_tid_mpi_partition_joink_loop);

/*  MPI_Barrier(comm);*/

  /* detect cycles and distribute them until finished */
  while (complete > 0)
  {
/*    if (rank == 0)
    {
      for (i = 0; i < size; i++)
      {
        for (j = 0; j < size; j++) printf("  %d", all_local_sizes[i * size + j]);
        printf("\n");
      }
      printf("%d here: complete = %d\n", rank, complete);
    }*/

    rti_tstart(rti_tid_mpi_partition_joink_loop_flow);

    sendrecv_size = 0;
    send_rank = recv_rank = -1;
    former_complete = complete;

    pepckeys_disjoint_cycle_reduce_back_edges(all_local_sizes, size, &complete, &sendrecv_size, &send_rank, &recv_rank, rank);
    
    rti_tstop(rti_tid_mpi_partition_joink_loop_flow);

    if (complete == former_complete)
    {
      printf("%d here: DEADLOCK, complete is unchanged %" sl_int_type_fmt "\n", rank, complete);
      break;
    }
    
/*    printf("%d here: complete = %d, sendrecv_size = %d, send_rank = %d, recv_rank = %d\n", rank, complete, sendrecv_size, send_rank ,recv_rank);
    printf("%d here: xchanging keys @ %p\n", rank, parts[send_rank].keys);*/

    if (sendrecv_size > 0)
    {
      rti_tstart(rti_tid_mpi_partition_joink_loop_dist);

#define xelem_call \
      MPI_Sendrecv_replace(xelem_buf(&parts[send_rank]), sendrecv_size, xelem_mpi_datatype, send_rank, rank, recv_rank, recv_rank, comm, &status);
#include "sl_xelem_call.h"

      elem_add(&parts[send_rank], sendrecv_size);

      rti_tstop(rti_tid_mpi_partition_joink_loop_dist);
    }

  }

  rti_tstop(rti_tid_mpi_partition_joink_loop);
  rti_tstop(rti_tid_mpi_partition_joink);

  return 0;
}
