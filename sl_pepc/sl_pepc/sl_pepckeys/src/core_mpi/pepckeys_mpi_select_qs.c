/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_select_qs.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* - solving "n" times the "selection problem" in parallel for non-distinct elements distributed on "size" processes
   - based on a quicksort like algorithm, called quickselect
   - expected linear time
   - uses a given function to detect pivot-elements
   - works with non-distinct elements
   
   - anstatt split3 mit split2 machen, aber wie dann den pivot extrahieren?
*/

#include "sl_common.h"


slint pepckeys_mpi_select_qs(elements_t *s, slint n, slint *iths, pivot_func pi, slint threshold, slkey_pure_t *e, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_select_qs */
{
  const slint first_pivot_root = 0;

  slint unfinished = n;

  slint i, j, p, nlt, nle;

  elements_t areas[n], former_area, t;

  slint local_partition_sizes[2 * n], global_partition_sizes[2 * n], all_local_partition_sizes[size * 2 * n];
  
/*  char outputname[128];
  FILE *output;
  
  sprintf(outputname, "out_%d", rank);
  output = fopen(outputname, "w");*/

  slkey_pure_t *pivots = e;
  slkey_pure_t pivot_proposals[n];
  slkey_pure_t all_pivot_proposals[size * n];
  
  if (n < 1) return 0;

  rti_treset(rti_tid_mpi_select_qs_part);
  rti_treset(rti_tid_mpi_select_qs_reduce_sizes);
  rti_treset(rti_tid_mpi_select_qs_area);
  rti_treset(rti_tid_mpi_select_qs_pivot_new);
  rti_treset(rti_tid_mpi_select_qs_pivot_gather);
  rti_treset(rti_tid_mpi_select_qs_pivot_detect);

  rti_tstart(rti_tid_mpi_select_qs);
  rti_tstart(rti_tid_mpi_select_qs_pre);

  pepckeys_elements_alloc(&t, 1, 1, 1);

  /* initialize the areas */
  for (i = 0; i < n; ++i) elem_assign(s, &areas[i]);

  /* compute & distribute the first pivot */
  if (rank == first_pivot_root) pivots[0] = key_purify(s->keys[(pi)(s)]);
  MPI_Bcast(&pivots[0], 1, pepckeys_key_pure_mpi_datatype, first_pivot_root, comm);
  for (i = 1; i < n; ++i) pivots[i] = pivots[0];
  
/*  pepckeys_elements_printf(s);
  pepckeys_elements_print_keys(s);*/

  rti_tstop(rti_tid_mpi_select_qs_pre);
  rti_tstart(rti_tid_mpi_select_qs_loop);

  while (1)
  {
    /* assume we have valid pivots */

    /* iths[i] < 0 signals that the requested element was found and resides in pivots[i] */

    rti_tstart(rti_tid_mpi_select_qs_part);
    
/*    if (rank == 0)
    for (i = 0; i < n; ++i)
    {
      printf("%d here: iths[%d] = %d, pivots[%d] = %d\n", rank, i, iths[i], i, pivots[i]);
    }*/
    
/*    fprintf(output, "round %d\ns = [", rounds);
    for (i = 0; i < s->size; ++i) fprintf(output, "  %d", s->keys[i]);
    fprintf(output, "]\n");*/

    /* partition the (max.) size-1 areas with the new pivots */
    elem_null(&former_area);
    for (i = 0; i < n; ++i)
    if (iths[i] >= 0) /* skip finished */
    {
      /* if the new area is unpartitioned, partition now */
      if (areas[i].keys != former_area.keys)
      {
/*        fprintf(output, "%d here: before partition with %d\n", rank, pivots[i]);*/
/*        if (rank == 0)
        {
          printf("%d here: before partition with %d\n", rank, pivots[i]);
          pepckeys_elements_print_keys(&areas[i]);
        }*/

        pepckeys_split3_lt_eq_gt(&areas[i], &pivots[i], &t, &nlt, &nle);

        elem_assign(&areas[i], &former_area);

/*        if (rank == 0)
        {
          printf("%d here: after partition retuned %d & %d\n", rank, nlt, nle);
          pepckeys_elements_print_keys(&areas[i]);
        }*/

/*        fprintf(output, "%d here: after partition retuned %d & %d\n", rank, nlt, nle);*/
      }

      local_partition_sizes[2 * i] = nlt;
      local_partition_sizes[2 * i + 1] = nle;
    }

/*    if (rank == 0)
    for (i = 0; i < size - 1; ++i) printf("%d here: local_partition_sizes[%d]: %d & %d\n", rank, i, local_partition_sizes[2 * i], local_partition_sizes[2 * i + 1]);*/

    rti_tstop(rti_tid_mpi_select_qs_part);
    rti_tstart(rti_tid_mpi_select_qs_reduce_sizes);

    /* distribute the sizes of the new partitions */

    /* distribution and sumation at once using reduction */
/*    MPI_Allreduce(local_partition_sizes, global_partition_sizes, 2 * n, sl_int_type_mpi, MPI_SUM, comm);*/

    /* distribution and local sumation */
    MPI_Allgather(local_partition_sizes, 2 * n, sl_int_type_mpi, all_local_partition_sizes, 2 * n, sl_int_type_mpi, comm);
    for (i = 0; i < n; i++)
    {
      global_partition_sizes[2 * i] = global_partition_sizes[2 * i + 1] = 0;

      for (j = 0; j < size; j++)
      {
        global_partition_sizes[2 * i] += all_local_partition_sizes[j * 2 * n + (2 * i)];
        global_partition_sizes[2 * i + 1] += all_local_partition_sizes[j * 2 * n + (2 * i) + 1];
      }
    }

/*    if (rank == 0)
    for (i = 0; i < n; ++i) printf("%d here: global_partition_sizes[%d]: %d & %d\n", rank, i, global_partition_sizes[2 * i], global_partition_sizes[2 * i + 1]);*/

    rti_tstop(rti_tid_mpi_select_qs_reduce_sizes);
    rti_tstart(rti_tid_mpi_select_qs_area);

    /* select new iths and areas */
    for (i = 0; i < n; ++i)
    if (iths[i] >= 0) /* skip finished */
    {
/*      if (rank == 0) printf("%d here: iths[%d] = %d - global_partition_size[%d] = %d & %d\n", rank, i, iths[i], i, global_partition_sizes[2 * i], global_partition_sizes[2 * i + 1]);*/

      if (iths[i] < global_partition_sizes[2 * i])
      {
/*        if (rank == 0) printf("%d here: %d < %d -> lower\n", rank, iths[i], global_partition_sizes[2 * i]);*/

        /* new partition is the lower one */
        elem_add(&areas[i], local_partition_sizes[2 * i + 1] - local_partition_sizes[2 * i]);
        areas[i].size = local_partition_sizes[2 * i];

      } else if (iths[i] >= global_partition_sizes[2 * i + 1])
      {
/*        if (rank == 0) printf("%d here: %d >= %d -> higher\n", rank, iths[i], global_partition_sizes[2 * i + 1]);*/

        /* new partition is the higher one */
        elem_add(&areas[i], local_partition_sizes[2 * i + 1]);
        areas[i].size -= local_partition_sizes[2 * i + 1];

        iths[i] -= global_partition_sizes[2 * i + 1];

      } else
      {
/*        if (rank == 0) printf("%d here: hit it\n", rank);*/

        /* pivots[i] is the requested one */
        iths[i] = -1;
        --unfinished;
      }
    }

    rti_tstop(rti_tid_mpi_select_qs_area);

/*    if (rank == 0) printf("%d: unfinsished = %d\n", rounds, unfinished);*/

    if (unfinished <= 0) break;

    rti_tstart(rti_tid_mpi_select_qs_pivot_new);

    /* detect pivot_proposals, use same pivot_proposals for same areas, ignore (= use former pivot) if iths[i] < 0 or areas[i].size <= 0 */
    for (i = 0; i < n; ++i)
    if (iths[i] >= 0 && areas[i].size > 0)
    {
      if (i > 0)
      if (iths[i - 1] >= 0 && areas[i - 1].size > 0 && areas[i].keys == areas[i - 1].keys)
      {
        pivot_proposals[i] = pivot_proposals[i - 1];
        continue;
      }

      /* chose a pivot proposal from the area */
      p = (pi)(&areas[i]);
      pivot_proposals[i] = key_purify(areas[i].keys[p]);
      
      /* position this pivot proposal element at the front of the area */
/*      elem_copy(&areas[i], &t);
      elem_copy_at(&areas[i], p, &areas[i], 0);
      elem_copy_at(&t, 0, &areas[i], p);*/

    } else pivot_proposals[i] = pivots[i]; /* chose the former pivot */

    rti_tstop(rti_tid_mpi_select_qs_pivot_new);
    rti_tstart(rti_tid_mpi_select_qs_pivot_gather);

    /* distribute pivot_proposals */
    MPI_Allgather(pivot_proposals, n, pepckeys_key_pure_mpi_datatype, all_pivot_proposals, n, pepckeys_key_pure_mpi_datatype, comm);

    rti_tstop(rti_tid_mpi_select_qs_pivot_gather);
    rti_tstart(rti_tid_mpi_select_qs_pivot_detect);

    /* select new pivots, pivot_proposals are equal to the former pivot, if
       (1) no more elements left to chose in a local partition
       (2) the same pivot was chosen accidently (if possible, using split3 this is impossible)
       -> if all proposals match the former pivot, at least one process chose this pivot */
    for (i = 0; i < n; ++i)
    if (iths[i] >= 0)
    {
      /* BETTER: chose median of pivot_proposals? */

      j = 0;
      while (j < size)
      if (key_pure_cmp_ne(pivots[i], all_pivot_proposals[j * n + i])) /* FORCE to chose a different pivot as before */
      {
        /* accept this proposal */
        pivots[i] = all_pivot_proposals[j * n + i];

        break;

      } else ++j;
    }

/*    if (rank == 0)
    for (i = 0; i < n; ++i)
    {
      printf("pivot_proposals[%d] = [", i);
      for (j = 0; j < size; ++j) printf("  %d", all_pivot_proposals[j * n + i]);
      printf("] my %d -> all %d\n", pivot_proposals[i], pivots[i]);
    }*/

    rti_tstop(rti_tid_mpi_select_qs_pivot_detect);
  }
  rti_tstop(rti_tid_mpi_select_qs_loop);

  pepckeys_elements_free(&t);

  rti_tstop(rti_tid_mpi_select_qs);
  
/*  fclose(output);*/
  
  return 0;
}
