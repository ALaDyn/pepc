/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_partition_radix.c
 *  timestamp: 2009-11-23 13:31:49 +0100
 *  
 */

/* wenn unsortiert und modifikation erlaubt -> bei Verfeinerung split ausführen und entsprechenden Bereich als neue area verwenden*/
/* wenn unsortiert und modifikation nicht erlaubt -> bei Verfeinerung immer den kompletten Bereich als neue area verwenden */
/* Abbruch an einer Grenze sofort oder erst wenn alle fertig
   -> kompliziert: min/max borders müssen eingeschränkt werden mit Auswirkung auf alle anderen borders
*/

/* final: was passiert wenn alle Bits abgearbeitet sind
   - 0 = genau in der Mitte von min/max aufteilen
   - 1 = "beste" Prozessesgrenze
   - 2 = "beste" letzte Grenze
*/

#include "sl_common.h"


#define radix_low                 0
#define radix_high                (sizeof(slkey_pure_t) * 8 - 1)
#define radix_key2class(k, x, y)  (((k) >> (x)) & y)

#define myabs(x)  (((x) >= 0)?(x):-(x))


#ifdef SL_DATA_IGNORE
# undef elem_weight
# define elem_weight(e, at)  (1)
#endif

/* sl_macro MPI_PARTITION_RADIX_REDUCEBCAST_THRESHOLD */
#define REDUCEBCAST_ROOT  0

/*#define SCOUNTS_SDISPLS*/
/*#define WEIGHT_STATS*/
/*#define TIMING_STATS*/

#define SYNC_ON_INIT
#define HAVENT_MPI_IN_PLACE


slint_t pepcparts_mpi_partition_radix(elements_t *s, partcond_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_partition_radix */
{
  slkey_pure_t max_nclasses;
  slkey_pure_t nclasses, bit_mask;
  slkey_pure_t k;

  const slint_t max_nareas = size - 1;
  slint_t nareas, nareas_new;
  elements_t areas0[max_nareas], areas1[max_nareas], *areas, *areas_new;

  double *locals, *globals;
  double *local_counts, *local_weights, *global_counts, *global_weights;

  const slint_t max_nparts = size - 1;
  slint_t nparts, nparts_removed, remove_this_part;
  slint_t parts[max_nparts], part_areas[max_nparts];
  
  double final_locals[2 * max_nparts];
  double *final_local_counts, *final_local_weights;

  slint_t i, j, l;
  slint_t current_width;

  double minmax_weights[4 * size];
  
  slint_t last_new_area, last_new_class;

#ifdef HAVENT_MPI_IN_PLACE
  double local_minmax_weights[4];
#endif

  slint_t lc, lcs;
/*  slint_t gc, gcs;*/
  double lw, gw, lws, gws;
  double d;

  elements_t xi, end;

  slint_t rounds = 0;

  slint_t nbsets, bsets[size], bsets_parent[size];
  double bsets_minmax[4 * size], bsets_minmax_new[4];
  double current_borders[4], current_minmax[2];
  double border_offsets[size], border_finalize[size];

#ifdef WEIGHT_STATS
  slint_t total_count = 0, partial_counts[size + 1];
  double total_weight = 0.0, partial_weights[size + 1];
  double vmin, vmax;
# ifdef HAVENT_MPI_IN_PLACE
  slint_t partial_counts2[size + 1];
  double partial_weights2[size + 1];
# endif
#endif

  rti_treset(rti_tid_mpi_partition_radix_while);
  rti_treset(rti_tid_mpi_partition_radix_while_count);
  rti_treset(rti_tid_mpi_partition_radix_while_allreduce);
  rti_treset(rti_tid_mpi_partition_radix_while_round1);
  rti_treset(rti_tid_mpi_partition_radix_while_round1_allgather);
  rti_treset(rti_tid_mpi_partition_radix_while_check);
  rti_treset(rti_tid_mpi_partition_radix_final);

  rti_tstart(rti_tid_mpi_partition_radix_sync);
#ifdef SYNC_ON_INIT
  MPI_Barrier(comm);
#endif
  rti_tstop(rti_tid_mpi_partition_radix_sync);

  rti_tstart(rti_tid_mpi_partition_radix);

  if (rhigh < 0) rhigh = radix_high;
  if (rlow < 0) rlow = radix_low;
  if (rwidth < 0) rwidth = sort_radix_width_default;
  
  max_nclasses = powof2_typed(rwidth, slkey_pure_t);

  locals = sl_alloc(2 * (max_nareas * max_nclasses + max_nareas), sizeof(double));
  globals = sl_alloc(2 * (max_nareas * max_nclasses + max_nareas), sizeof(double));

  final_local_counts = final_locals;
  if (pc->weighted) final_local_weights = final_locals + max_nparts;
  else final_local_weights = final_local_counts;

  areas = areas0;
  areas_new = areas1;

  /* init the first area (all elements) */
  nareas = 1;
  elem_assign(s, &areas[0]);

  /* init all parts */
  nparts = size - 1;
  for (i = 0; i < size - 1; ++i)
  {
    parts[i] = i;
    part_areas[i] = 0;
  }

  /* init sdispls */
  for (i = 0; i < size; ++i) sdispls[i] = 0;

  rti_tstart(rti_tid_mpi_partition_radix_while);

  while (nparts > 0 && rhigh >= 0)
  {
    ++rounds;

    /* setup bitmask */
    current_width = xmin(rwidth, rhigh - rlow + 1);
    rhigh -= (current_width > 0)?current_width - 1:rhigh;

    nclasses = (current_width > 0)?powof2_typed(current_width, slkey_pure_t):1;
    bit_mask = nclasses - 1;

    /* init counters */
    local_counts = locals;
    global_counts = globals;
    if (pc->weighted)
    {
      local_weights = locals + (nareas * nclasses) + nareas;
      global_weights = globals + (nareas * nclasses) + nareas;

    } else
    {
      local_weights = local_counts;
      global_weights = global_counts;
    }

    /* zero all counter */
    for (i = 0; i < nareas; ++i)
    for (k = 0; k < nclasses; ++k) local_counts[i * nclasses + k] = local_weights[i * nclasses + k] = 0.0;

    rti_tstart(rti_tid_mpi_partition_radix_while_count);

    /* for every area */
    for (i = 0; i < nareas; ++i)
    {
      elem_assign_at(&areas[i], areas[i].size, &end);

      if (nclasses > 1)
      {
        if (pc->weighted)
        {
          /* counts and weights in every class */
          for (elem_assign(&areas[i], &xi); xi.keys < end.keys; elem_inc(&xi))
          {
            k = radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);
            local_counts[i * nclasses + k] += 1;
            local_weights[i * nclasses + k] += elem_weight(&xi, 0);
          }

        } else
        {
          /* counts in every class */
          for (elem_assign(&areas[i], &xi); xi.keys < end.keys; elem_inc(&xi)) ++local_counts[i * nclasses + radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];
        }

      } else
      {
        /* total counts and weights */
        local_counts[i * nclasses + 0] = areas[i].size;

        if (pc->weighted)
        {
          for (elem_assign(&areas[i], &xi); xi.keys < end.keys; elem_inc(&xi)) local_weights[i * nclasses + 0] += elem_weight(&xi, 0);
        }
      }

      /* total counts and weights in this area */
      local_counts[nareas * nclasses + i] = areas[i].size;
      if (pc->weighted)
      {
        local_weights[nareas * nclasses + i] = 0.0;
        for (k = 0; k < nclasses; ++k) local_weights[nareas * nclasses + i] += local_weights[i * nclasses + k];
      }
    }

    rti_tstop(rti_tid_mpi_partition_radix_while_count);

    --rhigh;

    rti_tstart(rti_tid_mpi_partition_radix_while_allreduce);

    /* create global counts and weights */
#ifdef MPI_PARTITION_RADIX_REDUCEBCAST_THRESHOLD
    if (size >= MPI_PARTITION_RADIX_REDUCEBCAST_THRESHOLD)
    {
      MPI_Reduce(locals, globals, (1 + (pc->weighted != 0)) * (nareas * nclasses + nareas), MPI_DOUBLE, MPI_SUM, REDUCEBCAST_ROOT, comm);
      MPI_Bcast(globals, (1 + (pc->weighted != 0)) * (nareas * nclasses + nareas), MPI_DOUBLE, REDUCEBCAST_ROOT, comm);

    } else
#endif
      MPI_Allreduce(locals, globals, (1 + (pc->weighted != 0)) * (nareas * nclasses + nareas), MPI_DOUBLE, MPI_SUM, comm);

    rti_tstop(rti_tid_mpi_partition_radix_while_allreduce);

/*    if (rank == 0)
    {
      printf("%d: pepcparts_mpi_partition_radix: allreduce: %f (nareas: %" sl_int_type_fmt ", nclasses: %" sl_key_type_fmt " -> %" sl_int_type_fmt " doubles)\n",
        rank,
        rti_tlast(rti_tid_mpi_partition_radix_while_allreduce),
        nareas,
        nclasses,
        (1 + (pc->weighted != 0)) * (nareas * (slint_t) nclasses + nareas));
    }*/

/*    if (rank == 0)
    {
      printf("%d: locals\n", rank);
      for (i = 0; i < nareas; ++i)
      {
        printf("%d: %" sl_int_type_fmt ":", rank, i);
        for (k = 0; k < nclasses; ++k) printf("  %f", local_counts[i * nclasses + k]);
        printf(" = %f\n", local_counts[nareas * nclasses + i]);
        printf("%d: %" sl_int_type_fmt ":", rank, i);
        for (k = 0; k < nclasses; ++k) printf("  %f", local_weights[i * nclasses + k]);
        printf(" = %f\n", local_weights[nareas * nclasses + i]);
      }
      printf("%d: globals\n", rank);
      for (i = 0; i < nareas; ++i)
      {
        printf("%d: %" sl_int_type_fmt ":", rank, i);
        for (k = 0; k < nclasses; ++k) printf("  %f", global_counts[i * nclasses + k]);
        printf(" = %f\n", global_counts[nareas * nclasses + i]);
        printf("%d: %" sl_int_type_fmt ":", rank, i);
        for (k = 0; k < nclasses; ++k) printf("  %f", global_weights[i * nclasses + k]);
        printf(" = %f\n", global_weights[nareas * nclasses + i]);
      }
    }*/

    /* do some initializations */
    if (rounds == 1)
    {
      rti_tstart(rti_tid_mpi_partition_radix_while_round1);
    
      /* distribute min/max counts or weights */
      if (pc->weighted)
      {
        minmax_weights[rank * 4 + 0] = (pc->min_weight >= 0)?pc->min_weight:(-pc->min_weight * global_weights[nareas * nclasses + 0] / size);
        minmax_weights[rank * 4 + 1] = (pc->max_weight >= 0)?pc->max_weight:(-pc->max_weight * global_weights[nareas * nclasses + 0] / size);
        minmax_weights[rank * 4 + 2] = (pc->min_wpart >= 0)?pc->min_wpart:(-pc->min_wpart * global_weights[nareas * nclasses + 0]);
        minmax_weights[rank * 4 + 3] = (pc->max_wpart >= 0)?pc->max_wpart:(-pc->max_wpart * global_weights[nareas * nclasses + 0]);

      } else
      {
        minmax_weights[rank * 4 + 0] = (pc->min_count >= 0)?pc->min_count:(-pc->min_count * global_counts[nareas * nclasses + 0] / size);
        minmax_weights[rank * 4 + 1] = (pc->max_count >= 0)?pc->max_count:(-pc->max_count * global_counts[nareas * nclasses + 0] / size);
        minmax_weights[rank * 4 + 2] = (pc->min_cpart >= 0)?pc->min_cpart:(-pc->min_cpart * global_counts[nareas * nclasses + 0]);
        minmax_weights[rank * 4 + 3] = (pc->max_cpart >= 0)?pc->max_cpart:(-pc->max_cpart * global_counts[nareas * nclasses + 0]);
      }

      rti_tstart(rti_tid_mpi_partition_radix_while_round1_allgather);
#ifdef HAVENT_MPI_IN_PLACE
      local_minmax_weights[0] = minmax_weights[rank * 4 + 0];
      local_minmax_weights[1] = minmax_weights[rank * 4 + 1];
      local_minmax_weights[2] = minmax_weights[rank * 4 + 2];
      local_minmax_weights[3] = minmax_weights[rank * 4 + 3];
      MPI_Allgather(local_minmax_weights, 4, MPI_DOUBLE, minmax_weights, 4, MPI_DOUBLE, comm);
/*      MPI_Gather(local_minmax_weights, 4, MPI_DOUBLE, minmax_weights, 4, MPI_DOUBLE, 0, comm);
      MPI_Bcast(minmax_weights, 4 * size, MPI_DOUBLE, 0, comm);*/
#else
      MPI_Allgather(MPI_IN_PLACE, 4, MPI_DOUBLE, minmax_weights, 4, MPI_DOUBLE, comm);
#endif
      rti_tstop(rti_tid_mpi_partition_radix_while_round1_allgather);

      for (i = 0; i < size; ++i)
      {
        border_offsets[i] = 0.0;
        border_finalize[i] = 0.0;
      }

#ifdef WEIGHT_STATS
        total_count = global_counts[nareas * nclasses + 0];
        total_weight = (pc->weighted)?global_weights[nareas * nclasses + 0]:global_counts[nareas * nclasses + 0];
#endif

      nbsets = 1;
      bsets[0] = size - 1;
      
      for (i = 0; i < size; ++i)
      {
        bsets_minmax[4 * i + 0] = ((i == 0)?0:bsets_minmax[4 * (i - 1) + 0]) + minmax_weights[4 * i + 0];
        bsets_minmax[4 * i + 2] = ((i == 0)?0:bsets_minmax[4 * (i - 1) + 2]) + minmax_weights[4 * i + 1];

        j = size - i - 1;
        bsets_minmax[4 * j + 1] = (j >= size - 1)?0:(bsets_minmax[4 * (j + 1) + 1] + minmax_weights[4 * (j + 1) + 0]);
        bsets_minmax[4 * j + 3] = (j >= size - 1)?0:(bsets_minmax[4 * (j + 1) + 3] + minmax_weights[4 * (j + 1) + 1]);

        bsets_parent[i] = 0;
      }

      bsets_minmax[4 * (size - 1) + 0] = bsets_minmax[4 * (size - 1) + 2] = 0;
      bsets_minmax[4 * (size - 1) + 1] = bsets_minmax[4 * (size - 1) + 3] = global_weights[nareas * nclasses + 0];

      bsets_parent[size - 1] = 0;

      rti_tstop(rti_tid_mpi_partition_radix_while_round1);
    }

    nareas_new = 0;
    last_new_area = last_new_class = -1;

    rti_tstart(rti_tid_mpi_partition_radix_while_check);

    /* check all remaining parts */
    nparts_removed = 0;
    for (i = 0; i < nparts; ++i)
    {
      remove_this_part = 0;
    
      j = bsets[bsets_parent[parts[i]]];

/*      if (rank == 0) printf("%d: %" sl_int_type_fmt ": parts[]: %" sl_int_type_fmt ", bsets_parent[]: %" sl_int_type_fmt ", bsets[]: %" sl_int_type_fmt "\n", rank, i, parts[i], bsets_parent[parts[i]], j);*/
      
      current_borders[0] = bsets_minmax[4 * j + 0] + bsets_minmax[4 * parts[i] + 0];
      current_borders[2] = bsets_minmax[4 * j + 2] + bsets_minmax[4 * parts[i] + 2];
      current_borders[1] = bsets_minmax[4 * j + 1] - bsets_minmax[4 * parts[i] + 1];
      current_borders[3] = bsets_minmax[4 * j + 3] - bsets_minmax[4 * parts[i] + 3];

/*      if (rank == 0)
      {
        printf("%d: %d: %f + %f\n", rank, 0, bsets_minmax[4 * j + 0], bsets_minmax[4 * parts[i] + 0]);
        printf("%d: %d: %f + %f\n", rank, 2, bsets_minmax[4 * j + 2], bsets_minmax[4 * parts[i] + 2]);
        printf("%d: %d: %f - %f\n", rank, 1, bsets_minmax[4 * j + 1], bsets_minmax[4 * parts[i] + 1]);
        printf("%d: %d: %f - %f\n", rank, 3, bsets_minmax[4 * j + 3], bsets_minmax[4 * parts[i] + 3]);
      }*/

/*      if (rank == 0) printf("%d: %" sl_int_type_fmt ": current_borders: %f  %f  %f  %f\n", rank, parts[i], current_borders[0], current_borders[1], current_borders[2], current_borders[3]);*/

      SL_ASSERT(current_borders[0] <= current_borders[1]);
      SL_ASSERT(current_borders[2] >= current_borders[3]);
      SL_ASSERT(minmax_weights[(parts[i] + 1) * 4 + 2] <= minmax_weights[(parts[i] + 0) * 4 + 3]);

      current_minmax[0] = xmax3(current_borders[0], current_borders[3], minmax_weights[(parts[i] + 1) * 4 + 2]) - border_offsets[parts[i]];
      current_minmax[1] = xmin3(current_borders[2], current_borders[1], minmax_weights[(parts[i] + 0) * 4 + 3]) - border_offsets[parts[i]];

      SL_ASSERT(current_minmax[0] <= current_minmax[1]);

/*      if (rank == 0) printf("%d: %" sl_int_type_fmt ": border_offset: %f\n", rank, parts[i], border_offsets[parts[i]]);
      if (rank == 0) printf("%d: %" sl_int_type_fmt ": current_minmax: %f  %f\n", rank, parts[i], current_minmax[0], current_minmax[1]);*/

      lcs = /*gcs =*/ 0;
      lws = gws = 0;

      for (k = 0; k < nclasses; ++k)
      {
        lc = local_counts[part_areas[i] * nclasses + k];
/*        gc = global_counts[part_areas[i] * nclasses + k];*/
        lw = local_weights[part_areas[i] * nclasses + k];
        gw = global_weights[part_areas[i] * nclasses + k];

        current_minmax[0] -= gw;
        current_minmax[1] -= gw;

        /* stop if max border is skipped */
        if (current_minmax[1] < 0) break;

        lcs += lc;
/*        gcs += gc;*/
        lws += lw;
        gws += gw;

        sdispls[parts[i] + 1] += lc;

        /* if between min/max border */
        if (current_minmax[0] <= 0 && current_minmax[1] >= 0)
        {
          /* look ahead for a better stop */
          if (k + 1 < nclasses) 
          if (current_minmax[1] - global_weights[part_areas[i] * nclasses + k + 1] >= 0)
          if (myabs(current_minmax[0] + current_minmax[1]) > myabs(current_minmax[0] + current_minmax[1] - 2 * global_weights[part_areas[i] * nclasses + k + 1])) continue;

          /* stop */
          break;
        }
      }

      SL_ASSERT(k < nclasses);

      /* final partitioning required? */
      if (current_minmax[1] < 0 && rhigh < 0)
      {
        /* choose border in the middle (for final partitioning) */
        gw += (current_minmax[0] + current_minmax[1]) / 2;

        current_minmax[0] -= gw;
        current_minmax[1] -= gw;

/*        gcs += gc;*/
        gws += gw;

        /* sdispls update requires prefix sums of local counts */

        /* mark this border as 'final partitioning required!' */
        border_finalize[parts[i]] = gw;

        /* save local count/weight for the later prefix calculations */
        final_local_counts[i - nparts_removed] = lc;
        if (pc->weighted) final_local_weights[i - nparts_removed] = lw;
      }

      border_offsets[parts[i]] += gws;

      /* hit */
      if (current_minmax[0] <= 0 && current_minmax[1] >= 0)
      {
/*        if (rank == 0) printf("%d: hit %" sl_int_type_fmt " at %f\n", rank, parts[i], border_offsets[parts[i]]);*/

        j = bsets[bsets_parent[parts[i]]];
        
        /* split j at position part[i] */

        /* update min/max of the new subsets part[i] and j */
        bsets_minmax_new[0] = bsets_minmax[4 * j + 0];
        bsets_minmax_new[2] = bsets_minmax[4 * j + 2];
        bsets_minmax_new[1] = border_offsets[parts[i]] + bsets_minmax[4 * parts[i] + 1];
        bsets_minmax_new[3] = border_offsets[parts[i]] + bsets_minmax[4 * parts[i] + 3];

        bsets_minmax[4 * j + 0] = border_offsets[parts[i]] - bsets_minmax[4 * parts[i] + 0];
        bsets_minmax[4 * j + 2] = border_offsets[parts[i]] - bsets_minmax[4 * parts[i] + 2];
/*        bsets_minmax[4 * j + 1] = bsets_minmax[4 * j + 1];
        bsets_minmax[4 * j + 3] = bsets_minmax[4 * j + 3];*/

        bsets_minmax[4 * parts[i] + 0] = bsets_minmax_new[0];
        bsets_minmax[4 * parts[i] + 2] = bsets_minmax_new[2];
        bsets_minmax[4 * parts[i] + 1] = bsets_minmax_new[1];
        bsets_minmax[4 * parts[i] + 3] = bsets_minmax_new[3];

/*        if (rank == 0)
        {
          printf("%d: parts[i]: %f  %f  %f  %f\n", rank, bsets_minmax[4 * parts[i] + 0], bsets_minmax[4 * parts[i] + 1], bsets_minmax[4 * parts[i] + 2], bsets_minmax[4 * parts[i] + 3]);
          printf("%d:        j: %f  %f  %f  %f\n", rank, bsets_minmax[4 * j + 0], bsets_minmax[4 * j + 1], bsets_minmax[4 * j + 2], bsets_minmax[4 * j + 3]);
        }*/

        /* perform the split (the smaller subset becomes a new set) */
        if (parts[i] - bsets_parent[j] > j - parts[i] - 1)
        {
          for (l = parts[i] + 1; l < j; ++l) bsets_parent[l] = nbsets;
          bsets[nbsets] = j;
          bsets[bsets_parent[parts[i]]] = parts[i];

        } else
        {
          for (l = bsets_parent[j]; l < parts[i]; ++l) bsets_parent[l] = nbsets;
          bsets[nbsets] = parts[i];
          SL_ASSERT(bsets[bsets_parent[parts[i]]] == j);
        }
        ++nbsets;

        /* set the first border of the new subsets */
        bsets_parent[parts[i]] = bsets_parent[j];
        bsets_parent[j] = parts[i] + 1;

        /* remove this part if it is marked as 'done' */
        if (border_finalize[parts[i]] == 0.0) remove_this_part = 1;

      } else if (rhigh >= 0)
      {
        if (last_new_area == part_areas[i] && last_new_class == k) part_areas[i] = nareas_new - 1;
        else
        {
          /* updates last_new_... */
          last_new_area = part_areas[i];
          last_new_class = k;

          /* create new area */
          elem_assign_at(&areas[part_areas[i]], lcs, &areas_new[nareas_new]);
          areas_new[nareas_new].size = local_counts[part_areas[i] * nclasses + k];
          part_areas[i] = nareas_new;
          ++nareas_new;
        };
      }

      if (remove_this_part) ++nparts_removed;
      else
      {
        parts[i - nparts_removed] = parts[i];
        part_areas[i - nparts_removed] = part_areas[i];
      }
    }

    nparts -= nparts_removed;

/*    if (rank == 0) printf("nparts = %" sl_int_type_fmt " vs. nareas_new = %" sl_int_type_fmt "\n", nparts, nareas_new);*/

    rti_tstop(rti_tid_mpi_partition_radix_while_check);
    
    /* switch areas */
    nareas = nareas_new;
    if (areas == areas0)
    {
      areas = areas1;
      areas_new = areas0;
    } else
    {
      areas = areas0;
      areas_new = areas1;
    }
  }

  rti_tstop(rti_tid_mpi_partition_radix_while);

  SL_ASSERT(nparts == 0 || rhigh < 0);

/*  printf("%d: nparts: %" sl_int_type_fmt "\n", rank, nparts);*/

  rti_tstart(rti_tid_mpi_partition_radix_final);

  /* still parts left, but no more bits */
  if (nparts > 0)
  {
    MPI_Exscan(final_local_weights, local_weights, nparts, MPI_DOUBLE, MPI_SUM, comm);
    if (rank == 0) for (i = 0; i < nparts; ++i) local_weights[i] = 0;

    for (i = 0; i < nparts; ++i)
    {
      d = xmax(0, border_finalize[parts[i]] - local_weights[i]);
      
      if (d >= final_local_weights[i]) sdispls[1 + parts[i]] += final_local_counts[i];
      else
      {
        if (pc->weighted)
        {
          do
          {
            d -= elem_weight(s, sdispls[1 + parts[i]]);
            ++sdispls[1 + parts[i]];

          } while (d >= 0);

          --sdispls[1 + parts[i]];

        } else sdispls[1 + parts[i]] += d;
      }
    }
  }

  rti_tstop(rti_tid_mpi_partition_radix_final);

  /* create scounts */
  for (i = 0; i < size - 1; ++i) scounts[i] = sdispls[i + 1] - sdispls[i];
  scounts[size - 1] = s->size - sdispls[size - 1];

#ifdef SCOUNTS_SDISPLS
  printf("%d: scounts", rank);
  for (i = 0; i < size; ++i) printf("  %d", scounts[i]);
  printf("\n");
  printf("%d: sdispls", rank);
  for (i = 0; i < size; ++i) printf("  %d", sdispls[i]);
  printf("\n");
#endif

  sl_free(locals);
  sl_free(globals);

#ifdef WEIGHT_STATS
  partial_counts[size] = 0;
  partial_weights[size] = 0.0;
  for (i = 0; i < size; ++i)
  {
    partial_counts[i] = scounts[i];
    partial_weights[i] = 0.0;
    for (j = sdispls[i]; j < sdispls[i] + scounts[i]; ++j) partial_weights[i] += elem_weight(s, j);
    
    partial_counts[size] += partial_counts[i];
    partial_weights[size] += partial_weights[i];
  }

#ifdef HAVENT_MPI_IN_PLACE
  MPI_Reduce(partial_counts, partial_counts2, size + 1, sl_int_type_mpi, MPI_SUM, 0, comm);
  MPI_Reduce(partial_weights, partial_weights2, size + 1, MPI_DOUBLE, MPI_SUM, 0, comm);
# define partial_counts   partial_counts2
# define partial_weights  partial_weights2
#else
  /* recvbuf requires workaround for an in-place/aliased-buffer-check-bug in mpich2 (fixed with rev 5518) */
  MPI_Reduce((rank == 0)?MPI_IN_PLACE:partial_counts, (rank == 0)?partial_counts:NULL, size + 1, sl_int_type_mpi, MPI_SUM, 0, comm);
  MPI_Reduce((rank == 0)?MPI_IN_PLACE:partial_weights, (rank == 0)?partial_weights:NULL, size + 1, MPI_DOUBLE, MPI_SUM, 0, comm);
#endif

  if (rank == 0)
  {
    printf("%d: total_count: %" sl_int_type_fmt " vs. %" sl_int_type_fmt "\n", rank, total_count, partial_counts[size]);
    d = 0.0;
    vmin = 1.0;
    vmax = 0.0;
    for (i = 0; i < size; ++i)
    {
/*      printf("%d: %" sl_int_type_fmt " %" sl_int_type_fmt " / %f - %" sl_int_type_fmt " / %f\n", rank, i, partial_counts[i], (double) partial_counts[i] / partial_counts[size], (partial_counts[size] / size) - partial_counts[i], fabs(1.0 - ((double) partial_counts[i] * size / partial_counts[size])));*/
      d += fabs((partial_counts[size] / size) - partial_counts[i]);
      if (fabs(1.0 - ((double) partial_counts[i] * size / partial_counts[size])) < vmin) vmin = fabs(1.0 - ((double) partial_counts[i] * size / partial_counts[size]));
      if (fabs(1.0 - ((double) partial_counts[i] * size / partial_counts[size])) > vmax) vmax = fabs(1.0 - ((double) partial_counts[i] * size / partial_counts[size]));
    }
    printf("%d: min/max: %f / %f\n", rank, vmin, vmax);
    printf("%d: average_count: %" sl_int_type_fmt " - %f / %f\n", rank, partial_counts[size] / size, d / size, d / partial_counts[size]);

    printf("%d: total_weight: %f vs. %f\n", rank, total_weight, partial_weights[size]);
    d = 0.0;
    vmin = 1.0;
    vmax = 0.0;
    for (i = 0; i < size; ++i)
    {
/*      printf("%d: %" sl_int_type_fmt " %f / %f - %f / %f\n", rank, i, partial_weights[i], partial_weights[i] / partial_weights[size], (partial_weights[size] / size) - partial_weights[i], fabs(1.0 - (partial_weights[i] * size / partial_weights[size])));*/
      d += fabs((partial_weights[size] / size) - partial_weights[i]);
      if (fabs(1.0 - (partial_weights[i] * size / partial_weights[size])) < vmin) vmin = fabs(1.0 - (partial_weights[i] * size / partial_weights[size]));
      if (fabs(1.0 - (partial_weights[i] * size / partial_weights[size])) > vmax) vmax = fabs(1.0 - (partial_weights[i] * size / partial_weights[size]));
    }
    printf("%d: min/max: %f / %f\n", rank, vmin, vmax);
    printf("%d: average_weight: %f - %f / %f\n", rank, partial_weights[size] / size, d / size, d / partial_weights[size]);
  }
#endif

  rti_tstop(rti_tid_mpi_partition_radix);

#if defined(TIMING_STATS) && defined(SL_USE_RTI_TIM)
  if (rank == 0)
  {
    printf("%d: pepcparts_mpi_partition_radix: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix));
    printf("%d: pepcparts_mpi_partition_radix: sync: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_sync));
    printf("%d: pepcparts_mpi_partition_radix: while: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while));
    printf("%d: pepcparts_mpi_partition_radix:   count: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_count));
    printf("%d: pepcparts_mpi_partition_radix:   allreduce: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_allreduce));
    printf("%d: pepcparts_mpi_partition_radix:   round1: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_round1));
    printf("%d: pepcparts_mpi_partition_radix:     allgather: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_round1_allgather));
    printf("%d: pepcparts_mpi_partition_radix:   check: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_check));
    printf("%d: pepcparts_mpi_partition_radix: final: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_final));
  }
#endif

  return 0;
}
