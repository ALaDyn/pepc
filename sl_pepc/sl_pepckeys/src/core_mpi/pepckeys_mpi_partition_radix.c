/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_partition_radix.c
 *  timestamp: 2010-01-08 15:04:27 +0100
 *  
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

/*#define SCOUNTS_SDISPLS
#define RCOUNTS_RDISPLS
#define WEIGHT_STATS
#define TIMING_STATS*/

#define SYNC_ON_INIT
#define HAVENT_MPI_IN_PLACE


#define DEBUG_RANK  -1


slint_t pepckeys_mpi_partition_radix(elements_t *s, partcond_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_partition_radix */
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
  slint_t parts_low, parts_high, nparts_removed;
  slint_t parts[max_nparts], part_areas[max_nparts];

  double parts_range_[2 * 2 * (1 + max_nparts + 1)];
  double *parts_range = parts_range_ + (2 * 2);
  double parts_minmax_[2 * 4 * (1 + max_nparts + 1)];
  double *parts_minmax = parts_minmax_ + (2 * 4);
  slint_t parts_update_[1 + max_nparts + 1];
  slint_t *parts_update = parts_update_ + 1;

  double parts_minmax_new[2 * 4];
  double current_minmax[2 * 2];
  
  double final_locals[2 * max_nparts];

  slint_t i, j, jp1, jm1, l, lp1, lm1;
  slint_t current_width;

  double minmax[2 * 4 * size];
  
  slint_t last_new_area, last_new_class;

#ifdef HAVENT_MPI_IN_PLACE
  double local_minmax[2 * 4];
#endif

  slint_t lc, lcs, gc, gcs;
  double lw, gw, lws, gws;
  double d, m;

  elements_t xi, end;

  slint_t round = 0;
  slint_t direction = 1;

  slint_t refine, finalize;

#ifdef RCOUNTS_RDISPLS
  int *rcounts, *rdispls;
#endif

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
  rti_treset(rti_tid_mpi_partition_radix_while_exscan);
  rti_treset(rti_tid_mpi_partition_radix_while_check);
  rti_treset(rti_tid_mpi_partition_radix_while_check_pre);
  rti_treset(rti_tid_mpi_partition_radix_while_check_classes);
  rti_treset(rti_tid_mpi_partition_radix_while_check_final);
  rti_treset(rti_tid_mpi_partition_radix_while_check_post);

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

  areas = areas0;
  areas_new = areas1;

  /* init the first area (all elements) */
  nareas = 1;
  elem_assign(s, &areas[0]);

  /* init all parts */
  parts_low = 0;
  parts_high = max_nparts - 1;
  for (i = parts_low; i <= parts_high; ++i)
  {
    parts[i] = i;
    part_areas[i] = 0;
  }

  /* init sdispls */
  for (i = 0; i < size; ++i) sdispls[i] = 0;

  rti_tstart(rti_tid_mpi_partition_radix_while);

  while (parts_low <= parts_high)
  {
    ++round;

    /* setup bitmask */
    current_width = xmin(rwidth, rhigh - rlow + 1);
    rhigh -= (current_width > 0)?current_width - 1:rhigh;

    nclasses = (current_width > 0)?powof2_typed(current_width, slkey_pure_t):1;
    bit_mask = nclasses - 1;
    
    SL_TRACE_IF((rank == DEBUG_RANK), "ROUND: %" sl_int_type_fmt ", rhigh: %" sl_int_type_fmt ", current_width: %" sl_int_type_fmt ", nclasses: %" sl_key_pure_type_fmt, round, rhigh, current_width, nclasses);

    finalize = (current_width <= 0);

    if (!finalize || round == 1)
    {
      /* init counters */
      local_counts = locals;
      global_counts = globals;
      local_weights = locals + (nareas * nclasses) + nareas;
      global_weights = globals + (nareas * nclasses) + nareas;

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
          /* counts and weights in every class */
          for (elem_assign(&areas[i], &xi); xi.keys < end.keys; elem_inc(&xi))
          {
            k = radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);
            local_counts[i * nclasses + k] += 1;
            local_weights[i * nclasses + k] += elem_weight(&xi, 0);
          }

        } else
        {
          /* total counts and weights */
          local_counts[i * nclasses + 0] = areas[i].size;

          for (elem_assign(&areas[i], &xi); xi.keys < end.keys; elem_inc(&xi)) local_weights[i * nclasses + 0] += elem_weight(&xi, 0);
        }

        /* total counts and weights in this area */
        local_counts[nareas * nclasses + i] = areas[i].size;

        local_weights[nareas * nclasses + i] = 0.0;
        for (k = 0; k < nclasses; ++k) local_weights[nareas * nclasses + i] += local_weights[i * nclasses + k];
      }

      rti_tstop(rti_tid_mpi_partition_radix_while_count);

      --rhigh;

      rti_tstart(rti_tid_mpi_partition_radix_while_allreduce);

      /* create global counts and weights */
#ifdef MPI_PARTITION_RADIX_REDUCEBCAST_THRESHOLD
      if (size >= MPI_PARTITION_RADIX_REDUCEBCAST_THRESHOLD)
      {
        MPI_Reduce(locals, globals, (1 + 1) * (nareas * nclasses + nareas), MPI_DOUBLE, MPI_SUM, REDUCEBCAST_ROOT, comm);
        MPI_Bcast(globals, (1 + 1) * (nareas * nclasses + nareas), MPI_DOUBLE, REDUCEBCAST_ROOT, comm);

      } else
#endif
        MPI_Allreduce(locals, globals, (1 + 1) * (nareas * nclasses + nareas), MPI_DOUBLE, MPI_SUM, comm);

      rti_tstop(rti_tid_mpi_partition_radix_while_allreduce);
    }

#ifdef TIMING
    SL_TRACE_IF((rank == DEBUG_RANK), "allreduce: %f, nareas: %" sl_int_type_fmt ", nclasses: %" sl_key_type_fmt ", doubles: %" sl_int_type_fmt, rti_tlast(rti_tid_mpi_partition_radix_while_allreduce), nareas, nclasses, (1 + 1) * (nareas * nclasses + nareas));
#endif

/*    if (rank == DEBUG_RANK)
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
    if (round == 1)
    {
      rti_tstart(rti_tid_mpi_partition_radix_while_round1);
    
      /* distribute min/max counts and weights */
      minmax[rank * 2 * 4 + 0 + 0] = (pc->min_count >= 0)?pc->min_count:(-pc->min_count * global_counts[nareas * nclasses + 0] / size);
      minmax[rank * 2 * 4 + 0 + 1] = (pc->max_count >= 0)?pc->max_count:(-pc->max_count * global_counts[nareas * nclasses + 0] / size);
      minmax[rank * 2 * 4 + 0 + 2] = (pc->min_cpart >= 0)?pc->min_cpart:(-pc->min_cpart * global_counts[nareas * nclasses + 0]);
      minmax[rank * 2 * 4 + 0 + 3] = (pc->max_cpart >= 0)?pc->max_cpart:(-pc->max_cpart * global_counts[nareas * nclasses + 0]);

      minmax[rank * 2 * 4 + 4 + 0] = (pc->min_weight >= 0)?pc->min_weight:(-pc->min_weight * global_weights[nareas * nclasses + 0] / size);
      minmax[rank * 2 * 4 + 4 + 1] = (pc->max_weight >= 0)?pc->max_weight:(-pc->max_weight * global_weights[nareas * nclasses + 0] / size);
      minmax[rank * 2 * 4 + 4 + 2] = (pc->min_wpart >= 0)?pc->min_wpart:(-pc->min_wpart * global_weights[nareas * nclasses + 0]);
      minmax[rank * 2 * 4 + 4 + 3] = (pc->max_wpart >= 0)?pc->max_wpart:(-pc->max_wpart * global_weights[nareas * nclasses + 0]);

      rti_tstart(rti_tid_mpi_partition_radix_while_round1_allgather);
#ifdef HAVENT_MPI_IN_PLACE
      local_minmax[0 + 0] = minmax[rank * 2 * 4 + 0 + 0];
      local_minmax[0 + 1] = minmax[rank * 2 * 4 + 0 + 1];
      local_minmax[0 + 2] = minmax[rank * 2 * 4 + 0 + 2];
      local_minmax[0 + 3] = minmax[rank * 2 * 4 + 0 + 3];
      local_minmax[4 + 0] = minmax[rank * 2 * 4 + 4 + 0];
      local_minmax[4 + 1] = minmax[rank * 2 * 4 + 4 + 1];
      local_minmax[4 + 2] = minmax[rank * 2 * 4 + 4 + 2];
      local_minmax[4 + 3] = minmax[rank * 2 * 4 + 4 + 3];
      MPI_Allgather(local_minmax, 2 * 4, MPI_DOUBLE, minmax, 2 * 4, MPI_DOUBLE, comm);
/*      MPI_Gather(local_minmax_weights, 2 * 4, MPI_DOUBLE, minmax_weights, 2 * 4, MPI_DOUBLE, 0, comm);
      MPI_Bcast(minmax_weights, 2 * 4 * size, MPI_DOUBLE, 0, comm);*/
#else
      MPI_Allgather(MPI_IN_PLACE, 2 * 4, MPI_DOUBLE, minmax_weights, 2 * 4, MPI_DOUBLE, comm);
#endif
      rti_tstop(rti_tid_mpi_partition_radix_while_round1_allgather);

#ifdef WEIGHT_STATS
      total_count = global_counts[nareas * nclasses + 0];
      total_weight = global_weights[nareas * nclasses + 0];
#endif

      parts_minmax[2 * 4 * (parts_low - 1) + 0 + 0] = parts_minmax[2 * 4 * (parts_low - 1) + 0 + 2] = 0;
      parts_minmax[2 * 4 * (parts_low - 1) + 0 + 1] = parts_minmax[2 * 4 * (parts_low - 1) + 0 + 3] = 0;
      parts_minmax[2 * 4 * (parts_low - 1) + 4 + 0] = parts_minmax[2 * 4 * (parts_low - 1) + 4 + 2] = 0;
      parts_minmax[2 * 4 * (parts_low - 1) + 4 + 1] = parts_minmax[2 * 4 * (parts_low - 1) + 4 + 3] = 0;

      parts_minmax[2 * 4 * (parts_high + 1) + 0 + 0] = parts_minmax[2 * 4 * (parts_high + 1) + 0 + 2] = 0;
      parts_minmax[2 * 4 * (parts_high + 1) + 0 + 1] = parts_minmax[2 * 4 * (parts_high + 1) + 0 + 3] = global_counts[nareas * nclasses + 0];
      parts_minmax[2 * 4 * (parts_high + 1) + 4 + 0] = parts_minmax[2 * 4 * (parts_high + 1) + 4 + 2] = 0;
      parts_minmax[2 * 4 * (parts_high + 1) + 4 + 1] = parts_minmax[2 * 4 * (parts_high + 1) + 4 + 3] = global_weights[nareas * nclasses + 0];

      parts_range[2 * 2 * (parts_low - 1) + 0 + 0] = parts_range[2 * 2 * (parts_high + 1) + 0 + 0] = 0.0;
      parts_range[2 * 2 * (parts_low - 1) + 0 + 1] = parts_range[2 * 2 * (parts_high + 1) + 0 + 1] = global_counts[nareas * nclasses + 0];
      parts_range[2 * 2 * (parts_low - 1) + 2 + 0] = parts_range[2 * 2 * (parts_high + 1) + 2 + 0] = 0.0;
      parts_range[2 * 2 * (parts_low - 1) + 2 + 1] = parts_range[2 * 2 * (parts_high + 1) + 2 + 1] = global_weights[nareas * nclasses + 0];

      for (i = parts_high; i >= parts_low; --i)
      {
        parts_minmax[2 * 4 * parts[i] + 0 + 1] = parts_minmax[2 * 4 * (parts[i] + 1) + 0 + 1] - minmax[2 * 4 * parts[i] + 0 + 0];
        parts_minmax[2 * 4 * parts[i] + 0 + 3] = parts_minmax[2 * 4 * (parts[i] + 1) + 0 + 3] - minmax[2 * 4 * parts[i] + 0 + 1];
        parts_minmax[2 * 4 * parts[i] + 4 + 1] = parts_minmax[2 * 4 * (parts[i] + 1) + 4 + 1] - minmax[2 * 4 * parts[i] + 4 + 0];
        parts_minmax[2 * 4 * parts[i] + 4 + 3] = parts_minmax[2 * 4 * (parts[i] + 1) + 4 + 3] - minmax[2 * 4 * parts[i] + 4 + 1];
        
        parts_minmax[2 * 4 * parts[i] + 0 + 0] = parts_minmax[2 * 4 * parts[i] + 0 + 2] = parts_minmax[2 * 4 * parts[i] + 4 + 0] = parts_minmax[2 * 4 * parts[i] + 4 + 2] = -1;

        parts_range[2 * 2 * parts[i] + 0 + 0] = 0.0;
        parts_range[2 * 2 * parts[i] + 0 + 1] = global_counts[nareas * nclasses + 0];
        parts_range[2 * 2 * parts[i] + 2 + 0] = 0.0;
        parts_range[2 * 2 * parts[i] + 2 + 1] = global_weights[nareas * nclasses + 0];
/*        SL_ASSERT(minmax[2 * 4 * (parts[i] + 1) + 0 + 2] <= minmax[2 * 4 * (parts[i] + 0) + 0 + 3]);*/
/*        SL_ASSERT(minmax[2 * 4 * (parts[i] + 1) + 4 + 2] <= minmax[2 * 4 * (parts[i] + 0) + 4 + 3]);*/

        parts_update[parts[i]] = 1;

        if (finalize)
        {
          final_locals[2 * i + 0] = local_counts[nareas * nclasses + 0];
          final_locals[2 * i + 1] = local_weights[nareas * nclasses + 0];
        }
      }

      rti_tstop(rti_tid_mpi_partition_radix_while_round1);
    }

    if (finalize)
    {
      j = parts_high - parts_low + 1;
    
      SL_TRACE_IF((rank == DEBUG_RANK), "Exscan: finalizing %" sl_int_type_fmt " parts", j);

      rti_tstart(rti_tid_mpi_partition_radix_while_exscan);

      MPI_Exscan(&final_locals[2 * parts_low], &locals[2 * parts_low], 2 * j, MPI_DOUBLE, MPI_SUM, comm);
      if (rank == 0) for (i = parts_low; i <= parts_high; ++i) locals[2 * i + 0] = locals[2 * i + 1] = 0;

      rti_tstop(rti_tid_mpi_partition_radix_while_exscan);
    }

    nareas_new = 0;
    last_new_area = last_new_class = -1;

    /* check all remaining parts */

    SL_TRACE_IF((rank == DEBUG_RANK), "ROUND: %" sl_int_type_fmt ", %s", round, (direction > 0)?"forward":"backward");

    nparts_removed = 0;

    rti_tstart(rti_tid_mpi_partition_radix_while_check);

    i = (direction > 0)?parts_low:parts_high;
    while ((direction > 0)?(i <= parts_high):(i >= parts_low))
    {
      rti_tstart(rti_tid_mpi_partition_radix_while_check_pre);
    
      SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ": PART: %" sl_int_type_fmt ",%" sl_int_type_fmt, round, i, parts[i]);

      j = 2 * 4 * parts[i];
      jp1 = 2 * 4 * (parts[i] + 1);
      jm1 = 2 * 4 * (parts[i] - 1);
      l = 2 * 2 * parts[i];
      lp1 = 2 * 2 * (parts[i] + 1);
      lm1 = 2 * 2 * (parts[i] - 1);

      if (parts_update[parts[i]])
      {
        if (direction > 0)
        {
          parts_minmax_new[0 + 0] = parts_minmax[jm1 + 0 + 0] + minmax[j + 0 + 0];
          parts_minmax_new[0 + 2] = parts_minmax[jm1 + 0 + 2] + minmax[j + 0 + 1];
          parts_minmax_new[4 + 0] = parts_minmax[jm1 + 4 + 0] + minmax[j + 4 + 0];
          parts_minmax_new[4 + 2] = parts_minmax[jm1 + 4 + 2] + minmax[j + 4 + 1];

          SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": %f + %f, %f + %f  /  %f + %f, %f + %f", i, parts[i],
            parts_minmax[jm1 + 0 + 0], minmax[j + 0 + 0],
            parts_minmax[jm1 + 0 + 2], minmax[j + 0 + 1],
            parts_minmax[jm1 + 4 + 0], minmax[j + 4 + 0],
            parts_minmax[jm1 + 4 + 2], minmax[j + 4 + 1]);

          if (parts_minmax_new[0 + 0] < minmax[jp1 + 0 + 2]) parts_minmax_new[0 + 0] = minmax[jp1 + 0 + 2];
          if (parts_minmax_new[0 + 2] > minmax[j   + 0 + 3]) parts_minmax_new[0 + 2] = minmax[j   + 0 + 3];
          if (parts_minmax_new[4 + 0] < minmax[jp1 + 4 + 2]) parts_minmax_new[4 + 0] = minmax[jp1 + 4 + 2];
          if (parts_minmax_new[4 + 2] > minmax[j   + 4 + 3]) parts_minmax_new[4 + 2] = minmax[j   + 4 + 3];

          parts_minmax_new[0 + 1] = parts_minmax[j + 0 + 1];
          parts_minmax_new[0 + 3] = parts_minmax[j + 0 + 3];
          parts_minmax_new[4 + 1] = parts_minmax[j + 4 + 1];
          parts_minmax_new[4 + 3] = parts_minmax[j + 4 + 3];

        } else
        {
          parts_minmax_new[0 + 1] = parts_minmax[jp1 + 0 + 1] - minmax[j + 0 + 0];
          parts_minmax_new[0 + 3] = parts_minmax[jp1 + 0 + 3] - minmax[j + 0 + 1];
          parts_minmax_new[4 + 1] = parts_minmax[jp1 + 4 + 1] - minmax[j + 4 + 0];
          parts_minmax_new[4 + 3] = parts_minmax[jp1 + 4 + 3] - minmax[j + 4 + 1];

          SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": %f - %f, %f - %f  /  %f - %f, %f - %f", i, parts[i],
            parts_minmax[jp1 + 0 + 1], minmax[j + 0 + 0],
            parts_minmax[jp1 + 0 + 3], minmax[j + 0 + 1],
            parts_minmax[jp1 + 4 + 1], minmax[j + 4 + 0],
            parts_minmax[jp1 + 4 + 3], minmax[j + 4 + 1]);

          if (parts_minmax_new[0 + 3] < minmax[jp1 + 0 + 2]) parts_minmax_new[0 + 3] = minmax[jp1 + 0 + 2];
          if (parts_minmax_new[0 + 1] > minmax[j   + 0 + 3]) parts_minmax_new[0 + 1] = minmax[j   + 0 + 3];
          if (parts_minmax_new[4 + 3] < minmax[jp1 + 4 + 2]) parts_minmax_new[4 + 3] = minmax[jp1 + 4 + 2];
          if (parts_minmax_new[4 + 1] > minmax[j   + 4 + 3]) parts_minmax_new[4 + 1] = minmax[j   + 4 + 3];

          parts_minmax_new[0 + 0] = parts_minmax[j + 0 + 0];
          parts_minmax_new[0 + 2] = parts_minmax[j + 0 + 2];
          parts_minmax_new[4 + 0] = parts_minmax[j + 4 + 0];
          parts_minmax_new[4 + 2] = parts_minmax[j + 4 + 2];
        }

        SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ": 1. parts_minmax_new: %f  %f  %f  %f  /  %f  %f  %f  %f", parts[i], parts_minmax_new[0 + 0], parts_minmax_new[0 + 1], parts_minmax_new[0 + 2], parts_minmax_new[0 + 3], parts_minmax_new[4 + 0], parts_minmax_new[4 + 1], parts_minmax_new[4 + 2], parts_minmax_new[4 + 3]);
        SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ": minmax: %f  %f  /  %f  %f", parts[i], minmax[2 * 4 * (parts[i] + 1) + 0 + 2], minmax[2 * 4 * (parts[i] + 0) + 0 + 3], minmax[2 * 4 * (parts[i] + 1) + 4 + 2], minmax[2 * 4 * (parts[i] + 0) + 4 + 3]);

        if (parts_minmax_new[0 + 0] > parts_minmax_new[0 + 1]) parts_minmax_new[0 + 0] = parts_minmax_new[0 + 1] = (parts_minmax_new[0 + 0] + parts_minmax_new[0 + 1]) / 2;
        if (parts_minmax_new[0 + 2] < parts_minmax_new[0 + 3]) parts_minmax_new[0 + 2] = parts_minmax_new[0 + 3] = (parts_minmax_new[0 + 2] + parts_minmax_new[0 + 3]) / 2;

        if (parts_minmax_new[4 + 0] > parts_minmax_new[4 + 1]) parts_minmax_new[4 + 0] = parts_minmax_new[4 + 1] = (parts_minmax_new[4 + 0] + parts_minmax_new[4 + 1]) / 2;
        if (parts_minmax_new[4 + 2] < parts_minmax_new[4 + 3]) parts_minmax_new[4 + 2] = parts_minmax_new[4 + 3] = (parts_minmax_new[4 + 2] + parts_minmax_new[4 + 3]) / 2;

      } else
      {
        parts_minmax_new[0 + 0] = parts_minmax[j + 0 + 0];
        parts_minmax_new[0 + 1] = parts_minmax[j + 0 + 1];
        parts_minmax_new[0 + 2] = parts_minmax[j + 0 + 2];
        parts_minmax_new[0 + 3] = parts_minmax[j + 0 + 3];

        parts_minmax_new[4 + 0] = parts_minmax[j + 4 + 0];
        parts_minmax_new[4 + 1] = parts_minmax[j + 4 + 1];
        parts_minmax_new[4 + 2] = parts_minmax[j + 4 + 2];
        parts_minmax_new[4 + 3] = parts_minmax[j + 4 + 3];
      }

      SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": 2. parts_minmax_new: %f  %f  %f  %f  /  %f  %f  %f  %f", i, parts[i], parts_minmax_new[0 + 0], parts_minmax_new[0 + 1], parts_minmax_new[0 + 2], parts_minmax_new[0 + 3], parts_minmax_new[4 + 0], parts_minmax_new[4 + 1], parts_minmax_new[4 + 2], parts_minmax_new[4 + 3]);

      current_minmax[0 + 0] = xmax(parts_minmax_new[0 + 0], parts_minmax_new[0 + 3]) - parts_range[l + 0 + 0];
      current_minmax[0 + 1] = xmin(parts_minmax_new[0 + 2], parts_minmax_new[0 + 1]) - parts_range[l + 0 + 0];

      current_minmax[2 + 0] = xmax(parts_minmax_new[4 + 0], parts_minmax_new[4 + 3]) - parts_range[l + 2 + 0];
      current_minmax[2 + 1] = xmin(parts_minmax_new[4 + 2], parts_minmax_new[4 + 1]) - parts_range[l + 2 + 0];

      SL_ASSERT(current_minmax[0 + 0] <= current_minmax[0 + 1]);
      SL_ASSERT(current_minmax[2 + 0] <= current_minmax[2 + 1]);

      rti_tstop(rti_tid_mpi_partition_radix_while_check_pre);

      SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ": current_minmax: %f  %f / %f  %f", parts[i], current_minmax[0 + 0], current_minmax[0 + 1], current_minmax[2 + 0], current_minmax[2 + 1]);

      lcs = gcs = 0;
      lws = gws = 0;

      /* HIT is the default */
      refine = 0;

      if (!finalize)
      {
        rti_tstart(rti_tid_mpi_partition_radix_while_check_classes);
      
        for (k = 0; k < nclasses; ++k)
        {
          lc = local_counts[part_areas[i] * nclasses + k];
          gc = global_counts[part_areas[i] * nclasses + k];
          lw = local_weights[part_areas[i] * nclasses + k];
          gw = global_weights[part_areas[i] * nclasses + k];

          current_minmax[0 + 0] -= gc;
          current_minmax[0 + 1] -= gc;

          current_minmax[2 + 0] -= gw;
          current_minmax[2 + 1] -= gw;

          SL_TRACE_IF((rank == DEBUG_RANK), "k = %" sl_key_pure_type_fmt ", current_minmax: %f  %f  / %f  %f", k, current_minmax[0], current_minmax[1], current_minmax[2], current_minmax[3]);

          /* stop and refine if max count is skipped OR min count AND max weight is skipped */
          if ((current_minmax[0 + 1] < 0) || (current_minmax[0 + 0] < 0 && current_minmax[2 + 1] < 0))
          {
            refine = 1;
            break;
          }

          lcs += lc;
          gcs += gc;
          lws += lw;
          gws += gw;

          gc = gw = 0.0;

          /* if between min/max counts */
          if (current_minmax[0 + 0] <= 0 && current_minmax[0 + 1] >= 0)
          {
            /* go to next if max count not reached AND min weight not reached */
            if (current_minmax[0 + 1] > 0 && current_minmax[2 + 0] > 0) continue;

            /* look ahead for a better stop */
            if (k + 1 < nclasses && current_minmax[0 + 1] - global_counts[part_areas[i] * nclasses + k + 1] >= 0)
            {
              /* continue if weights will improve */
              if (myabs(current_minmax[2 + 0] + current_minmax[2 + 1]) > myabs(current_minmax[2 + 0] + current_minmax[2 + 1] - 2 * global_weights[part_areas[i] * nclasses + k + 1])) continue;
            }

            /* stop */
            break;
          }
        }

        SL_ASSERT(k < nclasses);

        SL_TRACE_IF((rank == DEBUG_RANK), "%s k = %" sl_key_pure_type_fmt, (refine)?"REFINE":"HIT", k);
      
        rti_tstop(rti_tid_mpi_partition_radix_while_check_classes);

      } else
      {
        rti_tstart(rti_tid_mpi_partition_radix_while_check_final);

        /* middle of min/max weight */
        m = (current_minmax[2 + 0] + current_minmax[2 + 1]) / 2;

        /* min. part of weight to contribute */
        d = xmax(0, m - locals[i * 2 + 1]);

        /* contribute all? */
        if (d >= final_locals[i * 2 + 1])
        {
          lc = final_locals[i * 2 + 0];
          lw = final_locals[i * 2 + 1];

        } else
        {
          /* contribute only a part */
          lc = 0;
          lw = 0; /* not required */

          do
          {
            d -= elem_weight(s, sdispls[1 + parts[i]] + lc);
            ++lc;

          } while (d >= 0 && lc < final_locals[i * 2 + 0]);

          --lc;
        
          /* if unweighted, then m = middle of min/max count, d = ..., lc = d */
        }

        /* check mc against min/max count borders */
        lc = xminmax(current_minmax[0 + 0] - locals[i * 2 + 0], lc, current_minmax[0 + 1] - locals[i * 2 + 0]);

        /* check agains 0 (don't step back!) and the local contribution */
        lc = xminmax(0, lc, final_locals[i * 2 + 0]);

        /* the exact global counts/weights are unknown (set gc/gw so that parts_range is not changed) */
        gc = 0;
        gw = 0;

        lcs += lc;
        gcs += gc;
        lws += lw;
        gws += gw;
        
        gc = (parts_range[2 * 2 * parts[i] + 0 + 1] - parts_range[2 * 2 * parts[i] + 0 + 0]);
        gw = (parts_range[2 * 2 * parts[i] + 2 + 1] - parts_range[2 * 2 * parts[i] + 2 + 0]);

        rti_tstop(rti_tid_mpi_partition_radix_while_check_final);
      }      

      rti_tstart(rti_tid_mpi_partition_radix_while_check_post);
      
      SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": sdispls[%" sl_int_type_fmt " + 1] = %d, lcs = %" sl_int_type_fmt, i, parts[i], parts[i], sdispls[parts[i] + 1], lcs);

      sdispls[parts[i] + 1] += lcs;

      if (gcs > 0 || gws > 0)
      {
        parts_range[l + 0 + 0] += gcs;
        parts_range[l + 0 + 1] = parts_range[l + 0 + 0] + gc;
        parts_range[l + 2 + 0] += gws;
        parts_range[l + 2 + 1] = parts_range[l + 2 + 0] + gw;

        SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": 3. part_minmax_new: %f  %f  %f  %f  /  %f  %f  %f  %f", i, parts[i], parts_minmax_new[0 + 0], parts_minmax_new[0 + 1], parts_minmax_new[0 + 2], parts_minmax_new[0 + 3], parts_minmax_new[4 + 0], parts_minmax_new[4 + 1], parts_minmax_new[4 + 2], parts_minmax_new[4 + 3]);
        SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": parts_range: %f  %f  /  %f  %f", i, parts[i], parts_range[2 * 2 * parts[i] + 0 + 0], parts_range[2 * 2 * parts[i] + 0 + 1], parts_range[2 * 2 * parts[i] + 2 + 0], parts_range[2 * 2 * parts[i] + 2 + 1]);

        parts_minmax_new[0 + 0] = xminmax(parts_range[l + 0 + 0], parts_minmax_new[0 + 0], parts_range[l + 0 + 1]);
        parts_minmax_new[0 + 2] = xminmax(parts_range[l + 0 + 0], parts_minmax_new[0 + 2], parts_range[l + 0 + 1]);
        parts_minmax_new[0 + 1] = xminmax(parts_range[l + 0 + 0], parts_minmax_new[0 + 1], parts_range[l + 0 + 1]);
        parts_minmax_new[0 + 3] = xminmax(parts_range[l + 0 + 0], parts_minmax_new[0 + 3], parts_range[l + 0 + 1]);
      
        parts_minmax_new[4 + 0] = xminmax(parts_range[l + 2 + 0], parts_minmax_new[4 + 0], parts_range[l + 2 + 1]);
        parts_minmax_new[4 + 2] = xminmax(parts_range[l + 2 + 0], parts_minmax_new[4 + 2], parts_range[l + 2 + 1]);
        parts_minmax_new[4 + 1] = xminmax(parts_range[l + 2 + 0], parts_minmax_new[4 + 1], parts_range[l + 2 + 1]);
        parts_minmax_new[4 + 3] = xminmax(parts_range[l + 2 + 0], parts_minmax_new[4 + 3], parts_range[l + 2 + 1]);
      }

      SL_TRACE_IF((rank == DEBUG_RANK), "%" sl_int_type_fmt ",%" sl_int_type_fmt ": 4. part_minmax_new: %f  %f  %f  %f  /  %f  %f  %f  %f", i, parts[i], parts_minmax_new[0 + 0], parts_minmax_new[0 + 1], parts_minmax_new[0 + 2], parts_minmax_new[0 + 3], parts_minmax_new[4 + 0], parts_minmax_new[4 + 1], parts_minmax_new[4 + 2], parts_minmax_new[4 + 3]);

      if (parts_minmax_new[0 + 0] != parts_minmax[j + 0 + 0] || parts_minmax_new[0 + 2] != parts_minmax[j + 0 + 2] || parts_minmax_new[4 + 0] != parts_minmax[j + 4 + 0] || parts_minmax_new[4 + 2] != parts_minmax[j + 4 + 2])
      {
        parts_minmax[j + 0 + 0] = parts_minmax_new[0 + 0];
        parts_minmax[j + 0 + 2] = parts_minmax_new[0 + 2];
        parts_minmax[j + 4 + 0] = parts_minmax_new[4 + 0];
        parts_minmax[j + 4 + 2] = parts_minmax_new[4 + 2];

        parts_update[parts[i] + 1] = 1;
      }

      if (parts_minmax_new[0 + 1] != parts_minmax[j + 0 + 1] || parts_minmax_new[0 + 3] != parts_minmax[j + 0 + 3] || parts_minmax_new[4 + 1] != parts_minmax[j + 4 + 1] || parts_minmax_new[4 + 3] != parts_minmax[j + 4 + 3])
      {
        parts_minmax[j + 0 + 1] = parts_minmax_new[0 + 1];
        parts_minmax[j + 0 + 3] = parts_minmax_new[0 + 3];
        parts_minmax[j + 4 + 1] = parts_minmax_new[4 + 1];
        parts_minmax[j + 4 + 3] = parts_minmax_new[4 + 3];

        parts_update[parts[i] - 1] = 1;
      }

      parts_update[parts[i]] = 0;

      /* refine or remove */
      if (refine)
      {
        /* bits left for partitioning? */
        if (rhigh >= rlow)
        {
          if (last_new_area == part_areas[i] && last_new_class == k) part_areas[i] = nareas_new - 1;
          else
          {
            /* update last_new_... */
            last_new_area = part_areas[i];
            last_new_class = k;

            /* create new area */
            elem_assign_at(&areas[part_areas[i]], lcs, &areas_new[nareas_new]);
            areas_new[nareas_new].size = local_counts[part_areas[i] * nclasses + k];
            part_areas[i] = nareas_new;
            ++nareas_new;
          }

        } else
        {
          /* save local count/weight for the later prefix calculations */
          final_locals[2 * (i - nparts_removed * direction) + 0] = lc;
          final_locals[2 * (i - nparts_removed * direction) + 1] = lw;
        }

        parts[i - nparts_removed * direction] = parts[i];
        part_areas[i - nparts_removed * direction] = part_areas[i];

      } else ++nparts_removed;

      rti_tstop(rti_tid_mpi_partition_radix_while_check_post);
      
      i += direction;
    }

    if (direction > 0) parts_high -= nparts_removed;
    else parts_low += nparts_removed;

    direction *= -1;

/*    SL_NOTICE_IF((rank == DEBUG_RANK), "nparts = %" sl_int_type_fmt " vs. nareas_new = %" sl_int_type_fmt, nparts, nareas_new);*/

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

  /* create scounts */
  for (i = 0; i < size - 1; ++i) scounts[i] = sdispls[i + 1] - sdispls[i];
  scounts[size - 1] = s->size - sdispls[size - 1];

#ifdef SCOUNTS_SDISPLS
  printf("%d: scounts", rank);
  for (i = 0, j = 0; i < size; ++i) { printf("  %d", scounts[i]); j += scounts[i]; }
  printf(" = %" sl_int_type_fmt "\n", j);
  printf("%d: sdispls", rank);
  for (i = 0; i < size; ++i) printf("  %d", sdispls[i]);
  printf("\n");
#endif

#ifdef RCOUNTS_RDISPLS
  rcounts = sl_alloc(size, sizeof(int));
  rdispls = sl_alloc(size, sizeof(int));

  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

  rdispls[0] = 0;
  for (i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];

  printf("%d: rcounts", rank);
  for (i = 0; i < size; ++i) printf("  %d", rcounts[i]);
  printf("\n");
  printf("%d: rdispls", rank);
  for (i = 0; i < size; ++i) printf("  %d", rdispls[i]);
  printf("\n");

  sl_free(rcounts);
  sl_free(rdispls);
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
    printf("%d: pepckeys_mpi_partition_radix: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix));
    printf("%d: pepckeys_mpi_partition_radix: sync: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_sync));
    printf("%d: pepckeys_mpi_partition_radix: while: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while));
    printf("%d: pepckeys_mpi_partition_radix:   count: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_count));
    printf("%d: pepckeys_mpi_partition_radix:   allreduce: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_allreduce));
    printf("%d: pepckeys_mpi_partition_radix:   round1: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_round1));
    printf("%d: pepckeys_mpi_partition_radix:     allgather: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_round1_allgather));
    printf("%d: pepckeys_mpi_partition_radix:   exscan: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while_exscan));
    printf("%d: pepckeys_mpi_partition_radix:   check: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_radix_while_check));
    printf("%d: pepckeys_mpi_partition_radix:     pre: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while_check_pre));
    printf("%d: pepckeys_mpi_partition_radix:     classes: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while_check_classes));
    printf("%d: pepckeys_mpi_partition_radix:     final: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while_check_final));
    printf("%d: pepckeys_mpi_partition_radix:     post: %f\n", rank, rti_tlast(rti_tid_mpi_partition_radix_while_check_post));
  }
#endif

  return 0;
}
