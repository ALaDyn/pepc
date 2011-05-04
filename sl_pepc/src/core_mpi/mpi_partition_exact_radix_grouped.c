/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_partition_exact_radix_grouped.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


/* sl_macro MPERG_TRACE_IF */


#include "sl_common.h"


/*#define PRINT_TIMINGS*/

#ifndef MPERG_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPERG_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPERG_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


#ifdef SL_INDEX

slint_t mpi_partition_exact_radix_ngroups(elements_t *s, partcond_t *pcond, slint_t ngroups, MPI_Comm *group_comms, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_radix_ngroups */
{
  slint_t g, i;

  slint_t lgcounts[2];
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  elements_t _sx;
  elements_t *es, ed0, ed1;
  
  int group_sizes[ngroups];
  int group_ranks[ngroups];
  
  MPI_Comm master_comms[ngroups];
  int master_sizes[ngroups];
  int master_ranks[ngroups];
  
  partcond_p group_pconds[ngroups];
  partcond_t merged_pconds[ngroups];

  int group_sdispls[size], group_scounts[size];
  int current_scounts[size], current_sdispls[size], current_rcounts[size], current_rdispls[size];

  int *_scounts, *_rcounts;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif


  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups);                  /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_pconds);           /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxin);            /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_up);               /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down);             /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_select);      /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall);    /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby);      /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_down_merge);       /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxout);           /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop);      /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall);  /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups);

  if (scounts == NULL) _scounts = z_alloca(size, sizeof(int)); else _scounts = scounts;
  if (rcounts == NULL) _rcounts = z_alloca(size, sizeof(int)); else _rcounts = rcounts;

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_pconds);

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute pconds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_pconds);

  if (sx == NULL || sx->size < 2 * pcond->count_max)
  {
    sx = &_sx;

    /* allocate scratch elements (need keys, indices and weights) */
    elements_alloc2(sx, 2 * pcond->count_max, 1, 1, 0, elem_has_weight);
  }

  /* use scratch for intermediate elements */  
  elem_assign(sx, &ed0);
  elem_assign_at(&ed0, (slint_t) pcond->count_max, &ed1);

  es = s;

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxin);

  /* assign source ranks */
  for (i = 0; i < es->size; ++i) es->indices[i] = rank;

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxin);


  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_up);

  /* collect and create merged pconds (bottom-up) */
  for (g = ngroups - 1; g >= 0; --g)
  {
    MPI_Comm_size(group_comms[g], &group_sizes[g]);
    MPI_Comm_rank(group_comms[g], &group_ranks[g]);
    
    if (g < ngroups - 1)
    {
      MPI_Comm_split(group_comms[g], (group_ranks[g + 1] == 0)?0:MPI_UNDEFINED, group_ranks[g], &master_comms[g]);
      if (master_comms[g] != MPI_COMM_NULL)
      {
        MPI_Comm_size(master_comms[g], &master_sizes[g]);
        MPI_Comm_rank(master_comms[g], &master_ranks[g]);

      } else master_ranks[g] = -1;

      MPI_Bcast(&master_sizes[g], 1, MPI_INT, 0, group_comms[g + 1]);

    } else
    {
      master_comms[g] = group_comms[g];
      master_sizes[g] = group_sizes[g];
      master_ranks[g] = group_ranks[g];
    }
    
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": group: %d of %d, master: %d of %d", g, group_ranks[g], group_sizes[g], master_ranks[g], master_sizes[g]);

    group_pconds[g] = z_alloca(master_sizes[g], sizeof(partcond_t));

    if (master_comms[g] != MPI_COMM_NULL) mpi_allgather_partconds((g < ngroups - 1)?&merged_pconds[g + 1]:pcond, group_pconds[g], master_sizes[g], master_ranks[g], master_comms[g]);

    if (g < ngroups - 1) mpi_bcast_partconds(master_sizes[g], group_pconds[g], 0, group_sizes[g + 1], group_ranks[g + 1], group_comms[g + 1]);

    merge_partconds(group_pconds[g], master_sizes[g], &merged_pconds[g]);
    
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": merged_pconds: %" slint_fmt "", g, merged_pconds[g].pcm);
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": merged_pconds: count: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].count_min, merged_pconds[g].count_max, merged_pconds[g].count_low, merged_pconds[g].count_high);
#ifdef elem_weight
    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": merged_pconds: weight: min/max: %" slweight_fmt "/%" slweight_fmt " - low/high: %" slweight_fmt "/%" slweight_fmt "", g, merged_pconds[g].weight_min, merged_pconds[g].weight_max, merged_pconds[g].weight_low, merged_pconds[g].weight_high);
#endif
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_up);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down);

  /* do selects (top-down) */
  for (g = 0; g < ngroups; ++g)
  {
    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_select);
    mpi_select_exact_radix(es, 1, master_sizes[g], group_pconds[g], rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, group_sizes[g], group_ranks[g], group_comms[g]);
    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_select);

    Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "%" slint_fmt ": group_sdispls = ", "%d  ", i, master_sizes[g], group_sdispls, g);

    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby);

    if (g < ngroups - 1)
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, group_scounts, es->size);

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "%" slint_fmt ": group_scounts = ", "%d  ", i, master_sizes[g], group_scounts, g);

      mpi_xcounts2ycounts_grouped(group_scounts, master_sizes[g], current_rcounts, group_comms[g + 1], master_comms[g], group_sizes[g], group_ranks[g], group_comms[g]);
    
      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "%" slint_fmt ": current_rcounts = ", "%d  ", i, group_sizes[g], current_rcounts, g);

/*#define SPARSE_SCOUNTS_FROM_RCOUNTS*/

      /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
      mpi_xcounts2ycounts_sparse(current_rcounts, current_scounts, es->size, group_sizes[g], group_ranks[g], group_comms[g]);
#else
      mpi_xcounts2ycounts_all2all(current_rcounts, current_scounts, group_sizes[g], group_ranks[g], group_comms[g]);
#endif

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "%" slint_fmt ": current_scounts = ", "%d  ", i, group_sizes[g], current_scounts, g);

    } else
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, current_scounts, es->size);

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "%" slint_fmt ": current_scounts = ", "%d  ", i, group_sizes[g], current_scounts, g);

      /* make rcounts from scounts */
      mpi_xcounts2ycounts_all2all(current_scounts, current_rcounts, group_sizes[g], group_ranks[g], group_comms[g]);

      Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "%" slint_fmt ": current_rcounts = ", "%d  ", i, group_sizes[g], current_rcounts, g);
    }

    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby);

    counts2displs(group_sizes[g], current_scounts, current_sdispls);
    counts2displs(group_sizes[g], current_rcounts, current_rdispls);

    Z_ASSERT(current_sdispls[group_sizes[g] - 1] + current_scounts[group_sizes[g] - 1] <= pcond->count_max);
    Z_ASSERT(current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1] <= pcond->count_max);

    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall);
#define xelem_data_not
#define xelem_data_weight
#define xelem_call \
    MPI_Alltoallv(xelem_buf(es), current_scounts, current_sdispls, xelem_mpi_datatype, xelem_buf(&ed0), current_rcounts, current_rdispls, xelem_mpi_datatype, group_comms[g]);
#include "sl_xelem_call.h"
    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall);

    ed0.size = current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1];

    Z_TRACE_IF(MPERG_TRACE_IF, "%" slint_fmt ": received: %" slint_fmt, g, ed0.size);

    rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_down_merge);

    /* merge received elements (but not in the last round) */
    if (g < ngroups - 1)
    {
      mergep_heap_int(&ed0, &ed1, group_sizes[g], current_rdispls, NULL);

      ed1.size = ed0.size;

      es = &ed1;
    
    } else es = &ed0;

    rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down_merge);

    z_freea(group_pconds[g]);
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_down);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxout);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop);

  /* create rcounts from source ranks (indexes) */
  for (i = 0; i < size; ++i) _rcounts[i] = 0;
  for (i = 0; i < es->size; ++i) ++_rcounts[es->indices[i]];

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "rcounts = ", "%d  ", i, size, _rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop);

  rti_tstart(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall);

  if (scounts)
  {
    /* make scounts from rcounts */
    mpi_xcounts2ycounts_all2all(_rcounts, _scounts, size, rank, comm);

    Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "scounts = ", "%d  ", i, size, _scounts);
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall);

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups_idxout);

  if (sx == &_sx) elements_free(sx);

  if (scounts == NULL) z_freea(_scounts);
  if (rcounts == NULL) z_freea(_rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_radix_ngroups);

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == 0)
  {
    printf("%d: mpi_partition_exact_radix_ngroups: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups));
    printf("%d: mpi_partition_exact_radix_ngroups: pconds: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_pconds));
    printf("%d: mpi_partition_exact_radix_ngroups: idxin: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxin));
    printf("%d: mpi_partition_exact_radix_ngroups: up: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_up));
    printf("%d: mpi_partition_exact_radix_ngroups: down: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down));
    printf("%d: mpi_partition_exact_radix_ngroups:  select: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_select));
    printf("%d: mpi_partition_exact_radix_ngroups:  alltoall: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_alltoall));
    printf("%d: mpi_partition_exact_radix_ngroups:  x2suby: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_x2suby));
    printf("%d: mpi_partition_exact_radix_ngroups:  merge: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_down_merge));
    printf("%d: mpi_partition_exact_radix_ngroups: idxout: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxout));
    printf("%d: mpi_partition_exact_radix_ngroups:  loop: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxout_loop));
    printf("%d: mpi_partition_exact_radix_ngroups:  alltoall: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_ngroups_idxout_alltoall));
  }
#endif

  return 0;
}

#endif /* SL_INDEX */


slint_t mpi_partition_exact_radix_2groups(elements_t *s, partcond_t *pcond, MPI_Comm group_comm, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_radix_2groups */
{
  slint_t i, j;

  slint_t lgcounts[2];
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  int group_size, group_rank;
  MPI_Comm master_comm;
  int master_size, master_rank;

  slint_t nparts;
  int group_sdispls[size], group_scounts[size];

  partcond_t *group_pconds, group_pcond;

  int sdispls[size], rdispls[size];
  slint_t rcounts_total, nsubelements;

  elements_t _sx;

  slint_t sub_elements_sources[size], sub_elements_sizes[size];
  elements_t sub_elements[size];

  int *sub_sdispls;

  int *_scounts, *_rcounts;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif


  rti_treset(rti_tid_mpi_partition_exact_radix_2groups);            /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_pconds);     /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_select1st);  /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_x2suby);     /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_alltoall);   /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_select2nd);  /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_radix_2groups_subx2y);     /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups);

  if (scounts == NULL) _scounts = z_alloca(size, sizeof(int)); else _scounts = scounts;
  if (rcounts == NULL) _rcounts = z_alloca(size, sizeof(int)); else _rcounts = rcounts;

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_pconds);

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute pconds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_pconds);

  /* init group */
  if (group_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(group_comm, &group_size);
    MPI_Comm_rank(group_comm, &group_rank);

  } else
  {
    group_size = 1;
    group_rank = 0;
  }
  
  /* init master */
  MPI_Comm_split(comm, (group_rank == 0)?0:MPI_UNDEFINED, rank, &master_comm);
  if (master_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(master_comm, &master_size);
    MPI_Comm_rank(master_comm, &master_rank);

  } else master_rank = -1;

  /* distribute num. of masters */
  if (group_comm != MPI_COMM_NULL) MPI_Bcast(&master_size, 1, MPI_INT, 0, group_comm);

  Z_TRACE_IF(MPERG_TRACE_IF, "%d: group: %d of %d, master: %d of %d", rank, group_rank, group_size, master_rank, master_size);

  nparts = master_size;

  /* create group partcond */
  group_pconds = z_alloca(group_size, sizeof(partcond_t));

  mpi_gather_partconds_grouped(pcond, group_comm, MPI_COMM_NULL, group_pconds, NULL, size, rank, comm);

  merge_partconds(group_pconds, group_size, &group_pcond);

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_select1st);

  /* perform 1st grouped select */  
  mpi_select_exact_radix_grouped(s, 1, &group_pcond, master_comm, group_comm, rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_select1st);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "group_sdispls = ", "%d  ", i, nparts, group_sdispls);

  /* create scounts from sdispls */
  displs2counts(nparts, group_sdispls, group_scounts, s->size);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "group_scounts = ", "%d  ", i, nparts, group_scounts);

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_x2suby);

  mpi_xcounts2ycounts_grouped(group_scounts, nparts, _rcounts, group_comm, master_comm, size, rank, comm);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "rcounts = ", "%d  ", i, size, _rcounts);

#define SPARSE_SCOUNTS_FROM_RCOUNTS

  /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
  mpi_xcounts2ycounts_sparse(_rcounts, _scounts, s->size, size, rank, comm);
#else
  mpi_xcounts2ycounts_all2all(_rcounts, _scounts, size, rank, comm);
#endif

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "scounts = ", "%d  ", i, size, _scounts);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_x2suby);

  /* create displs from counts */
  sdispls[0] = rdispls[0] = 0;
  nsubelements = (_rcounts[0] != 0);;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + _scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + _rcounts[i - 1];

    /* determine number of sub lists to receive */
    nsubelements += (_rcounts[i] != 0);
  }
  rcounts_total = rdispls[size - 1] + _rcounts[size - 1];

  Z_TRACE_IF(MPERG_TRACE_IF, "rcounts_total = %" slint_fmt, rcounts_total);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "sdispls = ", "%d  ", i, size, sdispls);
  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "rdispls = ", "%d  ", i, size, rdispls);

  if (sx == NULL || sx->size < rcounts_total)
  {
    sx = &_sx;

    /* allocate elements for 1st redistribution (need keys and weights) */
    elements_alloc2(sx, rcounts_total, 1, 0, 0, elem_has_weight);
  }

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_alltoall);

  /* 1st redistribution */
#define xelem_index_not
#define xelem_data_not
#define xelem_data_weight
#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), _scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), _rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_alltoall);

  /* create sub lists */
  j = 0;
  for (i = 0; i < size; ++i)
  if (_rcounts[i] != 0)
  {
    sub_elements_sources[j] = i;
    sub_elements_sizes[j] = _rcounts[i];

    elem_assign_at(sx, rdispls[i], &sub_elements[j]);
    sub_elements[j].size = _rcounts[i];

    ++j;
  }

  sub_sdispls = z_alloca(group_size * nsubelements, sizeof(int));

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_select2nd);

  /* perform 2nd select */
  mpi_select_exact_radix(sub_elements, nsubelements, group_size, group_pconds, rhigh, rlow, rwidth, 0, sub_sdispls, group_size, group_rank, group_comm);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_select2nd);

  rti_tstart(rti_tid_mpi_partition_exact_radix_2groups_subx2y);

  mpi_subxdispls2ycounts(nsubelements, sub_sdispls, sub_elements_sources, sub_elements_sizes, group_comm, group_size, _rcounts, size, rank, comm);

  Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "rcounts = ", "%d  ", i, size, _rcounts);

  if (scounts)
  {
    /* make scounts from rcounts */
    mpi_xcounts2ycounts_all2all(_rcounts, _scounts, size, rank, comm);

    Z_TRACE_ARRAY_IF(MPERG_TRACE_IF, "scounts = ", "%d  ", i, size, _scounts);
  }

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups_subx2y);

  z_freea(sub_sdispls);

  if (sx == &_sx) elements_free(sx);

  z_freea(group_pconds);

  if (master_comm != MPI_COMM_NULL) MPI_Comm_free(&master_comm);

  if (scounts == NULL) z_freea(_scounts);
  if (rcounts == NULL) z_freea(_rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_radix_2groups);

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == 0)
  {
    printf("%d: mpi_partition_exact_radix_2groups: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups));
    printf("%d: mpi_partition_exact_radix_2groups: pconds: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_pconds));
    printf("%d: mpi_partition_exact_radix_2groups: select1st: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_select1st));
    printf("%d: mpi_partition_exact_radix_2groups: x2suby: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_x2suby));
    printf("%d: mpi_partition_exact_radix_2groups: alltoall: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_alltoall));
    printf("%d: mpi_partition_exact_radix_2groups: select2nd: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_select2nd));
    printf("%d: mpi_partition_exact_radix_2groups: subx2y: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_radix_2groups_subx2y));
  }
#endif

  return 0;
}
