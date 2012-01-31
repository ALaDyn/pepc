/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_sort_partition.c
 *  
 */


/* sl_macro MSP_TRACE_IF */
/* sl_macro MSP_VERIFY */
/* sl_macro MSP_VERIFY_OUTPUT */


#include "sl_common.h"

/*#define MSP_VERIFY
#define MSP_VERIFY_OUTPUT*/


double msp_t[4];  /* sl_global, sl_var msp_t */

slint_t msp_sync = 0;  /* sl_global, sl_var msp_sync */

partcond_t *msp_r_pc = NULL;  /* sl_global, sl_var msp_r_pc */


#ifndef MSP_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSP_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSP_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_sort_partition(elements_t *s0, elements_t *s1, elements_t *xs, slint_t part_type, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition */
{
  const slint_t paex_rhigh = -1;
  const slint_t paex_rlow = -1;
  const slint_t paex_rwidth = 3;

  partcond_t pc;
  
  double imba = 0.01;
  
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];

  const slint_t max_nsubs = 4;
  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];


  if (msp_sync) MPI_Barrier(comm);
  msp_t[0] = z_time_get_s();
#ifdef key_integer
  sort_radix(s0, xs, -1, -1, -1);
#else
  sort_quick(s0, xs);
#endif
  msp_t[0] = z_time_get_s() - msp_t[0];

  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s();

  pc.pcm = SLPC_COUNTS_MM;
  pc.count_min = s0->size * (1.0 - imba);
  pc.count_max = s0->size * (1.0 + imba);

  if (part_type == 1)
  {
#ifdef key_integer
    mpi_partition_exact_radix(s0, &pc, paex_rhigh, paex_rlow, paex_rwidth, SL_SORTED_IN, scounts, rcounts, size, rank, comm);
#else
#endif
  
  } else if (part_type == 2)
  {
    nsubs = 2;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
#ifdef key_integer
    mpi_partition_exact_radix_2groups(s0, &pc, sub_comms[1], NULL, paex_rhigh, paex_rlow, paex_rwidth, scounts, rcounts, size, rank, comm);
#else
#endif
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  
  }
#ifdef SL_INDEX
    else if (part_type < 0)
  {
    nsubs = (-part_type <= max_nsubs)?-part_type:max_nsubs;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
#ifdef key_integer
    mpi_partition_exact_radix_ngroups(s0, &pc, nsubs, sub_comms, NULL, paex_rhigh, paex_rlow, paex_rwidth, scounts, rcounts, size, rank, comm);
#else
#endif
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  }
#endif

  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s() - msp_t[1];


  msp_t[2] = z_time_get_s();
  
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

#define xelem_index_not
#define xelem_call \
    MPI_Alltoallv(xelem_buf(s0), scounts, sdispls, xelem_mpi_datatype, xelem_buf(s1), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"
  
  msp_t[2] = z_time_get_s() - msp_t[2];

  s1->size = rdispls[size - 1] + rcounts[size - 1];

  if (msp_sync) MPI_Barrier(comm);
  msp_t[3] = z_time_get_s();

#ifdef key_integer
  sort_radix(s1, xs, -1, -1, -1);
#else
  sort_quick(s1, xs);
#endif

  msp_t[3] = z_time_get_s() - msp_t[3];

  return 0;
}


#ifdef key_integer

slint_t mpi_sort_partition_radix(elements_t *s0, elements_t *s1, elements_t *xs, slint_t part_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_radix */
{
  partcond_t pc;
  
  double imba = 0.01;
  
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];

  const slint_t max_nsubs = 4;
  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif


  /* local sort */
  if (msp_sync) MPI_Barrier(comm);
  msp_t[0] = z_time_get_s();
  sort_radix(s0, xs, rhigh, rlow, -1);
  msp_t[0] = z_time_get_s() - msp_t[0];

  /* partitioning */
  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s();

  if (msp_r_pc) pc = *msp_r_pc;
  else
  {
    pc.pcm = SLPC_COUNTS_MM;
    pc.count_min = -(1.0 - imba);
    pc.count_max = -(1.0 + imba);
  }

#ifdef MSEG_BORDER_UPDATE_REDUCTION
  mseg_border_update_count_reduction = 0.25;
#endif

  if (part_type == 1)
  {
    mpi_partition_exact_radix(s0, &pc, rhigh, rlow, rwidth, SL_SORTED_IN, scounts, rcounts, size, rank, comm);
  
  } else if (part_type == 2)
  {
    nsubs = 2;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
    mpi_partition_exact_radix_2groups(s0, &pc, sub_comms[1], NULL, rhigh, rlow, rwidth, scounts, rcounts, size, rank, comm);
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  
  }
#ifdef SL_INDEX
    else if (part_type < 0)
  {
    nsubs = (-part_type <= max_nsubs)?-part_type:max_nsubs;

    mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
    mpi_partition_exact_radix_ngroups(s0, &pc, nsubs, sub_comms, NULL, rhigh, rlow, rwidth, scounts, rcounts, size, rank, comm);
    mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);
  }
#endif

  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s() - msp_t[1];

  /* all-to-all */
  msp_t[2] = z_time_get_s();
  
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

#define xelem_index_not
#define xelem_call \
    MPI_Alltoallv(xelem_buf(s0), scounts, sdispls, xelem_mpi_datatype, xelem_buf(s1), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  /* local merge */
  if (msp_sync) MPI_Barrier(comm);
  msp_t[2] = z_time_get_s() - msp_t[2];

  s0->size = s1->size = rdispls[size - 1] + rcounts[size - 1];

#define HOW_TO_MERGEP  3

  if (msp_sync) MPI_Barrier(comm);
  msp_t[3] = z_time_get_s();
#if HOW_TO_MERGEP == 0
  sort_radix(s1, s0, rhigh, rlow, -1);
  elem_ncopy(s1, s0, s1->size);
#elif HOW_TO_MERGEP == 1
  mergep_heap_int(s1, s0, size, rdispls, rcounts);
#elif HOW_TO_MERGEP == 2
  mergep_2way_ip_int(s1, s0, size, rdispls, merge2_basic_straight_01_x);
  elem_ncopy(s1, s0, s1->size);
#else
  mergep_2way_ip_int_rec(s1, s0, size, rdispls, merge2_basic_straight_01_x);
  elem_ncopy(s1, s0, s1->size);
#endif
  msp_t[3] = z_time_get_s() - msp_t[3];

#undef HOW_TO_MERGEP

#ifdef MSP_VERIFY
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;
}


slint_t mpi_sort_partition_exact_radix(elements_t *s, elements_t *sx, partcond_t *pcond, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_exact_radix */
{
  elements_t _sx;

  int scounts[size], sdispls[size], rcounts[size], rdispls[size];

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif


  if (sx == NULL || sx->size < s->size)
  {
    elements_alloc2(&_sx, pcond->count_max, 1, 1, 1, 1);
    sx = &_sx;
  }

  if (msp_sync) MPI_Barrier(comm);
  msp_t[0] = z_time_get_s();
  sort_radix(s, sx, rhigh, rlow, -1);
  msp_t[0] = z_time_get_s() - msp_t[0];

  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s();
  mpi_partition_exact_radix(s, pcond, rhigh, rlow, rwidth, SL_SORTED_IN, scounts, rcounts, size, rank, comm);
  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s() - msp_t[1];

  msp_t[2] = z_time_get_s();
  
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

#define xelem_call \
    MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"
  
  msp_t[2] = z_time_get_s() - msp_t[2];

  s->size = sx->size = rdispls[size - 1] + rcounts[size - 1];

  if (msp_sync) MPI_Barrier(comm);
  msp_t[3] = z_time_get_s();
  mergep_heap_int(sx, s, size, rdispls, rcounts);
  msp_t[3] = z_time_get_s() - msp_t[3];

  if (sx == &_sx) elements_free(sx);

#ifdef MSP_VERIFY
  mpi_elements_check_order(s, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;
}


#ifdef SL_INDEX

slint_t mpi_sort_partition_exact_radix_ngroups(elements_t *s, elements_t *sx, partcond_t *pcond, slint_t ngroups, MPI_Comm *group_comms, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_exact_radix_ngroups */
{
  slint_t g, i;

  slint_t lgcounts[2];
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  elements_t _sx;
  
  int group_sizes[ngroups];
  int group_ranks[ngroups];
  
  MPI_Comm master_comms[ngroups];
  int master_sizes[ngroups];
  int master_ranks[ngroups];
  
  partcond_p group_pconds[ngroups];
  partcond_t merged_pconds[ngroups];

  int group_sdispls[size], group_scounts[size];
  int current_scounts[size], current_sdispls[size], current_rcounts[size], current_rdispls[size];

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif

  double _t;


  msp_t[0] = msp_t[1] = msp_t[2] = msp_t[3] = 0;

  if (sx == NULL || sx->size < s->size)
  {
    elements_alloc2(&_sx, pcond->count_max, 1, 1, 1, 1);
    sx = &_sx;
  }

  if (msp_sync) MPI_Barrier(comm);
  _t = z_time_get_s();
  sort_radix(s, sx, rhigh, rlow, -1);
  msp_t[0] += z_time_get_s() - _t;

  if (msp_sync) MPI_Barrier(comm);
  _t = z_time_get_s();

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
    
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": group: %d of %d, master: %d of %d", g, group_ranks[g], group_sizes[g], master_ranks[g], master_sizes[g]);

    group_pconds[g] = z_alloca(master_sizes[g], sizeof(partcond_t));

    if (master_comms[g] != MPI_COMM_NULL) mpi_allgather_partconds((g < ngroups - 1)?&merged_pconds[g + 1]:pcond, group_pconds[g], master_sizes[g], master_ranks[g], master_comms[g]);

    if (g < ngroups - 1) mpi_bcast_partconds(master_sizes[g], group_pconds[g], 0, group_sizes[g + 1], group_ranks[g + 1], group_comms[g + 1]);

    merge_partconds(group_pconds[g], master_sizes[g], &merged_pconds[g]);
    
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": merged_pconds: %" slint_fmt "", g, merged_pconds[g].pcm);
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": merged_pconds: count: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].count_min, merged_pconds[g].count_max, merged_pconds[g].count_low, merged_pconds[g].count_high);
#ifdef elem_weight
    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": merged_pconds: weight: min/max: %f/%f - low/high: %f/%f", g, merged_pconds[g].weight_min, merged_pconds[g].weight_max, merged_pconds[g].weight_low, merged_pconds[g].weight_high);
#endif
  }

  /* do selects (top-down) */
  for (g = 0; g < ngroups; ++g)
  {
    mpi_select_exact_radix(s, 1, master_sizes[g], group_pconds[g], rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, group_sizes[g], group_ranks[g], group_comms[g]);

/*    Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": group_sdispls = ", "%d  ", i, master_sizes[g], group_sdispls, g);*/

    if (g < ngroups - 1)
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, group_scounts, s->size);

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": group_scounts = ", "%d  ", i, master_sizes[g], group_scounts, g);*/

      mpi_xcounts2ycounts_grouped(group_scounts, master_sizes[g], current_rcounts, group_comms[g + 1], master_comms[g], group_sizes[g], group_ranks[g], group_comms[g]);
    
/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_rcounts = ", "%d  ", i, group_sizes[g], current_rcounts, g);*/

/*#define SPARSE_SCOUNTS_FROM_RCOUNTS*/

      /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
      mpi_xcounts2ycounts_sparse(current_rcounts, current_scounts, es->size, group_sizes[g], group_ranks[g], group_comms[g]);
#else
      mpi_xcounts2ycounts_all2all(current_rcounts, current_scounts, group_sizes[g], group_ranks[g], group_comms[g]);
#endif

#undef SPARSE_SCOUNTS_FROM_RCOUNTS

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_scounts = ", "%d  ", i, group_sizes[g], current_scounts, g);*/

    } else
    {
      /* create scounts from sdispls */
      displs2counts(master_sizes[g], group_sdispls, current_scounts, s->size);

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_scounts = ", "%d  ", i, group_sizes[g], current_scounts, g);*/

      /* make rcounts from scounts */
      mpi_xcounts2ycounts_all2all(current_scounts, current_rcounts, group_sizes[g], group_ranks[g], group_comms[g]);

/*      Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "%" slint_fmt ": current_rcounts = ", "%d  ", i, group_sizes[g], current_rcounts, g);*/
    }

    counts2displs(group_sizes[g], current_scounts, current_sdispls);
    counts2displs(group_sizes[g], current_rcounts, current_rdispls);

    Z_ASSERT(current_sdispls[group_sizes[g] - 1] + current_scounts[group_sizes[g] - 1] <= pcond->count_max);
    Z_ASSERT(current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1] <= pcond->count_max);

    msp_t[1] += z_time_get_s() - _t;

    _t = z_time_get_s();

#define xelem_call \
    MPI_Alltoallv(xelem_buf(s), current_scounts, current_sdispls, xelem_mpi_datatype, xelem_buf(sx), current_rcounts, current_rdispls, xelem_mpi_datatype, group_comms[g]);
#include "sl_xelem_call.h"

    sx->size = current_rdispls[group_sizes[g] - 1] + current_rcounts[group_sizes[g] - 1];

    Z_TRACE_IF(MSP_TRACE_IF, "%" slint_fmt ": received: %" slint_fmt, g, sx->size);

    msp_t[2] += z_time_get_s() - _t;

    _t = z_time_get_s();

    /* merge received elements (but not in the last round) */
    mergep_heap_int(sx, s, group_sizes[g], current_rdispls, NULL);

    s->size = sx->size;

    msp_t[3] += z_time_get_s() - _t;

    _t = z_time_get_s();

    z_freea(group_pconds[g]);
  }

  msp_t[1] += z_time_get_s() - _t;

  if (sx == &_sx) elements_free(sx);

#ifdef MSP_VERIFY
  mpi_elements_check_order(s, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;

#undef doweights
}

#endif


slint_t mpi_sort_partition_exact_radix_2groups(elements_t *s, elements_t *sx, partcond_t *pcond, MPI_Comm group_comm, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_partition_exact_radix_2groups */
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

  int scounts[size], rcounts[size];
  int sdispls[size], rdispls[size];
  slint_t rcounts_total, nsubelements;

  elements_t _sx;

  slint_t sub_elements_sources[size], sub_elements_sizes[size];
  elements_t sub_elements[size];

  int *sub_sdispls;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif

#ifdef MSP_VERIFY
  slint_t rorders[2];
#endif


  Z_TRACE_IF(MSP_TRACE_IF, "s: %" slint_fmt " (max %" slint_fmt "), sx: %" slint_fmt " (max: %" slint_fmt ")", s->size, s->max_size, (sx)?sx->size:0, (sx)?sx->max_size:0);

  if (msp_sync) MPI_Barrier(comm);
  msp_t[0] = z_time_get_s();
  sort_radix(s, sx, rhigh, rlow, -1);
  msp_t[0] = z_time_get_s() - msp_t[0];

  if (msp_sync) MPI_Barrier(comm);
  msp_t[1] = z_time_get_s();

#ifdef elem_weight
  doweights = ((pcond->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

  /* make absolute part.-conds */
#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds(1, pcond, size, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

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

  Z_TRACE_IF(MSP_TRACE_IF, "%d: group: %d of %d, master: %d of %d", rank, group_rank, group_size, master_rank, master_size);

  nparts = master_size;

  /* create group partcond */
  group_pconds = z_alloca(group_size, sizeof(partcond_t));

  mpi_gather_partconds_grouped(pcond, group_comm, MPI_COMM_NULL, group_pconds, NULL, size, rank, comm);

  merge_partconds(group_pconds, group_size, &group_pcond);

  Z_TRACE_IF(MSP_TRACE_IF, "1st: select grouped");

  /* perform 1st grouped select */  
  mpi_select_exact_radix_grouped(s, 1, &group_pcond, master_comm, group_comm, rhigh, rlow, rwidth, SL_SORTED_IN, group_sdispls, size, rank, comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "group_sdispls = ", "%d  ", i, nparts, group_sdispls);*/

  /* create scounts from sdispls */
  displs2counts(nparts, group_sdispls, group_scounts, s->size);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "group_scounts = ", "%d  ", i, nparts, group_scounts);*/

  mpi_xcounts2ycounts_grouped(group_scounts, nparts, rcounts, group_comm, master_comm, size, rank, comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rcounts = ", "%d  ", i, size, rcounts);*/

#define SPARSE_SCOUNTS_FROM_RCOUNTS

  /* make scounts from rcounts */
#ifdef SPARSE_SCOUNTS_FROM_RCOUNTS
  mpi_xcounts2ycounts_sparse(rcounts, scounts, s->size, size, rank, comm);
#else
  mpi_xcounts2ycounts_all2all(rcounts, scounts, size, rank, comm);
#endif

#undef SPARSE_SCOUNTS_FROM_RCOUNTS

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "scounts = ", "%d  ", i, size, scounts);*/

  /* create displs from counts */
  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "sdispls = ", "%d  ", i, size, sdispls);
  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rdispls = ", "%d  ", i, size, rdispls);*/
  
  /* determine number of sub lists to receive */
  nsubelements = 0;
  for (i = 0; i < size; ++i) nsubelements += (rcounts[i] != 0);

  Z_TRACE_IF(MSP_TRACE_IF, "nsubelements = %" slint_fmt, nsubelements);

  rcounts_total = rdispls[size - 1] + rcounts[size - 1];

  Z_TRACE_IF(MSP_TRACE_IF, "rcounts_total = %" slint_fmt, rcounts_total);

  if (sx == NULL || sx->size < rcounts_total)
  {
    sx = &_sx;

    /* allocate elements for 1st redistribution (need keys and weights) */
    elements_alloc2(sx, rcounts_total, 1, 0, 0, elem_has_weight);
  }

  /* 1st redistribution */
#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  /* create sub lists */
  j = 0;
  for (i = 0; i < size; ++i)
  if (rcounts[i] != 0)
  {
    sub_elements_sources[j] = i;
    sub_elements_sizes[j] = rcounts[i];

    elem_assign_at(sx, rdispls[i], &sub_elements[j]);
    sub_elements[j].size = rcounts[i];
    
    Z_TRACE_IF(MSP_TRACE_IF, "sub-elements: %" slint_fmt ": source: %" slint_fmt ", size: %" slint_fmt ", displs: %d", j, sub_elements_sources[j], sub_elements_sizes[j], rdispls[i]);

    ++j;
  }

  sub_sdispls = z_alloca(group_size * nsubelements, sizeof(int));

  Z_TRACE_IF(MSP_TRACE_IF, "2nd: select with %" slint_fmt " sub-elements", nsubelements);

  /* perform 2nd select */
  mpi_select_exact_radix(sub_elements, nsubelements, group_size, group_pconds, rhigh, rlow, rwidth, SL_SORTED_IN, sub_sdispls, group_size, group_rank, group_comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "sub_displs = ", "%d  ", i, nsubelements * group_size, sub_sdispls);*/

  mpi_subxdispls2ycounts(nsubelements, sub_sdispls, sub_elements_sources, sub_elements_sizes, group_comm, group_size, rcounts, size, rank, comm);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rcounts = ", "%d  ", i, size, rcounts);*/

  z_freea(sub_sdispls);
  z_freea(group_pconds);

  if (master_comm != MPI_COMM_NULL) MPI_Comm_free(&master_comm);

  mpi_xcounts2ycounts_all2all(rcounts, scounts, size, rank, comm);

  msp_t[1] = z_time_get_s() - msp_t[1];

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "scounts = ", "%d  ", i, size, scounts);*/

  if (msp_sync) MPI_Barrier(comm);
  msp_t[2] = z_time_get_s();

  counts2displs(size, scounts, sdispls);
  counts2displs(size, rcounts, rdispls);

/*  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "rdispls = ", "%d  ", i, size, rdispls);
  Z_TRACE_ARRAY_IF(MSP_TRACE_IF, "sdispls = ", "%d  ", i, size, sdispls);*/

#define xelem_call \
  MPI_Alltoallv(xelem_buf(s), scounts, sdispls, xelem_mpi_datatype, xelem_buf(sx), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  msp_t[2] = z_time_get_s() - msp_t[2];

  s->size = sx->size = rdispls[size - 1] + rcounts[size - 1];

  if (msp_sync) MPI_Barrier(comm);
  msp_t[3] = z_time_get_s();
  mergep_heap_int(sx, s, size, rdispls, rcounts);
  msp_t[3] = z_time_get_s() - msp_t[3];

  if (sx == &_sx) elements_free(sx);

#ifdef MSP_VERIFY
  mpi_elements_check_order(s, 1, rorders, size, rank, comm);

# ifndef MSP_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    SL_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;

#undef doweights
}

#endif
