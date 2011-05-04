/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_select_common.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


/* sl_macro MSC_TRACE_IF */


#include "sl_common.h"


#ifndef MSC_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSC_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSC_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t init_partconds(slint_t npconds, partcond_t *pconds, slint_t nparts, slint_t total_count, slweight_t total_weight) /* sl_proto, sl_func init_partconds */
{
  slint_t i;


  for (i = 0; i < npconds; ++i)
  {
    /* set default values and determine local (count/weight) limits */
    if (!(pconds[i].pcm & SLPC_COUNTS_MM)) { pconds[i].count_min = 0.0; pconds[i].count_max = total_count; }
    if (!(pconds[i].pcm & SLPC_COUNTS_LH)) { pconds[i].count_low = 0.0; pconds[i].count_high = total_count; }

    if (pconds[i].count_min < 0.0) pconds[i].count_min *= -total_count / nparts;
    if (pconds[i].count_max < 0.0) pconds[i].count_max *= -total_count / nparts;
    if (pconds[i].count_low < 0.0) pconds[i].count_low *= -total_count;
    if (pconds[i].count_high < 0.0) pconds[i].count_high *= -total_count;

/*    if (sl_mpi_rank == 0) printf("before: %e/%e - %e/%e\n", pconds[i].count_min, pconds[i].count_max, pconds[i].count_low, pconds[i].count_high);*/

    pconds[i].count_min = ceil(pconds[i].count_min);
    pconds[i].count_max = floor(pconds[i].count_max);
    pconds[i].count_low = ceil(pconds[i].count_low);
    pconds[i].count_high = floor(pconds[i].count_high);

/*    if (sl_mpi_rank == 0) printf("before: %e/%e - %e/%e\n", pconds[i].count_min, pconds[i].count_max, pconds[i].count_low, pconds[i].count_high);*/

#ifdef elem_weight
    if (!(pconds[i].pcm & SLPC_WEIGHTS_MM)) { pconds[i].weight_min = 0.0; pconds[i].weight_max = total_weight; }
    if (!(pconds[i].pcm & SLPC_WEIGHTS_LH)) { pconds[i].weight_low = 0.0; pconds[i].weight_high = total_weight; }

    if (pconds[i].weight_min < 0.0) pconds[i].weight_min *= -total_weight / nparts;
    if (pconds[i].weight_max < 0.0) pconds[i].weight_max *= -total_weight / nparts;
    if (pconds[i].weight_low < 0.0) pconds[i].weight_low *= -total_weight;
    if (pconds[i].weight_high < 0.0) pconds[i].weight_high *= -total_weight;

/*    if (sl_mpi_rank == 0) printf("before: %e/%e - %e/%e\n", pconds[i].weight_min, pconds[i].weight_max, pconds[i].weight_low, pconds[i].weight_high);*/
#endif
  }

  return 0;
}


slint_t init_partconds_intern(slint_t npconds, partcond_intern_t *pci, partcond_t *pc, slint_t nparts, slint_t total_count, slweight_t total_weight) /* sl_proto, sl_func init_partconds_intern */
{
  slint_t i;

  double avg_count = ((double) total_count) / nparts;
  slint_t sum_count[2] = { 0, 0 };
  slint_t tot_count[2] = { total_count, total_count };
#ifdef elem_weight
  slweight_t avg_weight = total_weight / nparts;
  slweight_t sum_weight[2] = { 0.0, 0.0 };
  slweight_t tot_weight[2] = { total_weight, total_weight };
#endif

  for (i = 0; i < npconds; ++i)
  {
    pci[i].pcm = pc[i].pcm;

    Z_TRACE_IF(MSC_TRACE_IF, "IN partcond %" slint_fmt ": %" slint_fmt, i, pc[i].pcm);

    Z_TRACE_IF(MSC_TRACE_IF, "IN partcond count %" slint_fmt ": %f / %f / %f / %f", i, pc[i].count_min, pc[i].count_max, pc[i].count_low, pc[i].count_high);

    if (pci[i].pcm & SLPC_COUNTS_MM)
    {
      pci[i].count_min = z_round((pc[i].count_min < 0.0)?(-pc[i].count_min * avg_count):pc[i].count_min);
      pci[i].count_max = z_round((pc[i].count_max < 0.0)?(-pc[i].count_max * avg_count):pc[i].count_max);

      /* check min/max consistency */
      if (pci[i].count_min > pci[i].count_max) pci[i].count_min = pci[i].count_max = (pci[i].count_min + pci[i].count_max) / 2;

      sum_count[0] += pci[i].count_min;
      sum_count[1] += pci[i].count_max;

    } else { pci[i].count_min = 0; pci[i].count_max = total_count; }

    if (pci[i].pcm & SLPC_COUNTS_LH)
    {
      pci[i].count_low  = z_round((pc[i].count_low  < 0.0)?(-pc[i].count_low  * total_count):pc[i].count_low);
      pci[i].count_high = z_round((pc[i].count_high < 0.0)?(-pc[i].count_high * total_count):pc[i].count_high);
      
      /* FIXME: low/high consistency not checked */

    } else { pci[i].count_low = 0; pci[i].count_high = total_count; }

    Z_TRACE_IF(MSC_TRACE_IF, "OUT partcond count %" slint_fmt ": %" slint_fmt " / %" slint_fmt " / %" slint_fmt " / %" slint_fmt, i, pci[i].count_min, pci[i].count_max, pci[i].count_low, pci[i].count_high);

#ifdef elem_weight
    Z_TRACE_IF(MSC_TRACE_IF, "IN partcond weight %" slint_fmt ": %f / %f / %f / %f", i, pc[i].weight_min, pc[i].weight_max, pc[i].weight_low, pc[i].weight_high);

    if (pci[i].pcm & SLPC_WEIGHTS_MM)
    {
      pci[i].weight_min = (pc[i].weight_min < 0.0)?(-pc[i].weight_min * avg_weight):pc[i].weight_min;
      pci[i].weight_max = (pc[i].weight_max < 0.0)?(-pc[i].weight_max * avg_weight):pc[i].weight_max;
      
      /* check min/max consistency */
      if (pci[i].weight_min > pci[i].weight_max) pci[i].weight_min = pci[i].weight_max = (pci[i].weight_min + pci[i].weight_max) / 2.0;

      sum_weight[0] += pci[i].weight_min;
      sum_weight[1] += pci[i].weight_max;

    } else { pci[i].weight_min = 0.0; pci[i].weight_max = total_weight; }

    if (pci[i].pcm & SLPC_WEIGHTS_LH)
    {
      pci[i].weight_low  = (pc[i].weight_low  < 0.0)?(-pc[i].weight_low  * total_weight):pc[i].weight_low;
      pci[i].weight_high = (pc[i].weight_high < 0.0)?(-pc[i].weight_high * total_weight):pc[i].weight_high;

      /* FIXME: low/high consistency not checked */

    } else { pci[i].weight_low = 0; pci[i].weight_high = total_weight; }

    Z_TRACE_IF(MSC_TRACE_IF, "OUT partcond weight %" slint_fmt ": %f / %f / %f / %f", i, pci[i].weight_min, pci[i].weight_max, pci[i].weight_low, pci[i].weight_high);
#endif
  }

  if (!(pci[0].pcm & SLPC_COUNTS_MM)) total_count = 0;

  Z_TRACE_IF(MSC_TRACE_IF, "total_count = %" slint_fmt ", sum_count = %" slint_fmt " / %" slint_fmt, total_count, sum_count[0], sum_count[1]);

#ifdef elem_weight
  if (!(pci[0].pcm & SLPC_WEIGHTS_MM)) total_weight = 0.0;

  Z_TRACE_IF(MSC_TRACE_IF, "total_weight = %" slweight_fmt ", sum_weight = %" slweight_fmt " / %" slweight_fmt, total_weight, sum_weight[0], sum_weight[1]);
#endif

  if (sum_count[0] <= total_count) sum_count[0] = -1;
  if (sum_count[1] >= total_count) sum_count[1] = -1;
#ifdef elem_weight
  if (sum_weight[0] <= total_weight) sum_weight[0] = -1;
  if (sum_weight[1] >= total_weight) sum_weight[1] = -1;
#endif

  if (sum_count[0] > 0 || sum_count[1] > 0
#ifdef elem_weight
   || sum_weight[0] > 0.0 || sum_weight[1] > 0.0
#endif
   )
  for (i = 0; i < npconds; ++i)
  {
    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt ": in: count_min/max: %" slint_fmt " / %" slint_fmt, i, pci[i].count_min, pci[i].count_max);

    if (sum_count[0] > 0)
    {
      sum_count[0] -= pci[i].count_min;
      pci[i].count_min = (tot_count[0] * pci[i].count_min) / (sum_count[0] + pci[i].count_min);
      tot_count[0] -= pci[i].count_min;
    }
    if (sum_count[1] > 0)
    {
      sum_count[1] -= pci[i].count_max;
      pci[i].count_max = (tot_count[1] * pci[i].count_max) / (sum_count[1] + pci[i].count_max);
      tot_count[1] -= pci[i].count_max;
    }

    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt ": out: count_min/max: %" slint_fmt " / %" slint_fmt, i, pci[i].count_min, pci[i].count_max);

#ifdef elem_weight
    if (sum_weight[0] > 0)
    {
      sum_weight[0] -= pci[i].weight_min;
      pci[i].weight_min = (tot_weight[0] * pci[i].weight_min) / (sum_weight[0] + pci[i].weight_min);
      tot_weight[0] -= pci[i].weight_min;
    }
    if (sum_weight[1] > 0)
    {
      sum_weight[1] -= pci[i].weight_max;
      pci[i].weight_max = (tot_weight[1] * pci[i].weight_max) / (sum_weight[1] + pci[i].weight_max);
      tot_weight[1] -= pci[i].weight_max;
    }
#endif
  }

  return 0;
}


slint_t merge_partconds(partcond_t *pconds_in, slint_t npconds_in, partcond_t *pcond_out)  /* sl_proto, sl_func merge_partconds */
{
  slint_t i;

  pcond_out->pcm = pconds_in[0].pcm;

  pcond_out->count_min = 0;
  pcond_out->count_max = 0;
  pcond_out->count_low = pconds_in[0].count_low;
  pcond_out->count_high = pconds_in[npconds_in - 1].count_high;
  for (i = 0; i < npconds_in; ++i)
  {
    pcond_out->count_min += pconds_in[i].count_min;
    pcond_out->count_max += pconds_in[i].count_max;
  }

#ifdef elem_weight  
  pcond_out->weight_min = 0;
  pcond_out->weight_max = 0;
  pcond_out->weight_low = pconds_in[0].weight_low;
  pcond_out->weight_high = pconds_in[npconds_in - 1].weight_high;
  for (i = 0; i < npconds_in; ++i)
  {
    pcond_out->weight_min += pconds_in[i].weight_min;
    pcond_out->weight_max += pconds_in[i].weight_max;
  }
#endif

  return 0;
}


slint_t mpi_gather_partconds_grouped(partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, partcond_t *pconds_out, slint_t *npconds_out, int size, int rank, MPI_Comm comm)  /* sl_proto, sl_func mpi_gather_partconds_grouped */
{
  int _npconds_out = 1;
#ifdef HAVENT_MPI_IN_PLACE
  partcond_t _pcond_in;
#endif

  if (npconds_out) _npconds_out = *npconds_out;

  if (_npconds_out < 0)
  {
    if (pcond_in_comm != MPI_COMM_NULL) MPI_Comm_size(pcond_in_comm, &_npconds_out);
    if (pconds_out_comm != MPI_COMM_NULL) MPI_Bcast(&_npconds_out, 1, MPI_INT, 0, pconds_out_comm);
  }
  
  if (pcond_in_comm != MPI_COMM_NULL)
  {
    if (pconds_out)
    {
      if (pcond_in) MPI_Allgather(pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, pcond_in_comm);
/*      else
      {
#ifdef HAVENT_MPI_IN_PLACE
        _pcond_in = pcond_out[?];
        MPI_Allgather(&_pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, pcond_in_comm);
#else
        MPI_Allgather(MPI_IN_PLACE, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, pcond_in_comm);
#endif
      }*/
    }

  }

  if (pconds_out_comm != MPI_COMM_NULL && pconds_out) MPI_Bcast(pconds_out, _npconds_out * sizeof(partcond_t), MPI_BYTE, 0, pconds_out_comm);

  if (npconds_out) *npconds_out = _npconds_out;

  return 0;
}


slint_t mpi_gather_partconds(partcond_t *pcond_in, partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_gather_partconds */
{
/*  if (comm != MPI_COMM_NULL)*/
  MPI_Gather(pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, root, comm);

  return 0;
}


slint_t mpi_allgather_partconds(partcond_t *pcond_in, partcond_t *pconds_out, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_allgather_partconds */
{
/*  if (comm != MPI_COMM_NULL)*/
  MPI_Allgather(pcond_in, sizeof(partcond_t), MPI_BYTE, pconds_out, sizeof(partcond_t), MPI_BYTE, comm);

  return 0;
}


slint_t mpi_bcast_partconds(slint_t npconds, partcond_t *pconds, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_bcast_partconds */
{
/*  if (comm != MPI_COMM_NULL)*/
  MPI_Bcast(pconds, npconds * sizeof(partcond_t), MPI_BYTE, root, comm);

  return 0;
}


slint_t mpi_post_check_partconds(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_post_check_partconds */
{
  slint_t i, j;
  int ssdispls[nparts], sscounts[nparts], rrcounts[nparts], ss;

  for (i = 0; i < nparts; ++i)
  {
    ssdispls[i] = 0;
    for (j = 0; j < nelements; ++j) ssdispls[i] += sdispls[i * nelements + j];
  }

  ss = 0;
  for (j = 0; j < nelements; ++j) ss += s[j].size;
  displs2counts(nparts, ssdispls, sscounts, ss);

  MPI_Reduce(sscounts, rrcounts, nparts,  MPI_INT, MPI_SUM, 0, comm);

  j = -1;
  if (rank == 0)
  for (i = 0; i < nparts; ++i)
  {
    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt " verifying %d against [%f  %f]", i, rrcounts[i], pconds[i].count_min, pconds[i].count_max);
  
    if (rrcounts[i] < pconds[i].count_min || rrcounts[i] > pconds[i].count_max) j = i;
  }

  MPI_Bcast(&j, 1, int_mpi_datatype, 0, comm);
  
  return j;
}


slint_t mpi_post_check_partconds_intern(elements_t *s, slint_t nelements, slint_t nparts, partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_post_check_partconds_intern */
{
  slint_t i, j;
  int ssdispls[nparts], sscounts[nparts], rrcounts[nparts], ss;

  for (i = 0; i < nparts; ++i)
  {
    ssdispls[i] = 0;
    for (j = 0; j < nelements; ++j) ssdispls[i] += sdispls[i * nelements + j];
  }

  ss = 0;
  for (j = 0; j < nelements; ++j) ss += s[j].size;
  displs2counts(nparts, ssdispls, sscounts, ss);

  MPI_Reduce(sscounts, rrcounts, nparts,  MPI_INT, MPI_SUM, 0, comm);

  j = -1;
  if (rank == 0)
  for (i = 0; i < nparts; ++i)
  {
    Z_TRACE_IF(MSC_TRACE_IF, "%" slint_fmt " verifying %d against [%" slint_fmt "  %" slint_fmt "]", i, rrcounts[i], pci[i].count_min, pci[i].count_max);
  
    if (rrcounts[i] < pci[i].count_min || rrcounts[i] > pci[i].count_max) j = i;
  }

  MPI_Bcast(&j, 1, int_mpi_datatype, 0, comm);
  
  return j;
}


slint_t mpi_select_stats(elements_t *s, slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_stats */
{
  slint_t i;
#ifdef elem_weight
  slint_t j;
#endif

  double v, vmin, vmax;

  slint_t partial_counts[nparts + 1];
#ifdef HAVENT_MPI_IN_PLACE
  slint_t partial_counts2[nparts + 1];
#endif

#ifdef elem_weight
  slweight_t partial_weights[nparts + 1];
# ifdef HAVENT_MPI_IN_PLACE
  slweight_t partial_weights2[nparts + 1];
# endif
#endif


  partial_counts[nparts] = 0;
#ifdef elem_weight
  partial_weights[nparts] = 0.0;
#endif

  for (i = 0; i < nparts; ++i)
  {
    partial_counts[i] = ((i < nparts - 1)?sdispls[i + 1]:s->size) - sdispls[i];
    partial_counts[nparts] += partial_counts[i];

#ifdef elem_weight
    partial_weights[i] = 0.0;
    for (j = sdispls[i]; j < ((i < nparts - 1)?sdispls[i + 1]:s->size); ++j) partial_weights[i] += elem_weight(s, j);
    partial_weights[nparts] += partial_weights[i];
#endif
  }

#ifdef HAVENT_MPI_IN_PLACE
  MPI_Reduce(partial_counts, partial_counts2, nparts + 1, int_mpi_datatype, MPI_SUM, 0, comm);
# define partial_counts   partial_counts2
# ifdef elem_weight
  MPI_Reduce(partial_weights, partial_weights2, nparts + 1, weight_mpi_datatype, MPI_SUM, 0, comm);
# define partial_weights  partial_weights2
# endif
#else
  /* recvbuf requires workaround for an in-place/aliased-buffer-check-bug in mpich2 (fixed with rev 5518) */
  MPI_Reduce((rank == 0)?MPI_IN_PLACE:partial_counts, (rank == 0)?partial_counts:NULL, nparts + 1, int_mpi_datatype, MPI_SUM, 0, comm);
# ifdef elem_weight
  MPI_Reduce((rank == 0)?MPI_IN_PLACE:partial_weights, (rank == 0)?partial_weights:NULL, nparts + 1, weight_mpi_datatype, MPI_SUM, 0, comm);
# endif
#endif

  if (rank == 0)
  {
    printf("%d: count total: %" sl_int_type_fmt "\n", rank, partial_counts[nparts]);
    v = 0.0;
    vmin = 1.0;
    vmax = 0.0;
    for (i = 0; i < nparts; ++i)
    {
/*      printf("%d: %" sl_int_type_fmt " %" sl_int_type_fmt " / %f - %" sl_int_type_fmt " / %f\n", rank, i, partial_counts[i], (double) partial_counts[i] / partial_counts[nparts], (partial_counts[nparts] / nparts) - partial_counts[i], fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts])));*/
      v += fabs((partial_counts[nparts] / nparts) - partial_counts[i]);
      if (fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts])) < vmin) vmin = fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts]));
      if (fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts])) > vmax) vmax = fabs(1.0 - ((double) partial_counts[i] * nparts / partial_counts[nparts]));
    }
    printf("%d: count min/max: %f / %f\n", rank, vmin, vmax);
    printf("%d: count average: %f - %f / %f\n", rank, (double) partial_counts[nparts] / nparts, v / nparts, v / partial_counts[nparts]);

#ifdef elem_weight
    printf("%d: weight total: %" slweight_fmt "\n", rank, partial_weights[nparts]);
    v = 0.0;
    vmin = 1.0;
    vmax = 0.0;
    for (i = 0; i < nparts; ++i)
    {
/*      printf("%d: %" sl_int_type_fmt " %f / %f - %f / %f\n", rank, i, partial_weights[i], partial_weights[i] / partial_weights[nparts], (partial_weights[nparts] / nparts) - partial_weights[i], fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts])));*/
      v += fabs((partial_weights[nparts] / nparts) - partial_weights[i]);
      if (fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts])) < vmin) vmin = fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts]));
      if (fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts])) > vmax) vmax = fabs(1.0 - (partial_weights[i] * nparts / partial_weights[nparts]));
    }
    printf("%d: weight min/max: %f / %f\n", rank, vmin, vmax);
    printf("%d: weight average: %f - %f / %f\n", rank, partial_weights[nparts] / nparts, v / nparts, v / partial_weights[nparts]);
#endif
  }
  
  return 0;
}
