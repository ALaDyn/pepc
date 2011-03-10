/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_partition_sample.c
 *  timestamp: 2011-01-30 00:28:58 +0100
 *  
 */


/* sl_macro MPS_TRACE_IF */


#include "sl_common.h"


/* config */
/*#define PRINT_SCOUNTS_RCOUNTS*/
/*#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/


#ifndef MPS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_partition_sample_regular(elements_t *s, partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_sample_regular */
{
  partcond_t pconds[size];

#ifdef PRINT_SCOUNTS_RCOUNTS
  slint_t i, j;
#endif
#ifdef PRINT_STATS
  int sdispls[size];
#endif


  rti_treset(rti_tid_mpi_partition_sample);          /* sl_tid */
  rti_treset(rti_tid_mpi_partition_sample_select);   /* sl_tid */
  rti_treset(rti_tid_mpi_partition_sample_rcounts);  /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_sample);

  rti_tstart(rti_tid_mpi_partition_sample_select);

  splitter_t sp;
  sp.displs = scounts;

  mpi_allgather_partconds(pcond, pconds, size, rank, comm);

  mpi_select_sample_regular(s, size, pconds, size - 1, &sp, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_sample_select);
  
  rti_tstart(rti_tid_mpi_partition_sample_rcounts);

  /* create scounts from sdispls */
  displs2counts(size, scounts, NULL, s->size);

  /* create rcounts if necessary */
  if (rcounts) mpi_xcounts2ycounts_all2all(scounts, rcounts, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_sample_rcounts);

  rti_tstop(rti_tid_mpi_partition_sample);

#ifdef PRINT_SCOUNTS_RCOUNTS
  printf("%d: scounts:", rank);
  for (i = 0, j = 0; i < size; ++i) { printf(" %d ", scounts[i]); j += scounts[i]; }
  printf(" = %" sl_int_type_fmt "\n", j);
  printf("%d: rcounts:", rank);
  if (rcounts)
  {
    for (i = 0, j = 0; i < size; ++i) { printf(" %d ", rcounts[i]); j += rcounts[i]; }
    printf(" = %" sl_int_type_fmt "\n", j);
  }
#endif

#ifdef PRINT_STATS
  counts2displs(size, scounts, sdispls);
  mpi_select_stats(s, size, sdispls, size, rank, comm);
#endif

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == PRINT_TIMINGS)
  {
    printf("%d: mpi_partition_sample: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_sample));
    printf("%d: mpi_partition_sample: select: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_sample_select));
    printf("%d: mpi_partition_sample: rcounts: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_sample_rcounts));
  }
#endif

  return 0;
}
