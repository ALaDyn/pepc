/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_partition_exact_generic.c
 *  timestamp: 2011-03-11 09:06:54 +0100
 *  
 */


/* sl_macro MPEG_TRACE_IF */


#include "sl_common.h"


/* config */
/*#define PRINT_SCOUNTS_RCOUNTS
#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/


#ifndef MPEG_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MPEG_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MPEG_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_partition_exact_generic(elements_t *s, partcond_t *pcond, binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_generic */
{
  splitter_t sp;
#ifdef PRINT_SCOUNTS_RCOUNTS
  slint_t i, j;
#endif
#ifdef PRINT_STATS
  int sdispls[size];
#endif


  rti_treset(rti_tid_mpi_partition_exact_generic);          /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_generic_select);   /* sl_tid */
  rti_treset(rti_tid_mpi_partition_exact_generic_rcounts);  /* sl_tid */


  rti_tstart(rti_tid_mpi_partition_exact_generic);

  rti_tstart(rti_tid_mpi_partition_exact_generic_select);

  sp.displs = scounts;
  mpi_select_exact_generic_grouped(s, 1, pcond, comm, MPI_COMM_NULL, bm, &sp, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_exact_generic_select);

  rti_tstart(rti_tid_mpi_partition_exact_generic_rcounts);

  /* create scounts from sdispls */
  displs2counts(size, scounts, NULL, s->size);

  /* create rcounts if necessary */
  if (rcounts) mpi_xcounts2ycounts_all2all(scounts, rcounts, size, rank, comm);

  rti_tstop(rti_tid_mpi_partition_exact_generic_rcounts);

  rti_tstop(rti_tid_mpi_partition_exact_generic);

#ifdef PRINT_SCOUNTS_RCOUNTS
  printf("%d: scounts:", rank);
  for (i = 0, j = 0; i < size; ++i) { printf(" %d ", scounts[i]); j += scounts[i]; }
  printf(" = %" sl_int_type_fmt "\n", j);
  printf("%d: rcounts:", rank);
  if (rcounts) for (i = 0, j = 0; i < size; ++i) { printf(" %d ", rcounts[i]); j += rcounts[i]; }
  printf("\n");
#endif

#ifdef PRINT_STATS
  counts2displs(size, scounts, sdispls);
  mpi_select_stats(s, size, sdispls, size, rank, comm);
#endif

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == PRINT_TIMINGS)
  {
    printf("%d: mpi_partition_exact_generic: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_generic));
    printf("%d: mpi_partition_exact_generic: select: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_generic_select));
    printf("%d: mpi_partition_exact_generic: rcounts: %f\n", rank, rti_tcumu(rti_tid_mpi_partition_exact_generic_rcounts));
  }
#endif

  return 0;
}


#if 0
slint_t mpi_partition_exact_generic2(elements_t *s, partcond_t *pcond, binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_partition_exact_generic2 */
#endif
