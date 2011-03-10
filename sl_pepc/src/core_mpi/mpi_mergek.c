/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_mergek.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


/* sl_macro MMK_TRACE_IF */
/* sl_macro MMK_TRACE_FAIL */
/* sl_macro MMK_EQUAL_SYNC */
/* sl_macro MMK_SYNC */
/* sl_macro MMK_PRINT_TIMINGS */


#include "sl_common.h"


#ifdef SLDEBUG
# define CHECK_ORDER
#endif

/*#define CHECK_ORDER_BREAK*/


#ifndef MMK_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MMK_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MMK_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_mergek_equal(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_mergek_equal */
{
  slint_t stage, other_rank, up, high_rank;

#ifdef CHECK_ORDER
  slint_t local_order, global_order;
#endif

#ifdef MMK_TRACE_FAIL
  elements_t sin;
  char trace_file[128];
#endif


  Z_TRACE_IF(MMK_TRACE_IF, "starting mpi_mergek_equal");

  /* sl_tid rti_tid_mpi_mergek_equal */

  rti_treset(rti_tid_mpi_mergek_equal_while);         /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_equal_while_merge2);  /* sl_tid */

  rti_tstart(rti_tid_mpi_mergek_equal);

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

#ifdef MMK_EQUAL_SYNC
    MPI_Barrier(comm);
#endif

  stage = 0;

  rti_tstart(rti_tid_mpi_mergek_equal_while);

  if (size > 1)
  while ((other_rank = (sn)(size, rank, stage, snd, &up)) >= 0)
  {
    Z_TRACE_IF(MMK_TRACE_IF, "stage: %" slint_fmt ", %d %s %" slint_fmt, stage, rank, ((rank == z_min(rank, other_rank))?(up?"<-":"->"):(up?"->":"<-")), other_rank);

#ifdef CHECK_ORDER
    local_order = elements_validate_order(s, 1);
    if (local_order) Z_ERROR("input order failed at %" slint_fmt "", local_order);

    local_order = (local_order)?1:0;
    MPI_Allreduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_SUM, comm);
    if (global_order)
    {
# ifdef CHECK_ORDER_BREAK
      Z_ERROR("break: input order failed on %" slint_fmt " process(es)", global_order);
      break;
# endif
    }
#endif

#ifdef MMK_TRACE_FAIL
    elements_alloc(&sin, s->size, SLCM_ALL);
    elem_ncopy(s, &sin, s->size);
#endif

    rti_tstart(rti_tid_mpi_mergek_equal_while_merge2);
    high_rank = (up)?(z_min(rank, other_rank)):(z_max(rank, other_rank));
    mpi_merge2(s, other_rank, high_rank, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_equal_while_merge2);

#ifdef CHECK_ORDER
    local_order = elements_validate_order(s, 1);
    if (local_order) Z_ERROR("output order failed at %" slint_fmt "", local_order);

    local_order = (local_order)?1:0;
    MPI_Allreduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_SUM, comm);
    if (global_order)
    {
#ifdef MMK_TRACE_FAIL
      sprintf(trace_file, "mergek_trace_rank%d_stage%" slint_fmt "_other%" slint_fmt "_high%" slint_fmt ".in", rank, stage, other_rank, high_rank);
      elements_save_keys_to_file(&sin, trace_file);
      sprintf(trace_file, "mergek_trace_rank%d_stage%" slint_fmt "_other%" slint_fmt "_high%" slint_fmt ".out", rank, stage, other_rank, high_rank);
      elements_save_keys_to_file(s, trace_file);

# ifdef CHECK_ORDER_BREAK
      elements_free(&sin);
# endif
#endif

# ifdef CHECK_ORDER_BREAK
      Z_ERROR("break: output order failed on %" slint_fmt " process(es)", global_order);
      break;
# endif
    }
#endif

#ifdef MMK_TRACE_FAIL
    elements_free(&sin);
#endif

#ifdef MMK_EQUAL_SYNC
    MPI_Barrier(comm);
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek_equal_while);

  rti_tstop(rti_tid_mpi_mergek_equal);

#if defined(MMK_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MMK_PRINT_TIMINGS)
  {
    printf("%d: mpi_mergek_equal: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_equal));
    printf("%d: mpi_mergek_equal: while: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_equal_while));
    printf("%d: mpi_mergek_equal:  merge2: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_equal_while_merge2));
    printf("%d: mpi_mergek_equal: stages: %" slint_fmt "\n", rank, stage);
  }
#endif

  return stage;
}


slint_t mpi_mergek(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_mergek */
{
  slint_t stage, done;

  /* FIXME: use something more efficient that works for all kind of (non-equal) distributions? */


  Z_TRACE_IF(MMK_TRACE_IF, "starting mpi_mergek");

  /* sl_tid rti_tid_mpi_mergek */

  rti_treset(rti_tid_mpi_mergek_equalike);       /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_while);          /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_while_check);    /* sl_tid */
  rti_treset(rti_tid_mpi_mergek_while_oddeven);  /* sl_tid */

  rti_tstart(rti_tid_mpi_mergek);

  /* assume equal distribution */
  rti_tstart(rti_tid_mpi_mergek_equalike);
  stage = mpi_mergek_equal(s, sn, snd, m2x, xs, size, rank, comm);
  rti_tstop(rti_tid_mpi_mergek_equalike);

/*  return finalize;*/

#ifdef MMK_SYNC
  MPI_Barrier(comm);
#endif

  rti_tstart(rti_tid_mpi_mergek_while);

  /* use alternating odd/even rounds (simliar to odd-even-tranpose) until sorted */
  if (size > 1)
  while (1)
  {
    rti_tstart(rti_tid_mpi_mergek_while_check);
    done = mpi_check_global_order(key_purify(s->keys[0]), key_purify(s->keys[s->size - 1]), -1, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_while_check);

    if (done) break;

    Z_TRACE_IF(MMK_TRACE_IF, "stage: %" slint_fmt " order failed, do %s", stage, ((stage % 2)?"odd":"even"));

    rti_tstart(rti_tid_mpi_mergek_while_oddeven);
    mpi_mergek_equal(s, (stage % 2)?sn_odd:sn_even, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_while_oddeven);

#ifdef MMK_SYNC
    MPI_Barrier(comm);
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek_while);

  rti_tstop(rti_tid_mpi_mergek);

#if defined(MMK_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MMK_PRINT_TIMINGS)
  {
    printf("%d: mpi_mergek: %f\n", rank, rti_tlast(rti_tid_mpi_mergek));
    printf("%d: mpi_mergek: equal: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_equalike));
    printf("%d: mpi_mergek: while: %f\n", rank, rti_tlast(rti_tid_mpi_mergek_while));
    printf("%d: mpi_mergek:  check: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_while_check));
    printf("%d: mpi_mergek:  oddeven: %f\n", rank, rti_tcumu(rti_tid_mpi_mergek_while_oddeven));
    printf("%d: mpi_mergek: stages: %" slint_fmt "\n", rank, stage);
  }
#endif

  return stage;
}


slint_t mpi_mergek_equal2(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms) /* sl_proto, sl_func mpi_mergek_equal2 */
{
  slint_t stage = 0, other_ranks[2], ups[2], high_ranks[2];

  rti_tstart(rti_tid_mpi_mergek);

  if ((comms[0] != MPI_COMM_NULL && sizes[0] > 1) || (comms[1] != MPI_COMM_NULL && sizes[1] > 1))
  while (1)
  {
    other_ranks[0] = sn(sizes[0], ranks[0], stage, snd, &ups[0]);
    other_ranks[1] = sn(sizes[1], ranks[1], stage, snd, &ups[1]);
    
/*    printf("%d: doing merge2 %d vs. %d and %d vs. %d\n", sl_mpi_rank, ranks[0], other_ranks[0], ranks[1], other_ranks[1]);*/

    if ((other_ranks[0] < 0) && (other_ranks[1] < 0)) break;

    if (other_ranks[0] >= 0)
    {
      high_ranks[0] = (ups[0])?(z_min(ranks[0], other_ranks[0])):(z_max(ranks[0], other_ranks[0]));
      Z_TRACE_IF(MMK_TRACE_IF, "#0 merge2 with %" slint_fmt ", high: %" slint_fmt ", size = %d, rank = %d", other_ranks[0], high_ranks[0], sizes[0], ranks[0]);
      mpi_merge2(&s[0], other_ranks[0], high_ranks[0], NULL, m2x, xs, sizes[0], ranks[0], comms[0]);
    }

    if (other_ranks[1] >= 0)
    {
      high_ranks[1] = (ups[1])?(z_min(ranks[1], other_ranks[1])):(z_max(ranks[1], other_ranks[1]));
      Z_TRACE_IF(MMK_TRACE_IF, "#1 merge2 with %" slint_fmt ", high: %" slint_fmt ", size = %d, rank = %d", other_ranks[1], high_ranks[1], sizes[1], ranks[1]);
      mpi_merge2(&s[1], other_ranks[1], high_ranks[1], NULL, m2x, xs, sizes[1], ranks[1], comms[1]);
    }

    rti_tstop(rti_tid_mpi_mergek_merge2);

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek);

  return stage;
}
