/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_merge2.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


/* sl_macro MM2_ELEMENTS_SENDRECV_REPLACE */
/* sl_macro MM2_TRACE_IF */
/* sl_macro MM2_PRINT_TIMINGS */

#define MM2_ELEMENTS_SENDRECV_REPLACE

#include "sl_common.h"


#ifdef SLDEBUG
# define CHECK_ORDER
#endif

#ifndef MM2_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MM2_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MM2_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_merge2(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, merge2x_f m2, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_merge2 */
{
  const int tag = 1;

  slint_t ex_start, ex_sizes[2], nx_move, ex_size;
  elements_t s0, s1;

  MPI_Status status;

#ifdef CHECK_ORDER
  slint_t check_order;
#endif


  Z_TRACE_IF(MM2_TRACE_IF, "starting mpi_merge2");

  /* sl_tid rti_tid_mpi_merge2 */

  rti_treset(rti_tid_mpi_merge2_find);       /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_moveright);  /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_exchange);   /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_moveleft);   /* sl_tid */
  rti_treset(rti_tid_mpi_merge2_local);     /* sl_tid */

  rti_tclear(rti_tid_mpi_merge2);

  if (other_rank < 0 || other_rank >= size) return -1;

  if (rank == other_rank) return 0;

  rti_tstart(rti_tid_mpi_merge2);

#ifdef CHECK_ORDER
  check_order = elements_validate_order(s, 1);
  if (check_order) Z_ERROR("input order failed at %" slint_fmt "", check_order);
#endif

  Z_TRACE_IF(MM2_TRACE_IF, "find_exact: s->size = %" slint_fmt ", other_rank / high_rank = %" slint_fmt " / %" slint_fmt, s->size, other_rank, high_rank);

  rti_tstart(rti_tid_mpi_merge2_find);
  mpi_find_exact(s, other_rank, high_rank, dst_size, &ex_start, ex_sizes, &nx_move, size, rank, comm);
  rti_tstop(rti_tid_mpi_merge2_find);

  Z_TRACE_IF(MM2_TRACE_IF, "find_exact: ex_start = %" slint_fmt ", ex_sizes = { %" slint_fmt ", %" slint_fmt " }, nx_move = %" slint_fmt, ex_start, ex_sizes[0], ex_sizes[1], nx_move);

  /* move the nx-block to the right (before exchange) */
  rti_tstart(rti_tid_mpi_merge2_moveright);

  if (nx_move > 0 && s->size - ex_sizes[0] > 0)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "moving right %" slint_fmt "", nx_move);

    if (rank != high_rank) elem_nmove_at(s, 0, s, nx_move, s->size - ex_sizes[0]);
    else elem_nmove_at(s, ex_sizes[0], s, ex_sizes[0] + nx_move, s->size - ex_sizes[0]);
  }

  rti_tstop(rti_tid_mpi_merge2_moveright);

  /* exchange elements */
  rti_tstart(rti_tid_mpi_merge2_exchange);

  elem_assign_at(s, ex_start, &s0);
  ex_size = z_min(ex_sizes[0], ex_sizes[1]);

  if (ex_size > 0)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "exchanging %" slint_fmt " elements at %" slint_fmt "", ex_size, ex_start);

#ifdef MM2_ELEMENTS_SENDRECV_REPLACE
    mpi_elements_sendrecv_replace(&s0, ex_size, other_rank, tag, other_rank, tag, size, rank, comm);
#else
#define xelem_call \
    MPI_Sendrecv_replace(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, other_rank, tag, comm, &status);
#include "sl_xelem_call.h"
#endif
  }

  elem_add(&s0, ex_size);

  if (ex_size < ex_sizes[0])
  {
    ex_size = ex_sizes[0] - ex_size;
    
    Z_TRACE_IF(MM2_TRACE_IF, "sending %" slint_fmt " at %" slint_fmt "", ex_size, (slint_t) (s0.keys - s->keys));

#define xelem_call \
    MPI_Send(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, comm);
#include "sl_xelem_call.h"

  } else if (ex_size < ex_sizes[1])
  {
    ex_size = ex_sizes[1] - ex_size;

    Z_TRACE_IF(MM2_TRACE_IF, "receiving %" slint_fmt " at %" slint_fmt "", ex_size, (slint_t) (s0.keys - s->keys));

#define xelem_call \
    MPI_Recv(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, comm, &status);
#include "sl_xelem_call.h"
  }

  rti_tstop(rti_tid_mpi_merge2_exchange);

  /* move the nx-block to the left (after exchange) */
  rti_tstart(rti_tid_mpi_merge2_moveleft);

  if (nx_move < 0 && s->size - ex_sizes[0] > 0)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "moving left %" slint_fmt "", nx_move);

    if (rank != high_rank) elem_nmove_at(s, 0, s, nx_move, s->size - ex_sizes[0]);
    else elem_nmove_at(s, ex_sizes[0], s, ex_sizes[0] + nx_move, s->size - ex_sizes[0]);
  }

  rti_tstop(rti_tid_mpi_merge2_moveleft);

  /* prepare the local merge2 */
  if (rank != high_rank)
  {
    elem_assign_at(s, 0, &s0);
    s0.size = s->size - ex_sizes[0];
    
    elem_assign_at(s, s0.size, &s1);
    s1.size = ex_sizes[1];

  } else
  {
    elem_assign_at(s, 0, &s0);
    s0.size = ex_sizes[1];
    
    elem_assign_at(s, s0.size, &s1);
    s1.size = s->size - ex_sizes[0];
  }

#ifdef CHECK_ORDER
  check_order = elements_validate_order(&s0, 1);
  if (check_order) Z_ERROR("intermediate lower order failed at %" slint_fmt "", check_order);
  check_order = elements_validate_order(&s1, 1);
  if (check_order) Z_ERROR("intermediate higher order failed at %" slint_fmt "", check_order);
#endif

  s->size = s0.size + s1.size;

  /* local merge */
  rti_tstart(rti_tid_mpi_merge2_local);

  if (s0.size > 0 && s1.size > 0 && m2 != NULL)
  {
    Z_TRACE_IF(MM2_TRACE_IF, "local merge2 %" slint_fmt " with %" slint_fmt "", s0.size, s1.size);

    m2(&s0, &s1, xs);
  }

  rti_tstop(rti_tid_mpi_merge2_local);

#ifdef CHECK_ORDER
  check_order = elements_validate_order(s, 1);
  if (check_order) Z_ERROR("output order failed at %" slint_fmt "", check_order);
#endif

  rti_tstop(rti_tid_mpi_merge2);

#if defined(MM2_PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (MM2_PRINT_TIMINGS)
  {
    printf("%d: mpi_merge2: %f\n", rank, rti_tlast(rti_tid_mpi_merge2));
    printf("%d: mpi_merge2: find: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_find));
    printf("%d: mpi_merge2: move-right: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_moveright));
    printf("%d: mpi_merge2: exchange: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_exchange));
    printf("%d: mpi_merge2: move-left: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_moveleft));
    printf("%d: mpi_merge2: local: %f\n", rank, rti_tlast(rti_tid_mpi_merge2_local));
  }
#endif

  return 0;
}
