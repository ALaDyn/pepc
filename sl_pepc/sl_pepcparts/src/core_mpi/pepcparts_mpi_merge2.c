/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_merge2.c
 *  timestamp: 2009-11-12 08:41:36 +0100
 *  
 */


#include "sl_common.h"


#ifdef DEBUG
# define CHECK_ORDER
#endif


slint_t pepcparts_mpi_merge2(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, merge2x_f m2, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_merge2 */
{
  const int tag = 1;

  slint_t ex_start, ex_sizes[2], nx_move, ex_size;
  elements_t s0, s1;

  MPI_Status status;

#ifdef CHECK_ORDER
  slint_t check_order;
#endif


  rti_tclear(rti_tid_mpi_merge2);

  if (other_rank < 0 || other_rank >= size) return -1;

  if (rank == other_rank) return 0;

  rti_tstart(rti_tid_mpi_merge2);

#ifdef CHECK_ORDER
  check_order = pepcparts_elements_validate_order(s, 1);
  if (check_order) fprintf(stderr, "%d: pepcparts_mpi_merge2: warning: input order failed (at %" sl_int_type_fmt ")\n", rank, check_order);
#endif

  rti_tstart(rti_tid_mpi_merge2_fe);
  pepcparts_mpi_find_exact(s, other_rank, high_rank, dst_size, &ex_start, ex_sizes, &nx_move, size, rank, comm);
  rti_tstop(rti_tid_mpi_merge2_fe);

/*  printf("%d here: ex_start = %ld, ex_sizes = { %ld, %ld }, nx_move = %ld\n", rank, ex_start, ex_sizes[0], ex_sizes[1], nx_move);*/

  /* move the nx-block to the right (before exchange) */
  if (nx_move > 0 && s->size - ex_sizes[0] > 0)
  {
/*    printf("%d here: nx_move = %ld\n", rank, nx_move);*/
    if (rank != high_rank) elem_nmove_at(s, 0, s, nx_move, s->size - ex_sizes[0]);
    else elem_nmove_at(s, ex_sizes[0], s, ex_sizes[0] + nx_move, s->size - ex_sizes[0]);
  }


  /* exchange elements */
  rti_tstart(rti_tid_mpi_merge2_xchg);

  elem_assign_at(s, ex_start, &s0);
  ex_size = xmin(ex_sizes[0], ex_sizes[1]);

  if (ex_size > 0)
  {
/*    printf("%d here: exchanging %ld elements at %ld\n", rank, ex_size, ex_start);*/

#define xelem_call \
    MPI_Sendrecv_replace(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, other_rank, tag, comm, &status);
#include "sl_xelem_call.h"
  }

  elem_add(&s0, ex_size);

  if (ex_size < ex_sizes[0])
  {
    ex_size = ex_sizes[0] - ex_size;
    
/*    printf("%d here: sending %ld at %d\n", rank, ex_size, s0.keys - s->keys);*/

#define xelem_call \
    MPI_Send(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, comm);
#include "sl_xelem_call.h"

  } else if (ex_size < ex_sizes[1])
  {
    ex_size = ex_sizes[1] - ex_size;

/*    printf("%d here: receiving %ld at %d\n", rank, ex_size, s0.keys - s->keys);*/

#define xelem_call \
    MPI_Recv(xelem_buf(&s0), ex_size, xelem_mpi_datatype, other_rank, tag, comm, &status);
#include "sl_xelem_call.h"
  }

  rti_tstop(rti_tid_mpi_merge2_xchg);


  /* move the nx-block to the right (before exchange) */
  if (nx_move < 0 && s->size - ex_sizes[0] > 0)
  {
/*    printf("%d here: nx_move = %ld\n", rank, nx_move);*/

    if (rank != high_rank) elem_nmove_at(s, 0, s, nx_move, s->size - ex_sizes[0]);
    else elem_nmove_at(s, ex_sizes[0], s, ex_sizes[0] + nx_move, s->size - ex_sizes[0]);
  }


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
  check_order = pepcparts_elements_validate_order(&s0, 1);
  if (check_order) fprintf(stderr, "%d: pepcparts_mpi_merge2: warning: intermediate lower order failed (at %" sl_int_type_fmt ")\n", rank, check_order);
  check_order = pepcparts_elements_validate_order(&s1, 1);
  if (check_order) fprintf(stderr, "%d: pepcparts_mpi_merge2: warning: intermediate higher order failed (at %" sl_int_type_fmt ")\n", rank, check_order);
#endif

  s->size = s0.size + s1.size;

  /* local merge */
  if (s0.size > 0 && s1.size > 0 && m2 != NULL)
  {
    rti_tstart(rti_tid_mpi_merge2_local);
    (m2)(&s0, &s1, xs);
    rti_tstop(rti_tid_mpi_merge2_local);
  }

#ifdef CHECK_ORDER
  check_order = pepcparts_elements_validate_order(s, 1);
  if (check_order) fprintf(stderr, "%d: pepcparts_mpi_merge2: warning: output order failed (at %" sl_int_type_fmt ")\n", rank, check_order);
#endif

  rti_tstop(rti_tid_mpi_merge2);

  return 0;
}
