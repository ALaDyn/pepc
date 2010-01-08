/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_merge2_old.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


/*#define VALIDATE_ORDER
#define TRACE_MERGE2_INOUT*/

#ifndef VALIDATE_ORDER
 #undef TRACE_MERGE2_INOUT
#endif


slint pepcparts_mpi_merge2_old(elements_t *s, slint counterpart, slint high, m2x_func m2, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_merge2_old */
{
  slint xa_start, xa_end, xa;
  elements_t send, s0, s1;
  MPI_Status status;

#ifdef VALIDATE_ORDER
  static slint round = 0;
  slint o0, o1, of;
#ifdef TRACE_MERGE2_INOUT
  elements_t backup;
  char buffer[256];
#endif
#endif

  rti_tclear(rti_tid_mpi_merge2);

  if ((rank == counterpart) || (rank >= size) || (counterpart >= size)) return 0;

#ifdef VALIDATE_ORDER
  round++;
  if (o0 = pepcparts_elements_validate_order(s, 1)) printf("%d here: round = %d - input is unsorted (%d)\n", rank, round, o0);
#endif

  rti_tstart(rti_tid_mpi_merge2);

  rti_tstart(rti_tid_mpi_merge2_fe);
  xa = pepcparts_mpi_find_exact_old(s, counterpart, high, xs, &xa_start, &xa_end, size, rank, comm);
  rti_tstop(rti_tid_mpi_merge2_fe);

  /* exit, we have something to exchange */
  if (xa > 0)
  {
/*    printf("sending %d from %d to %d\n", xa, rank, counterpart);*/

    elem_assign_at(s, xa_start, &send);

    /* exchange elements */
    rti_tstart(rti_tid_mpi_merge2_xchg);

#define xelem_call \
    MPI_Sendrecv_replace(xelem_buf(&send), xa, xelem_mpi_datatype, counterpart, 0, counterpart, 0, comm, &status);
#include "sl_xelem_call.h"

/*    printf("sending %d from %d to %d DONE (%d == %d?)\n", xa, rank, counterpart, mpi_error, MPI_SUCCESS);*/

    rti_tstop(rti_tid_mpi_merge2_xchg);

    /* exit, if not all elemenents were exchanged */
    if (xa < s->size)
    {
      /* prepare the two sequences */
      if (rank == high) pepcparts_elements_extract(s, xa, &s0, &s1);
      else pepcparts_elements_extract(s, s->size - xa, &s0, &s1);

#ifdef VALIDATE_ORDER
      if (o0 = pepcparts_elements_validate_order(&s0, 1)) printf("%d here: round = %d - merge2_input_s0 is unsorted (%d)\n", rank, round, o0);
      if (o1 = pepcparts_elements_validate_order(&s1, 1)) printf("%d here: round = %d - merge2_input_s1 is unsorted (%d)\n", rank, round, o1);
#endif

#ifdef TRACE_MERGE2_INOUT
      pepcparts_elements_alloc(&backup, s->size, 1, 1);
      elem_ncopy(s, &backup, s->size);
#endif

      /* merge the 2 sequences local */
      rti_tstart(rti_tid_mpi_merge2_local);
      (m2)(&s0, &s1, xs);
      rti_tstop(rti_tid_mpi_merge2_local);

#ifdef VALIDATE_ORDER
      if (of = pepcparts_elements_validate_order(&s0, 1))
      {
        printf("%d here: round = %d - participants = %d & %d (%d) - merge2 (%d & %d) - output unsorted (%d & %d -> %d)\n", rank, round, rank, counterpart, high, s0.size, s1.size, o0, o1, of);

 #ifdef TRACE_MERGE2_INOUT

        sprintf(buffer, "failed_%d_%d_%d_%d_vs_%d_in", rank, round, s->size, s0.size, s1.size);
        elements_save_to_file(&backup, buffer);

        sprintf(buffer, "failed_%d_%d_%d_%d_vs_%d_out", rank, round, s->size, s0.size, s1.size);
        elements_save_to_file(s, buffer);

 #endif /* TRACE_MERGE2_INOUT */

      }

      pepcparts_elements_free(&backup);

#endif /* VALIDATE_ORDER */

    }
  }

  rti_tstop(rti_tid_mpi_merge2);

  return xa;
}
