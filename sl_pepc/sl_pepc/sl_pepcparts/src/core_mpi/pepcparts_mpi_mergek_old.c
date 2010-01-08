/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_mergek_old.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


/*#define VALIDATE_ORDER_BREAK*/


slint pepcparts_mpi_mergek_old(elements_t *s, sn_func sn, void *snp, m2x_func m2, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_mergek_old */
{
  slint stage = 0, cp, up;

#ifdef VALIDATE_ORDER_BREAK
  slint local_error, global_error;
#endif

  rti_tstart(rti_tid_mpi_mergek);

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

  if (size > 1)
  while ((cp = (sn)(size, rank, stage++, snp, &up)) >= 0)
  {

#ifdef VALIDATE_ORDER_BREAK

    local_error = pepcparts_elements_validate_order(s, 1);
    MPI_Allreduce(&local_error, &global_error, 1, sl_int_type_mpi, MPI_SUM, comm);
    if (global_error) { printf("%d here: unsorted input, me is %d\n", rank, local_error); break; }

#endif /* VALIDATE_ORDER_BREAK */

    rti_tstart(rti_tid_mpi_mergek_merge2);
    pepcparts_mpi_merge2_old(s, cp, (up)?(xmin(rank, cp)):(xmax(rank,cp)), m2, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_merge2);

#ifdef VALIDATE_ORDER_BREAK

    local_error = pepcparts_elements_validate_order(s, 1);
    MPI_Allreduce(&local_error, &global_error, 1, sl_int_type_mpi, MPI_SUM, comm);
    if (global_error) { printf("%d here: unsorted output, me is %d\n", rank, local_error); break; }

#endif /* VALIDATE_ORDER_BREAK */
  }

  rti_tstop(rti_tid_mpi_mergek);

  return stage - 1;
}
