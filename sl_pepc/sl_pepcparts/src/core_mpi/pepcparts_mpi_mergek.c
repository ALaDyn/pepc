/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_mergek.c
 *  timestamp: 2009-11-26 17:36:08 +0100
 *  
 */


#include "sl_common.h"


#ifdef SLDEBUG
# define CHECK_ORDER
#endif


slint_t pepcparts_mpi_mergek_equal(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_mergek_equal */
{
  slint_t stage = 0, other_rank, up, high_rank;

#ifdef CHECK_ORDER
  slint_t local_order, global_order;
#endif

  rti_tstart(rti_tid_mpi_mergek);

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

  if (size > 1)
  while ((other_rank = (sn)(size, rank, stage, snd, &up)) >= 0)
  {
#ifdef CHECK_ORDER
    local_order = pepcparts_elements_validate_order(s, 1);
    MPI_Allreduce(&local_order, &global_order, 1, sl_int_type_mpi, MPI_SUM, comm);
    if (global_order)
    {
      fprintf(stderr, "%d: pepcparts_mpi_mergek_equal: warning: input order failed (at %" sl_int_type_fmt ")\n", rank, local_order);
# ifdef CHECK_ORDER_BREAK
      break;
# endif
    }
#endif

    rti_tstart(rti_tid_mpi_mergek_merge2);
    high_rank = (up)?(xmin(rank, other_rank)):(xmax(rank, other_rank));
    pepcparts_mpi_merge2(s, other_rank, high_rank, NULL, m2x, xs, size, rank, comm);
    rti_tstop(rti_tid_mpi_mergek_merge2);

#ifdef CHECK_ORDER
    local_order = pepcparts_elements_validate_order(s, 1);
    MPI_Allreduce(&local_order, &global_order, 1, sl_int_type_mpi, MPI_SUM, comm);
    if (global_order)
    {
      fprintf(stderr, "%d: pepcparts_mpi_mergek_equal: warning: output order failed (at %" sl_int_type_fmt ")\n", rank, local_order);
# ifdef CHECK_ORDER_BREAK
      break;
# endif
    }
#endif

    ++stage;
  }

  rti_tstop(rti_tid_mpi_mergek);

  return stage - 1;
}


slint_t pepcparts_mpi_mergek(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_mergek */
{
  int i;
  slint_t *vsizes, *vranks, vsize;
  slint_t size_total, size_max;

  vsizes = sl_alloc(2 * size, sizeof(slint_t));
  vranks = vsizes + size;


  MPI_Allgather(&s->size, 1, sl_int_type_mpi, vsizes, 1, sl_int_type_mpi, comm);

  size_total = size_max = vsizes[0];
  for (i = 1; i < size; ++i)
  {
    size_total += vsizes[i];
    size_max = xmin(size_max, vsizes[i]);
  }

  vsize = 0;
  vranks[0] = 0;
  for (i = 0; i < size; ++i)
  {
    vsizes[i] = (vsizes[i] / size_max) + (vsizes[i] % size_max > 0);
    vsize += vsizes[i];
    if (i > 0) vranks[i] = vranks[i - 1] + vsizes[i - 1];
  }

/*  printf("%d here: size_total = %" sl_int_type_fmt ", size_max = %" sl_int_type_fmt "\n", rank, size_total, size_max);

  printf("%d here: vranks:vsizes = [ ", rank);
  for (i = 0; i < size; ++i) printf(" %" sl_int_type_fmt ":%" sl_int_type_fmt " ", vranks[i], vsizes[i]);
  printf("]\n");*/


  free(vsizes);

  return 0;
}
