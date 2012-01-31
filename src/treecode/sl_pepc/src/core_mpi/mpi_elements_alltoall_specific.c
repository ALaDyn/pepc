/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_elements_alltoall_specific.c
 *  
 */


/* sl_macro MEAS_TRACE_IF */

#include "sl_common.h"


double meas_t[2];  /* sl_global, sl_var meas_t */


#ifndef MEAS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MEAS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MEAS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_elements_alltoall_specific(elements_t *s0, elements_t *s1, elements_t *xs, tproc_f tproc, void *data, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoall_specific */
{
  slint_t i, j;
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];
  tproc_data_t *tproc_data = data;


  /* local rearrange */
  meas_t[0] = z_time_get_s();

  for (i = 0; i < size; ++i) scounts[i] = 0;

#define TPROC_LOOP \
  for (i = 0; i < s0->size; ++i) \
  { \
    j = TPROC; \
    ++scounts[j]; \
  }

  if (tproc == SL_TPROC_INT)
  {
#define TPROC  tproc_data->int_tprocs[i]
    TPROC_LOOP
#undef TPROC

  } else if (tproc == SL_TPROC_INT_MASK)
  {
#define TPROC  (tproc_data->int_tprocs[i] & tproc_data->int_mask)
    TPROC_LOOP
#undef TPROC

  } else
  {
#define TPROC  tproc(s0, i, data)
    TPROC_LOOP
#undef TPROC
  }

#undef TPROC_LOOP

  scounts[rank] = 0;

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];

  sdispls[rank] = 0;

#define TPROC_LOOP \
  for (i = 0; i < s0->size; ++i) \
  { \
    j = TPROC; \
    if (j == rank) elem_copy_at(s0, i, s0, sdispls[j]); \
    else elem_copy_at(s0, i, xs, sdispls[j]); \
    ++sdispls[j]; \
  }

  if (tproc == SL_TPROC_INT)
  {
#define TPROC  tproc_data->int_tprocs[i]
    TPROC_LOOP
#undef TPROC

  } else if (tproc == SL_TPROC_INT_MASK)
  {
#define TPROC  (tproc_data->int_tprocs[i] & tproc_data->int_mask)
    TPROC_LOOP
#undef TPROC

  } else
  {
#define TPROC  tproc(s0, i, data)
    TPROC_LOOP
#undef TPROC
  }

#undef TPROC_LOOP

  meas_t[0] = z_time_get_s() - meas_t[0];

  /* all-to-all */
  meas_t[1] = z_time_get_s();

  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

  rdispls[0] = sdispls[rank];
  sdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  mpi_elements_alltoallv_db(xs, scounts, sdispls, s0, rcounts, rdispls, size, rank, comm);
  
  s0->size = rdispls[size - 1] + rcounts[size - 1];

  meas_t[1] = z_time_get_s() - meas_t[1];

  return 0;
}
