/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_find_exact.c
 *  timestamp: 2011-02-08 21:32:54 +0100
 *  
 */

/* - basierend auf "bitonic merge-exchange", gefunden in [taxonomy] mit Verweis auf [Alekseyev 1969] und [Knuth 1973]
   - binäre Suche durch [Tridgell,Brent] Übergang zur Verwendung von "bisection" beschrieben in [Zhou,Tridgell]
*/

#include "sl_common.h"


#define DEBUG_OR_NOT  (sl_mpi_rank == -1)


/*#define CHECK_EQUAL*/

slint_t mpi_find_exact_equal(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *ex_start, slint_t *ex_size, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_find_exact_equal */
{
  const int tag = 1;

  elements_t interval;

  slint_t low_cmp_i, high_cmp_i, *my_cmp_i, *other_cmp_i;
  slkey_pure_t low_cmp_key, high_cmp_key, *my_cmp_key, *other_cmp_key;
  
  slint_t local_ex_start, local_ex_size;

  MPI_Status status;


  if (ex_start == NULL) ex_start = &local_ex_start;
  if (ex_size == NULL) ex_size = &local_ex_size;

  if (rank != high_rank)
  {
    my_cmp_i = &low_cmp_i;
    other_cmp_i = &high_cmp_i;

    my_cmp_key = &low_cmp_key;
    other_cmp_key = &high_cmp_key;

  } else
  {
    my_cmp_i = &high_cmp_i;
    other_cmp_i = &low_cmp_i;

    my_cmp_key = &high_cmp_key;
    other_cmp_key = &low_cmp_key;
  }

#ifdef CHECK_EQUAL
  *my_cmp_i = s->size;
  MPI_Sendrecv(my_cmp_i, 1, int_mpi_datatype, other_rank, tag, other_cmp_i, 1, int_mpi_datatype, other_rank, tag, comm, &status);
  if (*other_cmp_i != s->size)
  {
    fprintf(stderr, "mpi_find_exact_equal: error: element lists not equal (%" sl_int_type_fmt " vs. %" sl_int_type_fmt ")\n", *my_cmp_i, *other_cmp_i);
    return -1;
  }
#endif

  elem_assign(s, &interval);

  while (interval.size > 0)
  {
    low_cmp_i = interval.size / 2;
    high_cmp_i = interval.size - low_cmp_i - 1;
    
    *my_cmp_key = key_purify(*elem_key_at(&interval, *my_cmp_i));

    MPI_Sendrecv(my_cmp_key, key_pure_size_mpi, key_pure_type_mpi, other_rank, tag, other_cmp_key, key_pure_size_mpi, key_pure_type_mpi, other_rank, tag, comm, &status);

    if (key_pure_cmp_le(low_cmp_key, high_cmp_key))
    {
      if (rank != high_rank) elem_add(&interval, low_cmp_i + 1);
      
      interval.size = high_cmp_i;

    } else
    {
      if (rank == high_rank) elem_add(&interval, high_cmp_i + 1);
      
      interval.size = low_cmp_i;
    }
  }

  if (rank != high_rank)
  {
    *ex_start = interval.keys - s->keys;
    *ex_size = s->size - *ex_start;

  } else
  {
    *ex_start = 0;
    *ex_size = interval.keys - s->keys;
  }

  return *ex_size;
}


slint_t mpi_find_exact(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, slint_t *ex_start, slint_t *ex_sizes, slint_t *nx_move, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_find_exact */
{
  const int tag = 1;

  slint_t low_sizes[2], high_sizes[2], *my_sizes, *other_sizes;
  slint_t low_diff, high_diff;
  slint_t low_int[2], high_int[2], int_size;
  slint_t local_ex_start, local_ex_sizes[2], local_nx_move;
  
  elements_t interval;

  MPI_Status status;


  if (ex_start == NULL) ex_start = &local_ex_start;
  if (ex_sizes == NULL) ex_sizes = local_ex_sizes;
  if (nx_move == NULL) nx_move = &local_nx_move;

  if (rank != high_rank)
  {
    my_sizes = low_sizes;
    other_sizes = high_sizes;

  } else
  {
    my_sizes = high_sizes;
    other_sizes = low_sizes;
  }
  
  my_sizes[0] = s->size;
  if (dst_size) my_sizes[1] = *dst_size; else my_sizes[1] = s->size;

  MPI_Sendrecv(my_sizes, 2, int_mpi_datatype, other_rank, tag, other_sizes, 2, int_mpi_datatype, other_rank, tag, comm, &status);

  low_diff = low_sizes[1] - low_sizes[0];
  high_diff = high_sizes[1] - high_sizes[0];

  /* correct diffs */  
  if (low_diff < -low_sizes[0])
  {
    if (high_diff < -high_sizes[0])
    {
      low_diff = 0;
      high_diff = 0;

    } else
    {
      low_diff = -z_min(high_diff, low_sizes[0]);
      high_diff = -low_diff;
    }

  } else
  {
    if (high_diff < -high_sizes[0] || low_diff + high_diff >= 0)
    {
      high_diff = -z_min(low_diff, high_sizes[0]);
      low_diff = -high_diff;

    } else if (low_diff + high_diff < 0)
    {
      fprintf(stderr, "mpi_find_exact2: error: destination sizes too small (%" sl_int_type_fmt " + %" sl_int_type_fmt " < 0!)\n", low_diff, high_diff);
      return -1;
    }
  }

  low_sizes[1] = low_diff + low_sizes[0];
  high_sizes[1] = high_diff + high_sizes[0];

  /* FIXME: APP_ASSERT(low_diff == -high_diff) */
  
  low_int[0] = 0;
  low_int[1] = low_sizes[0] - 1;
  
  high_int[0] = 0;
  high_int[1] = high_sizes[0] - 1;
  
  if (low_diff > 0) high_int[0] += low_diff;
  if (high_diff > 0) low_int[1] -= high_diff;
  
  int_size = z_min(low_int[1] - low_int[0] + 1, high_int[1] - high_int[0] + 1);
  
  low_int[0] = low_int[1] - int_size + 1;
  high_int[1] = high_int[0] + int_size - 1;
  
  if (rank != high_rank) elem_assign_at(s, low_int[0], &interval);
  else elem_assign_at(s, high_int[0], &interval);
  
  interval.size = int_size;
  
  mpi_find_exact_equal(&interval, other_rank, high_rank, ex_start, ex_sizes, size, rank, comm);

  low_int[0] = low_int[1] - *ex_sizes + 1;
  high_int[1] = high_int[0] + *ex_sizes - 1;

  if (low_diff > 0) high_int[0] -= low_diff;
  if (high_diff > 0) low_int[1] += high_diff;

  if (rank != high_rank)
  {
    if (dst_size) *dst_size = low_sizes[1];

    *ex_start = low_int[0];

    ex_sizes[0] = low_int[1] - low_int[0] + 1;
    ex_sizes[1] = high_int[1] - high_int[0] + 1;

    *nx_move = 0;

  } else
  {
    if (dst_size) *dst_size = high_sizes[1];

    *ex_start = high_int[0];

    ex_sizes[0] = high_int[1] - high_int[0] + 1;
    ex_sizes[1] = low_int[1] - low_int[0] + 1;

    *nx_move = high_diff;
  }
  
  return 0;
}
