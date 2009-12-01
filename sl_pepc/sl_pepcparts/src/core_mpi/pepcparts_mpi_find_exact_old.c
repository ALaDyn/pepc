/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepcparts_mpi_find_exact_old.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* - basierend auf "bitonic merge-exchange", gefunden in [taxonomy] mit Verweis auf [Alekseyev 1969] und [Knuth 1973]
   - binäre Suche durch [Tridgell,Brent] Übergang zur Verwendung von "bisection" beschrieben in [Zhou,Tridgell]
*/

#include "sl_common.h"


slint_t pepcparts_mpi_find_exact_old(elements_t *s, slint_t counterpart, slint_t high, elements_t *xs, slint_t *start, slint_t *end, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_find_exact_old */
{
  elements_t x;
  
  elements_t intervall;
  
  slint_t n0, n1, n;
  slint_t tstart, tend;

  slint_t high_i, low_i;
  slkey_t high_key, low_key, *foreign_key, *my_key;
  
  MPI_Status status;
  

  if (xs == NULL)
  {
    xs = &x;
    pepcparts_elements_alloc(xs, 1, 1, 0);
  }

  if (rank != high) n0 = n = s->size;
  else n1 = n = s->size;
  
  MPI_Sendrecv_replace(&n, 1, sl_int_type_mpi, counterpart, 1, counterpart, 1, comm, &status);
  
  if (rank != high) n1 = n;
  else n0 = n;
  
  n = xmin(n0, n1);

/*  printf("%d here: n0 %d - n1 %d - min %d\n", rank, n0, n1, n);*/

  if (rank != high) elem_assign_at(s, s->size - n, &intervall);
  else elem_assign(s, &intervall);
  
  intervall.size = n;
    
  /* doit */

  if (rank != high)
  {
    my_key = &low_key;
    foreign_key = &high_key;
  } else
  {
    my_key = &high_key;
    foreign_key = &low_key;
  }

  while (intervall.size > 0)
  {
    high_i = intervall.size / 2;
    low_i = intervall.size - high_i - 1;
    
    if (rank != high)
    {
      key_copy_at(intervall.keys, low_i, my_key, 0);
    } else
    {
      key_copy_at(intervall.keys, high_i, my_key, 0);
    }

    MPI_Sendrecv(my_key, key_size_mpi, key_type_mpi, counterpart, 1, foreign_key, key_size_mpi, key_type_mpi, counterpart, 1, comm, &status);

/*    printf("%d here: intervall", rank); pepcparts_elements_printf(&intervall);
    printf("%d here: low_i %d - high_i %d\n", rank, low_i, high_i);
    printf("%d here: low_key %d - high_key %d\n", rank, low_key, high_key);
    printf("%d here: my_key %d - foreign_key %d\n", rank, *my_key, *foreign_key);*/
    
    if (key_cmp_gt(low_key, high_key))
    {
      if (rank == high) elem_add(&intervall, high_i + 1);

      intervall.size = low_i;

    } else
    {
      if (rank != high) elem_add(&intervall, low_i + 1);

      intervall.size = high_i;
    }
  }

/*  printf("%d here: finished with ", rank); pepcparts_elements_printf(&intervall);*/
  
  if (start == NULL) start = &tstart;
  if (end == NULL) end = &tend;

  if (rank != high)
  {
    *start = intervall.keys - s->keys;
    *end = s->size - 1;

  } else
  {
    *start = 0;
    *end = intervall.keys - s->keys - 1;
  }

  if (xs == &x) pepcparts_elements_free(xs);

  return *end - *start + 1;
}
