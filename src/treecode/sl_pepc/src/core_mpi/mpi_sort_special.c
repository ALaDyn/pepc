/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_sort_special.c
 *  
 */


/* sl_macro MSS_TRACE_IF */
/* sl_macro MSS_VERIFY */

#include "sl_common.h"

/*#define MSS_VERIFY*/

double mss_i_t[3];  /* sl_global, sl_var mss_i_t */
double mss_p_t[3];  /* sl_global, sl_var mss_p_t */
double mss_b_t[3];  /* sl_global, sl_var mss_b_t */

slint_t mss_sync = 0;  /* sl_global, sl_var mss_sync */
slint_t mss_i_sync = 0;  /* sl_global, sl_var mss_i_sync */
slint_t mss_p_sync = 0;  /* sl_global, sl_var mss_p_sync */
slint_t mss_b_sync = 0;  /* sl_global, sl_var mss_b_sync */


#ifndef MSS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


static inline slint_t key2proc(slpkey_t k, slint_t nmm, slpkey_t *mm)
{
  slint_t l = 0;
  slint_t h = nmm - 1;
  slint_t m = (l + h) / 2;

  while (l < h)
  {
    if (sl_key_pure_cmp_le(k, mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }

  return m;
}


#ifdef key_integer

slint_t mpi_sort_insert_radix(elements_t *s0, elements_t *s1, elements_t *xs, slpkey_t *mmkeys, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_insert_radix */
{
  slint_t i, j;

  slpkey_t all_mmkeys[2 * size];
  
  int scounts[size], sdispls[size], rcounts[size], rdispls[size];


  if (mss_sync || mss_i_sync) MPI_Barrier(comm);

  /* local rearrange */
  mss_i_t[0] = z_time_get_s();

  MPI_Allgather(mmkeys, 2, pkey_mpi_datatype, all_mmkeys, 2, pkey_mpi_datatype, comm);
  
  for (i = 0; i < size; ++i) scounts[i] = 0;

  for (i = 0; i < s0->size; ++i)
  {
    j = key2proc(s0->keys[i], size, all_mmkeys);
    ++scounts[j];
  }

  sdispls[0] = 0;
  for (i = 1; i < size; ++i) sdispls[i] = sdispls[i - 1] + scounts[i - 1];
  
  for (i = 0; i < s0->size; ++i)
  {
    j = key2proc(s0->keys[i], size, all_mmkeys);
    elem_copy_at(s0, i, s1, sdispls[j]);
    ++sdispls[j];
  }

  if (mss_sync || mss_i_sync) MPI_Barrier(comm);
  mss_i_t[0] = z_time_get_s() - mss_i_t[0];


  /* all-to-all */
  mss_i_t[1] = z_time_get_s();

  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

  sdispls[0] = rdispls[0] = 0;
  for (i = 1; i < size; ++i)
  {
    sdispls[i] = sdispls[i - 1] + scounts[i - 1];
    rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  }

  s0->size = s1->size = rdispls[size - 1] + rcounts[size - 1];

  mpi_elements_alltoallv_db(s1, scounts, sdispls, s0, rcounts, rdispls, size, rank, comm);

  if (mss_sync || mss_i_sync) MPI_Barrier(comm);
  mss_i_t[1] = z_time_get_s() - mss_i_t[1];


  /* local sort */
  mss_i_t[2] = z_time_get_s();

  sort_radix(s0, s1, rhigh, rlow, -1);

  if (mss_sync || mss_i_sync) MPI_Barrier(comm);
  mss_i_t[2] = z_time_get_s() - mss_i_t[2];

  return 0;
}


slint_t mpi_sort_presorted_radix(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_presorted_radix */
{
  merge2x_f m2x;


  if (mss_sync || mss_p_sync) MPI_Barrier(comm);

  /* local sort */
  mss_p_t[0] = z_time_get_s();

  sort_radix(s0, xs, rhigh, rlow, -1);

  if (mss_sync || mss_p_sync) MPI_Barrier(comm);
  mss_p_t[0] = z_time_get_s() - mss_p_t[0];


  /* merge sorted */
  mss_p_t[1] = z_time_get_s();

  switch (merge_type)
  {
    case 0:
      m2x = merge2_basic_straight_01_x;
      break;
    case 1:
      m2x = merge2_compo_tridgell;
      break;
    case 2:
      m2x = merge2_compo_hula;
      break;
    default:
      m2x = merge2_memory_adaptive;
      break;
  }

  mpi_mergek_sorted(s0, m2x, xs, size, rank, comm);

  if (mss_sync || mss_p_sync) MPI_Barrier(comm);
  mss_p_t[1] = z_time_get_s() - mss_p_t[1];

  return 0;
}

#endif


typedef struct _tproc_sort_back_data
{
  slint_t nmm;
  slpkey_t *mm;

} tproc_sort_back_data;


static inline int tproc_sort_back(elements_t *s, slint_t x, void *data)
{
  tproc_sort_back_data *d = data;

  int l = 0;
  int h = d->nmm - 1;
  int m = (l + h) / 2;

  while (l < h)
  {
    if (sl_key_pure_cmp_le(s->keys[x], d->mm[2 * m + 1])) h = m;
    else l = m + 1;

    m = (l + h) / 2;
  }

  return m;
}


slint_t mpi_sort_back(elements_t *s0, elements_t *s1, elements_t *xs, slpkey_t *lh, slint_t ntotal, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_back */
{
  slint_t i, j;

  slpkey_t lhs[2 * size];
  
  tproc_sort_back_data data;


  if (mss_sync || mss_b_sync) MPI_Barrier(comm);

  mss_b_t[0] = z_time_get_s();

  MPI_Allgather(lh, 2, pkey_mpi_datatype, lhs, 2, pkey_mpi_datatype, comm);

  mss_b_t[0] = z_time_get_s() - mss_b_t[0];
  

  /* all-to-all specific */
  mss_b_t[1] = z_time_get_s();

  data.nmm = size;
  data.mm = lhs;

  mpi_elements_alltoall_specific(s0, NULL, s1, tproc_sort_back, &data, size, rank, comm);

  mss_b_t[1] = z_time_get_s() - mss_b_t[1];


  /* local permute back */
  mss_b_t[2] = z_time_get_s();

  for (i = 0; i < s0->size; ++i)
  {
    j = s0->keys[i] - lh[0];
    elem_copy_at(s0, i, s1, j);
  }
  elem_ncopy(s1, s0, s0->size);

  if (mss_sync || mss_b_sync) MPI_Barrier(comm);
  mss_b_t[2] = z_time_get_s() - mss_b_t[2];

  return 0;
}
