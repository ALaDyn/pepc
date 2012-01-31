/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_sort_merge.c
 *  
 */


/* sl_macro MSM_TRACE_IF */
/* sl_macro MSM_VERIFY */
/* sl_macro MSM_VERIFY_OUTPUT */

#include "sl_common.h"

/*#define MSM_VERIFY
#define MSM_VERIFY_OUTPUT*/

double msm_t[4];  /* sl_global, sl_var msm_t */

slint_t msm_sync = 0;  /* sl_global, sl_var msm_sync */


#ifndef MSM_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSM_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSM_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


#ifdef SLDEBUG
# define CHECK_ORDER
#endif


slint_t mpi_sort_merge(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_merge */
{
  merge2x_f m2x = merge2_memory_adaptive;
  sortnet_f sn = sn_batcher;

#ifdef CHECK_ORDER
  slint_t stages;
  slint_t rorders[2];
#endif

#ifdef key_integer
  sort_radix(s0, xs, -1, -1, -1);
#else
  sort_quick(s0, xs);
#endif

#ifdef CHECK_ORDER
  stages =
#endif
    mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

#ifdef CHECK_ORDER
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

  if (rank == 0) printf("%d: %s (%" slint_fmt ", %" slint_fmt ") with %" slint_fmt " stages\n", rank, (rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1], stages);
#endif

  return 0;
}


slint_t mpi_sort_merge2(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t sort_type, double *times, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_merge2 */
{
  double _times[2];

  merge2x_f m2x;
  sortnet_f sn;

#ifdef CHECK_ORDER
  slint_t stages;
  slint_t rorders[2];
#endif


  if (times == NULL) times = _times;

  MPI_Barrier(comm);

  times[0] = z_time_get_s();
  sort_radix(s0, xs, -1, -1, -1);

  MPI_Barrier(comm);
  times[0] = z_time_get_s() - times[0];

  switch (merge_type)
  {
    case 0:
      m2x = merge2_basic_auto_01_x;
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

  switch (sort_type)
  {
    case 0:
      sn = sn_batcher;
      break;
    case 1:
      sn = sn_bitonic;
      break;
    case 2:
      sn = sn_odd_even_trans;
      break;
  }

  times[1] = z_time_get_s();

#ifdef CHECK_ORDER
  stages =
#endif
    mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

  MPI_Barrier(comm);
  times[1] = z_time_get_s() - times[1];

#ifdef CHECK_ORDER
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

  if (rank == 0) printf("%d: %s (%" slint_fmt ", %" slint_fmt ") with %" slint_fmt " stages\n", rank, (rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1], stages);
#endif

  return 0;
}


#ifdef key_integer

slint_t mpi_sort_merge_radix(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t sort_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_sort_merge_radix */
{
  merge2x_f m2x;
  sortnet_f sn;

#ifdef MSM_VERIFY
  slint_t rorders[2];
#endif


  if (msm_sync) MPI_Barrier(comm);
  msm_t[0] = z_time_get_s();
  sort_radix(s0, xs, -1, -1, -1);
  msm_t[0] = z_time_get_s() - msm_t[0];

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

  switch (sort_type)
  {
    case 0:
      sn = sn_batcher;
      break;
    case 1:
      sn = sn_bitonic;
      break;
    case 2:
      sn = sn_odd_even_trans;
      break;
  }

  if (msm_sync) MPI_Barrier(comm);
  msm_t[1] = z_time_get_s();

  mpi_mergek(s0, sn, NULL, m2x, xs, size, rank, comm);

  if (msm_sync) MPI_Barrier(comm);
  msm_t[1] = z_time_get_s() - msm_t[1];

#ifdef MSM_VERIFY
  mpi_elements_check_order(s0, 1, rorders, size, rank, comm);

# ifndef MSM_VERIFY_OUTPUT
  if (!rorders[0] || !rorders[1])
# endif
    Z_NOTICE_IF((rank == 0), "%s (%" slint_fmt ", %" slint_fmt ")", (rorders[0] && rorders[1])?"SUCCESS":"FAILED", rorders[0], rorders[1]);
#endif

  return 0;
}

#endif


#undef CHECK_ORDER
