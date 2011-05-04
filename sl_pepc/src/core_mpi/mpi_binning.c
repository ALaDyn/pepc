/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_binning.c
 *  timestamp: 2011-02-24 21:37:35 +0100
 *  
 */


/* sl_macro MB_REDUCEBCAST_THRESHOLD */
/* sl_macro MB_TRACE_IF */


#include "sl_common.h"


#define REDUCEBCAST_ROOT  0


#if !defined(MB_REDUCEBCAST_THRESHOLD) && defined(GLOBAL_REDUCEBCAST_THRESHOLD)
# define MB_REDUCEBCAST_THRESHOLD  GLOBAL_REDUCEBCAST_THRESHOLD
#endif

#ifndef MB_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MB_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MB_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_binning_create(global_bins_t *gb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t doweights, binning_t *bm, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_create */
{
  gb->bm = bm;

  binning_create(&gb->lb, max_nbins, max_nbinnings, s, nelements, doweights, bm);

#ifdef elem_weight
# define MY_WEIGHT_FACTOR  gb->lb.weight_factor
#else
# define MY_WEIGHT_FACTOR  1
#endif

  gb->cws = z_alloc(gb->lb.max_nbinnings * MY_WEIGHT_FACTOR * gb->bm->max_nbins, sizeof(slweight_t));

  gb->bcws = gb->lb.bcws;

  return 0;
}


slint_t mpi_binning_destroy(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_destroy */
{
  binning_destroy(&gb->lb);

  z_free(gb->cws);

  return 0;
}


slint_t mpi_binning_pre(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_pre */
{
  binning_pre(&gb->lb);

  return 0;
}


slint_t mpi_binning_exec_reset(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_exec_reset */
{
  binning_exec_reset(&gb->lb);

  return 0;
}


slint_t mpi_binning_exec_local(global_bins_t *gb, slint_t b, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_exec_local */
{
  return binning_exec(&gb->lb, b);
}


slint_t mpi_binning_exec_global(global_bins_t *gb, slint_t root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_exec_global */
{
  slint_t b, k;


#ifdef elem_weight  
# define MY_WEIGHT_FACTOR  gb->lb.weight_factor
#else
# define MY_WEIGHT_FACTOR  1
#endif

  Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " doubles", ((root < 0)?"all-":""), (gb->lb.nbinnings * MY_WEIGHT_FACTOR * gb->bm->nbins));

  if (root < 0)
  {
#ifdef MB_REDUCEBCAST_THRESHOLD
    if (size >= MB_REDUCEBCAST_THRESHOLD)
    {
      Z_TRACE_IF(MB_TRACE_IF, "allreduce = reduce + bcast");

      MPI_Reduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * MY_WEIGHT_FACTOR * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, REDUCEBCAST_ROOT, comm);
      MPI_Bcast(gb->cws, gb->lb.nbinnings * MY_WEIGHT_FACTOR * gb->bm->nbins, weight_mpi_datatype, REDUCEBCAST_ROOT, comm);

    } else
#endif
      MPI_Allreduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * MY_WEIGHT_FACTOR * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, comm);

  } else MPI_Reduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * MY_WEIGHT_FACTOR * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, root, comm);

  if (root < 0 || root == rank)
  {
    for (b = 0; b < gb->lb.nbinnings; ++b)
    {
      Z_TRACE_ARRAY_IF(MB_TRACE_IF, "%" slint_fmt ": counts =", " %" slweight_fmt, k, gb->bm->nbins, gb_counts(gb, b, 0), b);
#ifdef elem_weight
      if (gb->lb.doweights)
        Z_TRACE_ARRAY_IF(MB_TRACE_IF, "%" slint_fmt ": weights =", " %" slweight_fmt, k, gb->bm->nbins, gb_weights(gb, b, 0), b);
#endif
    }
  }

  return 0;
}


slint_t mpi_binning_refine(global_bins_t *gb, slint_t b, slint_t k, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_refine */
{
  return binning_refine(&gb->lb, b, k, sp, s);
}


slint_t mpi_binning_hit(global_bins_t *gb, slint_t b, slint_t k, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_hit */
{
  binning_hit(&gb->lb, b, k, sp, s);
  
  return 0;
}


slint_t mpi_binning_finalize(global_bins_t *gb, slint_t b, slweight_t dcw, slint_t lc_min, slint_t lc_max, slweight_t *lcw, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_finalize */
{
  binning_finalize(&gb->lb, b, dcw, lc_min, lc_max, lcw, sp, s);
  
  return 0;
}


slint_t mpi_binning_post(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_post */
{
  binning_post(&gb->lb);

  return 0;
}
