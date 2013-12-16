/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_binning.c
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


slint_t mpi_binning_create(global_bins_t *gb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t docounts, slint_t doweights, binning_t *bm, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_create */
{
  gb->bm = bm;

  binning_create(&gb->lb, max_nbins, max_nbinnings, s, nelements, docounts, doweights, bm);

  gb->bcws = gb->lb.bcws;

#if defined(elem_weight) && defined(sl_weight_intequiv)
  gb->cws = z_alloc(gb->lb.max_nbinnings * gb->lb.cw_factor * gb->bm->max_nbins, sizeof(slweight_t));
#else
  gb->cs = z_alloc(gb->lb.max_nbinnings * 1 * gb->bm->max_nbins, sizeof(slint_t));
# ifdef elem_weight
  if (gb->lb.doweights)
    gb->ws = z_alloc(gb->lb.max_nbinnings * 1 * gb->bm->max_nbins, sizeof(slweight_t));
  else gb->ws = NULL;
# endif
#endif

  return 0;
}


slint_t mpi_binning_destroy(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_destroy */
{
  binning_destroy(&gb->lb);

#if defined(elem_weight) && defined(sl_weight_intequiv)
  z_free(gb->cws);
#else
  z_free(gb->cs);
# ifdef elem_weight
  if (gb->lb.doweights)
    z_free(gb->ws);
# endif
#endif

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


#if defined(elem_weight) && defined(sl_weight_intequiv)
  Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " weights", ((root < 0)?"all-":""), (gb->lb.nbinnings * gb->lb.cw_factor * gb->bm->nbins));
#else
  Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " ints", ((root < 0)?"all-":""), (gb->lb.nbinnings * 1 * gb->bm->nbins));
# ifdef elem_weight
  if (gb->lb.doweights)
    Z_TRACE_IF(MB_TRACE_IF, "%sreducing %" slint_fmt " weights", ((root < 0)?"all-":""), (gb->lb.nbinnings * 1 * gb->bm->nbins));
# endif
#endif

  if (root < 0)
  {
#if defined(elem_weight) && defined(sl_weight_intequiv)
    sl_MPI_Allreduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * gb->lb.cw_factor * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, comm, size, rank);
#else
    sl_MPI_Allreduce(gb->lb.cs, gb->cs, gb->lb.nbinnings * 1 * gb->bm->nbins, int_mpi_datatype, MPI_SUM, comm, size, rank);
# ifdef elem_weight
    if (gb->lb.doweights)
      sl_MPI_Allreduce(gb->lb.ws, gb->ws, gb->lb.nbinnings * 1 * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, comm, size, rank);
# endif
#endif

  } else
  {
#if defined(elem_weight) && defined(sl_weight_intequiv)
    MPI_Reduce(gb->lb.cws, gb->cws, gb->lb.nbinnings * gb->lb.cw_factor * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, root, comm);
#else
    MPI_Reduce(gb->lb.cs, gb->cs, gb->lb.nbinnings * 1 * gb->bm->nbins, int_mpi_datatype, MPI_SUM, root, comm);
# ifdef elem_weight
    if (gb->lb.doweights)
      MPI_Reduce(gb->lb.ws, gb->ws, gb->lb.nbinnings * 1 * gb->bm->nbins, weight_mpi_datatype, MPI_SUM, root, comm);
# endif
#endif
  }

  if (root < 0 || root == rank)
  {
    for (b = 0; b < gb->lb.nbinnings; ++b)
    {
      if (gb->lb.docounts)
        Z_TRACE_ARRAY_IF(MB_TRACE_IF, k, gb->bm->nbins, " %" slcount_fmt, gb_counts(gb, b, 0)[k], "%" slint_fmt ": counts =", b);
#ifdef elem_weight
      if (gb->lb.doweights)
        Z_TRACE_ARRAY_IF(MB_TRACE_IF, k, gb->bm->nbins, " %" slweight_fmt, gb_weights(gb, b, 0)[k], "%" slint_fmt ": weights =", b);
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


slint_t mpi_binning_finalize(global_bins_t *gb, slint_t b, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_finalize */
{
  binning_finalize(&gb->lb, b, dc, dw, lc_min, lc_max, lcs, lws, sp, s);
  
  return 0;
}


slint_t mpi_binning_post(global_bins_t *gb, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_binning_post */
{
  binning_post(&gb->lb);

  return 0;
}


#undef REDUCEBCAST_ROOT
