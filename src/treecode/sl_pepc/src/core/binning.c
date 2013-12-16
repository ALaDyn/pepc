/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/binning.c
 *  
 */


/* sl_macro B_TRACE_IF */


#include "sl_common.h"


#ifndef B_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define B_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define B_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t binning_create(local_bins_t *lb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t docounts, slint_t doweights, binning_t *bm) /* sl_proto, sl_func binning_create */
{
  slint_t j;

  bm->docounts = docounts;
#ifdef elem_weight
  bm->doweights = doweights;
#endif
  lb->bm = bm;

  lb->nbins = 0;
  lb->max_nbins = max_nbins;

  lb->nbinnings = 0;
  lb->max_nbinnings = (max_nbinnings > 0)?z_min(max_nbinnings,max_nbins):max_nbins;

  lb->nelements = nelements;

  lb->docounts = docounts;
#ifdef elem_weight
  lb->doweights = doweights;
#endif

  lb->bins0 = z_alloc(lb->max_nbins * lb->nelements, sizeof(bin_t));
  lb->bins1 = z_alloc(lb->max_nbins * lb->nelements, sizeof(bin_t));

  lb->bcws = z_alloc(lb->max_nbins, sizeof(slint_t));

#if defined(elem_weight) && defined(sl_weight_intequiv)
  lb->cw_factor = (lb->docounts != 0) + (lb->doweights != 0);
  lb->w_index = (lb->docounts != 0);
  lb->cws = z_alloc(lb->max_nbinnings * lb->cw_factor * lb->bm->max_nbins, sizeof(slweight_t));
  lb->bin_cw_factor = 1 + (lb->doweights != 0);
  lb->bin_cws = z_alloc(lb->max_nbinnings * lb->bin_cw_factor * lb->nelements * bm->max_nbins, sizeof(slweight_t));
#else
  lb->cs = z_alloc(lb->max_nbinnings * 1 * lb->bm->max_nbins, sizeof(slint_t));
  lb->bin_cs = z_alloc(lb->max_nbinnings * 1 * lb->nelements * bm->max_nbins, sizeof(slint_t));
# ifdef elem_weight
  if (lb->doweights)
  {
    lb->ws = z_alloc(lb->max_nbinnings * 1 * lb->bm->max_nbins, sizeof(slweight_t));
    lb->bin_ws = z_alloc(lb->max_nbinnings * 1 * lb->nelements * bm->max_nbins, sizeof(slweight_t));

  } else lb->ws = lb->bin_ws = NULL;
# endif
#endif

  lb->last_exec_b = -1;

  lb->bins = lb->bins0;
  lb->bins_new = lb->bins1;
  
  if (lb->max_nbins > 0)
  {
    lb->nbins = 1;
    for (j = 0; j < lb->nelements; ++j) elem_assign(&s[j], &lb->bins[0 * lb->nelements + j].s);
  }

  return 0;
}


slint_t binning_destroy(local_bins_t *lb) /* sl_proto, sl_func binning_destroy */
{
  z_free(lb->bins0);
  z_free(lb->bins1);
  
  z_free(lb->bcws);

#if defined(elem_weight) && defined(sl_weight_intequiv)
  z_free(lb->cws);
  z_free(lb->bin_cws);
#else
  z_free(lb->cs);
  z_free(lb->bin_cs);
# ifdef elem_weight
  if (lb->doweights)
  {
    z_free(lb->ws);
    z_free(lb->bin_ws);
  }
# endif
#endif

  return 0;
}


slint_t binning_pre(local_bins_t *lb) /* sl_proto, sl_func binning_pre */
{
  lb->bm->pre(lb->bm);

  lb->nbins_new = 0;
  lb->last_new_b = lb->last_new_k = -1;

  return 0;
}


slint_t binning_exec_reset(local_bins_t *lb) /* sl_proto, sl_func binning_exec_reset */
{
  lb->nbinnings = 0;

  return 0;
}


slint_t binning_exec(local_bins_t *lb, slint_t b) /* sl_proto, sl_func binning_exec */
{
  slint_t j;
  slkey_pure_t k;

  slcount_t *counts, *bin_counts;
#ifdef elem_weight
  slweight_t *weights, *bin_weights;
#endif


  Z_TRACE_IF(B_TRACE_IF, "b: %" slint_fmt ", last_exec_b: %" slint_fmt ", nbinnings: %" slint_fmt " of %" slint_fmt, b, lb->last_exec_b, lb->nbinnings, lb->max_nbinnings);

  if (lb->last_exec_b == b) return 0;

  if (lb->nbinnings >= lb->max_nbinnings) return -1;
  
  lb->bcws[b] = lb->nbinnings;
  
  ++lb->nbinnings;
  
  lb->last_exec_b = b;

  counts = lb_counts(lb, b, 0);
  bin_counts = lb_bin_counts(lb, b, 0, 0);

#ifdef elem_weight
  weights = lb_weights(lb, b, 0);
  bin_weights = lb_bin_weights(lb, b, 0, 0);
#endif

#ifdef elem_weight
  if (lb->doweights)
  {
    for (k = 0; k < lb->bm->nbins; ++k) counts[k] = weights[k] = 0.0;
  
    for (j = 0; j < lb->nelements; ++j)
    for (k = 0; k < lb->bm->nbins; ++k) bin_counts[j * lb->bm->nbins + k] = bin_weights[j * lb->bm->nbins + k] = 0.0;

    /* for every list of elements */
    for (j = 0; j < lb->nelements; ++j)
    {
      Z_TRACE_IF(B_TRACE_IF, "bin %" slint_fmt ",%" slint_fmt ": size = %" slint_fmt, b, j, lb->bins[b * lb->nelements + j].s.size);
    
      lb->bm->exec(lb->bm, &lb->bins[b * lb->nelements + j], bin_counts, elem_weight_ifelse(bin_weights, NULL));

      lb->bins[b * lb->nelements + j].weight = 0;

      for (k = 0; k < lb->bm->nbins; ++k)
      {
        if (lb->docounts) counts[k] += bin_counts[k];

        weights[k] += bin_weights[k];
        lb->bins[b * lb->nelements + j].weight += bin_weights[k];
      }

      Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slcount_fmt, bin_counts[k], "%" slint_fmt ",%" slint_fmt ": bin_counts =", b, j);
      Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slweight_fmt, bin_weights[k], "%" slint_fmt ",%" slint_fmt ": bin_weights =", b, j);

      bin_counts += lb->bm->nbins;
      bin_weights += lb->bm->nbins;
    }

  } else
#endif
  {
    for (k = 0; k < lb->bm->nbins; ++k) counts[k] = 0;
  
    for (j = 0; j < lb->nelements; ++j)
    for (k = 0; k < lb->bm->nbins; ++k) bin_counts[j * lb->bm->nbins + k] = 0.0;

    /* for every list of elements */
    for (j = 0; j < lb->nelements; ++j)
    {
      Z_TRACE_IF(B_TRACE_IF, "bin %" slint_fmt ",%" slint_fmt ": size = %" slint_fmt, b, j, lb->bins[b * lb->nelements + j].s.size);
    
      lb->bm->exec(lb->bm, &lb->bins[b * lb->nelements + j], bin_counts, elem_weight_ifelse(bin_weights, NULL));

      for (k = 0; k < lb->bm->nbins; ++k)
      {
        counts[k] += bin_counts[k];
      }
    
      Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slcount_fmt, bin_counts[k], "%" slint_fmt ",%" slint_fmt ": bin_counts =", b, j);

      bin_counts += lb->bm->nbins;
    }
  }

  if (lb->docounts)
    Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slcount_fmt, counts[k], "%" slint_fmt ": counts =", b);
#ifdef elem_weight
  if (lb->doweights)
    Z_TRACE_ARRAY_IF(B_TRACE_IF, k, lb->bm->nbins, " %" slweight_fmt, weights[k], "%" slint_fmt ": weights =", b);
#endif

  return 0;
}


slint_t binning_refine(local_bins_t *lb, slint_t b, slint_t k, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_refine */
{
  slint_t j;
  bin_t *new_bin = NULL;

  if (lb->last_new_b != b || lb->last_new_k != k)
  {
    /* update last_new_... */
    lb->last_new_b = b;
    lb->last_new_k = k;
    
    new_bin = &lb->bins_new[lb->nbins_new * lb->nelements];
    
    ++lb->nbins_new;

  } else Z_TRACE_IF(B_TRACE_IF, "no new bin, with b = %" slint_fmt " and k = %" slint_fmt, b, k);

  /* create new bin */
  for (j = 0; j < lb->nelements; ++j)
  {
    lb->bm->refine(lb->bm, &lb->bins[b * lb->nelements + j], k, lb_bin_counts(lb, b, j, 0), lb_bin_weights(lb, b, j, 0), sp, s * lb->nelements + j, new_bin);

    if (new_bin)
    {
      Z_TRACE_IF(B_TRACE_IF, "new bin %td count: %" slint_fmt, new_bin - lb->bins_new, new_bin->s.size);
#ifdef elem_weight
      if (lb->doweights)
        Z_TRACE_IF(B_TRACE_IF, "new bin %td count: %" slweight_fmt, new_bin - lb->bins_new, new_bin->weight);
#endif
      ++new_bin;

    }
  }

/*  Z_TRACE_IF(B_TRACE_IF, "b: %" slint_fmt ", k: %" slint_fmt ", returning %" slint_fmt, b, k, lb->nbins_new - 1);*/

  return lb->nbins_new - 1;
}


slint_t binning_hit(local_bins_t *lb, slint_t b, slint_t k, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_hit */
{
  slint_t j;

  for (j = 0; j < lb->nelements; ++j)
  {
    lb->bm->hit(lb->bm, &lb->bins[b * lb->nelements + j], k, lb_bin_counts(lb, b, j, 0), sp, s * lb->nelements + j);
  }

  return 0;
}


slint_t binning_finalize(local_bins_t *lb, slint_t b, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_finalize */
{
  slint_t j, dc_left;
#ifdef elem_weight
  slweight_t dw_left;
#else
# define dw_left  0
#endif


  Z_TRACE_IF(B_TRACE_IF, "b: %" slint_fmt ", dc: %" slint_fmt ", dw: %" slweight_fmt ", %" slint_fmt ", %" slint_fmt, b, dc, dw, lc_min, lc_max);

  dc_left = dc;
  *lcs = 0.0;
#ifdef elem_weight
  dw_left = dw;
  *lws = 0.0;
#endif

  for (j = 0; j < lb->nelements; ++j)
  {
#ifdef elem_weight
    if (lb->doweights)
      dw_left = dw - *lws;
    else
#endif
      dc_left = dc - *lcs;

    if (lb->bm->finalize(lb->bm, &lb->bins[b * lb->nelements + j], dc_left, dw_left, lc_min - *lcs, lc_max - *lcs, lcs, lws, sp, s * lb->nelements + j)) break;
  }

  return 0;
}


slint_t binning_post(local_bins_t *lb) /* sl_proto, sl_func binning_post */
{
  lb->bm->post(lb->bm);
  
  lb->last_exec_b = -1;

  lb->nbins = lb->nbins_new;

  if (lb->bins == lb->bins0)
  {
    lb->bins = lb->bins1;
    lb->bins_new = lb->bins0;

  } else
  {
    lb->bins = lb->bins0;
    lb->bins_new = lb->bins1;
  }

  return 0;
}
