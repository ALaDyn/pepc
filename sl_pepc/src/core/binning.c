/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/binning.c
 *  timestamp: 2011-02-10 21:20:49 +0100
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


slint_t binning_create(local_bins_t *lb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t doweights, binning_t *bm) /* sl_proto, sl_func binning_create */
{
  slint_t j;

#ifdef elem_weight  
  bm->doweights = doweights;
#endif
  lb->bm = bm;

  lb->nbins = 0;
  lb->max_nbins = max_nbins;

  lb->nbinnings = 0;
  lb->max_nbinnings = (max_nbinnings > 0)?z_min(max_nbinnings,max_nbins):max_nbins;

  lb->nelements = nelements;

#ifdef elem_weight
  lb->doweights = doweights;
  lb->weight_factor = 1 + (lb->doweights != 0);
# define MY_WEIGHT_FACTOR  lb->weight_factor
#else
# define MY_WEIGHT_FACTOR  1
#endif

  lb->bins0 = z_alloc(lb->max_nbins * lb->nelements, sizeof(bin_t));
  lb->bins1 = z_alloc(lb->max_nbins * lb->nelements, sizeof(bin_t));

  lb->bcws = z_alloc(lb->max_nbins, sizeof(slint_t));

  lb->cws = z_alloc(lb->max_nbinnings * MY_WEIGHT_FACTOR * lb->bm->max_nbins, sizeof(slweight_t));

  lb->bin_cws = z_alloc(lb->max_nbinnings * MY_WEIGHT_FACTOR * lb->nelements * bm->max_nbins, sizeof(slweight_t));

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

  z_free(lb->cws);
  z_free(lb->bin_cws);

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

  slweight_t *counts, *bin_counts;
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
        counts[k] += bin_counts[k];

        weights[k] += bin_weights[k];
        lb->bins[b * lb->nelements + j].weight += bin_weights[k];
      }
    
      Z_TRACE_ARRAY_IF(B_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": bin_counts =", " %" slweight_fmt, k, lb->bm->nbins, bin_counts, b, j);
      Z_TRACE_ARRAY_IF(B_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": bin_weights =", " %" slweight_fmt, k, lb->bm->nbins, bin_weights, b, j);

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
    
      Z_TRACE_ARRAY_IF(B_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": bin_counts =", " %" slweight_fmt, k, lb->bm->nbins, bin_counts, b, j);

      bin_counts += lb->bm->nbins;
    }
  }

  Z_TRACE_ARRAY_IF(B_TRACE_IF, "%" slint_fmt ": counts =", " %" slweight_fmt, k, lb->bm->nbins, counts, b);
#ifdef elem_weight
  if (lb->doweights)
    Z_TRACE_ARRAY_IF(B_TRACE_IF, "%" slint_fmt ": weights =", " %" slweight_fmt, k, lb->bm->nbins, weights, b);
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


slint_t binning_finalize(local_bins_t *lb, slint_t b, slweight_t dcw, slint_t lc_min, slint_t lc_max, slweight_t *lcw, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_finalize */
{
  slint_t j;
  slweight_t dcw_left;


  Z_TRACE_IF(B_TRACE_IF, "b: %" slint_fmt ", dcw: %" slweight_fmt ", %" slint_fmt ", %" slint_fmt, b, dcw, lc_min, lc_max);

  lcw[0] =
#ifdef elem_weight
    lcw[1] =
#endif
    0.0;

  dcw_left = dcw;

  for (j = 0; j < lb->nelements; ++j)
  {
#ifdef elem_weight
    if (lb->doweights)
      dcw_left = dcw - lcw[1];
    else
#endif
      dcw_left = dcw - lcw[0];

    if (lb->bm->finalize(lb->bm, &lb->bins[b * lb->nelements + j], dcw_left, lc_min - lcw[0], lc_max - lcw[0], lcw, sp, s * lb->nelements + j)) break;
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
