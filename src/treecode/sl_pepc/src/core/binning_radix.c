/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/binning_radix.c
 *  
 */


/* sl_macro BR_TRACE_IF */


#include "sl_common.h"


#ifndef BR_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define BR_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define BR_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t binning_radix_create(binning_t *bm, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted) /* sl_proto, sl_func binning_radix_create */
{
  bm->pre = binning_radix_pre;
  bm->exec = binning_radix_exec;
  bm->refine = binning_radix_refine;
  bm->hit = binning_radix_hit;
  bm->finalize = binning_radix_finalize;
  bm->post = binning_radix_post;

  bm->sorted = sorted;

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth < 0) rwidth = sort_radix_width_default;

  bm->nbins = 0;
  bm->max_nbins = z_powof2_typed(rwidth, slkey_pure_t);

  bm->bd.radix.rhigh = rhigh;
  bm->bd.radix.rlow = rlow;
  bm->bd.radix.rwidth = rwidth;

  elements_alloc(&bm->bd.radix.sx, 1, SLCM_ALL);

  return 0;
}


slint_t binning_radix_destroy(binning_t *bm) /* sl_proto, sl_func binning_radix_destroy */
{
  elements_free(&bm->bd.radix.sx);

  return 0;
}


slint_t binning_radix_pre(binning_t *bm) /* sl_proto, sl_func binning_radix_pre */
{
  bm->bd.radix.rcurrent = z_min(bm->bd.radix.rwidth, bm->bd.radix.rhigh - bm->bd.radix.rlow + 1);
  bm->bd.radix.rhigh -= (bm->bd.radix.rcurrent > 0)?bm->bd.radix.rcurrent - 1:bm->bd.radix.rhigh;

  bm->nbins = (bm->bd.radix.rcurrent > 0)?z_powof2(bm->bd.radix.rcurrent):1;
  bm->bd.radix.bit_mask = bm->nbins - 1;
  
  return 0;
}


slint_t binning_radix_exec(binning_t *bm, bin_t *bin, slcount_t *counts, slweight_t *weights) /* sl_proto, sl_func binning_radix_exec */
{
  elements_t xi, end;
  slkey_pure_t k;
  slint_t i, *c, special;


  elem_assign_at(&bin->s, bin->s.size, &end);

  if (bm->nbins > 1)
  {
#ifdef elem_weight
    if (bm->doweights)
    {
      /* counts and weights in every class */
      for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi))
      {
        k = key_radix_key2class(key_purify(*xi.keys), bm->bd.radix.rhigh, bm->bd.radix.bit_mask);
        counts[k] += 1;
        weights[k] += elem_weight(&xi, 0);
      }

    } else
#endif
    {
      /* counts in every class */
      if (bm->sorted & SL_SORTED_IN)
      {
        special = ((bm->bd.radix.bit_mask << bm->bd.radix.rhigh) < 0);
      
        Z_TRACE_IF(BR_TRACE_IF, "bitmask = %" key_pure_type_fmt " -> %s", (bm->bd.radix.bit_mask << bm->bd.radix.rhigh), special?"special signedness handling":"normal");

/*        if (BR_TRACE_IF) elements_print_keys(&bin->s);*/

        if (special)
        {
          special = sl_search_binary_sign_switch(&bin->s);
        
          Z_TRACE_IF(BR_TRACE_IF, "sign switch @ %" slint_fmt, special);
          
          /* make bit mask 1 bit smaller (and erase the sign bit) */
          bm->bd.radix.bit_mask = (bm->bd.radix.bit_mask >> 1) & ~(~((slkey_pure_t) 0) << (sizeof(slkey_pure_t)*8-1));

          elem_assign(&bin->s, &xi);
          xi.size = special;
          
          for (i = 0; i < bm->bd.radix.bit_mask; ++i)
          {
            k = ((slkey_pure_t) i + 1) << bm->bd.radix.rhigh;
            counts[i] = sl_search_binary_lt_bmask(&xi, k, bm->bd.radix.bit_mask << bm->bd.radix.rhigh);
            Z_TRACE_IF(BR_TRACE_IF, "i = %" slint_fmt " of %" key_pure_type_fmt " << %" slint_fmt ", k = %" key_pure_type_fmt ", searching in %" slint_fmt " -> %" slcount_fmt, i, bm->bd.radix.bit_mask, bm->bd.radix.rhigh, k, xi.size, counts[i]);
            elem_add(&xi, (slint_t) counts[i]);
            xi.size -= counts[i];
          }
          counts[i] = xi.size;
          
          elem_assign_at(&bin->s, special, &xi);
          xi.size = bin->s.size - special;

          for (i = 0; i < bm->bd.radix.bit_mask; ++i)
          {
            k = ((slkey_pure_t) i + 1) << bm->bd.radix.rhigh;
            counts[i + bm->bd.radix.bit_mask + 1] = sl_search_binary_lt_bmask(&xi, k, bm->bd.radix.bit_mask << bm->bd.radix.rhigh);
            Z_TRACE_IF(BR_TRACE_IF, "i = %" slint_fmt " of %" key_pure_type_fmt " << %" slint_fmt ", k = %" key_pure_type_fmt ", searching in %" slint_fmt " -> %" slcount_fmt, i, bm->bd.radix.bit_mask, bm->bd.radix.rhigh, k, xi.size, counts[i + bm->bd.radix.bit_mask + 1]);
            elem_add(&xi, (slint_t) counts[i + bm->bd.radix.bit_mask + 1]);
            xi.size -= counts[i + bm->bd.radix.bit_mask + 1];
          }
          counts[i + bm->bd.radix.bit_mask + 1] = xi.size;

        } else
        {
          elem_assign(&bin->s, &xi);
          for (i = 0; i < bm->bd.radix.bit_mask; ++i)
          {
            k = ((slkey_pure_t) i + 1) << bm->bd.radix.rhigh;
            counts[i] = sl_search_binary_lt_bmask(&xi, k, bm->bd.radix.bit_mask << bm->bd.radix.rhigh);
            Z_TRACE_IF(BR_TRACE_IF, "i = %" slint_fmt " of %" key_pure_type_fmt " << %" slint_fmt ", k = %" key_pure_type_fmt ", searching in %" slint_fmt " -> %" slcount_fmt, i, bm->bd.radix.bit_mask, bm->bd.radix.rhigh, k, xi.size, counts[i]);
            elem_add(&xi, (slint_t) counts[i]);
            xi.size -= counts[i];
          }
          counts[i] = xi.size;
        }
        
      } else
      {
        for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi))
        {
          ++counts[key_radix_key2class(key_purify(*xi.keys), bm->bd.radix.rhigh, bm->bd.radix.bit_mask)];
        }
      }
    }

    if (!(bm->sorted & SL_SORTED_IN) && (bm->sorted & SL_SORTED_OUT))
    {
      c = z_alloca(bm->bd.radix.bit_mask + 1, sizeof(slint_t));
      for (i = 0; i < bm->bd.radix.bit_mask + 1; ++i) c[i] = counts[i];
      splitx_radix(&bin->s, &bm->bd.radix.sx, bm->bd.radix.bit_mask + 1, bm->bd.radix.rhigh, c);
      z_freea(c);
    }

  } else
  {
    /* total counts */
    counts[0] += bin->s.size;

#ifdef elem_weight
    /* total weights */
    if (bm->doweights)
    {
      for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi)) weights[0] += elem_weight(&xi, 0);
    }
#endif
  }

  Z_TRACE_ARRAY_IF(BR_TRACE_IF, k, bm->nbins, " %" slcount_fmt, counts[k], "counts of %" slint_fmt " @ %p:", bin->s.size, bin->s.keys);
#ifdef elem_weight
  if (bm->doweights)
    Z_TRACE_ARRAY_IF(BR_TRACE_IF, k, bm->nbins, " %" slweight_fmt, weights[k], "weights of %" slint_fmt " @ %p:", bin->s.size, bin->s.keys);
#endif

  return 0;
}


slint_t binning_radix_refine(binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, slweight_t *weights, splitter_t *sp, slint_t s, bin_t *new_bin) /* sl_proto, sl_func binning_radix_refine */
{
  slint_t l, lcs;

  lcs = 0;
  for (l = 0; l < k; ++l) lcs += counts[l];

  if (new_bin)
  {
    elem_assign_at(&bin->s, lcs, &new_bin->s);
    new_bin->s.size = counts[k];
  
    Z_TRACE_IF(BR_TRACE_IF, "new bin count: %" slint_fmt, new_bin->s.size);

#ifdef elem_weight
    if (bm->doweights)
    {
      new_bin->weight = weights[k];

      Z_TRACE_IF(BR_TRACE_IF, "new bin weight: %" slweight_fmt, new_bin->weight);
    }
#endif
  }

  sp->displs[s] += lcs;

  Z_TRACE_IF(BR_TRACE_IF, "displs[%" slint_fmt "] += %" slint_fmt " = %d", s, lcs, sp->displs[s]);

  return 0;
}


slint_t binning_radix_hit(binning_t *bm, bin_t *bin, slint_t k, slcount_t *counts, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_radix_hit */
{
  slint_t l;

  for (l = 0; l < k; ++l) sp->displs[s] += counts[l];

  Z_TRACE_IF(BR_TRACE_IF, "displs[%" slint_fmt "] += ... = %d", s, sp->displs[s]);

  return 0;
}


slint_t binning_radix_finalize(binning_t *bm, bin_t *bin, slint_t dc, slweight_t dw, slint_t lc_min, slint_t lc_max, slcount_t *lcs, slweight_t *lws, splitter_t *sp, slint_t s) /* sl_proto, sl_func binning_radix_finalize */
{
  slint_t lc, r;
#ifdef elem_weight
  elements_t xi, end;
  slweight_t lw = 0.0;
#endif


  Z_TRACE_IF(BR_TRACE_IF, "bin size: %" slint_fmt ", dc = %" slint_fmt ", lc: %" slint_fmt " - %" slint_fmt ", *lcs = %" slcount_fmt, bin->s.size, dc, lc_min, lc_max, *lcs);
#ifdef elem_weight
  if (bm->doweights)
    Z_TRACE_IF(BR_TRACE_IF, "bin weight: %" slweight_fmt ", dw = %" slweight_fmt ", lc: %" slint_fmt " - %" slint_fmt ", *lws = %" slweight_fmt, bin->weight, dw, lc_min, lc_max, *lws);
#endif

  r = 0;

#ifdef elem_weight
  if (bm->doweights)
  {
    lc = 0;
    lw = 0.0;

    if (bin->s.size <= lc_min || (dw >= bin->weight && bin->s.size <= lc_max))
    {
      lc = bin->s.size;
      lw = bin->weight;

    } else
    {
      if (0 < lc_max)
      {
        elem_assign_at(&bin->s, bin->s.size, &end);

        lw = dw;

        for (elem_assign(&bin->s, &xi); xi.keys < end.keys; elem_inc(&xi))
        {
          ++lc;
          lw -= elem_weight(&xi, 0);
        
          if (lc <= lc_min) continue;

          if (lw < 0.0 || lc > lc_max)
          {
            lw += elem_weight(&xi, 0);
            --lc;
            break;
          }
        }
      
        lw = dw - lw;
      }

      r = 1;
    }

  } else
#endif
  {
    lc = z_min(dc, bin->s.size);
    
    r = (lc >= dc);
  }

  *lcs += lc;
  Z_TRACE_IF(BR_TRACE_IF, "*lcs = %" slcount_fmt " + %" slint_fmt " = %" slcount_fmt, (slcount_t) (*lcs - lc), lc, *lcs);
#ifdef elem_weight
  if (bm->doweights)
  {
    *lws += lw;
    Z_TRACE_IF(BR_TRACE_IF, "*lws = %" slweight_fmt " + %" slweight_fmt " = %" slweight_fmt, *lws - lw, lw, *lws);
  }
#endif

  sp->displs[s] += lc;

  Z_TRACE_IF(BR_TRACE_IF, "displs[%" slint_fmt "] += %" slint_fmt " = %d", s, lc, sp->displs[s]);

  return r;
}


slint_t binning_radix_post(binning_t *bm) /* sl_proto, sl_func binning_radix_post */
{
  --bm->bd.radix.rhigh;

  return 0;
}
