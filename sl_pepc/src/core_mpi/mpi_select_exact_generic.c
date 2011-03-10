/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_select_exact_generic.c
 *  timestamp: 2011-03-03 22:37:02 +0100
 *  
 */


/* sl_macro MSEG_ROOT */
/* sl_macro MSEG_REDUCEBCAST_THRESHOLD */
/* sl_macro MSEG_BORDER_UPDATE_REDUCTION */
/* sl_macro MSEG_BORDER_UPDATE_FULL */
/* sl_macro MSEG_INFO */
/* sl_macro MSEG_TRACE_IF */

  
#include "sl_common.h"


/* config */
/*#define SYNC_ON_INIT
#define SYNC_ON_EXIT*/
#define HAVENT_MPI_IN_PLACE

/*#define PRINT_SDISPLS*/
/*#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/

/*#define VERIFY*/

/*#define OLD*/

typedef struct _border_info_t {
#ifdef OLD
#define LO_LE  0
#define HI_LE  1
#define LO_RI  2
#define HI_RI  3
#ifdef MSEG_BORDER_UPDATE_FULL
  slint_t done;
#endif
  slint_t update;
  slint_t crange[2], cmmlr[4];
#ifdef elem_weight
  slweight_t wrange[2], wmmlr[4];
#endif
#else
  slint_t crange[2], ccurrent[2];
#ifdef elem_weight
  slweight_t wrange[2], wcurrent[2];
#endif
#endif

} border_info_t;


#define LO  0
#define HI  1

#define CNT_LO  0
#define CNT_HI  1
#define WHT_LO  2
#define WHT_HI  3


#ifdef MSEG_ROOT
int mseg_root = -1;  /* sl_global, sl_var mseg_root */
#endif

#define REDUCEBCAST_ROOT  0

#if !defined(MSEG_REDUCEBCAST_THRESHOLD) && defined(GLOBAL_REDUCEBCAST_THRESHOLD)
# define MSEG_REDUCEBCAST_THRESHOLD  GLOBAL_REDUCEBCAST_THRESHOLD
#endif

#ifdef MSEG_BORDER_UPDATE_REDUCTION
double mseg_border_update_count_reduction = 0.0;  /* sl_global, sl_var mseg_border_update_count_reduction */
# ifdef elem_weight
double mseg_border_update_weight_reduction = 0.0;  /* sl_global, sl_var mseg_border_update_weight_reduction */
# endif
#endif

#ifdef MSEG_BORDER_UPDATE_FULL
slint_t mseg_border_update_full = 0;  /* sl_global, sl_var mseg_border_update_full */
#endif

#ifdef MSEG_INFO
slint_t mseg_info_rounds = 0;             /* sl_global, sl_var mseg_info_rounds */
slint_t *mseg_info_finish_rounds = NULL;  /* sl_global, sl_var mseg_info_finish_rounds */
double mseg_info_finish_rounds_avg = 0;   /* sl_global, sl_var mseg_info_finish_rounds_avg */
#endif

slint_t mseg_binnings = -1;  /* sl_global, sl_var mseg_binnings */

slint_t mseg_finalize_mode = SL_MSEG_FM_EXACT;  /* sl_global, sl_var mseg_finalize_mode */

#ifndef MSEG_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSEG_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSEG_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif

#define MSEG_ASSERT_IF  (rank == 0)


#ifdef OLD
void border_update_old(slint_t doweights, border_info_t *bi, partcond_intern_t *pc, slint_t dir, slint_t check_inconsistence, slint_t reduction) /* sl_func border_update_old */
{
#ifdef MSEG_BORDER_UPDATE_REDUCTION
  slint_t count_reduction;
# ifdef elem_weight
  slweight_t weight_reduction;
# endif
#endif

  /* forward */
  if (dir > 0)
  {
    /* init from range */
    bi[0].cmmlr[LO_LE] = bi[0].crange[0];
    bi[0].cmmlr[HI_LE] = bi[0].crange[1];

    /* init from min/max */
    if (pc[0].pcm & SLPC_COUNTS_MM)
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction)
      {
        count_reduction = z_round((bi[-1].cmmlr[HI_LE] - bi[-1].cmmlr[LO_LE]) * 0.5 * mseg_border_update_count_reduction);
        Z_ASSERT(count_reduction >= 0);

        Z_TRACE_IF(MSEG_TRACE_IF, "forward: count_reduction: %" slint_fmt, count_reduction);

        bi[0].cmmlr[LO_LE] = z_min(bi[-1].cmmlr[LO_LE] + count_reduction + pc[0].count_min, bi[0].cmmlr[HI_RI]);
        bi[0].cmmlr[HI_LE] = z_max(bi[-1].cmmlr[HI_LE] - count_reduction + pc[0].count_max, bi[0].cmmlr[LO_RI]);

        Z_TRACE_IF(MSEG_TRACE_IF, "forward: count[min/max-left]: min(%" slint_fmt " + %" slint_fmt " + %" slint_fmt ", %" slint_fmt "), max(%" slint_fmt " - % "slint_fmt " + %" slint_fmt ", %" slint_fmt ")",
          bi[-1].cmmlr[LO_LE], count_reduction, pc[0].count_min, bi[0].cmmlr[HI_RI], bi[-1].cmmlr[HI_LE], count_reduction, pc[0].count_max, bi[0].cmmlr[LO_RI]);

      } else
#endif
      {
        bi[0].cmmlr[LO_LE] = bi[-1].cmmlr[LO_LE] + pc[0].count_min;
        bi[0].cmmlr[HI_LE] = bi[-1].cmmlr[HI_LE] + pc[0].count_max;

        Z_TRACE_IF(MSEG_TRACE_IF, "forward: count[min/max-left]: %" slint_fmt " + %" slint_fmt ", %" slint_fmt " + %" slint_fmt,
          bi[-1].cmmlr[LO_LE], pc[0].count_min, bi[-1].cmmlr[HI_LE], pc[0].count_max);
      }
    }

    /* check against low/high */
    if (pc[0].pcm & SLPC_COUNTS_LH)
    {
      if (bi[0].cmmlr[LO_LE] < pc[1].count_low)  bi[0].cmmlr[LO_LE] = pc[1].count_low;
      if (bi[0].cmmlr[HI_LE] > pc[0].count_high) bi[0].cmmlr[HI_LE] = pc[0].count_high;
    }

    /* fit to range */
    bi[0].cmmlr[LO_LE] = z_minmax(bi[0].crange[0], bi[0].cmmlr[LO_LE], bi[0].crange[1]);
    bi[0].cmmlr[HI_LE] = z_minmax(bi[0].crange[0], bi[0].cmmlr[HI_LE], bi[0].crange[1]);

#ifdef elem_weight
    if (doweights)
    {
      /* init from range */
      bi[0].wmmlr[LO_LE] = bi[0].wrange[0];
      bi[0].wmmlr[HI_LE] = bi[0].wrange[1];

      /* init from min/max */
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
# ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          weight_reduction = (bi[-1].wmmlr[HI_LE] - bi[-1].wmmlr[LO_LE]) * 0.5 * mseg_border_update_weight_reduction;
          Z_ASSERT(weight_reduction >= 0.0);

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: weight_reduction: %" slweight_fmt, weight_reduction);
  
          bi[0].wmmlr[LO_LE] = z_min(bi[-1].wmmlr[LO_LE] + weight_reduction + pc[0].weight_min, bi[0].wmmlr[HI_RI]);
          bi[0].wmmlr[HI_LE] = z_max(bi[-1].wmmlr[HI_LE] - weight_reduction + pc[0].weight_max, bi[0].wmmlr[LO_RI]);

          Z_TRACE_IF(MSEG_TRACE_IF, "weight_reduction: weight[min/max-left]: min(%" slweight_fmt " + %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt "), max(%" slweight_fmt " - %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt ")",
            bi[-1].wmmlr[LO_LE], weight_reduction, pc[0].weight_min, bi[0].wmmlr[HI_RI], bi[-1].wmmlr[HI_LE], weight_reduction, pc[0].weight_max, bi[0].wmmlr[LO_RI]);

        } else
# endif
        {
          bi[0].wmmlr[LO_LE] = bi[-1].wmmlr[LO_LE] + pc[0].weight_min;
          bi[0].wmmlr[HI_LE] = bi[-1].wmmlr[HI_LE] + pc[0].weight_max;

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: weight[min/max-left]: %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt " + %" slweight_fmt,
            bi[-1].wmmlr[LO_LE], pc[0].weight_min, bi[-1].wmmlr[HI_LE], pc[0].weight_max);
        }
      }

      /* check against low/high (on demand) */
      if (pc[0].pcm & SLPC_WEIGHTS_LH)
      {
        if (bi[0].wmmlr[LO_LE] < pc[1].weight_low)  bi[0].wmmlr[LO_LE] = pc[1].weight_low;
        if (bi[0].wmmlr[HI_LE] > pc[0].weight_high) bi[0].wmmlr[HI_LE] = pc[0].weight_high;
      }

      /* fit to range */
      bi[0].wmmlr[LO_LE] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[LO_LE], bi[0].wrange[1]);
      bi[0].wmmlr[HI_LE] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[HI_LE], bi[0].wrange[1]);

    } else
#endif
    { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }
      
  } else /* backward */
  {
    /* init from range */
    bi[0].cmmlr[LO_RI] = bi[0].crange[0];
    bi[0].cmmlr[HI_RI] = bi[0].crange[1];

    /* init from min/max */
    if (pc[0].pcm & SLPC_COUNTS_MM)
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction)
      {
        count_reduction = z_round((bi[1].cmmlr[HI_RI] - bi[1].cmmlr[LO_RI]) * 0.5 * mseg_border_update_count_reduction);
        Z_ASSERT(count_reduction >= 0);

        Z_TRACE_IF(MSEG_TRACE_IF, "backward: count_reduction: %" slint_fmt, count_reduction);

        bi[0].cmmlr[HI_RI] = z_max(bi[1].cmmlr[HI_RI] - count_reduction - pc[1].count_min, bi[0].cmmlr[LO_LE]);
        bi[0].cmmlr[LO_RI] = z_min(bi[1].cmmlr[LO_RI] + count_reduction - pc[1].count_max, bi[0].cmmlr[HI_LE]);

        Z_TRACE_IF(MSEG_TRACE_IF, "backward: count[min/max-right]: max(%" slint_fmt " - %" slint_fmt " + %" slint_fmt ", %" slint_fmt "), min(%" slint_fmt " + %" slint_fmt " - %" slint_fmt ", %" slint_fmt ")",
          bi[1].cmmlr[HI_RI], count_reduction, pc[1].count_min, bi[0].cmmlr[LO_LE], bi[1].cmmlr[LO_RI], count_reduction, pc[1].count_max, bi[0].cmmlr[HI_LE]);

      } else
#endif
      {
        bi[0].cmmlr[HI_RI] = bi[1].cmmlr[HI_RI] - pc[1].count_min;
        bi[0].cmmlr[LO_RI] = bi[1].cmmlr[LO_RI] - pc[1].count_max;

        Z_TRACE_IF(MSEG_TRACE_IF, "backward: count[min/max-right]: %" slint_fmt " - %" slint_fmt ", %" slint_fmt " - %" slint_fmt "",
          bi[1].cmmlr[HI_RI], pc[1].count_min, bi[1].cmmlr[LO_RI], pc[1].count_max);
      }
    }

    /* check against low/high (on demand) */
    if (pc[0].pcm & SLPC_COUNTS_LH)
    {
      if (bi[0].cmmlr[LO_RI] < pc[1].count_low)  bi[0].cmmlr[LO_RI] = pc[1].count_low;
      if (bi[0].cmmlr[HI_RI] > pc[0].count_high) bi[0].cmmlr[HI_RI] = pc[0].count_high;
    }

    /* fit to range */
    bi[0].cmmlr[HI_RI] = z_minmax(bi[0].crange[0], bi[0].cmmlr[HI_RI], bi[0].crange[1]);
    bi[0].cmmlr[LO_RI] = z_minmax(bi[0].crange[0], bi[0].cmmlr[LO_RI], bi[0].crange[1]);

#ifdef elem_weight
    if (doweights)
    {
      /* init from range */
      bi[0].wmmlr[LO_RI] = bi[0].wrange[0];
      bi[0].wmmlr[HI_RI] = bi[0].wrange[1];

      /* init from min/max */
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          weight_reduction = (bi[1].wmmlr[HI_RI] - bi[1].wmmlr[LO_RI]) * 0.5 * mseg_border_update_weight_reduction;
          Z_ASSERT(weight_reduction >= 0.0);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight_reduction: %" slweight_fmt, weight_reduction);

          bi[0].wmmlr[HI_RI] = z_max(bi[1].wmmlr[HI_RI] - weight_reduction - pc[1].weight_min, bi[0].wmmlr[LO_LE]);
          bi[0].wmmlr[LO_RI] = z_min(bi[1].wmmlr[LO_RI] + weight_reduction - pc[1].weight_max, bi[0].wmmlr[HI_LE]);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight[min/max-right]: max(%" slweight_fmt " - %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt "), min(%" slweight_fmt " + %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt ")",
            bi[1].wmmlr[HI_RI], weight_reduction, pc[1].weight_min, bi[0].wmmlr[LO_LE], bi[1].wmmlr[LO_RI], weight_reduction, pc[1].weight_max, bi[0].wmmlr[HI_LE]);

        } else
#endif
        {
          bi[0].wmmlr[HI_RI] = bi[1].wmmlr[HI_RI] - pc[1].weight_min;
          bi[0].wmmlr[LO_RI] = bi[1].wmmlr[LO_RI] - pc[1].weight_max;

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight[min/max-right]: %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt " - %" slweight_fmt,
            bi[1].wmmlr[HI_RI], pc[1].weight_min, bi[1].wmmlr[LO_RI], pc[1].weight_max);
        }
      }

      /* check against low/high (on demand) */
      if (pc[0].pcm & SLPC_WEIGHTS_LH)
      {
        if (bi[0].wmmlr[LO_RI] < pc[1].weight_low)  bi[0].wmmlr[LO_RI] = pc[1].weight_low;
        if (bi[0].wmmlr[HI_RI] > pc[0].weight_high) bi[0].wmmlr[HI_RI] = pc[0].weight_high;
      }

      /* fit to range */
      bi[0].wmmlr[HI_RI] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[HI_RI], bi[0].wrange[1]);
      bi[0].wmmlr[LO_RI] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[LO_RI], bi[0].wrange[1]);

    } else
#endif
    { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }
  }

  Z_TRACE_IF(MSEG_TRACE_IF, "count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
    bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
  if (doweights)
    Z_TRACE_IF(MSEG_TRACE_IF, "weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
      bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);
  else
#endif
    Z_TRACE_IF(MSEG_TRACE_IF, "");

  if (check_inconsistence)
  {
    /* check against inconsistence */
    if (bi[0].cmmlr[LO_LE] > bi[0].cmmlr[HI_RI]) bi[0].cmmlr[LO_LE] = bi[0].cmmlr[HI_RI] = (bi[0].cmmlr[LO_LE] + bi[0].cmmlr[HI_RI]) / 2;
    if (bi[0].cmmlr[HI_LE] < bi[0].cmmlr[LO_RI]) bi[0].cmmlr[HI_LE] = bi[0].cmmlr[LO_RI] = (bi[0].cmmlr[HI_LE] + bi[0].cmmlr[LO_RI]) / 2;

    Z_TRACE_IF(MSEG_TRACE_IF, "consistence count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
      bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
    if (doweights)
    {
      if (bi[0].wmmlr[LO_LE] > bi[0].wmmlr[HI_RI]) bi[0].wmmlr[LO_LE] = bi[0].wmmlr[HI_RI] = (bi[0].wmmlr[LO_LE] + bi[0].wmmlr[HI_RI]) / 2;
      if (bi[0].wmmlr[HI_LE] < bi[0].wmmlr[LO_RI]) bi[0].wmmlr[HI_LE] = bi[0].wmmlr[LO_RI] = (bi[0].wmmlr[HI_LE] + bi[0].wmmlr[LO_RI]) / 2;

      Z_TRACE_IF(MSEG_TRACE_IF, "consistence weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
        bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);

    } else
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");
  }
}


#ifdef MSEG_BORDER_UPDATE_FULL
void border_update_full_old(slint_t doweights, border_info_t *bi, partcond_intern_t *pc, slint_t check_inconsistence) /* sl_func border_update_full_old */
{
  slint_t lo = -1, hi = -1;

  bi[0].cmmlr[LO_LE] = bi[0].cmmlr[LO_RI] = bi[0].crange[0];
  bi[0].cmmlr[HI_LE] = bi[0].cmmlr[HI_RI] = bi[0].crange[1];

#ifdef elem_weight
  if (doweights)
  {
    bi[0].wmmlr[LO_LE] = bi[0].wmmlr[LO_RI] = bi[0].wrange[0];
    bi[0].wmmlr[HI_LE] = bi[0].wmmlr[HI_RI] = bi[0].wrange[1];
  }
#endif

  /* init from min/max */
  if (pc[0].pcm & (SLPC_COUNTS_MM|SLPC_WEIGHTS_MM))
  {
    /* seach backward */
    for (lo = 0; (!bi[lo - 1].done); --lo);

    /* seach forward */
    for (hi = 0; (!bi[hi + 1].done); ++hi);

    Z_TRACE_IF(MSEG_TRACE_IF, "found: lo: %" slint_fmt ", hi: %" slint_fmt "", lo, hi);

    if (pc[0].pcm & SLPC_COUNTS_MM)
    {
      bi[0].cmmlr[LO_LE] = bi[lo - 1].cmmlr[LO_LE];
      bi[0].cmmlr[HI_LE] = bi[lo - 1].cmmlr[HI_LE];
      Z_ASSERT(bi[0].cmmlr[LO_LE] == bi[0].cmmlr[HI_LE]);

      bi[0].cmmlr[HI_RI] = bi[hi + 1].cmmlr[HI_RI]; 
      bi[0].cmmlr[LO_RI] = bi[hi + 1].cmmlr[LO_RI]; 
      Z_ASSERT(bi[0].cmmlr[HI_RI] == bi[0].cmmlr[LO_RI]);
    }

    Z_TRACE_IF(MSEG_TRACE_IF, "start: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
      bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
    if (doweights)
    {
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
        bi[0].wmmlr[LO_LE] = bi[lo - 1].wmmlr[LO_LE];
        bi[0].wmmlr[HI_LE] = bi[lo - 1].wmmlr[HI_LE];
        Z_ASSERT(bi[0].wmmlr[LO_LE] == bi[0].wmmlr[HI_LE]);

        bi[0].wmmlr[HI_RI] = bi[hi + 1].wmmlr[HI_RI];
        bi[0].wmmlr[LO_RI] = bi[hi + 1].wmmlr[LO_RI];
        Z_ASSERT(bi[0].wmmlr[HI_RI] == bi[0].wmmlr[LO_RI]);
      }
    }
#endif
  
    /* go backward */
    for (; lo <= 0; ++lo)
    {
      if (pc[0].pcm & SLPC_COUNTS_MM)
      {
        Z_TRACE_IF(MSEG_TRACE_IF, "- count[min/max-left]: %" slint_fmt " + %" slint_fmt ", %" slint_fmt " + %" slint_fmt "",
          bi[0].cmmlr[LO_LE], pc[0].count_min, bi[0].cmmlr[HI_LE], pc[0].count_max);

        bi[0].cmmlr[LO_LE] = bi[0].cmmlr[LO_LE] + pc[lo].count_min;
        bi[0].cmmlr[HI_LE] = bi[0].cmmlr[HI_LE] + pc[lo].count_max;
      }

#ifdef elem_weight
      if (doweights && pc[0].pcm & SLPC_WEIGHTS_MM)
      {
        bi[0].wmmlr[LO_LE] = bi[0].wmmlr[LO_LE] + pc[lo].weight_min;
        bi[0].wmmlr[HI_LE] = bi[0].wmmlr[HI_LE] + pc[lo].weight_max;
      }
#endif
    }

    /* go forward */
    for (; 0 <= hi; --hi)
    {
      if (pc[0].pcm & SLPC_COUNTS_MM)
      {
        Z_TRACE_IF(MSEG_TRACE_IF, "+ count[min/max-left]: %" slint_fmt " - %" slint_fmt ", %" slint_fmt " - %" slint_fmt "",
          bi[0].cmmlr[HI_RI], pc[hi + 1].count_min, bi[0].cmmlr[LO_RI], pc[hi + 1].count_max);

        bi[0].cmmlr[HI_RI] = bi[0].cmmlr[HI_RI] - pc[hi + 1].count_min;
        bi[0].cmmlr[LO_RI] = bi[0].cmmlr[LO_RI] - pc[hi + 1].count_max;
      }

#ifdef elem_weight
      if (doweights && pc[0].pcm & SLPC_WEIGHTS_MM)
      {
        bi[0].wmmlr[HI_RI] = bi[0].wmmlr[HI_RI] + pc[hi + 1].weight_min;
        bi[0].wmmlr[LO_RI] = bi[0].wmmlr[LO_RI] + pc[hi + 1].weight_max;
      }
#endif
    }

  } else Z_TRACE_IF(MSEG_TRACE_IF, "no search");

  /* check against low/high */
  if (pc[0].pcm & SLPC_COUNTS_LH)
  {
    if (bi[0].cmmlr[LO_LE] < pc[1].count_low)  bi[0].cmmlr[LO_LE] = pc[1].count_low;
    if (bi[0].cmmlr[HI_LE] > pc[0].count_high) bi[0].cmmlr[HI_LE] = pc[0].count_high;
    if (bi[0].cmmlr[LO_RI] < pc[1].count_low)  bi[0].cmmlr[LO_RI] = pc[1].count_low;
    if (bi[0].cmmlr[HI_RI] > pc[0].count_high) bi[0].cmmlr[HI_RI] = pc[0].count_high;
  }

#ifdef elem_weight
  if (doweights && pc[0].pcm & SLPC_WEIGHTS_MM)
  {
    if (bi[0].wmmlr[LO_LE] < pc[1].weight_low)  bi[0].wmmlr[LO_LE] = pc[1].weight_low;
    if (bi[0].wmmlr[HI_LE] > pc[0].weight_high) bi[0].wmmlr[HI_LE] = pc[0].weight_high;
    if (bi[0].wmmlr[LO_RI] < pc[1].weight_low)  bi[0].wmmlr[LO_RI] = pc[1].weight_low;
    if (bi[0].wmmlr[HI_RI] > pc[0].weight_high) bi[0].wmmlr[HI_RI] = pc[0].weight_high;
  }
#endif

  /* fit to range */
  bi[0].cmmlr[LO_LE] = z_minmax(bi[0].crange[0], bi[0].cmmlr[LO_LE], bi[0].crange[1]);
  bi[0].cmmlr[HI_LE] = z_minmax(bi[0].crange[0], bi[0].cmmlr[HI_LE], bi[0].crange[1]);
  bi[0].cmmlr[HI_RI] = z_minmax(bi[0].crange[0], bi[0].cmmlr[HI_RI], bi[0].crange[1]);
  bi[0].cmmlr[LO_RI] = z_minmax(bi[0].crange[0], bi[0].cmmlr[LO_RI], bi[0].crange[1]);

#ifdef elem_weight
  if (doweights)
  {
    bi[0].wmmlr[LO_LE] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[LO_LE], bi[0].wrange[1]);
    bi[0].wmmlr[HI_LE] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[HI_LE], bi[0].wrange[1]);
    bi[0].wmmlr[HI_RI] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[HI_RI], bi[0].wrange[1]);
    bi[0].wmmlr[LO_RI] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[LO_RI], bi[0].wrange[1]);
  }
#endif

  Z_TRACE_IF(MSEG_TRACE_IF, "count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
    bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
  if (doweights)
    Z_TRACE_IF(MSEG_TRACE_IF, "weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
      bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);
#endif

  if (check_inconsistence)
  {
    /* check against inconsistence */
    if (bi[0].cmmlr[LO_LE] > bi[0].cmmlr[HI_RI]) bi[0].cmmlr[LO_LE] = bi[0].cmmlr[HI_RI] = (bi[0].cmmlr[LO_LE] + bi[0].cmmlr[HI_RI]) / 2;
    if (bi[0].cmmlr[HI_LE] < bi[0].cmmlr[LO_RI]) bi[0].cmmlr[HI_LE] = bi[0].cmmlr[LO_RI] = (bi[0].cmmlr[HI_LE] + bi[0].cmmlr[LO_RI]) / 2;

    Z_TRACE_IF(MSEG_TRACE_IF, "consistence count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
      bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
    if (doweights)
    {
      if (bi[0].wmmlr[LO_LE] > bi[0].wmmlr[HI_RI]) bi[0].wmmlr[LO_LE] = bi[0].wmmlr[HI_RI] = (bi[0].wmmlr[LO_LE] + bi[0].wmmlr[HI_RI]) / 2;
      if (bi[0].wmmlr[HI_LE] < bi[0].wmmlr[LO_RI]) bi[0].wmmlr[HI_LE] = bi[0].wmmlr[LO_RI] = (bi[0].wmmlr[HI_LE] + bi[0].wmmlr[LO_RI]) / 2;

      Z_TRACE_IF(MSEG_TRACE_IF, "consistence weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
        bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);
    }
#endif
  }
}
#endif


void border_change_old(slint_t doweights, border_info_t *bi, border_info_t *bi_old, slint_t gcs, slint_t gc, slweight_t gws, slweight_t gw, slint_t dir) /* sl_func border_change_old */
{
  Z_TRACE_IF(MSEG_TRACE_IF, "change: gcs = %" slint_fmt ", gc = %" slint_fmt "", gcs, gc);

  bi[0].crange[0] += gcs;
  bi[0].crange[1] = bi[0].crange[0] + gc;

  Z_TRACE_IF(MSEG_TRACE_IF, "counts_range: %" slint_fmt "  %" slint_fmt "", bi[0].crange[0], bi[0].crange[1]);

#ifdef elem_weight
  if (doweights)
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "change: gws = %" slweight_fmt ", gc = %" slweight_fmt, gws, gw);

    bi[0].wrange[0] += gws;
    bi[0].wrange[1] = bi[0].wrange[0] + gw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weights_range: %" slweight_fmt "  %" slweight_fmt, bi[0].wrange[0], bi[0].wrange[1]);

  } else
#endif
  { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }

  /* forward or hit */
  if (dir >= 0)
  {
    bi[0].cmmlr[LO_LE] = z_minmax(bi[0].crange[0], bi[0].cmmlr[LO_LE], bi[0].crange[1]);
    bi[0].cmmlr[HI_LE] = z_minmax(bi[0].crange[0], bi[0].cmmlr[HI_LE], bi[0].crange[1]);

    Z_TRACE_IF(MSEG_TRACE_IF, "count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
      bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
    if (doweights)
    {
      bi[0].wmmlr[LO_LE] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[LO_LE], bi[0].wrange[1]);
      bi[0].wmmlr[HI_LE] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[HI_LE], bi[0].wrange[1]);

      Z_TRACE_IF(MSEG_TRACE_IF, "weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
        bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);

    } else
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");


#ifndef MSEG_BORDER_UPDATE_FULL
    if (bi[0].cmmlr[LO_LE] != bi_old[0].cmmlr[LO_LE] || bi[0].cmmlr[HI_LE] != bi_old[0].cmmlr[HI_LE]
# ifdef elem_weight
     || (doweights && (bi[0].wmmlr[LO_LE] != bi_old[0].wmmlr[LO_LE] || bi[0].wmmlr[HI_LE] != bi_old[0].wmmlr[HI_LE]))
# endif
     ) bi[1].update = 1;
#endif
  }
  
  /* backward or hit */
  if (dir <= 0)
  {
    bi[0].cmmlr[HI_RI] = z_minmax(bi[0].crange[0], bi[0].cmmlr[HI_RI], bi[0].crange[1]);
    bi[0].cmmlr[LO_RI] = z_minmax(bi[0].crange[0], bi[0].cmmlr[LO_RI], bi[0].crange[1]);

    Z_TRACE_IF(MSEG_TRACE_IF, "count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
      bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
    if (doweights)
    {
      bi[0].wmmlr[HI_RI] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[HI_RI], bi[0].wrange[1]);
      bi[0].wmmlr[LO_RI] = z_minmax(bi[0].wrange[0], bi[0].wmmlr[LO_RI], bi[0].wrange[1]);

      Z_TRACE_IF(MSEG_TRACE_IF, "weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
        bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);

    } else
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");

#ifndef MSEG_BORDER_UPDATE_FULL
    if (bi[0].cmmlr[HI_RI] != bi_old[0].cmmlr[HI_RI] || bi[0].cmmlr[LO_RI] != bi_old[0].cmmlr[LO_RI]
# ifdef elem_weight
     || (doweights && (bi[0].wmmlr[HI_RI] != bi_old[0].wmmlr[HI_RI] || bi[0].wmmlr[LO_RI] != bi_old[0].wmmlr[LO_RI]))
# endif
     ) bi[-1].update = 1;
#endif
  }

#ifndef MSEG_BORDER_UPDATE_FULL
  bi[0].update = 0;
#endif

  Z_TRACE_IF(MSEG_TRACE_IF, "range diff 0: %" slint_fmt "-%" slint_fmt " | %" slint_fmt "-%" slint_fmt,
    bi[0].crange[0] - bi[-1].crange[1], bi[0].crange[0] - bi[-1].crange[0],
    bi[1].crange[0] - bi[ 0].crange[0], bi[1].crange[1] - bi[ 0].crange[0]);
  Z_TRACE_IF(MSEG_TRACE_IF, "range diff 1: %" slint_fmt "-%" slint_fmt " | %" slint_fmt "-%" slint_fmt,
    bi[0].crange[1] - bi[-1].crange[1], bi[0].crange[1] - bi[-1].crange[0],
    bi[1].crange[0] - bi[ 0].crange[1], bi[1].crange[1] - bi[ 0].crange[1]);
}


void border_init_old(slint_t doweights, border_info_t *bi, partcond_intern_t *pc, slint_t tc, slweight_t tw) /* sl_func border_init_old */
{
#ifdef MSEG_BORDER_UPDATE_FULL
  bi[0].done = (pc == NULL);
#endif
  bi[0].update = (pc != NULL);

  bi[0].crange[0] = 0;
  bi[0].crange[1] = tc;

  Z_TRACE_IF(MSEG_TRACE_IF, "count range: %" slint_fmt " - %" slint_fmt "",
    bi[0].crange[0], bi[0].crange[1]);

  bi[0].cmmlr[LO_LE] = bi[0].cmmlr[HI_LE] = bi[0].cmmlr[HI_RI] = bi[0].cmmlr[LO_RI] = -1;

#ifdef elem_weight
  if (doweights)
  {
    bi[0].wrange[0] = 0.0;
    bi[0].wrange[1] = tw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight range: %" slweight_fmt " - %" slweight_fmt,
      bi[0].wrange[0], bi[0].wrange[1]);

    bi[0].wmmlr[LO_LE] = bi[0].wmmlr[HI_LE] = bi[0].wmmlr[HI_RI] = bi[0].wmmlr[LO_RI] = -1.0;

  } else
#endif
    Z_TRACE_IF(MSEG_TRACE_IF, "");

  if (!pc)
  {
    bi[0].cmmlr[LO_LE] = bi[0].cmmlr[HI_LE] = 0;
    bi[0].cmmlr[HI_RI] = bi[0].cmmlr[LO_RI] = tc;

#ifdef elem_weight
    if (doweights)
    {
      bi[0].wmmlr[LO_LE] = bi[0].wmmlr[HI_LE] = 0.0;
      bi[0].wmmlr[HI_RI] = bi[0].wmmlr[LO_RI] = tw;
    }
#endif
  }

/*  Z_TRACE_IF(MSEG_TRACE_IF, "count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
    bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

#ifdef elem_weight
  if (doweights)
    Z_TRACE_IF(MSEG_TRACE_IF, "weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
      bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);
#endif*/
}


inline void border_currents_old(slint_t doweights, border_info_t *bi, slweight_t *currents) /* sl_func border_currents_old */
{
  Z_TRACE_IF(MSEG_TRACE_IF, "count[min/max-left/right]: %" slint_fmt " / %" slint_fmt " - %" slint_fmt " / %" slint_fmt "",
    bi[0].cmmlr[LO_LE], bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI], bi[0].cmmlr[LO_RI]);

  Z_TRACE_IF(MSEG_TRACE_IF, "crange: %" slint_fmt " - %" slint_fmt "",
    bi[0].crange[0], bi[0].crange[1]);

  /* select max. of low and min. of high */
  currents[CNT_LO] = z_max(bi[0].cmmlr[LO_LE], bi[0].cmmlr[LO_RI]) - bi[0].crange[0];
  currents[CNT_HI] = z_min(bi[0].cmmlr[HI_LE], bi[0].cmmlr[HI_RI]) - bi[0].crange[0];

  Z_TRACE_IF(MSEG_TRACE_IF, "currents count: %" slweight_fmt " - %" slweight_fmt, currents[CNT_LO], currents[CNT_HI]);

#ifdef elem_weight
  if (doweights)
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "weight[min/max-left/right]: %" slweight_fmt " / %" slweight_fmt " - %" slweight_fmt " / %" slweight_fmt,
      bi[0].wmmlr[LO_LE], bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI], bi[0].wmmlr[LO_RI]);

    Z_TRACE_IF(MSEG_TRACE_IF, "wrange: %" slweight_fmt " - %" slweight_fmt, bi[0].wrange[0], bi[0].wrange[1]);

    /* select highest min and lowest max */
    currents[WHT_LO] = z_max(bi[0].wmmlr[LO_LE], bi[0].wmmlr[LO_RI]) - bi[0].wrange[0];
    currents[WHT_HI] = z_min(bi[0].wmmlr[HI_LE], bi[0].wmmlr[HI_RI]) - bi[0].wrange[0];

    Z_TRACE_IF(MSEG_TRACE_IF, "currents weight: %" slweight_fmt " - %" slweight_fmt, currents[WHT_LO], currents[WHT_HI]);

  } else
#endif
  { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }
}

#else /* OLD */

void border_init(slint_t doweights, border_info_t *bi, slint_t current, slint_t tc, slweight_t tw) /* sl_func border_init */
{
  bi[0].crange[0] = 0;
  bi[0].crange[1] = tc;

  Z_TRACE_IF(MSEG_TRACE_IF, "count range: %" slint_fmt " - %" slint_fmt "",
    bi[0].crange[0], bi[0].crange[1]);

  bi[0].ccurrent[LO] = (current == 0 || current != 1)?0:tc;
  bi[0].ccurrent[HI] = (current == 1 || current != 0)?tc:0;

  Z_TRACE_IF(MSEG_TRACE_IF, "count[low/high]: %" slint_fmt " / %" slint_fmt,
    bi[0].ccurrent[LO], bi[0].ccurrent[HI]);

#ifdef elem_weight
  if (doweights)
  {
    bi[0].wrange[0] = 0.0;
    bi[0].wrange[1] = tw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight range: %" slweight_fmt " - %" slweight_fmt,
      bi[0].wrange[0], bi[0].wrange[1]);

    bi[0].wcurrent[LO] = (current == 0 || current != 1)?0:tw;
    bi[0].wcurrent[HI] = (current == 1 || current != 0)?tw:0;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight[low/high]: %" slweight_fmt " / %" slweight_fmt,
      bi[0].wcurrent[LO], bi[0].wcurrent[HI]);

  } else
#endif
    Z_TRACE_IF(MSEG_TRACE_IF, "");
}


void border_update(slint_t doweights, border_info_t *bi, partcond_intern_t *pc, slint_t dir, slint_t reduction) /* sl_func border_update */
{
  slint_t ccurrent_new[2];
#ifdef elem_weight
  slweight_t wcurrent_new[2];
#endif

#ifdef MSEG_BORDER_UPDATE_REDUCTION
  slint_t count_reduction;
# ifdef elem_weight
  slweight_t weight_reduction;
# endif
#endif

  /* init from range */
  ccurrent_new[LO] = bi[0].crange[0];
  ccurrent_new[HI] = bi[0].crange[1];
#ifdef elem_weight
  if (doweights)
  {
    wcurrent_new[LO] = bi[0].wrange[0];
    wcurrent_new[HI] = bi[0].wrange[1];
  }
#endif

  /* forward */
  if (dir > 0)
  {
    /* init from min/max */
    if (pc[0].pcm & SLPC_COUNTS_MM)
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction)
      {
        count_reduction = z_round((bi[-1].ccurrent[HI] - bi[-1].ccurrent[LO]) * 0.5 * mseg_border_update_count_reduction);
        Z_ASSERT(count_reduction >= 0);

        Z_TRACE_IF(MSEG_TRACE_IF, "forward: count_reduction: %" slint_fmt, count_reduction);

        ccurrent_new[0] = z_min(bi[-1].ccurrent[LO] + count_reduction + pc[0].count_min, bi[0].ccurrent[HI]);
        ccurrent_new[1] = z_max(bi[-1].ccurrent[HI] - count_reduction + pc[0].count_max, bi[0].ccurrent[LO]);

        Z_TRACE_IF(MSEG_TRACE_IF, "forward: count[low/high]: min(%" slint_fmt " + %" slint_fmt " + %" slint_fmt ", %" slint_fmt "), max(%" slint_fmt " - % "slint_fmt " + %" slint_fmt ", %" slint_fmt ")",
          bi[-1].ccurrent[LO], count_reduction, pc[0].count_min, bi[0].ccurrent[HI], bi[-1].ccurrent[HI], count_reduction, pc[0].count_max, bi[0].ccurrent[LO]);

      } else
#endif
      {
        ccurrent_new[LO] = bi[-1].ccurrent[LO] + pc[0].count_min;
        ccurrent_new[HI] = bi[-1].ccurrent[HI] + pc[0].count_max;

        Z_TRACE_IF(MSEG_TRACE_IF, "forward: count[low/high]: %" slint_fmt " + %" slint_fmt ", %" slint_fmt " + %" slint_fmt,
          bi[-1].ccurrent[LO], pc[0].count_min, bi[-1].ccurrent[HI], pc[0].count_max);
      }
    }

#ifdef elem_weight
    if (doweights)
    {
      /* init from min/max */
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
# ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          weight_reduction = (bi[-1].wcurrent[HI] - bi[-1].wcurrent[LO]) * 0.5 * mseg_border_update_weight_reduction;
          Z_ASSERT(weight_reduction >= 0.0);

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: weight_reduction: %" slweight_fmt, weight_reduction);
  
          wcurrent_new[LO] = z_min(bi[-1].wcurrent[LO] + weight_reduction + pc[0].weight_min, bi[0].wcurrent[HI]);
          wcurrent_new[HI] = z_max(bi[-1].wcurrent[HI] - weight_reduction + pc[0].weight_max, bi[0].wcurrent[LO]);

          Z_TRACE_IF(MSEG_TRACE_IF, "weight_reduction: weight[low/high]: min(%" slweight_fmt " + %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt "), max(%" slweight_fmt " - %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt ")",
            bi[-1].wcurrent[LO], weight_reduction, pc[0].weight_min, bi[0].wcurrent[HI], bi[-1].wcurrent[HI], weight_reduction, pc[0].weight_max, bi[0].wcurrent[LO]);

        } else
# endif
        {
          wcurrent_new[LO] = bi[-1].wcurrent[LO] + pc[0].weight_min;
          wcurrent_new[HI] = bi[-1].wcurrent[HI] + pc[0].weight_max;

          Z_TRACE_IF(MSEG_TRACE_IF, "forward: weight[low/high]: %" slweight_fmt " + %" slweight_fmt ", %" slweight_fmt " + %" slweight_fmt,
            bi[-1].wcurrent[LO], pc[0].weight_min, bi[-1].wcurrent[HI], pc[0].weight_max);
        }
      }

    } else
#endif
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction) Z_TRACE_IF(MSEG_TRACE_IF, "");
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");
    }

  } else /* backward */
  {
    /* init from min/max */
    if (pc[0].pcm & SLPC_COUNTS_MM)
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction)
      {
        count_reduction = z_round((bi[1].ccurrent[HI] - bi[1].ccurrent[LO]) * 0.5 * mseg_border_update_count_reduction);
        Z_ASSERT(count_reduction >= 0);

        Z_TRACE_IF(MSEG_TRACE_IF, "backward: count_reduction: %" slint_fmt, count_reduction);

        ccurrent_new[LO] = z_min(bi[1].ccurrent[LO] + count_reduction - pc[1].count_max, bi[0].ccurrent[HI]);
        ccurrent_new[HI] = z_max(bi[1].ccurrent[HI] - count_reduction - pc[1].count_min, bi[0].ccurrent[LO]);

        Z_TRACE_IF(MSEG_TRACE_IF, "backward: count[low/high]: min(%" slint_fmt " + %" slint_fmt " - %" slint_fmt ", %" slint_fmt "), max(%" slint_fmt " - %" slint_fmt " + %" slint_fmt ", %" slint_fmt ")",
          bi[1].ccurrent[LO], count_reduction, pc[1].count_max, bi[0].ccurrent[HI], bi[1].ccurrent[HI], count_reduction, pc[1].count_min, bi[0].ccurrent[LO]);

      } else
#endif
      {
        ccurrent_new[LO] = bi[1].ccurrent[LO] - pc[1].count_max;
        ccurrent_new[HI] = bi[1].ccurrent[HI] - pc[1].count_min;

        Z_TRACE_IF(MSEG_TRACE_IF, "backward: count[low/high]: %" slint_fmt " - %" slint_fmt ", %" slint_fmt " - %" slint_fmt "",
          bi[1].ccurrent[LO], pc[1].count_max, bi[1].ccurrent[HI], pc[1].count_min);
      }
    }

#ifdef elem_weight
    if (doweights)
    {
      /* init from min/max */
      if (pc[0].pcm & SLPC_WEIGHTS_MM)
      {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
        if (reduction)
        {
          weight_reduction = (bi[1].wcurrent[1] - bi[1].wcurrent[LO]) * 0.5 * mseg_border_update_weight_reduction;
          Z_ASSERT(weight_reduction >= 0.0);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight_reduction: %" slweight_fmt, weight_reduction);

          wcurrent_new[LO] = z_min(bi[1].wcurrent[LO] + weight_reduction - pc[1].weight_max, bi[0].wcurrent[HI]);
          wcurrent_new[HI] = z_max(bi[1].wcurrent[HI] - weight_reduction - pc[1].weight_min, bi[0].wcurrent[LO]);

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight[low/high]: min(%" slweight_fmt " + %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt "), max(%" slweight_fmt " - %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt ")",
            bi[1].wcurrent[LO], weight_reduction, pc[1].weight_max, bi[0].wcurrent[HI], bi[1].wcurrent[HI], weight_reduction, pc[1].weight_min, bi[0].wcurrent[LO]);

        } else
#endif
        {
          wcurrent_new[LO] = bi[1].wcurrent[LO] - pc[1].weight_max;
          wcurrent_new[HI] = bi[1].wcurrent[HI] - pc[1].weight_min;

          Z_TRACE_IF(MSEG_TRACE_IF, "backward: weight[low/high]: %" slweight_fmt " - %" slweight_fmt ", %" slweight_fmt " - %" slweight_fmt,
            bi[1].wcurrent[LO], pc[1].weight_max, bi[1].wcurrent[HI], pc[1].weight_min);
        }
      }

    } else
#endif
    {
#ifdef MSEG_BORDER_UPDATE_REDUCTION
      if (reduction) Z_TRACE_IF(MSEG_TRACE_IF, "");
#endif
      Z_TRACE_IF(MSEG_TRACE_IF, "");
    }
  }

  /* check against low/high */
  if (pc[0].pcm & SLPC_COUNTS_LH)
  {
    if (ccurrent_new[LO] < pc[1].count_low)  ccurrent_new[LO] = pc[1].count_low;
    if (ccurrent_new[HI] > pc[0].count_high) ccurrent_new[HI] = pc[0].count_high;
  }

  /* fit to range (NOT REQUIRED!?) */
/*  ccurrent_new[LO] = z_minmax(bi[0].crange[0], ccurrent_new[LO], bi[0].crange[1]);
  ccurrent_new[HI] = z_minmax(bi[0].crange[0], ccurrent_new[HI], bi[0].crange[1]);*/

  bi[0].ccurrent[LO] = z_max(bi[0].ccurrent[LO], ccurrent_new[LO]);
  bi[0].ccurrent[HI] = z_min(bi[0].ccurrent[HI], ccurrent_new[HI]);

  Z_TRACE_IF(MSEG_TRACE_IF, "current count[low/high]: %" slint_fmt " / %" slint_fmt, bi[0].ccurrent[LO], bi[0].ccurrent[HI]);

#ifdef elem_weight
  if (doweights)
  {
    if (pc[0].pcm & SLPC_WEIGHTS_LH)
    {
      if (wcurrent_new[LO] < pc[1].weight_low)  wcurrent_new[LO] = pc[1].weight_low;
      if (wcurrent_new[HI] > pc[0].weight_high) wcurrent_new[HI] = pc[0].weight_high;
    }

    /* fit to range (NOT REQUIRED!?) */
/*    wcurrent_new[LO] = z_minmax(bi[0].wrange[0], wcurrent_new[LO], bi[0].wrange[1]);
    wcurrent_new[HI] = z_minmax(bi[0].wrange[0], wcurrent_new[HI], bi[0].wrange[1]);*/

    bi[0].wcurrent[LO] = z_max(bi[0].wcurrent[LO], wcurrent_new[LO]);
    bi[0].wcurrent[HI] = z_min(bi[0].wcurrent[HI], wcurrent_new[HI]);

    Z_TRACE_IF(MSEG_TRACE_IF, "current weight[low/high]: %" slweight_fmt " / %" slweight_fmt, bi[0].wcurrent[LO], bi[0].wcurrent[HI]);
  }
#endif
}


void border_change(slint_t doweights, border_info_t *bi, slint_t gcs, slint_t gc, slweight_t gws, slweight_t gw) /* sl_func border_change */
{
  Z_TRACE_IF(MSEG_TRACE_IF, "change: gcs = %" slint_fmt ", gc = %" slint_fmt "", gcs, gc);

  bi[0].crange[0] += gcs;
  bi[0].crange[1] = bi[0].crange[0] + gc;

  Z_TRACE_IF(MSEG_TRACE_IF, "counts_range: %" slint_fmt "  %" slint_fmt "", bi[0].crange[0], bi[0].crange[1]);

  bi[0].ccurrent[LO] = z_minmax(bi[0].crange[0], bi[0].ccurrent[LO], bi[0].crange[1]);
  bi[0].ccurrent[HI] = z_minmax(bi[0].crange[0], bi[0].ccurrent[HI], bi[0].crange[1]);

  Z_TRACE_IF(MSEG_TRACE_IF, "count[low/high]: %" slint_fmt " / %" slint_fmt,
    bi[0].ccurrent[LO], bi[0].ccurrent[HI]);

#ifdef elem_weight
  if (doweights)
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "change: gws = %" slweight_fmt ", gc = %" slweight_fmt, gws, gw);

    bi[0].wrange[0] += gws;
    bi[0].wrange[1] = bi[0].wrange[0] + gw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weights_range: %" slweight_fmt "  %" slweight_fmt, bi[0].wrange[0], bi[0].wrange[1]);

    bi[0].wcurrent[LO] = z_minmax(bi[0].wrange[0], bi[0].wcurrent[LO], bi[0].wrange[1]);
    bi[0].wcurrent[HI] = z_minmax(bi[0].wrange[0], bi[0].wcurrent[HI], bi[0].wrange[1]);

    Z_TRACE_IF(MSEG_TRACE_IF, "weight[low/high]: %" slweight_fmt " / %" slweight_fmt,
      bi[0].wcurrent[LO], bi[0].wcurrent[HI]);

  } else
#endif
  { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }

#if 0
  Z_TRACE_IF(MSEG_TRACE_IF, "range diff 0: %" slint_fmt "-%" slint_fmt " | %" slint_fmt "-%" slint_fmt,
    bi[0].crange[0] - bi[-1].crange[1], bi[0].crange[0] - bi[-1].crange[0],
    bi[1].crange[0] - bi[ 0].crange[0], bi[1].crange[1] - bi[ 0].crange[0]);
  Z_TRACE_IF(MSEG_TRACE_IF, "range diff 1: %" slint_fmt "-%" slint_fmt " | %" slint_fmt "-%" slint_fmt,
    bi[0].crange[1] - bi[-1].crange[1], bi[0].crange[1] - bi[-1].crange[0],
    bi[1].crange[0] - bi[ 0].crange[1], bi[1].crange[1] - bi[ 0].crange[1]);
#endif
}

#endif


slint_t mpi_select_exact_generic_bulk(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic_bulk */
{
  const slint_t max_nborders = nparts - 1;
  slint_t border_lo, border_hi, nborders_removed;
  slint_t borders[max_nborders], border_bins[max_nborders];

  border_info_t border_infos_[1 + max_nborders + 1], *border_infos = border_infos_ + 1;
#ifdef OLD
  border_info_t  border_info_old;
#endif

  slint_t total_counts;
#ifdef elem_weight
  slweight_t total_weights;
#endif

  partcond_intern_t pci[nparts];

  slweight_t currents[2 * WEIGHT_FACTOR];

  slweight_t final_locals[WEIGHT_FACTOR], final_globals[WEIGHT_FACTOR];

  slint_t round, direction, refine, finalize;

  slint_t nothing, lc_min, lc_max, lcw2gcw;

  slweight_t mcw, dcw, lcw[WEIGHT_FACTOR], gcw[WEIGHT_FACTOR];

  slint_t gc, gcs;
#ifdef elem_weight
  slweight_t gw, gws;
#endif

  slint_t i, j, k, ix;

#ifdef elem_weight
  slint_t doweights, weight_factor;
#else
# define doweights  0
#endif

#ifdef VERIFY
  slint_t v;
#endif

  global_bins_t gb;


  Z_TRACE_IF(MSEG_TRACE_IF, "starting mpi_select_exact_generic");

  /* sl_tid rti_tid_mpi_select_exact_generic rti_tid_mpi_select_exact_generic_sync_init rti_tid_mpi_select_exact_generic_sync_exit */

  rti_treset(rti_tid_mpi_select_exact_generic_while);                    /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check);              /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_bins);         /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_bins_local);   /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_bins_global);  /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_round1);       /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_pre);          /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_part);         /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_part_root);    /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_final);        /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_final_root);   /* sl_tid */
  rti_treset(rti_tid_mpi_select_exact_generic_while_check_post);         /* sl_tid */

  rti_tstart(rti_tid_mpi_select_exact_generic);

  rti_tstart(rti_tid_mpi_select_exact_generic_sync_init);
#ifdef SYNC_ON_INIT
  MPI_Barrier(comm);
#endif
  rti_tstop(rti_tid_mpi_select_exact_generic_sync_init);

#ifdef VERIFY
  v = elements_validate_order(s, 1);
  
  Z_TRACE_IF(MSEG_TRACE_IF, "elements order: %s (%" slint_fmt ")", (v > 0)?"FAILED":"SUCCESS", v);
#endif

#ifdef elem_weight
  doweights = ((pconds->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

#ifdef elem_weight
  weight_factor = 1 + (doweights != 0);
# define MY_WEIGHT_FACTOR  weight_factor
#else
# define MY_WEIGHT_FACTOR  1
#endif

  mpi_binning_create(&gb, max_nborders, mseg_binnings, s, nelements, doweights, bm, size, rank, comm);

  /* init parts */
  border_lo = 0;
  border_hi = max_nborders - 1;
  for (i = border_lo; i <= border_hi; ++i)
  {
    borders[i] = i;
    border_bins[i] = 0;
  }

  /* reset splitter */
  sp->n = nparts * nelements;
  splitter_reset(sp);

#ifdef MSEG_INFO
  mseg_info_finish_rounds_avg = 0;
#endif

  rti_tstart(rti_tid_mpi_select_exact_generic_while);

  direction = 1;

  round = 0;
  while (border_lo <= border_hi)
  {
    ++round;

    Z_TRACE_IF(MSEG_TRACE_IF, "ROUND: %" slint_fmt ", %s, %" slint_fmt " border(s)", round, (direction > 0)?"forward":"backward", border_hi - border_lo + 1);

    nborders_removed = 0;

    mpi_binning_pre(&gb, size, rank, comm);

    Z_TRACE_IF(MSEG_TRACE_IF, "ROUND: %" slint_fmt ", bm_nbins: %" slint_fmt, round, gb.bm->nbins);

    finalize = (gb.bm->nbins <= 1);

    rti_tstart(rti_tid_mpi_select_exact_generic_while_check);

    i = (direction > 0)?border_lo:border_hi;
    while ((direction > 0)?(i <= border_hi):(i >= border_lo))
    {
      Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ": PART: %" slint_fmt ",%" slint_fmt ": %s", round, i, borders[i], ((direction > 0)?"forward":"backward"));

      ix = i;

      rti_tstart(rti_tid_mpi_select_exact_generic_while_check_bins);

      if (!finalize || (finalize && round == 1 && i == border_lo))
      {
        mpi_binning_exec_reset(&gb, size, rank, comm);

        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_bins_local);

        while ((direction > 0)?(ix <= border_hi):(ix >= border_lo))
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": bin %" slint_fmt, ix, borders[ix], border_bins[ix]);

          if (mpi_binning_exec_local(&gb, border_bins[ix], size, rank, comm) < 0)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "break");
            break;
          }

          ix += direction;
        }

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_bins_local);

        Z_TRACE_IF(MSEG_TRACE_IF, "global %d (i = %" slint_fmt ", ix = %" slint_fmt ")", abs(ix - i), i, ix);

        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_bins_global);

        mpi_binning_exec_global(&gb,
#ifdef MSEG_ROOT
          mseg_root,
#else
          -1,
#endif
          size, rank, comm);

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_bins_global);

      } else ix += direction;

      rti_tstop(rti_tid_mpi_select_exact_generic_while_check_bins);

      if (round == 1 && i == border_lo)
      {
        /* do initialization */
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_round1);

#ifdef MSEG_ROOT
        if (mseg_root < 0 || mseg_root == rank)
#endif
        {
          total_counts = 0;
#ifdef elem_weight
          total_weights = 0.0;
#endif
          for (j = 0; j < gb.bm->nbins; ++j)
          {
            total_counts += *gb_counts(&gb, border_bins[i], j);
#ifdef elem_weight
            if (doweights)
              total_weights += *gb_weights(&gb, border_bins[i], j);
#endif
          }

          Z_TRACE_IF(MSEG_TRACE_IF, "total_counts = %" slint_fmt, total_counts);
#ifdef elem_weight
          if (doweights)
            Z_TRACE_IF(MSEG_TRACE_IF, "total_weights = %" slweight_fmt , total_weights);
          else
#endif
            Z_TRACE_IF(MSEG_TRACE_IF, "");

          init_partconds_intern(nparts, pci, pconds, nparts, total_counts, elem_weight_ifelse(total_weights, 0));

          /* init lowest and highest part (sentinels) */
#ifdef OLD
          Z_TRACE_IF(MSEG_TRACE_IF, "init lowest border:");
          border_init_old(doweights, &border_infos[border_lo - 1], NULL, total_counts, elem_weight_ifelse(total_weights, 0));
          Z_TRACE_IF(MSEG_TRACE_IF, "init highest border:");
          border_init_old(doweights, &border_infos[border_hi + 1], NULL, total_counts, elem_weight_ifelse(total_weights, 0));
#else
          Z_TRACE_IF(MSEG_TRACE_IF, "init lowest border:");
          border_init(doweights, &border_infos[border_lo - 1], 0, total_counts, elem_weight_ifelse(total_weights, 0));
          Z_TRACE_IF(MSEG_TRACE_IF, "init highest border:");
          border_init(doweights, &border_infos[border_hi + 1], 1, total_counts, elem_weight_ifelse(total_weights, 0));
#endif

#ifdef MSEG_BORDER_UPDATE_REDUCTION
          /* init+update forwards */
          for (j = border_lo; j <= border_hi; ++j)
          {
#ifdef OLD
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init_old(doweights, &border_infos[borders[j]], &pci[borders[j]], total_counts, elem_weight_ifelse(total_weights, 0));

            Z_TRACE_IF(MSEG_TRACE_IF, "update update %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update_old(doweights, &border_infos[borders[j]], &pci[borders[j]], 1, 0, 0);
#else
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init(doweights, &border_infos[borders[j]], -1, total_counts, elem_weight_ifelse(total_weights, 0));

            Z_TRACE_IF(MSEG_TRACE_IF, "update update %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update(doweights, &border_infos[borders[j]], &pci[borders[j]], 1, 0);
#endif
          }
#endif

          /* [init+]update backwards */
          for (j = border_hi; j >= border_lo; --j)
          {
#ifdef OLD
#ifndef MSEG_BORDER_UPDATE_REDUCTION
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init_old(doweights, &border_infos[borders[j]], &pci[borders[j]], total_counts, elem_weight_ifelse(total_weights, 0));
#endif
            Z_TRACE_IF(MSEG_TRACE_IF, "update border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update_old(doweights, &border_infos[borders[j]], &pci[borders[j]], -1, 0, 1);
#else
#ifndef MSEG_BORDER_UPDATE_REDUCTION
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init(doweights, &border_infos[borders[j]], -1, total_counts, elem_weight_ifelse(total_weights, 0));
#endif
            Z_TRACE_IF(MSEG_TRACE_IF, "update border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update(doweights, &border_infos[borders[j]], &pci[borders[j]], -1, 0);
#endif
          }
        }
        
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_round1);
      }

do_partitioning:

      rti_tstart(rti_tid_mpi_select_exact_generic_while_check_pre);

#ifdef MSEG_ROOT
      if (mseg_root < 0 || mseg_root == rank)
#endif
      {
        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": check", i, borders[i]);

        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": crange: %" slint_fmt " - %" slint_fmt, i, borders[i], border_infos[borders[i]].crange[0], border_infos[borders[i]].crange[1]);
#ifdef elem_weight
        if (doweights)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": wrange: %" slweight_fmt " - %" slweight_fmt, i, borders[i], border_infos[borders[i]].wrange[0], border_infos[borders[i]].wrange[1]);
#endif

#ifdef OLD
        /* save old limits */
        border_info_old = border_infos[borders[i]];

#ifdef MSEG_BORDER_UPDATE_FULL
        if (mseg_border_update_full)
          border_update_full_old(doweights, &border_infos[borders[i]], &pci[borders[i]], 1);
        else
#endif
        {
          /* is an update required? */
          if (border_infos[borders[i]].update) border_update_old(doweights, &border_infos[borders[i]], &pci[borders[i]], direction, 1, 1);
          else { Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); Z_TRACE_IF(MSEG_TRACE_IF, ""); }
        }

        /* get currents */
        border_currents_old(doweights, &border_infos[borders[i]], currents);
#else
        border_update(doweights, &border_infos[borders[i]], &pci[borders[i]], direction, 1);
        currents[CNT_LO] = border_infos[borders[i]].ccurrent[LO] - border_infos[borders[i]].crange[0];
        currents[CNT_HI] = border_infos[borders[i]].ccurrent[HI] - border_infos[borders[i]].crange[0];
#ifdef elem_weight
        if (doweights)
        {
          currents[WHT_LO] = border_infos[borders[i]].wcurrent[LO] - border_infos[borders[i]].wrange[0];
          currents[WHT_HI] = border_infos[borders[i]].wcurrent[HI] - border_infos[borders[i]].wrange[0];
        }
#endif
#endif

        Z_ASSERT_IF(MSEG_ASSERT_IF, currents[CNT_LO] <= currents[CNT_HI]);
        Z_ASSERT_IF(MSEG_ASSERT_IF, 0 <= currents[CNT_LO]);
        
        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": currents count: %" slweight_fmt " - %" slweight_fmt " (range: %" slweight_fmt ")",
          i, borders[i], currents[CNT_LO], currents[CNT_HI], currents[CNT_HI] - currents[CNT_LO]);

#ifdef elem_weight
        if (doweights)
        {
          Z_ASSERT_IF(MSEG_ASSERT_IF, currents[WHT_LO] <= currents[WHT_HI]);

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": currents weight: %" slweight_fmt " - %" slweight_fmt " (range: %" slweight_fmt ")",
            i, borders[i], currents[WHT_LO], currents[WHT_HI], currents[WHT_HI] - currents[WHT_LO]);

        } else
#endif
          Z_TRACE_IF(MSEG_TRACE_IF, "");
      }

      rti_tstop(rti_tid_mpi_select_exact_generic_while_check_pre);

      refine = 0;

      if (!finalize)
      {
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_part);

#ifdef MSEG_ROOT
        if (mseg_root < 0 || mseg_root == rank)
#endif
        {
          gcs = 0;
#ifdef elem_weight
          gws = 0.0;
#endif

          for (k = 0; k < gb.bm->nbins; ++k)
          {
            gc = *gb_counts(&gb, border_bins[i], k);

            currents[CNT_LO] -= gc;
            currents[CNT_HI] -= gc;

            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": k = %" slint_fmt ", currents count: %" slweight_fmt " - %" slweight_fmt ", gc = %" slint_fmt ", gcs = %" slint_fmt,
              i, borders[i], k, currents[CNT_LO], currents[CNT_HI], gc, gcs);

#ifdef elem_weight
            if (doweights)
            {
              gw = *gb_weights(&gb, border_bins[i], k);

              currents[WHT_LO] -= gw;
              currents[WHT_HI] -= gw;

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": k = %" slint_fmt ", currents weight: %" slweight_fmt " - %" slweight_fmt ", gw = %" slweight_fmt ", gws = %" slweight_fmt,
                i, borders[i], k, currents[WHT_LO], currents[WHT_HI], gw, gws);

            } else
#endif
            {
#ifdef elem_weight
              gw = 0;
#endif
              Z_TRACE_IF(MSEG_TRACE_IF, "");
            }

            /* stop and refine if max count is skipped OR (min count AND max weight is skipped) */
            if ((currents[CNT_HI] < 0)
#ifdef elem_weight
              || (doweights && currents[CNT_LO] < 0 && currents[WHT_HI] < 0.0)
#endif
              )
            {
              /* stop for REFINE only if there are more than one elements to refine */
              if (gc > 1)
              {
                refine = 1;
                break;
              }
            }

            gcs += gc;
            gc = 0;

#ifdef elem_weight
            gws += gw;
            gw = 0.0;
#endif

            /* if between min/max counts */
            if (currents[CNT_LO] <= 0 && currents[CNT_HI] >= 0)
            {
#ifdef elem_weight
              if (doweights)
              {
                Z_TRACE_IF(MSEG_TRACE_IF, "go to next: %d && %d", (currents[CNT_HI] > 0), (currents[WHT_LO] > 0));

                /* go to next if max count not reached AND min weight not reached */
                if (currents[CNT_HI] > 0 && currents[WHT_LO] > 0) continue;
              }
#endif

              /* look ahead for a better stop */
              if (k + 1 < gb.bm->nbins && currents[CNT_HI] - *gb_counts(&gb, border_bins[i], k + 1) >= 0)
              {
#ifdef elem_weight
                if (doweights)
                {
                  /* continue if weights will improve */
                  if (z_abs(currents[WHT_LO] + currents[WHT_HI]) > z_abs(currents[WHT_LO] + currents[WHT_HI] - 2 * *gb_weights(&gb, border_bins[i], k + 1))) continue;

                } else
#endif
                {
                  /* continue if counts will improve */
                  if (z_abs(currents[CNT_LO] + currents[CNT_HI]) > z_abs(currents[CNT_LO] + currents[CNT_HI] - 2 * *gb_counts(&gb, border_bins[i], k + 1))) continue;
                }
              }

              /* stop */
              break;
            }
          }

          Z_ASSERT_IF(MSEG_ASSERT_IF, k < gb.bm->nbins);

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s k = %" slint_fmt, i, borders[i], (refine)?"REFINE":"HIT", k);

          /* make sure k is safe (it is used as index later) */
          if (k >= gb.bm->nbins) k = gb.bm->nbins - 1;

          k = (k + 1) * ((refine)?-1:1);
        }

#ifdef MSEG_ROOT
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_part_root);
        if (mseg_root >= 0) MPI_Bcast(&k, 1, int_mpi_datatype, mseg_root, comm);
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_part_root);
#endif

        refine = (k < 0);
        if (k < 0) k = -k;
        --k;

        /* refine or hit */
        if (refine)
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": refine bin: %" slint_fmt " @ k = %" slint_fmt, i, borders[i], border_bins[i], k);

          border_bins[i] = mpi_binning_refine(&gb, border_bins[i], k, sp, borders[i] + 1, size, rank, comm);
          
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": new bin: %" slint_fmt, i, borders[i], border_bins[i]);

        } else
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": hit bin: %" slint_fmt " @ k = %" slint_fmt, i, borders[i], border_bins[i], k);

          mpi_binning_hit(&gb, border_bins[i], k, sp, borders[i] + 1, size, rank, comm);
        }

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_part);

      } else
      {
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_final);

#ifdef MSEG_ROOT
        rti_tstart(rti_tid_mpi_select_exact_generic_while_check_final_root);
        if (mseg_root >= 0) MPI_Bcast(currents, 2 * MY_WEIGHT_FACTOR, weight_mpi_datatype, mseg_root, comm);
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_final_root);
#endif

        switch (mseg_finalize_mode)
        {
          case SL_MSEG_FM_ALLORNOTHING:
            Z_TRACE_IF(MSEG_TRACE_IF, "finalize mode: all or nothing");
#ifdef elem_weight
            if (doweights)
            {
              nothing = (currents[WHT_LO] < ((border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]) - currents[WHT_HI]));
              Z_TRACE_IF(MSEG_TRACE_IF, "weight: %" slweight_fmt " vs. %" slweight_fmt " -> %s", currents[WHT_LO], (border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]) - currents[WHT_HI], ((nothing)?"NOTHING":"ALL"));

            } else
#endif
            {
              nothing = (currents[CNT_LO] < ((border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]) - currents[CNT_HI]));
              Z_TRACE_IF(MSEG_TRACE_IF, "count: %" slint_fmt " vs. %" slint_fmt " -> %s", (slint_t) currents[CNT_LO], (slint_t) ((border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]) - currents[CNT_HI]), ((nothing)?"NOTHING":"ALL"));
            }

            if (nothing)
            {
              mcw = dcw = 0.0;
              lc_min = lc_max = 0;

            } else
            {
#ifdef elem_weight
              if (doweights)
              {
                mcw = dcw = border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0];

              } else
#endif
              {
                mcw = dcw = border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0];
              }
              lc_min = lc_max = border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0];
            }

            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s: mcw = %" slweight_fmt ", dcw = %" slweight_fmt, i, borders[i], ((doweights)?"weights":"counts"), mcw, dcw);

            gcw[0] = (nothing)?0.0:(border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]);
#ifdef elem_weight
            if (doweights)
              gcw[1] = (nothing)?0.0:(border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]);
#endif
            lcw2gcw = 0;
            break;
          
          case SL_MSEG_FM_MIDDLE:
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": finalize mode: middle (CURRENTLY MISSING!)", i, borders[i]);
            break;
          
          case SL_MSEG_FM_EXACT:
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": finalize mode: exact", i, borders[i]);

            final_locals[0] =
#ifdef elem_weight
              final_locals[1] =
#endif
              0.0;

            for (j = 0; j < nelements; ++j)
            {
              final_locals[0] += lb_bin_count(&gb.lb, border_bins[i], j);
#ifdef elem_weight
              final_locals[1] += lb_bin_weight(&gb.lb, border_bins[i], j);
#endif
            }

            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_locals[0]: %" slweight_fmt, i, borders[i], final_locals[0]);
#ifdef elem_weight
            if (doweights)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_locals[1]: %" slweight_fmt, i, borders[i], final_locals[1]);
#endif

            MPI_Exscan(final_locals, final_globals, MY_WEIGHT_FACTOR, weight_mpi_datatype, MPI_SUM, comm);
            if (rank == 0) final_globals[0] =
#ifdef elem_weight
              final_globals[1] =
#endif
              0.0;

            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_globals[0]: %" slweight_fmt, i, borders[i], final_globals[0]);
#ifdef elem_weight
            if (doweights)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_globals[1]: %" slweight_fmt, i, borders[i], final_globals[1]);
#endif

#ifdef elem_weight
            if (doweights)
            {
              /* middle of min/max weight */
              mcw = (currents[WHT_LO] + currents[WHT_HI]) / 2.0;

              /* min. part of weight to contribute */
              dcw = z_max(0, mcw - final_globals[1]);

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": weights: mcw = %" slweight_fmt ", dcw = %" slweight_fmt, i, borders[i], mcw, dcw);

            } else
#endif
            {
              /* middle of min/max count */
              mcw = (currents[CNT_LO] + currents[CNT_HI]) / 2;

              /* min. part of count to contribute */
              dcw = z_max(0, mcw - final_globals[0]);

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": counts: mcw = %" slweight_fmt ", dcw = %" slweight_fmt, i, borders[i], mcw, dcw);
            }

            lc_min = currents[CNT_LO] - final_globals[0];
            lc_max = currents[CNT_HI] - final_globals[0];

            lcw2gcw = 1;
            break;
        }

        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lc_min = %" slint_fmt ", lc_max = %" slint_fmt, i, borders[i], lc_min, lc_max);

        mpi_binning_finalize(&gb, border_bins[i], dcw, lc_min, lc_max, lcw, sp, borders[i] + 1, size, rank, comm);
        
        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lcw[0] = %" slweight_fmt, i, borders[i], lcw[0]);
#ifdef elem_weight
        if (doweights)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lcw[1] = %" slweight_fmt, i, borders[i], lcw[1]);
#endif

        gcs = gc = 0;
#ifdef elem_weight
        gws = gw = 0.0;
#endif

        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": next border: %" slint_fmt " <= %" slint_fmt " + %" slint_fmt " <= %" slint_fmt,
          i, borders[i], border_lo, i, direction, border_hi);

        /* if the next open border is really the _next_ border */
        if (border_lo <= i + direction && i + direction <= border_hi && borders[i + direction] == borders[i] + direction)
        {
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": next border: %" slint_fmt " == %" slint_fmt " + %" slint_fmt,
            i, borders[i], borders[i + direction], borders[i], direction);

#ifdef elem_weight
          if (doweights)
          {
            if (lcw2gcw)
            {
              /* need to determine the exact global counts/weights from the local counts/weights */
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": get gcw from lcw", i, borders[i]);
              
# ifdef MSEG_ROOT
              if (mseg_root >= 0) MPI_Reduce(lcw, gcw, MY_WEIGHT_FACTOR, weight_mpi_datatype, MPI_SUM, mseg_root, comm);
              else
# endif
              {
# ifdef MSEG_REDUCEBCAST_THRESHOLD
                if (size >= MSEG_REDUCEBCAST_THRESHOLD)
                {
                  Z_TRACE_IF(MSEG_TRACE_IF, "%d >= %d: allreduce = reduce + bcast", size, (int) MSEG_REDUCEBCAST_THRESHOLD);

                  MPI_Reduce(lcw, gcw, MY_WEIGHT_FACTOR, weight_mpi_datatype, MPI_SUM, REDUCEBCAST_ROOT, comm);
                  MPI_Bcast(gcw, MY_WEIGHT_FACTOR, weight_mpi_datatype, REDUCEBCAST_ROOT, comm);

                } else
# endif
                  MPI_Allreduce(lcw, gcw, MY_WEIGHT_FACTOR, weight_mpi_datatype, MPI_SUM, comm);
              }
            }

          } else
#endif
          {
            /* global counts is just what we selected above */
            gcw[0] = mcw;
          }

#ifdef MSEG_ROOT
          if (mseg_root < 0 || mseg_root == rank)
#endif
          {
            gcs = gcw[0];
#ifdef elem_weight
            gws = gcw[1];
#endif
          }

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": gcs = %" slint_fmt, i, borders[i], gcs);
/*          Z_ASSERT_IF(MSEG_ASSERT_IF, currents[CNT_LO] <= gcs && gcs <= currents[CNT_HI]);*/ /* FIXME: only if exact */
#ifdef elem_weight
          if (doweights)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": gws = %" slweight_fmt, i, borders[i], gws);
/*            Z_ASSERT_IF(MSEG_ASSERT_IF, currents[WHT_LO] <= gws && gws <= currents[WHT_HI]);*/  /* FIXME: only if exact */
          }
#endif
        }

        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_final);
      }

      rti_tstart(rti_tid_mpi_select_exact_generic_while_check_post);

#ifdef MSEG_ROOT
      if (mseg_root < 0 || mseg_root == rank)
#endif
      {
#ifdef OLD
        border_change_old(doweights, &border_infos[borders[i]], &border_info_old, gcs, gc, elem_weight_ifelse(gws, 0.0), elem_weight_ifelse(gw, 0.0), (refine)?direction:0);
#else
        border_change(doweights, &border_infos[borders[i]], gcs, gc, elem_weight_ifelse(gws, 0.0), elem_weight_ifelse(gw, 0.0));
#endif
      }
      
      Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s", i, borders[i], (refine)?"REFINE":"REMOVE");

      if (refine)
      {
        borders[i - nborders_removed * direction] = borders[i];
        border_bins[i - nborders_removed * direction] = border_bins[i];

      } else
      {
        ++nborders_removed;
#ifdef OLD
#ifdef MSEG_BORDER_UPDATE_FULL
        border_infos[borders[i]].done = 1;
#endif
#endif

#ifdef MSEG_INFO
        if (mseg_info_finish_rounds) mseg_info_finish_rounds[borders[i]] = round;
        mseg_info_finish_rounds_avg += round;
#endif
      }

      rti_tstop(rti_tid_mpi_select_exact_generic_while_check_post);

      i += direction;
      
      Z_TRACE_IF(MSEG_TRACE_IF, "do partitioning: %" slint_fmt " vs. %" slint_fmt, i, ix);
      
      if (i != ix) goto do_partitioning;
    }
    rti_tstop(rti_tid_mpi_select_exact_generic_while_check);

    mpi_binning_post(&gb, size, rank, comm);

    /* restrict the parts */
    if (direction > 0) border_hi -= nborders_removed;
    else border_lo += nborders_removed;

    Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ": remove %" slint_fmt", lo: %" slint_fmt ", hi: %" slint_fmt "", round, nborders_removed, border_lo, border_hi);

    /* change direction */
    direction *= -1;
  }

#ifdef MSEG_INFO  
  mseg_info_rounds = round;
  mseg_info_finish_rounds_avg /= (size - 1);
#endif

  rti_tstop(rti_tid_mpi_select_exact_generic_while);

  mpi_binning_destroy(&gb, size, rank, comm);

  rti_tstop(rti_tid_mpi_select_exact_generic);

  rti_tstart(rti_tid_mpi_select_exact_generic_sync_exit);
#ifdef SYNC_ON_EXIT
  MPI_Barrier(comm);
#endif
  rti_tstop(rti_tid_mpi_select_exact_generic_sync_exit);

#ifdef VERIFY
  if (size == 1) v = -1;
  else v = mpi_post_check_partconds_intern(s, nelements, nparts, pci, sp->displs, size, rank, comm);
  
  Z_ASSERT_IF(MSEG_ASSERT_IF, v < 0);
  
  Z_NOTICE_IF(rank == 0, "post_check_partconds: %s (%" slint_fmt ")", (v >= 0)?"FAILED":"SUCCESS", v);
#endif

#ifdef PRINT_SDISPLS
  printf("%d: sdispls:", rank);
  for (i = 0; i < nparts; ++i) printf(" %d ", sp->displs[i]);
  printf("\n");
#endif

#ifdef PRINT_STATS
  mpi_select_stats(s, nparts, sp->displs, size, rank, comm);
#endif

#if defined(PRINT_TIMINGS) && defined(SL_USE_RTI_TIM)
  if (rank == PRINT_TIMINGS)
  {
    printf("%d: mpi_select_exact_generic: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic));
    printf("%d: mpi_select_exact_generic: sync init: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic_sync_init));
    printf("%d: mpi_select_exact_generic: while: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic_while));
    printf("%d: mpi_select_exact_generic:  check: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check));
    printf("%d: mpi_select_exact_generic:   bins: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_bins));
    printf("%d: mpi_select_exact_generic:    local: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_bins_local));
    printf("%d: mpi_select_exact_generic:    global: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_bins_global));
    printf("%d: mpi_select_exact_generic:   round1: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_round1));
    printf("%d: mpi_select_exact_generic:   pre: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_pre));
    printf("%d: mpi_select_exact_generic:   part: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_part));
    printf("%d: mpi_select_exact_generic:    root: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_part_root));
    printf("%d: mpi_select_exact_generic:   final: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_final));
    printf("%d: mpi_select_exact_generic:    root: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_final_root));
    printf("%d: mpi_select_exact_generic:   post: %f\n", rank, rti_tcumu(rti_tid_mpi_select_exact_generic_while_check_post));
    printf("%d: mpi_select_exact_generic: sync exit: %f\n", rank, rti_tlast(rti_tid_mpi_select_exact_generic_sync_exit));
    printf("%d: mpi_select_exact_generic: rounds: %" slint_fmt "\n", rank, round);
  }
#endif

  return 0;
}


slint_t mpi_select_exact_generic_grouped(elements_t *s, slint_t nelements, partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic_grouped */
{
  slint_t npconds = -1;
  partcond_t *pconds;


  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, NULL, &npconds, size, rank, comm);

  pconds = z_alloca(npconds, sizeof(partcond_t));
  
  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, pconds, &npconds, size, rank, comm);

  mpi_select_exact_generic_bulk(s, nelements, npconds, pconds, bm, sp, size, rank, comm);
  
  z_freea(pconds);
  
  return 0;
}


#if 0
slint_t mpi_select_exact_generic(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic */
{
  const slint_t max_nborders = nparts - 1;
  slint_t border_lo, border_hi, nborders_removed;
/*  slint_t borders[max_nborders], border_bins[max_nborders];

  border_info_t border_infos_[1 + max_nborders + 1], *border_infos = border_infos_ + 1, border_info_old;*/

  slint_t round, direction;

  slint_t i;
/*  slint_t i, j, k, l;*/

  slint_t curr_l, curr_g, curr_p;
  slint_t next_l, next_g, next_p;
  
  slint_t binning_at_once = 1;
  

  border_lo = 0;
  border_hi = max_nborders - 1;

  direction = 1;

  round = 0;
  while (border_lo <= border_hi)
  {
    ++round;

    nborders_removed = 0;

    i = (direction > 0)?border_lo:border_hi;

    next_l = i;
    next_g = -1;
    next_p = -1;
    
    while ((direction > 0)?(i <= border_hi):(i >= border_lo))
    {
      curr_l = next_l;
      curr_g = next_g;
      curr_p = next_p;


      if (border_lo <= curr_g && curr_g <= border_hi)
      {
        /* init global binning at curr_g */
      
        next_p = curr_g;
      }

      if (border_lo <= curr_p && curr_p <= border_hi)
      {
        /* wait global binning at curr_p */
        
      
        /* partitioning at curr_p */
        
        i += binning_at_once * direction;
      }
      
      if (border_lo <= curr_l && curr_l <= border_hi)
      {
        /* local binning at curr_l */

        next_l += binning_at_once * direction;
        next_g = curr_l;
      }
    }

    /* restrict the parts */
    if (direction > 0) border_hi -= nborders_removed;
    else border_lo += nborders_removed;

    /* change direction */
    direction *= -1;
  }
  
  return 0;
}
#endif
