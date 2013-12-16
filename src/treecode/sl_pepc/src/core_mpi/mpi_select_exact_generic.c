/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_select_exact_generic.c
 *  
 */


/* sl_macro MSEG_ROOT */
/* sl_macro MSEG_BORDER_UPDATE_REDUCTION */
/* sl_macro MSEG_DISABLE_BEST_CHOICE */
/* sl_macro MSEG_DISABLE_MINMAX */
/* sl_macro MSEG_ENABLE_OPTIMZED_LOWHIGH */
/* sl_macro MSEG_FORWARD_ONLY */
/* sl_macro MSEG_INFO */
/* sl_macro MSEG_TRACE_IF */


#include "sl_common.h"


/* config */
/*#define SYNC_ON_INIT
#define SYNC_ON_EXIT*/

/*#define PRINT_SDISPLS*/
/*#define PRINT_STATS*/
/*#define PRINT_TIMINGS  0*/

/*#define VERIFY*/

typedef struct _border_info_t {
  slint_t crange[2], ccurrent[2];
#ifdef elem_weight
  slweight_t wrange[2], wcurrent[2];
#endif

} border_info_t;


#define LO  0
#define HI  1


#ifdef MSEG_ROOT
int mseg_root = -1;  /* sl_global, sl_var mseg_root */
#endif

#ifdef MSEG_BORDER_UPDATE_REDUCTION
double mseg_border_update_count_reduction = 0.0;  /* sl_global, sl_var mseg_border_update_count_reduction */
# ifdef elem_weight
double mseg_border_update_weight_reduction = 0.0;  /* sl_global, sl_var mseg_border_update_weight_reduction */
# endif
#endif

#ifdef MSEG_FORWARD_ONLY
slint_t mseg_forward_only = 0;  /* sl_global, sl_var mseg_forward_only */
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


void border_init(slint_t docounts, slint_t doweights, border_info_t *bi, slint_t current, slint_t tc, slweight_t tw) /* sl_func border_init */
{
  if (docounts)
  {
    bi[0].crange[0] = 0;
    bi[0].crange[1] = tc;

    Z_TRACE_IF(MSEG_TRACE_IF, "count range: %" slint_fmt " - %" slint_fmt "",
      bi[0].crange[0], bi[0].crange[1]);

    bi[0].ccurrent[LO] = (current == 0 || current != 1)?0:tc;
    bi[0].ccurrent[HI] = (current == 1 || current != 0)?tc:0;

    Z_TRACE_IF(MSEG_TRACE_IF, "count[low/high]: %" slint_fmt " / %" slint_fmt,
      bi[0].ccurrent[LO], bi[0].ccurrent[HI]);
  }

#ifdef elem_weight
  if (doweights)
  {
    bi[0].wrange[0] = 0;
    bi[0].wrange[1] = tw;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight range: %" slweight_fmt " - %" slweight_fmt,
      bi[0].wrange[0], bi[0].wrange[1]);

    bi[0].wcurrent[LO] = (current == 0 || current != 1)?0:tw;
    bi[0].wcurrent[HI] = (current == 1 || current != 0)?tw:0;

    Z_TRACE_IF(MSEG_TRACE_IF, "weight[low/high]: %" slweight_fmt " / %" slweight_fmt,
      bi[0].wcurrent[LO], bi[0].wcurrent[HI]);

  }
#endif
}


void border_update(slint_t docounts, slint_t doweights, border_info_t *bi, partcond_intern_t *pc, slint_t dir, slint_t reduction) /* sl_func border_update */
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


#if 1
  /* if counts are used and the current interval is already as small as possible, then skip the update */
  if (docounts && bi[0].ccurrent[LO] == bi[0].ccurrent[HI])
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "skip border_update");
    return;
  }
#endif

  /* init from range */
  if (docounts)
  {
    ccurrent_new[LO] = bi[0].crange[0];
    ccurrent_new[HI] = bi[0].crange[1];
  }
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
    if (docounts)
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

          ccurrent_new[LO] = z_min(bi[-1].ccurrent[LO] + count_reduction + pc[0].count_min, bi[0].ccurrent[HI]);
          ccurrent_new[HI] = z_max(bi[-1].ccurrent[HI] - count_reduction + pc[0].count_max, bi[0].ccurrent[LO]);

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
    if (docounts)
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
  if (docounts)
  {
    if (pc[0].pcm & SLPC_COUNTS_LH)
    {
      if (ccurrent_new[LO] < pc[1].count_low)  ccurrent_new[LO] = pc[1].count_low;
      if (ccurrent_new[HI] > pc[0].count_high) ccurrent_new[HI] = pc[0].count_high;
    }

    /* fit to range (NOT REQUIRED!?) */
/*    ccurrent_new[LO] = z_minmax(bi[0].crange[0], ccurrent_new[LO], bi[0].crange[1]);
    ccurrent_new[HI] = z_minmax(bi[0].crange[0], ccurrent_new[HI], bi[0].crange[1]);*/

    bi[0].ccurrent[LO] = z_max(bi[0].ccurrent[LO], ccurrent_new[LO]);
    bi[0].ccurrent[HI] = z_min(bi[0].ccurrent[HI], ccurrent_new[HI]);

    Z_TRACE_IF(MSEG_TRACE_IF, "current count[low/high]: %" slint_fmt " / %" slint_fmt, bi[0].ccurrent[LO], bi[0].ccurrent[HI]);
  }

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


void border_change(slint_t docounts, slint_t doweights, border_info_t *bi, slint_t gcs, slint_t gc, slweight_t gws, slweight_t gw) /* sl_func border_change */
{
  if (docounts)
  {
    Z_TRACE_IF(MSEG_TRACE_IF, "change: gcs = %" slint_fmt ", gc = %" slint_fmt "", gcs, gc);

    bi[0].crange[0] += gcs;
    bi[0].crange[1] = bi[0].crange[0] + gc;

    Z_TRACE_IF(MSEG_TRACE_IF, "counts_range: %" slint_fmt "  %" slint_fmt "", bi[0].crange[0], bi[0].crange[1]);

    bi[0].ccurrent[LO] = z_minmax(bi[0].crange[0], bi[0].ccurrent[LO], bi[0].crange[1]);
    bi[0].ccurrent[HI] = z_minmax(bi[0].crange[0], bi[0].ccurrent[HI], bi[0].crange[1]);

    Z_TRACE_IF(MSEG_TRACE_IF, "count[low/high]: %" slint_fmt " / %" slint_fmt,
      bi[0].ccurrent[LO], bi[0].ccurrent[HI]);
  }

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


#ifdef MSEG_ENABLE_OPTIMZED_LOWHIGH
# define border_update_update(_dc_, _dw_, _bi_, _pc_, _dir_, _red_)        Z_NOP()
# ifdef elem_weight
# define border_change_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)  Z_MOP( \
  if (_dc_) { (_bi_)[0].crange[0] += (_gcs_); (_bi_)[0].crange[1] = (_bi_)[0].crange[0] + (_gc_); } \
  if (_dw_) { (_bi_)[0].wrange[0] += (_gws_); (_bi_)[0].wrange[1] = (_bi_)[0].wrange[0] + (_gw_); })
# else
# define border_change_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)  Z_MOP( \
  if (_dc_) { (_bi_)[0].crange[0] += (_gcs_); (_bi_)[0].crange[1] = (_bi_)[0].crange[0] + (_gc_); })
# endif
#else
# define border_update_update(_dc_, _dw_, _bi_, _pc_, _dir_, _red_)        border_update(_dc_, _dw_, _bi_, _pc_, _dir_, _red_)
# define border_change_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)  border_change(_dc_, _dw_, _bi_, _gcs_, _gc_, _gws_, _gw_)
#endif


slint_t mpi_select_exact_generic_bulk(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_generic_bulk */
{
  const slint_t max_nborders = nparts - 1;
  slint_t border_lo, border_hi, nborders_removed;
  slint_t borders[max_nborders], border_bins[max_nborders];

  border_info_t border_infos_[1 + max_nborders + 1], *border_infos = border_infos_ + 1;

  slint_t total_counts;
#ifdef elem_weight
  slweight_t total_weights;
#endif

  slint_t pcm;
  partcond_intern_t pci[nparts];

#if defined(elem_weight) && defined(sl_weight_intequiv)
  slweight_t current_cw[4];
# define current_clo  current_cw[0]
# define current_chi  current_cw[1]
# define current_wlo  current_cw[2]
# define current_whi  current_cw[3]
#else
  slint_t current_c[2];
# define current_clo  current_c[0]
# define current_chi  current_c[1]
# ifdef elem_weight
  slweight_t current_w[2];
#  define current_wlo  current_w[0]
#  define current_whi  current_w[1]
# endif
#endif

  slint_t round, direction, refine, finalize;

  slint_t nothing, lc_min, lc_max, lcw2gcw;

#if defined(elem_weight) && defined(sl_weight_intequiv)
  slweight_t final_lcw[2], final_gcw[2];
# define final_lc   final_lcw[0]
# define final_gc   final_gcw[0]
# define final_lw   final_lcw[1]
# define final_gw   final_gcw[1]
  slweight_t final_lcws[2], final_gcws[2];
# define final_lcs  final_lcws[0]
# define final_gcs  final_gcws[0]
# define final_lws  final_lcws[1]
# define final_gws  final_gcws[1]
#else
  slint_t final_lc, final_gc;
  slint_t final_lcs, final_gcs;
# ifdef elem_weight
  slweight_t final_lw, final_gw;
  slweight_t final_lws, final_gws;
# endif
#endif

  slint_t gc, gcs;
  slint_t final_mc, final_dc;
#ifdef elem_weight
  slweight_t gw, gws;
  slweight_t final_mw, final_dw;
#endif

  slint_t i, j, k, ix;

  slint_t docounts;
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

  pcm = pconds->pcm;
#ifdef MSEG_DISABLE_MINMAX
  pcm &= ~(SLPC_COUNTS_MM|SLPC_WEIGHTS_MM);
#endif

  docounts = ((pcm & (SLPC_COUNTS_MM|SLPC_COUNTS_LH)) != 0);
#ifdef elem_weight
  doweights = ((pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
  weight_factor = 1 + (doweights != 0);
#endif

  mpi_binning_create(&gb, max_nborders, mseg_binnings, s, nelements, docounts, doweights, bm, size, rank, comm);

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
  while (border_lo <= border_hi && (docounts || doweights))
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
          total_weights = 0;
#endif
          for (j = 0; j < gb.bm->nbins; ++j)
          {
            if (docounts) total_counts += *gb_counts(&gb, border_bins[i], j);
#ifdef elem_weight
            if (doweights) total_weights += *gb_weights(&gb, border_bins[i], j);
#endif
          }

          if (docounts) Z_TRACE_IF(MSEG_TRACE_IF, "total_counts = %" slint_fmt, total_counts);
#ifdef elem_weight
          if (doweights) Z_TRACE_IF(MSEG_TRACE_IF, "total_weights = %" slweight_fmt , total_weights);
#endif

          init_partconds_intern(nparts, pci, pconds, nparts, total_counts, elem_weight_ifelse(total_weights, 0));

          /* init lowest and highest part (sentinels) */
          Z_TRACE_IF(MSEG_TRACE_IF, "init lowest border:");
          border_init(docounts, doweights, &border_infos[border_lo - 1], 0, total_counts, elem_weight_ifelse(total_weights, 0));
          Z_TRACE_IF(MSEG_TRACE_IF, "init highest border:");
          border_init(docounts, doweights, &border_infos[border_hi + 1], 1, total_counts, elem_weight_ifelse(total_weights, 0));

#ifdef MSEG_BORDER_UPDATE_REDUCTION
          /* init+update forwards */
          for (j = border_lo; j <= border_hi; ++j)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init(docounts, doweights, &border_infos[borders[j]], -1, total_counts, elem_weight_ifelse(total_weights, 0));

            Z_TRACE_IF(MSEG_TRACE_IF, "update update %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update(docounts, doweights, &border_infos[borders[j]], &pci[borders[j]], 1, 0);
          }
#endif

          /* [init+]update backwards */
          for (j = border_hi; j >= border_lo; --j)
          {
#ifndef MSEG_BORDER_UPDATE_REDUCTION
            Z_TRACE_IF(MSEG_TRACE_IF, "init border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_init(docounts, doweights, &border_infos[borders[j]], -1, total_counts, elem_weight_ifelse(total_weights, 0));
#endif
            Z_TRACE_IF(MSEG_TRACE_IF, "update border %" slint_fmt ",%" slint_fmt ":", j, borders[j]);
            border_update(docounts, doweights, &border_infos[borders[j]], &pci[borders[j]], -1, 1);
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

        if (docounts)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": crange: %" slint_fmt " - %" slint_fmt, i, borders[i], border_infos[borders[i]].crange[0], border_infos[borders[i]].crange[1]);
#ifdef elem_weight
        if (doweights)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": wrange: %" slweight_fmt " - %" slweight_fmt, i, borders[i], border_infos[borders[i]].wrange[0], border_infos[borders[i]].wrange[1]);
#endif

        border_update_update(docounts, doweights, &border_infos[borders[i]], &pci[borders[i]], direction, 1);

        if (docounts)
        {
          current_clo = border_infos[borders[i]].ccurrent[LO] - border_infos[borders[i]].crange[0];
          current_chi = border_infos[borders[i]].ccurrent[HI] - border_infos[borders[i]].crange[0];

          Z_ASSERT_IF(MSEG_ASSERT_IF, current_clo <= current_chi);
          Z_ASSERT_IF(MSEG_ASSERT_IF, 0 <= current_clo);
        
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": currents count: %" slcount_fmt " - %" slcount_fmt " (range: %" slcount_fmt ")",
            i, borders[i], current_clo, current_chi, current_chi - current_clo);

        } else { current_clo = -1; current_chi = 1; }

#ifdef elem_weight
        if (doweights)
        {
          current_wlo = border_infos[borders[i]].wcurrent[LO] - border_infos[borders[i]].wrange[0];
          current_whi = border_infos[borders[i]].wcurrent[HI] - border_infos[borders[i]].wrange[0];

          Z_ASSERT_IF(MSEG_ASSERT_IF, current_wlo <= current_whi);

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": currents weight: %" slweight_fmt " - %" slweight_fmt " (range: %" slweight_fmt ")",
            i, borders[i], current_wlo, current_whi, current_whi - current_wlo);

        } else { current_wlo = -1; current_whi = 1; }
#endif
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
          gws = 0;
#endif

          k = 0;
          
          while (1)
          {
            /* check for HIT */

            /* HIT if max count already skipped */
            if (k == 0 && current_chi < 0) break;

            /* if between min/max counts */
            if (current_clo <= 0 && current_chi >= 0)
            {
#ifdef elem_weight
              if (doweights)
              {
                Z_TRACE_IF(MSEG_TRACE_IF, "go to next: %d && %d", (current_chi > 0), (current_wlo > 0));

                /* go to next if max count not reached AND min weight not reached */
                if (current_chi > 0 && current_wlo > 0) goto donthit;
              }
#endif

#ifndef MSEG_DISABLE_BEST_CHOICE
              /* look ahead for a better stop */
              if (k < gb.bm->nbins && (!docounts || current_chi - *gb_counts(&gb, border_bins[i], k) >= 0))
              {
#ifdef elem_weight
                if (doweights)
                {
                  /* continue if weights will improve */
                  if (z_abs(current_wlo + current_whi) > z_abs(current_wlo + current_whi - 2 * *gb_weights(&gb, border_bins[i], k))) goto donthit;

                } else
#endif
                {
                  /* continue if counts will improve */
                  if (z_abs(current_clo + current_chi) > z_abs(current_clo + current_chi - 2 * *gb_counts(&gb, border_bins[i], k))) goto donthit;
                }
              }
#endif

              /* HIT if there is no better stop */
              break;
            }

donthit:

            /* HIT in the worst case */
            if (k >= gb.bm->nbins) break;

            /* skip k-th bin */
            
            if (docounts)
            {
              gc = *gb_counts(&gb, border_bins[i], k);

              current_clo -= gc;
              current_chi -= gc;

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": k = %" slint_fmt ", currents count: %" slcount_fmt " - %" slcount_fmt ", gc = %" slint_fmt ", gcs = %" slint_fmt,
                i, borders[i], k, current_clo, current_chi, gc, gcs);

            } else gc = 0;

#ifdef elem_weight
            if (doweights)
            {
              gw = *gb_weights(&gb, border_bins[i], k);

              current_wlo -= gw;
              current_whi -= gw;

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": k = %" slint_fmt ", currents weight: %" slweight_fmt " - %" slweight_fmt ", gw = %" slweight_fmt ", gws = %" slweight_fmt,
                i, borders[i], k, current_wlo, current_whi, gw, gws);

            } else gw = 0;
#endif

            /* check for REFINE */

            /* stop and refine if max count is skipped OR (min count AND max weight is skipped) */
            if (current_chi < 0
#ifdef elem_weight
              /* '(!docounts || current_clo < 0)' is omitted, because if !docounts then current_clo = -1 */
              /* 'doweights &&' is omitted, because if !doweights then current_whi = 1 */
              || (current_clo < 0 && current_whi < 0)
#endif
              )
            {
              /* stop for REFINE if we do not know the counts, or */
              /* if counts are known and there are more than one elements to refine (otherwise a HIT follows next) */
              if (!docounts || gc > 1)
              {
                refine = 1;
                break;
              }
            }

            gcs += gc;
            gc = 0;

#ifdef elem_weight
            gws += gw;
            gw = 0;
#endif
            
            ++k;
          }

          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s k = %" slint_fmt, i, borders[i], (refine)?"REFINE":"HIT", k);

          Z_ASSERT_IF(MSEG_ASSERT_IF, ((refine && k < gb.bm->nbins) || (!refine && k <= gb.bm->nbins)));

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
        if (mseg_root >= 0)
        {
#if defined(elem_weight) && defined(sl_weight_intequiv)
          MPI_Bcast(current_cw, 4, weight_mpi_datatype, mseg_root, comm);
#else
          if (docounts)
            MPI_Bcast(current_c, 2, int_mpi_datatype, mseg_root, comm);
# ifdef elem_weight
          if (doweights)
            MPI_Bcast(current_w, 2, weight_mpi_datatype, mseg_root, comm);
# endif
#endif
        }
        rti_tstop(rti_tid_mpi_select_exact_generic_while_check_final_root);
#endif

        switch (mseg_finalize_mode)
        {
          case SL_MSEG_FM_ALLORNOTHING:
            Z_TRACE_IF(MSEG_TRACE_IF, "finalize mode: all or nothing");
#ifdef elem_weight
            if (doweights)
            {
              nothing = (current_wlo < ((border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]) - current_whi));
              Z_TRACE_IF(MSEG_TRACE_IF, "weight: %" slweight_fmt " vs. %" slweight_fmt " -> %s", current_wlo, (border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]) - current_whi, ((nothing)?"NOTHING":"ALL"));

            } else
#endif
            {
              nothing = (current_clo < ((border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]) - current_chi));
              Z_TRACE_IF(MSEG_TRACE_IF, "count: %" slint_fmt " vs. %" slint_fmt " -> %s", (slint_t) current_clo, (slint_t) ((border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]) - current_chi), ((nothing)?"NOTHING":"ALL"));
            }

            if (nothing)
            {
              final_mc = final_dc = 0;
#ifdef elem_weight
              final_mw = final_dw = 0;
#endif
              lc_min = lc_max = 0;

            } else
            {
#ifdef elem_weight
              if (doweights)
              {
                final_mw = final_dw = border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0];

              } else
#endif
              {
                final_mc = final_dc = border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0];
              }
              lc_min = lc_max = border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0];
            }

#ifdef elem_weight
            if (doweights)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": weights: final_mw = %" slweight_fmt ", final_dw = %" slweight_fmt, i, borders[i], final_mw, final_dw);
            else
#endif
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": counts: final_mc = %" slint_fmt ", final_dc = %" slint_fmt, i, borders[i], final_mc, final_dc);

            if (docounts)
              final_gcs = (nothing)?0:(border_infos[borders[i]].crange[1] - border_infos[borders[i]].crange[0]);
#ifdef elem_weight
            if (doweights)
              final_gws = (nothing)?0:(border_infos[borders[i]].wrange[1] - border_infos[borders[i]].wrange[0]);
#endif
            lcw2gcw = 0;
            break;
          
          case SL_MSEG_FM_MIDDLE:
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": finalize mode: middle (CURRENTLY MISSING!)", i, borders[i]);
            break;
          
          case SL_MSEG_FM_EXACT:
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": finalize mode: exact", i, borders[i]);

            final_lc = 0;
#ifdef elem_weight
            final_lw = 0;
#endif

            for (j = 0; j < nelements; ++j)
            {
              final_lc += lb_bin_count(&gb.lb, border_bins[i], j);
#ifdef elem_weight
              final_lw += lb_bin_weight(&gb.lb, border_bins[i], j);
#endif
            }

            if (docounts)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_lc: %" slcount_fmt, i, borders[i], final_lc);
#ifdef elem_weight
            if (doweights)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_lw: %" slweight_fmt, i, borders[i], final_lw);
#endif

#if defined(elem_weight) && defined(sl_weight_intequiv)
            MPI_Exscan(final_lcw, final_gcw, weight_factor, weight_mpi_datatype, MPI_SUM, comm);
            if (rank == 0) final_gcw[0] = final_gcw[1] = 0;
#else
            if (docounts)
            {
              MPI_Exscan(&final_lc, &final_gc, 1, int_mpi_datatype, MPI_SUM, comm);
              if (rank == 0) final_gc = 0;
            }
# ifdef elem_weight
            if (doweights)
            {
              MPI_Exscan(&final_lw, &final_gw, 1, weight_mpi_datatype, MPI_SUM, comm);
              if (rank == 0) final_gw = 0;
            }
# endif
#endif

            if (docounts)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_gc: %" slcount_fmt, i, borders[i], final_gc);
#ifdef elem_weight
            if (doweights)
              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": final_gw: %" slweight_fmt, i, borders[i], final_gw);
#endif

#ifdef elem_weight
            if (doweights)
            {
              /* middle of min/max weight */
              final_mw = (current_wlo + current_whi) / 2.0;

              /* min. part of weight to contribute */
              final_dw = z_max(0, final_mw - final_gw);

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": weights: final_mw = %" slweight_fmt ", final_dw = %" slweight_fmt, i, borders[i], final_mw, final_dw);

            } else
#endif
            {
              /* middle of min/max count */
              final_mc = (current_clo + current_chi) / 2;

              /* min. part of count to contribute */
              final_dc = z_max(0, final_mc - final_gc);

              Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": counts: final_mc = %" slint_fmt ", final_dc = %" slint_fmt, i, borders[i], final_mc, final_dc);
            }

            lc_min = current_clo - final_gc;
            lc_max = current_chi - final_gc;

            lcw2gcw = 1;
            break;
        }

        Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lc_min = %" slint_fmt ", lc_max = %" slint_fmt, i, borders[i], lc_min, lc_max);

        mpi_binning_finalize(&gb, border_bins[i], elem_weight_ifelse(0, final_dc), elem_weight_ifelse(final_dw, 0), lc_min, lc_max, &final_lcs, elem_weight_ifelse(&final_lws, NULL), sp, borders[i] + 1, size, rank, comm);

        if (docounts)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lcs_final = %" slcount_fmt, i, borders[i], final_lcs);
#ifdef elem_weight
        if (doweights)
          Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": lws_final = %" slweight_fmt, i, borders[i], final_lws);
#endif

        gcs = gc = 0;
#ifdef elem_weight
        gws = gw = 0;
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
              if (mseg_root >= 0)
              {
#  ifdef sl_weight_intequiv
                MPI_Reduce(final_lcws, final_gcws, weight_factor, weight_mpi_datatype, MPI_SUM, mseg_root, comm);
#  else
                MPI_Reduce(&final_lcs, &final_gcs, 1, int_mpi_datatype, MPI_SUM, mseg_root, comm);
                MPI_Reduce(&final_lws, &final_gws, 1, weight_mpi_datatype, MPI_SUM, mseg_root, comm);
#  endif
              } else
# endif
              {
# ifdef sl_weight_intequiv
                sl_MPI_Allreduce(final_lcws, final_gcws, weight_factor, weight_mpi_datatype, MPI_SUM, comm, size, rank);
# else
                sl_MPI_Allreduce(&final_lcs, &final_gcs, 1, int_mpi_datatype, MPI_SUM, comm, size, rank);
                sl_MPI_Allreduce(&final_lws, &final_gws, 1, weight_mpi_datatype, MPI_SUM, comm, size, rank);
# endif
              }
            }

          } else
#endif
          {
            /* global counts is just what we selected above */
            final_gcs = final_mc;
          }

#ifdef MSEG_ROOT
          if (mseg_root < 0 || mseg_root == rank)
#endif
          {
            gcs = final_gcs;
#ifdef elem_weight
            gws = final_gws;
#endif
          }

          if (docounts)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": gcs = %" slint_fmt, i, borders[i], gcs);
/*            Z_ASSERT_IF(MSEG_ASSERT_IF, current_clo <= gcs && gcs <= current_chi);*/ /* FIXME: only if exact */
          }
#ifdef elem_weight
          if (doweights)
          {
            Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": gws = %" slweight_fmt, i, borders[i], gws);
/*            Z_ASSERT_IF(MSEG_ASSERT_IF, current_wlo <= gws && gws <= current_whi);*/  /* FIXME: only if exact */
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
        border_change_change(docounts, doweights, &border_infos[borders[i]], gcs, gc, elem_weight_ifelse(gws, 0), elem_weight_ifelse(gw, 0));
      }
      
      Z_TRACE_IF(MSEG_TRACE_IF, "%" slint_fmt ",%" slint_fmt ": %s", i, borders[i], (refine)?"REFINE":"REMOVE");

      if (refine)
      {
        borders[i - nborders_removed * direction] = borders[i];
        border_bins[i - nborders_removed * direction] = border_bins[i];

      } else
      {
        ++nborders_removed;

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

#ifdef MSEG_FORWARD_ONLY
    if (mseg_forward_only)
    {
      /* do not change direction, but if there are min-max bounds, then perform a separate backward pass to update the (remaining) borders */
      if (pci->pcm & (SLPC_COUNTS_MM|SLPC_WEIGHTS_MM))
      {
        for (i = border_hi; i >= border_lo; --i) border_update_update(docounts, doweights, &border_infos[borders[i]], &pci[borders[i]], direction, -1);
      }

    } else
#endif
    {
      /* change direction */
      direction *= -1;
    }
  }

#ifdef MSEG_INFO  
  mseg_info_rounds = round;
  if (size > 1) mseg_info_finish_rounds_avg /= (size - 1); else mseg_info_finish_rounds_avg = 0.0;
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

#undef doweights
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


#undef SYNC_ON_INIT
#undef SYNC_ON_EXIT
#undef PRINT_SDISPLS
#undef PRINT_STATS
#undef PRINT_TIMINGS
#undef VERIFY
#undef LO
#undef HI
