/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_rti.h
 *  timestamp: 2009-11-04 13:37:31 +0100
 *  
 */


#ifndef __SL_RTI_H__
#define __SL_RTI_H__


typedef struct
{
  int cmp, movek, moved;  /* DONE: kein slint mehr in konfigurationsweitem sl_rti.h*/
} rti_cmc;

#define my_rti_ccmp(m)           m.cmc.cmp
#define my_rti_cmovek(m)         m.cmc.movek
#define my_rti_cmoved(m)         m.cmc.moved

enum rti_tid
{
  rti_tid_all,
  rti_tid_sort_insert,
  rti_tid_sort_quick,
  rti_tid_sort_radix,
  rti_tid_sort_radix_iter,
  rti_tid_sort_permute_forward,
  rti_tid_sort_permute_backward,

  rti_tid_mpi_all,
  rti_tid_mpi_merge2,
  rti_tid_mpi_merge2_fe,
  rti_tid_mpi_merge2_xchg,
  rti_tid_mpi_merge2_local,
  rti_tid_mpi_mergek,
  rti_tid_mpi_mergek_merge2,

  rti_tid_mpi_splitk_exact,
  rti_tid_mpi_splitk_exact_init,
  rti_tid_mpi_splitk_exact_loop,
  rti_tid_mpi_splitk_exact_loop_walk,
  rti_tid_mpi_splitk_exact_loop_flow,
  rti_tid_mpi_splitk_exact_loop_flow_gather,
  rti_tid_mpi_splitk_exact_loop_flow_create,
  rti_tid_mpi_splitk_exact_loop_flow_reduce,
  rti_tid_mpi_splitk_exact_loop_flow_unbalance,
  rti_tid_mpi_splitk_exact_loop_dist,
  rti_tid_mpi_splitk_exact_loop_dist_pre,
  rti_tid_mpi_splitk_exact_loop_dist_a2av,

  rti_tid_mpi_splitk_dummy,
  rti_tid_mpi_splitk_dummy_init,
  rti_tid_mpi_splitk_dummy_loop,

  rti_tid_mpi_partition_joink,
  rti_tid_mpi_partition_joink_init,
  rti_tid_mpi_partition_joink_loop,
  rti_tid_mpi_partition_joink_loop_flow,
  rti_tid_mpi_partition_joink_loop_dist,
  
  rti_tid_mpi_partition_radix,
  rti_tid_mpi_partition_radix_sync,
  rti_tid_mpi_partition_radix_while,
  rti_tid_mpi_partition_radix_while_count,
  rti_tid_mpi_partition_radix_while_allreduce,
  rti_tid_mpi_partition_radix_while_round1,
  rti_tid_mpi_partition_radix_while_round1_allgather,
  rti_tid_mpi_partition_radix_while_check,
  rti_tid_mpi_partition_radix_final,

  rti_tid_mpi_sample_complete,
  rti_tid_mpi_sample_complete_gather,
  rti_tid_mpi_sample_complete_detect,
  rti_tid_mpi_sample_complete_bcast,

  rti_tid_mpi_select_qs,
  rti_tid_mpi_select_qs_pre,
  rti_tid_mpi_select_qs_loop,
  rti_tid_mpi_select_qs_part,
  rti_tid_mpi_select_qs_reduce_sizes,
  rti_tid_mpi_select_qs_area,
  rti_tid_mpi_select_qs_pivot_new,
  rti_tid_mpi_select_qs_pivot_gather,
  rti_tid_mpi_select_qs_pivot_detect,

  rti_tid_mpi_sample_select_qs,
  rti_tid_mpi_sample_select_qs_pre,
  rti_tid_mpi_sample_select_qs_select,

  rti_tid_mpi_sample_precise,
  rti_tid_mpi_sample_precise_llec,
  rti_tid_mpi_sample_precise_gather,
  rti_tid_mpi_sample_precise_detect,
  
  rti_tid_mpi_sample_permutation,

  rti_tid_mpi_sm_simple,
  rti_tid_mpi_sm_simple_sort,
  rti_tid_mpi_sm_simple_merge,

  rti_tids
};

typedef struct
{
  double start, stop;
  double last, cumu;
} rti_tim[rti_tids];

#define my_rti_tlast(m, t)       m.tim[t].last
#define my_rti_tcumu(m, t)       m.tim[t].cumu

/* DONE: kein slint mehr in konfigurationsweitem sl_rti.h*/
typedef struct
{
  int nalloc, nfree;
  int cur, max;

  /* tracing the size of an allocated piece of memory for bookkeeping and using by rti_mfree */
  int *alloc_sizes_hash_table;
} rti_mem;

typedef struct _rti
{
  /* compare-move-counter */
  rti_cmc cmc;
  /* timer */
  rti_tim tim;
  /* memory */
  rti_mem mem;
} rti;

#define my_rti_reset(m)          memset((void *) &m, 0, sizeof(m))


#endif /* __SL_RTI_H__ */
