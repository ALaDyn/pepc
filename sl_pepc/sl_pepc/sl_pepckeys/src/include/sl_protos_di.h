/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos.h
 *  timestamp: 2010-01-05 17:56:41 +0100
 *  
 */


#ifndef __SL_PROTOS_DI_H__
#define __SL_PROTOS_DI_H__

/* src/core/checksum_crc.c */
unsigned short pepckeys_cs_crc16_di(elements_t *s, slint n, slint keys, slint data);
unsigned int pepckeys_cs_crc32_di(elements_t *s, slint n, slint keys, slint data);

/* src/core/elements.c */
slint_t pepckeys_elements_alloc_di(elements_t *s, slint_t nelements, slint_t keys, slint_t data);
slint_t pepckeys_elements_free_di(elements_t *s);
slint_t pepckeys_elements_alloc_from_block_di(elements_t *s, void *block, slint_t blocksize, slint_t alignment);
slint_t pepckeys_elements_copy_di(elements_t *s, elements_t *d);
slint_t pepckeys_elements_copy_at_di(elements_t *s, slint_t sat, elements_t *d, slint_t dat);
slint_t pepckeys_elements_ncopy_di(elements_t *s, elements_t *d, slint_t n);
slint_t pepckeys_elements_nmove_di(elements_t *s, elements_t *d, slint_t n);
slint_t pepckeys_elements_printf_di(elements_t *s);
slint_t pepckeys_elements_extract_di(elements_t *src, slint_t nelements, elements_t *dst0, elements_t *dst1);
slint_t pepckeys_elements_touch_di(elements_t *s);
slint_t pepckeys_elements_random_exchange_di(elements_t *s, slint_t rounds, elements_t *xs);
slint_t pepckeys_elements_init_keys_di(elements_t *s, slint_t dtype, slint_t _min, slint_t _max);
slint_t pepckeys_elements_init_keys_from_file_di(elements_t *s, slint_t data, char *filename, slint_t from, slint_t to, slint_t const_bytes_per_line);
slint_t pepckeys_elements_save_keys_to_file_di(elements_t *s, char *filename);
slint_t pepckeys_elements_validate_order_di(elements_t *s, slint_t n);
slint_t pepckeys_elements_validate_order_bmask_di(elements_t *s, slint_t n, slkey_pure_t bmask);
slint_t pepckeys_elements_validate_order_weight_di(elements_t *s, slint_t n, slkey_pure_t weight);
slint_t pepckeys_elements_print_keys_di(elements_t *s);

/* src/core/elements_packed.c */
slint_t pepckeys_elements_alloc_packed_di(packed_elements_t *s, slint_t nelements);
slint_t pepckeys_elements_free_packed_di(packed_elements_t *s);
slint_t pepckeys_elements_pack_indexed_di(elements_t *s, packed_elements_t *d, slindex_t *rindx, slindex_t *windx);
slint_t pepckeys_elements_pack_di(elements_t *s, packed_elements_t *d);
slint_t pepckeys_elements_unpack_indexed_di(packed_elements_t *s, elements_t *d, slindex_t *rindx, slindex_t *windx);
slint_t pepckeys_elements_unpack_di(packed_elements_t *s, elements_t *d);
slint_t pepckeys_elements_unpack_keys_di(packed_elements_t *s, slkey_t *k);

/* src/core/key2class.c */
slint pepckeys_key2class_equal_di(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_split_di(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_split_keys_di(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_random_di(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_ci_nocounts_di(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_ci_di(slkey_t *k, slint i, void *ci);

/* src/core/merge2_basic.c */
slint pepckeys_merge2_basic_01_x_di(elements_t *s0, elements_t *s1, elements_t *sx, m2x_func _x0_1, m2x_func _0x_1);
slint pepckeys_merge2_basic_01_X_di(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t, m2X_func _X0_1, m2X_func _0X_1);

/* src/core/merge2_basic_auto.c */
slint pepckeys_merge2_basic_auto_01_x_di(elements_t *s0, elements_t *s1, elements_t *sx);

/* src/core/merge2_basic_sbin.c */
slint pepckeys_merge2_basic_sbin_x0_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sbin_0x_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sbin_01_x_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sbin_01_di(elements_t *s0, elements_t *s1, elements_t *t);

/* src/core/merge2_basic_shyb.c */
slint pepckeys_merge2_basic_shyb_x0_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_shyb_0x_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_shyb_01_x_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_shyb_01_di(elements_t *s0, elements_t *s1, elements_t *t);

/* src/core/merge2_basic_sseq.c */
slint pepckeys_merge2_basic_sseq_x0_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sseq_0x_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sseq_01_x_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sseq_01_di(elements_t *s0, elements_t *s1, elements_t *t);

/* src/core/merge2_basic_straight.c */
slint pepckeys_merge2_basic_straight_x0_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_0x_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_01_x_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_x_0_1_di(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_X0_1_di(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint pepckeys_merge2_basic_straight_0X_1_di(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint pepckeys_merge2_basic_straight_01_X_di(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint pepckeys_merge2_basic_straight_X0_1u_di(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);

/* src/core/pepckeys_merge2_compo_hula_di.c */
slint pepckeys_merge2_compo_hula_di(elements_t *s0, elements_t *s1, elements_t *xs);

/* src/core/pepckeys_merge2_compo_tridgell_di.c */
slint pepckeys_merge2_compo_tridgell_di(elements_t *s0, elements_t *s1, elements_t *sx);

/* src/core/pepckeys_merge2_memory_adaptive_di.c */
slint pepckeys_merge2_memory_adaptive_di(elements_t *s0, elements_t *s1, elements_t *sx);

/* src/core/merge2_simplify.c */
slint pepckeys_merge2_simplify_s1_di(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements);

/* src/core/pepckeys_mergep_heap_di.c */
slint pepckeys_mergep_heap_di(elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts);
slint pepckeys_mergep_heap_unpack_di(packed_elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts);

/* src/core/search.c */
slint pepckeys_sl_search_sequential_lt_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_sequential_le_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_sequential_gt_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_sequential_ge_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_lt_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_le_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_gt_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_ge_di(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_hybrid_lt_di(elements_t *s, slkey_t *k, slint t);
slint pepckeys_sl_search_hybrid_le_di(elements_t *s, slkey_t *k, slint t);
slint pepckeys_sl_search_hybrid_gt_di(elements_t *s, slkey_t *k, slint t);
slint pepckeys_sl_search_hybrid_ge_di(elements_t *s, slkey_t *k, slint t);

/* src/core/sl_common.c */
slint pepckeys_ilog2c_di(slint x);
slint pepckeys_ilog2f_di(slint x);
slint pepckeys_print_bits_di(slint v);
slint pepckeys_pivot_random_di(elements_t *s);

/* src/core/sl_elem.c */
slint_t pepckeys_elem_set_data_di(elements_t *e, ...);
slint_t pepckeys_elem_reverse_di(elements_t *e, elements_t *t);
slint_t pepckeys_elem_nxchange_at_di(elements_t *e0, slint_t at0, elements_t *e1, slint_t at1, slint_t n, elements_t *t);
slint_t pepckeys_elem_nxchange_di(elements_t *e0, elements_t *e1, slint_t n, elements_t *t);
slint_t pepckeys_elem_nxchange_ro0_di(elements_t *e0, elements_t *e1, slint_t n, elements_t *t);
slint_t pepckeys_elem_rotate_di(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t pepckeys_elem_rotate_ro0_di(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t pepckeys_elem_rotate_ro1_di(elements_t *e, slint_t m, slint_t n, elements_t *t);

/* src/core/pepckeys_sort_counting_di.c */
slint_t pepckeys_sort_counting_use_displs_di(elements_t *s, elements_t *d, slint_t ndispls, slint_t *displs);
slint_t pepckeys_sort_counting_use_counts_di(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts);
slint_t pepckeys_sort_counting_get_counts_di(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts);
slint_t pepckeys_sort_counting_di(elements_t *s, elements_t *d, slint_t ncounts);

/* src/core/pepckeys_sort_heap_di.c */
slint pepckeys_sort_heap_di(elements_t *s, elements_t *xs);

/* src/core/pepckeys_sort_insert_di.c */
slint pepckeys_sort_insert_di(elements_t *s, elements_t *sx);

/* src/core/sort_permute.c */
slint pepckeys_sort_permute_forward_di(elements_t *s, elements_t *sx, slint *perm, slint offset, slint mask_bit);
slint pepckeys_sort_permute_backward_di(elements_t *s, elements_t *sx, slint *perm, slint offset, slint mask_bit);

/* src/core/pepckeys_sort_quick_di.c */
slint pepckeys_sort_quick_di(elements_t *s, elements_t *xs);

/* src/core/pepckeys_sort_radix_di.c */
slint pepckeys_sort_radix_di(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth);

/* src/core/pepckeys_sort_radix_1bit_di.c */
slint pepckeys_sort_radix_1bit_di(elements_t *s, elements_t *sx, slint rhigh, slint rlow);

/* src/core/pepckeys_sort_radix_af_di.c */
slint pepckeys_sort_radix_af_di(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth);

/* src/core/pepckeys_sort_radix_iter_di.c */
slint pepckeys_sort_radix_iter_di(elements_t *s, elements_t *sx, slint presorted, slint rhigh, slint rlow, slint rwidth);

/* src/core/sortnet.c */
slint pepckeys_sn_hypercube_lh_di(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_hypercube_hl_di(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_odd_even_trans_di(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_batcher_di(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_bitonic_di(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_connected_di(slint size, slint rank, slint stage, void *snp, slint *up);

/* src/core/splitx.c */
slint pepckeys_split2_lt_ge_di(elements_t *s, slkey_pure_t *k, elements_t *t);
slint pepckeys_split2_le_gt_di(elements_t *s, slkey_pure_t *k, elements_t *t);
slint pepckeys_split3_lt_eq_gt_di(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle);
slint pepckeys_split3_lt_eq_gt_old_di(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle);
slint pepckeys_split2_b_di(elements_t *s, elements_t *sx, slkey_pure_t bmask);
slint pepckeys_splitk_k2c_af_di(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data);
slint pepckeys_splitk_k2c_di(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data);
slint pepckeys_splitk_k2c_count_di(elements_t *s, slint k, slint *c, k2c_func k2c, void *k2c_data);


#endif /* __SL_PROTOS_DI_H__ */
