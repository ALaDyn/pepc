/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos.h
 *  timestamp: 2009-11-23 22:06:40 +0100
 *  
 */


#ifndef __SL_PROTOS_H__
#define __SL_PROTOS_H__

/* src/core/checksum_crc.c */
unsigned short pepckeys_cs_crc16(elements_t *s, slint n, slint keys, slint data);
unsigned int pepckeys_cs_crc32(elements_t *s, slint n, slint keys, slint data);

/* src/core/elements.c */
slint_t pepckeys_elements_alloc(elements_t *s, slint_t nelements, slint_t keys, slint_t data);
slint_t pepckeys_elements_free(elements_t *s);
slint_t pepckeys_elements_alloc_from_block(elements_t *s, void *block, slint_t blocksize, slint_t alignment);
slint_t pepckeys_elements_copy(elements_t *s, elements_t *d);
slint_t pepckeys_elements_copy_at(elements_t *s, slint_t sat, elements_t *d, slint_t dat);
slint_t pepckeys_elements_ncopy(elements_t *s, elements_t *d, slint_t n);
slint_t pepckeys_elements_nmove(elements_t *s, elements_t *d, slint_t n);
slint_t pepckeys_elements_printf(elements_t *s);
slint_t pepckeys_elements_extract(elements_t *src, slint_t nelements, elements_t *dst0, elements_t *dst1);
slint_t pepckeys_elements_touch(elements_t *s);
slint_t pepckeys_elements_random_exchange(elements_t *s, slint_t rounds, elements_t *xs);
slint_t pepckeys_elements_init_keys(elements_t *s, slint_t dtype, slint_t _min, slint_t _max);
slint_t pepckeys_elements_init_keys_from_file(elements_t *s, slint_t data, char *filename, slint_t from, slint_t to, slint_t const_bytes_per_line);
slint_t pepckeys_elements_save_keys_to_file(elements_t *s, char *filename);
slint_t pepckeys_elements_validate_order(elements_t *s, slint_t n);
slint_t pepckeys_elements_validate_order_bmask(elements_t *s, slint_t n, slkey_pure_t bmask);
slint_t pepckeys_elements_validate_order_weight(elements_t *s, slint_t n, slkey_pure_t weight);
slint_t pepckeys_elements_print_keys(elements_t *s);

/* src/core/elements_packed.c */
slint_t pepckeys_elements_alloc_packed(packed_elements_t *s, slint_t nelements);
slint_t pepckeys_elements_free_packed(packed_elements_t *s);
slint_t pepckeys_elements_pack_indexed(elements_t *s, packed_elements_t *d, slindex_t *rindx, slindex_t *windx);
slint_t pepckeys_elements_pack(elements_t *s, packed_elements_t *d);
slint_t pepckeys_elements_unpack_indexed(packed_elements_t *s, elements_t *d, slindex_t *rindx, slindex_t *windx);
slint_t pepckeys_elements_unpack(packed_elements_t *s, elements_t *d);
slint_t pepckeys_elements_unpack_keys(packed_elements_t *s, slkey_t *k);

/* src/core/key2class.c */
slint pepckeys_key2class_equal(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_split(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_split_keys(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_random(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_ci_nocounts(slkey_t *k, slint i, void *ci);
slint pepckeys_key2class_ci(slkey_t *k, slint i, void *ci);

/* src/core/merge2_basic.c */
slint pepckeys_merge2_basic_01_x(elements_t *s0, elements_t *s1, elements_t *sx, m2x_func _x0_1, m2x_func _0x_1);
slint pepckeys_merge2_basic_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t, m2X_func _X0_1, m2X_func _0X_1);

/* src/core/merge2_basic_auto.c */
slint pepckeys_merge2_basic_auto_01_x(elements_t *s0, elements_t *s1, elements_t *sx);

/* src/core/merge2_basic_sbin.c */
slint pepckeys_merge2_basic_sbin_x0_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sbin_0x_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sbin_01_x(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sbin_01(elements_t *s0, elements_t *s1, elements_t *t);

/* src/core/merge2_basic_shyb.c */
slint pepckeys_merge2_basic_shyb_x0_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_shyb_0x_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_shyb_01_x(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_shyb_01(elements_t *s0, elements_t *s1, elements_t *t);

/* src/core/merge2_basic_sseq.c */
slint pepckeys_merge2_basic_sseq_x0_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sseq_0x_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sseq_01_x(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_sseq_01(elements_t *s0, elements_t *s1, elements_t *t);

/* src/core/merge2_basic_straight.c */
slint pepckeys_merge2_basic_straight_x0_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_0x_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_01_x(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_x_0_1(elements_t *s0, elements_t *s1, elements_t *sx);
slint pepckeys_merge2_basic_straight_X0_1(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint pepckeys_merge2_basic_straight_0X_1(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint pepckeys_merge2_basic_straight_01_X(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);
slint pepckeys_merge2_basic_straight_X0_1u(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);

/* src/core/pepckeys_merge2_compo_hula.c */
slint pepckeys_merge2_compo_hula(elements_t *s0, elements_t *s1, elements_t *xs);

/* src/core/pepckeys_merge2_compo_tridgell.c */
slint pepckeys_merge2_compo_tridgell(elements_t *s0, elements_t *s1, elements_t *sx);

/* src/core/pepckeys_merge2_memory_adaptive.c */
slint pepckeys_merge2_memory_adaptive(elements_t *s0, elements_t *s1, elements_t *sx);

/* src/core/merge2_simplify.c */
slint pepckeys_merge2_simplify_s1(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements);

/* src/core/search.c */
slint pepckeys_sl_search_sequential_lt(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_sequential_le(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_sequential_gt(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_sequential_ge(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_lt(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_le(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_gt(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_binary_ge(elements_t *s, slkey_t *k);
slint pepckeys_sl_search_hybrid_lt(elements_t *s, slkey_t *k, slint t);
slint pepckeys_sl_search_hybrid_le(elements_t *s, slkey_t *k, slint t);
slint pepckeys_sl_search_hybrid_gt(elements_t *s, slkey_t *k, slint t);
slint pepckeys_sl_search_hybrid_ge(elements_t *s, slkey_t *k, slint t);

/* src/core/sl_common.c */
slint pepckeys_ilog2c(slint x);
slint pepckeys_ilog2f(slint x);
slint pepckeys_print_bits(slint v);
slint pepckeys_pivot_random(elements_t *s);

/* src/core/sl_elem.c */
slint_t pepckeys_elem_set_data(elements_t *e, ...);
slint_t pepckeys_elem_reverse(elements_t *e, elements_t *t);
slint_t pepckeys_elem_nxchange_at(elements_t *e0, slint_t at0, elements_t *e1, slint_t at1, slint_t n, elements_t *t);
slint_t pepckeys_elem_nxchange(elements_t *e0, elements_t *e1, slint_t n, elements_t *t);
slint_t pepckeys_elem_nxchange_ro0(elements_t *e0, elements_t *e1, slint_t n, elements_t *t);
slint_t pepckeys_elem_rotate(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t pepckeys_elem_rotate_ro0(elements_t *e, slint_t m, slint_t n, elements_t *t);
slint_t pepckeys_elem_rotate_ro1(elements_t *e, slint_t m, slint_t n, elements_t *t);

/* src/core/pepckeys_sort_counting.c */
slint_t pepckeys_sort_counting_use_displs(elements_t *s, elements_t *d, slint_t ndispls, slint_t *displs);
slint_t pepckeys_sort_counting_use_counts(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts);
slint_t pepckeys_sort_counting_get_counts(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts);
slint_t pepckeys_sort_counting(elements_t *s, elements_t *d, slint_t ncounts);

/* src/core/pepckeys_sort_heap.c */
slint pepckeys_sort_heap(elements_t *s, elements_t *xs);

/* src/core/pepckeys_sort_insert.c */
slint pepckeys_sort_insert(elements_t *s, elements_t *sx);

/* src/core/sort_permute.c */
slint pepckeys_sort_permute_forward(elements_t *s, elements_t *sx, slint *perm, slint offset, slint mask_bit);
slint pepckeys_sort_permute_backward(elements_t *s, elements_t *sx, slint *perm, slint offset, slint mask_bit);

/* src/core/pepckeys_sort_quick.c */
slint pepckeys_sort_quick(elements_t *s, elements_t *xs);

/* src/core/pepckeys_sort_radix.c */
slint pepckeys_sort_radix(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth);

/* src/core/pepckeys_sort_radix_1bit.c */
slint pepckeys_sort_radix_1bit(elements_t *s, elements_t *sx, slint rhigh, slint rlow);

/* src/core/pepckeys_sort_radix_af.c */
slint pepckeys_sort_radix_af(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth);

/* src/core/pepckeys_sort_radix_iter.c */
slint pepckeys_sort_radix_iter(elements_t *s, elements_t *sx, slint presorted, slint rhigh, slint rlow, slint rwidth);

/* src/core/sortnet.c */
slint pepckeys_sn_hypercube_lh(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_hypercube_hl(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_odd_even_trans(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_batcher(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_bitonic(slint size, slint rank, slint stage, void *snp, slint *up);
slint pepckeys_sn_connected(slint size, slint rank, slint stage, void *snp, slint *up);

/* src/core/splitx.c */
slint pepckeys_split2_lt_ge(elements_t *s, slkey_pure_t *k, elements_t *t);
slint pepckeys_split2_le_gt(elements_t *s, slkey_pure_t *k, elements_t *t);
slint pepckeys_split3_lt_eq_gt(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle);
slint pepckeys_split3_lt_eq_gt_old(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle);
slint pepckeys_split2_b(elements_t *s, elements_t *sx, slkey_pure_t bmask);
slint pepckeys_splitk_k2c_af(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data);
slint pepckeys_splitk_k2c(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data);
slint pepckeys_splitk_k2c_count(elements_t *s, slint k, slint *c, k2c_func k2c, void *k2c_data);


#endif /* __SL_PROTOS_H__ */
