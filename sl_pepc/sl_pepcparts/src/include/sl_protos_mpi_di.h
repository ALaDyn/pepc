/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  timestamp: 2010-01-08 15:04:49 +0100
 *  
 */


#ifndef __SL_PROTOS_MPI_DI_H__
#define __SL_PROTOS_MPI_DI_H__

/* src/core_mpi/mpi_common.c */
slint_t pepcparts_mpi_datatypes_init_di();
slint_t pepcparts_mpi_datatypes_release_di();

/* src/core_mpi/mpi_elements.c */
slint pepcparts_mpi_elements_init_keys_from_file_di(elements_t *s, char *filename, slint from, slint to, slint const_bytes_per_line, slint root, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_elements_validate_order_di(elements_t *s, slint n, int size, int rank, MPI_Comm comm);
unsigned short pepcparts_mpi_cs16_di(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm);
unsigned int pepcparts_mpi_cs32_di(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
slint_t pepcparts_mpi_elements_packed_datatype_create_di(MPI_Datatype *pdt, slint_t structured);
slint_t pepcparts_mpi_elements_packed_datatype_destroy_di(MPI_Datatype *pdt);

/* src/core_mpi/pepcparts_mpi_find_exact_di.c */
slint_t pepcparts_mpi_find_exact_equal_di(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *ex_start, slint_t *ex_size, int size, int rank, MPI_Comm comm);
slint_t pepcparts_mpi_find_exact_di(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, slint_t *ex_start, slint_t *ex_sizes, slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_find_exact_old_di.c */
slint_t pepcparts_mpi_find_exact_old_di(elements_t *s, slint_t counterpart, slint_t high, elements_t *xs, slint_t *start, slint_t *end, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2_di.c */
slint_t pepcparts_mpi_merge2_di(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, merge2x_f m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2_old_di.c */
slint pepcparts_mpi_merge2_old_di(elements_t *s, slint counterpart, slint high, m2x_func m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek_di.c */
slint_t pepcparts_mpi_mergek_equal_di(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm);
slint_t pepcparts_mpi_mergek_di(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek_old_di.c */
slint pepcparts_mpi_mergek_old_di(elements_t *s, sn_func sn, void *snp, m2x_func m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_joink_di.c */
slint pepcparts_mpi_partition_joink_di(elements_t *s, slint *sizes, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_radix_di.c */
slint_t pepcparts_mpi_partition_radix_di(elements_t *s, partcond_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_radix_old_di.c */
slint_t pepcparts_mpi_partition_radix_old_di(elements_t *s, partcond_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_rebalance_di.c */
slint_t pepcparts_mpi_rebalance_di(elements_t *s0, elements_t *s1, slint_t stable, slint_t *dst_size, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sample_ci.c */
slint pepcparts_mpi_sample_ci_init_di(classification_info *ci, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_sample_ci_duplicate_di(classification_info *ci_src, classification_info *ci_dup, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_sample_ci_free_di(classification_info *ci, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_sample_ci_print_di(classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_complete_di.c */
slint pepcparts_mpi_sample_complete_di(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_permutation_di.c */
slint pepcparts_mpi_sample_permutation_di(elements_t *s, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_precise_counts_di.c */
slint pepcparts_mpi_sample_precise_counts_di(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_select_qs_di.c */
slint pepcparts_mpi_sample_select_qs_di(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_select_qs_di.c */
slint pepcparts_mpi_select_qs_di(elements_t *s, slint n, slint *iths, pivot_func pi, slint threshold, slkey_pure_t *e, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sort_merge_di.c */
slint_t pepcparts_mpi_sort_merge_di(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_splitk_di.c */
slint pepcparts_mpi_splitk_di(elements_t *s, k2c_func k2c, void *ci, elements_t *sx, elements_t *sa, slint *nne, slint *nue, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_splitk_dummy_di.c */
slint pepcparts_mpi_splitk_dummy_di(elements_t *s, k2c_func k2c, void *ci, elements_t *sx, slint *send_stats, int size, int rank, MPI_Comm comm);


#endif /* __SL_PROTOS_MPI_DI_H__ */
