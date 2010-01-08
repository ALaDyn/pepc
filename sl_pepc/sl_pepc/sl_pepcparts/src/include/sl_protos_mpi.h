/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  timestamp: 2010-01-05 17:56:41 +0100
 *  
 */


#ifndef __SL_PROTOS_MPI_H__
#define __SL_PROTOS_MPI_H__

/* src/core_mpi/mpi_common.c */
slint_t pepcparts_mpi_datatypes_init();
slint_t pepcparts_mpi_datatypes_release();

/* src/core_mpi/mpi_elements.c */
slint pepcparts_mpi_elements_init_keys_from_file(elements_t *s, char *filename, slint from, slint to, slint const_bytes_per_line, slint root, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_elements_validate_order(elements_t *s, slint n, int size, int rank, MPI_Comm comm);
unsigned short pepcparts_mpi_cs16(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm);
unsigned int pepcparts_mpi_cs32(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
slint_t pepcparts_mpi_elements_packed_datatype_create(MPI_Datatype *pdt, slint_t structured);
slint_t pepcparts_mpi_elements_packed_datatype_destroy(MPI_Datatype *pdt);

/* src/core_mpi/pepcparts_mpi_find_exact.c */
slint_t pepcparts_mpi_find_exact_equal(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *ex_start, slint_t *ex_size, int size, int rank, MPI_Comm comm);
slint_t pepcparts_mpi_find_exact(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, slint_t *ex_start, slint_t *ex_sizes, slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_find_exact_old.c */
slint_t pepcparts_mpi_find_exact_old(elements_t *s, slint_t counterpart, slint_t high, elements_t *xs, slint_t *start, slint_t *end, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2.c */
slint_t pepcparts_mpi_merge2(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, merge2x_f m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2_old.c */
slint pepcparts_mpi_merge2_old(elements_t *s, slint counterpart, slint high, m2x_func m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek.c */
slint_t pepcparts_mpi_mergek_equal(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm);
slint_t pepcparts_mpi_mergek(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek_old.c */
slint pepcparts_mpi_mergek_old(elements_t *s, sn_func sn, void *snp, m2x_func m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_joink.c */
slint pepcparts_mpi_partition_joink(elements_t *s, slint *sizes, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_radix.c */
slint_t pepcparts_mpi_partition_radix(elements_t *s, partcond_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_radix_old.c */
slint_t pepcparts_mpi_partition_radix_old(elements_t *s, partcond_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_rebalance.c */
slint_t pepcparts_mpi_rebalance(elements_t *s0, elements_t *s1, slint_t stable, slint_t *dst_size, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sample_ci.c */
slint pepcparts_mpi_sample_ci_init(classification_info *ci, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_sample_ci_duplicate(classification_info *ci_src, classification_info *ci_dup, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_sample_ci_free(classification_info *ci, int size, int rank, MPI_Comm comm);
slint pepcparts_mpi_sample_ci_print(classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_complete.c */
slint pepcparts_mpi_sample_complete(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_permutation.c */
slint pepcparts_mpi_sample_permutation(elements_t *s, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_precise_counts.c */
slint pepcparts_mpi_sample_precise_counts(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_select_qs.c */
slint pepcparts_mpi_sample_select_qs(elements_t *s, slint threshold, classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_select_qs.c */
slint pepcparts_mpi_select_qs(elements_t *s, slint n, slint *iths, pivot_func pi, slint threshold, slkey_pure_t *e, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sort_merge.c */
slint_t pepcparts_mpi_sort_merge(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_splitk.c */
slint pepcparts_mpi_splitk(elements_t *s, k2c_func k2c, void *ci, elements_t *sx, elements_t *sa, slint *nne, slint *nue, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_splitk_dummy.c */
slint pepcparts_mpi_splitk_dummy(elements_t *s, k2c_func k2c, void *ci, elements_t *sx, slint *send_stats, int size, int rank, MPI_Comm comm);


#endif /* __SL_PROTOS_MPI_H__ */
