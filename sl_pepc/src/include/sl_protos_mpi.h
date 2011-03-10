/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  timestamp: 2011-03-06 21:59:33 +0100
 *  
 */


#ifndef __SL_PROTOS_MPI_H__
#define __SL_PROTOS_MPI_H__


/* src/core_mpi/mpi_binning.c */
slint_t SL_PROTO(mpi_binning_create)(global_bins_t *gb, slint_t max_nbins, slint_t max_nbinnings, elements_t *s, slint_t nelements, slint_t doweights, binning_t *bm, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_destroy)(global_bins_t *gb, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_pre)(global_bins_t *gb, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_exec_reset)(global_bins_t *gb, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_exec_local)(global_bins_t *gb, slint_t b, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_exec_global)(global_bins_t *gb, slint_t root, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_refine)(global_bins_t *gb, slint_t b, slint_t k, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_hit)(global_bins_t *gb, slint_t b, slint_t k, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_finalize)(global_bins_t *gb, slint_t b, slweight_t dcw, slint_t lc_min, slint_t lc_max, slweight_t *lcw, splitter_t *sp, slint_t s, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_binning_post)(global_bins_t *gb, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_common.c */
slint_t SL_PROTO(mpi_datatypes_init)();
slint_t SL_PROTO(mpi_datatypes_release)();
slint_t SL_PROTO(mpi_get_grid_properties)(slint_t ndims, slint_t *dims, slint_t *pos, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_get_grid)(slint_t ndims, slint_t *dims, slint_t *pos, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_subgroups_create)(slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_subgroups_delete)(slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements.c */
slint SL_PROTO(mpi_elements_init_keys_from_file)(elements_t *s, char *filename, slint from, slint to, slint const_bytes_per_line, slint root, int size, int rank, MPI_Comm comm);
slint SL_PROTO(mpi_elements_validate_order)(elements_t *s, slint n, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_linear_exchange_pure_keys)(slkey_pure_t *in, slkey_pure_t *out, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_check_order)(elements_t *s, slint_t nelements, slint_t *orders, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_check_global_order)(slkey_pure_t local_min, slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_digest_sum)(elements_t *s, slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(mpi_cs32)(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_digest_hash)(elements_t *s, slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_get_counts)(elements_t *s, slint_t *clocal, slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
slweight_t SL_PROTO(mpi_elements_get_weights)(elements_t *s, slweight_t *wlocal, slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_get_counts_and_weights)(elements_t *s, slint_t nelements, slint_t *counts, slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_sendrecv_replace)(elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(mpi_elements_crc32)(elements_t *s, slint_t n, slint_t keys, slint_t data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_alltoallv.c */
slint_t SL_PROTO(mpi_elements_alltoallv_db)(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_elements_alltoallv_ip)(elements_t *sbuf, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
slint_t SL_PROTO(mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, slint_t structured);
slint_t SL_PROTO(mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);

/* src/core_mpi/mpi_find_exact.c */
slint_t SL_PROTO(mpi_find_exact_equal)(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *ex_start, slint_t *ex_size, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_find_exact)(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, slint_t *ex_start, slint_t *ex_sizes, slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_linsplit.c */
slint_t SL_PROTO(mpi_linsplit)(MPI_Comm comm_in, slkey_pure_t *keys_in, MPI_Comm *comms_out, slint_t *parity, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_linsplit2)(MPI_Comm comm_in, slkey_pure_t *keys_in, MPI_Comm *comms_out, slint_t *parity, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_merge2.c */
slint_t SL_PROTO(mpi_merge2)(elements_t *s, slint_t other_rank, slint_t high_rank, slint_t *dst_size, merge2x_f m2, elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_mergek.c */
slint_t SL_PROTO(mpi_mergek_equal)(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_mergek)(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_mergek_equal2)(elements_t *s, sortnet_f sn, sortnet_data_t snd, merge2x_f m2x, elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);

/* src/core_mpi/mpi_partition_exact_generic.c */
slint_t SL_PROTO(mpi_partition_exact_generic)(elements_t *s, partcond_t *pcond, binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_partition_exact_generic2)(elements_t *s, partcond_t *pcond, binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_exact_radix.c */
slint_t SL_PROTO(mpi_partition_exact_radix)(elements_t *s, partcond_t *pcond, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_exact_radix_grouped.c */
slint_t SL_PROTO(mpi_partition_exact_radix_ngroups)(elements_t *s, partcond_t *pcond, slint_t ngroups, MPI_Comm *group_comms, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_partition_exact_radix_2groups)(elements_t *s, partcond_t *pcond, MPI_Comm group_comm, elements_t *sx, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_radix2.c */
slint_t SL_PROTO(mpi_partition_radix2)(elements_t *s, partcond2_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_radix_old.c */
slint_t SL_PROTO(mpi_partition_radix_old)(elements_t *s, partcond2_t *pc, slint_t rhigh, slint_t rlow, slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_sample.c */
slint_t SL_PROTO(mpi_partition_sample_regular)(elements_t *s, partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_rebalance.c */
slint_t SL_PROTO(mpi_rebalance)(elements_t *s0, elements_t *s1, slint_t stable, slint_t *dst_size, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_rebalance_alltoallv)(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_common.c */
slint_t SL_PROTO(init_partconds)(slint_t npconds, partcond_t *pconds, slint_t nparts, slint_t total_count, slweight_t total_weight);
slint_t SL_PROTO(init_partconds_intern)(slint_t npconds, partcond_intern_t *pci, partcond_t *pc, slint_t nparts, slint_t total_count, slweight_t total_weight);
slint_t SL_PROTO(merge_partconds)(partcond_t *pconds_in, slint_t npconds_in, partcond_t *pcond_out);
slint_t SL_PROTO(mpi_gather_partconds_grouped)(partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, partcond_t *pconds_out, slint_t *npconds_out, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_gather_partconds)(partcond_t *pcond_in, partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_allgather_partconds)(partcond_t *pcond_in, partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_bcast_partconds)(slint_t npconds, partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_post_check_partconds)(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_post_check_partconds_intern)(elements_t *s, slint_t nelements, slint_t nparts, partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_select_stats)(elements_t *s, slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_exact_generic.c */
slint_t SL_PROTO(mpi_select_exact_generic_bulk)(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_select_exact_generic_grouped)(elements_t *s, slint_t nelements, partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_select_exact_generic)(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, binning_t *bm, splitter_t *sp, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_exact_radix.c */
slint_t SL_PROTO(mpi_select_exact_radix)(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_select_exact_radix_grouped)(elements_t *s, slint_t nelements, partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_sample.c */
slint_t SL_PROTO(mpi_select_sample_regular)(elements_t *s, slint_t nparts, partcond_t *pconds, slint_t nsamples, splitter_t *sp, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sort_merge.c */
slint_t SL_PROTO(mpi_sort_merge)(elements_t *s0, elements_t *s1, elements_t *xs, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_sort_merge2)(elements_t *s0, elements_t *s1, elements_t *xs, slint_t merge_type, slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sort_partition.c */
slint_t SL_PROTO(mpi_sort_partition)(elements_t *s0, elements_t *s1, elements_t *xs, slint_t part_type, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_sort_partition_radix)(elements_t *s0, elements_t *s1, elements_t *xs, slint_t part_type, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_sort_partition_exact_radix)(elements_t *s, elements_t *sx, partcond_t *pcond, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_sort_partition_exact_radix_ngroups)(elements_t *s, elements_t *sx, partcond_t *pcond, slint_t ngroups, MPI_Comm *group_comms, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_sort_partition_exact_radix_2groups)(elements_t *s, elements_t *sx, partcond_t *pcond, MPI_Comm group_comm, slint_t rhigh, slint_t rlow, slint_t rwidth, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_xcounts2ycounts.c */
slint_t SL_PROTO(mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, slint_t ytotal, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_xcounts2ycounts_grouped)(int *xcounts, slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
slint_t SL_PROTO(mpi_subxdispls2ycounts)(slint_t nsubs, int *sub_xdispls, slint_t *sub_sources, slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* __SL_PROTOS_MPI_H__ */
