/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_globals.h
 *  
 */


#ifndef __SL_GLOBALS_H__
#define __SL_GLOBALS_H__


/* src/core/sl_common.c */
extern rti rti_env;
extern int sl_mpi_rank_dummy;

/* src/core/sort_radix.c */
extern slint_t sr_ip_threshold;
extern slint_t sr_db_threshold;
extern slint_t sr_ma_threshold;

/* src/core_mpi/mpi_common.c */
#ifdef SL_USE_MPI
extern MPI_Datatype int_mpi_datatype;
extern MPI_Datatype key_mpi_datatype;
extern MPI_Datatype pkey_mpi_datatype;
extern MPI_Datatype pwkey_mpi_datatype;
extern MPI_Datatype index_mpi_datatype;
extern MPI_Datatype weight_mpi_datatype;
extern MPI_Datatype data_mpi_datatype[];
#endif
#ifdef SL_USE_MPI
extern int sl_mpi_rank;
#endif

/* src/core_mpi/mpi_elements.c */
extern void *me_sendrecv_replace_mem;
extern slint_t me_sendrecv_replace_memsize;
extern slint_t me_sendrecv_replace_mpi_maxsize;

/* src/core_mpi/mpi_elements_alltoall_specific.c */
extern double meas_t[];

/* src/core_mpi/mpi_select_exact_generic.c */
extern int mseg_root;
extern double mseg_border_update_count_reduction;
extern double mseg_border_update_weight_reduction;
extern slint_t mseg_forward_only;
extern slint_t mseg_info_rounds;
extern slint_t *mseg_info_finish_rounds;
extern double mseg_info_finish_rounds_avg;
extern slint_t mseg_binnings;
extern slint_t mseg_finalize_mode;

/* src/core_mpi/mpi_select_sample.c */
extern int mss_root;

/* src/core_mpi/mpi_sort_merge.c */
extern double msm_t[];
extern slint_t msm_sync;

/* src/core_mpi/mpi_sort_partition.c */
extern double msp_t[];
extern slint_t msp_sync;
extern partcond_t *msp_r_pc;

/* src/core_mpi/mpi_sort_special.c */
extern double mss_i_t[];
extern double mss_p_t[];
extern double mss_b_t[];
extern slint_t mss_sync;
extern slint_t mss_i_sync;
extern slint_t mss_p_sync;
extern slint_t mss_b_sync;


#endif /* __SL_GLOBALS_H__ */
