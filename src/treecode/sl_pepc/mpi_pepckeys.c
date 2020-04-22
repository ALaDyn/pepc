
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "sl_pepckeys.h"


#include "fortran2c_types.h"


typedef pepckeys_slint_t slint_t;
#define slint_fmt pepckeys_sl_int_type_fmt

/*#define MAX_IMBALANCE  0.01*/

#define PART_MINMAX

/*#define MPI_PARTITION_RADIX_2GROUPS*/
/*#define MPI_PARTITION_RADIX_NGROUPS  2*/
/*#define MPI_PARTITION_SAMPLE*/

/*#define SORT_INSTEAD_OF_MERGE*/

/*#define VERBOSE*/
/*#define VALIDATE*/
/*#define TIMING*/


#define SORT_RHIGH   -1
#define SORT_RLOW    -1
#define SORT_RWIDTH  -1

#define PART_RHIGH   62
#define PART_RLOW    -1
#define PART_RWIDTH  3


#ifdef VERBOSE
# define VERBOSE_MOP(_mop_)  do { _mop_ ; } while (0)
#else
# define VERBOSE_MOP(_mop_)  do { } while (0)
#endif

#ifdef TIMING
# define TSTART(tid)  tid = MPI_Wtime()
# define TSTOP(tid)   tid = MPI_Wtime() - tid
#else
# define TSTART(tid)  do { } while (0)
# define TSTOP(tid)   do { } while (0)
#endif


void receive_stats(pepckeys_sl_globalcount_t nmax, int *scounts, int *sdispls, int *rcounts, int *rdispls, pepckeys_sldata0_t *work, int size, int rank, MPI_Comm comm);
void border_stats(slint_t nkeys, pepckeys_slkey_t *keys, int size, int rank, MPI_Comm comm);


void slsort_keys(pepckeys_sl_globalcount_t *nin,                             /* IN */
                 pepckeys_sl_globalcount_t *nmax,                            /* IN */
                 pepckeys_slkey_t *keys,                                     /* INOUT */
                 pepckeys_sldata0_t *work,                                   /* INOUT */
                 pepckeys_slindex_t *balance_weight,                         /* IN */
                 pepckeys_sldata0_t *max_imbalance,                          /* IN */
                 pepckeys_sl_globalcount_t *nout,                            /* OUT */
                 pepckeys_slindex_t *indxl, pepckeys_slindex_t *irnkl,       /* OUT */
                 pepckeys_slindex_t *fscounts, pepckeys_slindex_t *frcounts, /* OUT */
                 pepckeys_slindex_t *fsdispls, pepckeys_slindex_t *frdispls, /* OUT */
                 pepckeys_slkey_t *keys2,                                    /* SCRATCH */
                 pepckeys_slindex_t *irnkl2,                                 /* SCRATCH */
                 MPI_Fint *fsize, MPI_Fint *frank, MPI_Fint *fcomm)          /* IN */
{
  int size = *fsize;
  int rank = *frank;
  // see http://www.open-mpi.org/community/lists/users/2006/09/1812.php for conversion of fortran communicator to C
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);

  slint_t i;

  pepckeys_elements_t s0, s1;
  pepckeys_partcond_t pc;

#ifdef MAX_IMBALANCE
  double imba = MAX_IMBALANCE;
#else
  double imba = *max_imbalance;
#endif

  int scounts[size], rcounts[size], sdispls[size], rdispls[size];

#if defined(MPI_PARTITION_RADIX_2GROUPS) || defined(MPI_PARTITION_RADIX_NGROUPS)
  const slint_t max_nsubs = 4;

  slint_t nsubs;
  MPI_Comm sub_comms[max_nsubs];
  int sub_sizes[max_nsubs], sub_ranks[max_nsubs];
#endif

#ifdef VALIDATE
  slint_t o, l;
#endif

#ifdef TIMING
  double ttotal, tinitpre, tpresort, tpartition, talltoall, talltoallv, tinitpost, tpostmerge;
#endif

#if defined(VERBOSE) || defined(TIMING)
  slint_t ntotal = *n * size;
  const pepcparts_slint_t ndims = 4;
  pepcparts_slint_t dims[ndims], pos[ndims];
#endif


#if defined(VERBOSE) || defined(TIMING)
  pepcparts_mpi_get_grid_properties(ndims, dims, pos, size, rank, comm);
  if (rank == 0)
    printf("# np: %d, grid: %dx%dx%dx%d, n: %" finteger_fmt ", nmax: %" finteger_fmt ", ntotal: %" slint_fmt ", partitioning: %s, minmax: %s, weighted: %" finteger_fmt ", imba: %f\n",
      size, (int) dims[3], (int) dims[2], (int) dims[1], (int) dims[0], *n, *nmax, ntotal,
# ifdef MPI_PARTITION_SAMPLE
      "rs",
# else
      "ep",
# endif
# ifdef PART_MINMAX
      "yes",
# else
      "no",
# endif
      *balance_weight, imba);
#endif

  TSTART(ttotal);

#ifdef VERBOSE
  printf("%d: slsort_keys with %d processes\n", rank, size);
  printf("%d:  nin: %" finteger_fmt "\n", rank, *nin);
  printf("%d:  nmax: %" finteger_fmt "\n", rank, *nmax);
  printf("%d:  balance_weight: %" finteger_fmt "\n", rank, *balance_weight);
  printf("%d:  max_imbalance: %f\n", rank, *max_imbalance);

  printf("%d:  sizeof(integer) = %d\n", rank, (int) sizeof(FINT_TYPE_C));
  printf("%d:  sizeof(integer*8) = %d\n", rank, (int) sizeof(FINT8_TYPE_C));
  printf("%d:  sizeof(pepckeys_key) = %d\n", rank, (int) sizeof(pepckeys_slkey_t));
  printf("%d:  sizeof(pepckeys_index) = %d\n", rank, (int) sizeof(pepckeys_slindex_t));
#endif

  pepckeys_sl_mpi_rank = rank;

  pepckeys_mpi_datatypes_init();

  TSTART(tinitpre);

  /* init pre indexes */
  for (i = 0; i < *nin; ++i) indxl[i] = i + 1;

  pepckeys_elem_set_size(&s0, *nin);
  pepckeys_elem_set_max_size(&s0, *nmax);
  pepckeys_elem_set_keys(&s0, keys);
  pepckeys_elem_set_indices(&s0, indxl);
  pepckeys_elem_set_data(&s0, work);

  TSTOP(tinitpre);

  /* pre sort local */
  VERBOSE_MOP(printf("%d: slsort_keys: 1. sort keys\n", rank));

  TSTART(tpresort);
  pepckeys_sort_radix(&s0, NULL, SORT_RHIGH, SORT_RLOW, SORT_RWIDTH);
  TSTOP(tpresort);

  VERBOSE_MOP(printf("%d: slsort_keys: 1. sort keys done\n", rank));


#ifdef VALIDATE
  l = pepckeys_elements_validate_order(&s0, 1);
#endif


  /* partitioning */
  VERBOSE_MOP(printf("%d: slsort_keys: 2. partitioning\n", rank));

#define REDUCTION  0.25

#if defined(PART_MINMAX) && !defined(MPI_PARTITION_SAMPLE)
  if (*balance_weight == 0)
  {
    pc.pcm = SLPC_COUNTS_MM;
    pc.count_min = -(1.0 - imba);
    pc.count_max = -(1.0 + imba);

    pepckeys_mseg_border_update_count_reduction = REDUCTION;

  } else
  {
    pc.pcm = SLPC_COUNTS_MM|SLPC_WEIGHTS_MM;
    pc.count_min = 0;
    pc.count_max = *nmax;
    pc.weight_min = -(1.0 - imba);
    pc.weight_max = -(1.0 + imba);

    pepckeys_mseg_border_update_weight_reduction = REDUCTION;
  }
#else
  if (*balance_weight == 0)
  {
    pc.pcm = SLPC_COUNTS_LH;
    pc.count_low = ((double) (rank + 0) / (double) size) - (0.5 * imba / size);
    pc.count_low = -((pc.count_low > 0)?pc.count_low:0);
    pc.count_high = ((double) (rank + 1) / (double) size) + (0.5 * imba / size);
    pc.count_high = -((pc.count_high < 1)?pc.count_high:1);

  } else
  {
    pc.pcm = SLPC_WEIGHTS_LH;
    pc.weight_low = ((double) (rank + 0) / (double) size) - (0.5 * imba / size);
    pc.weight_low = -((pc.weight_low > 0)?pc.weight_low:0);
    pc.weight_high = ((double) (rank + 1) / (double) size) + (0.5 * imba / size);
    pc.weight_high = -((pc.weight_high < 1)?pc.weight_high:1);

    pc.pcm |= SLPC_COUNTS_MM;
    pc.count_min = 0;
    pc.count_max = *nmax;
  }
#endif

/*  pepckeys_mseg_finalize_mode = SL_MSEG_FM_ALLORNOTHING;*/

  TSTART(tpartition);

#if defined(MPI_PARTITION_SAMPLE)

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: regular sampling\n", rank));

/*  pepckeys_mss_root = -1;*/

  pepckeys_mpi_partition_sample_regular(&s0, &pc, scounts, NULL, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_2GROUPS)

  nsubs = 2;

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: 2groups\n", rank));

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_2groups(&s0, &pc, sub_comms[1], NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_NGROUPS)

  nsubs = (MPI_PARTITION_RADIX_NGROUPS <= max_nsubs)?MPI_PARTITION_RADIX_NGROUPS:max_nsubs;

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: ngroups (%" slint_fmt ")\n", rank, nsubs));

  pepckeys_elem_set_indices(&k0, indxl2); /* ..._ngroups requires indices that can be modified (indxl2 is not used otherwise) */

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_ngroups(&s0, &pc, nsubs, sub_comms, NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#else

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: direct\n", rank));

  pepckeys_mpi_partition_exact_radix(&s0, &pc, PART_RHIGH, PART_RLOW, PART_RWIDTH, SL_SORTED_IN, scounts, NULL, size, rank, comm);

  VERBOSE_MOP(printf("%d: average finish round: %f\n", rank, pepckeys_mseg_info_finish_rounds_avg));

#endif

  TSTOP(tpartition);

  pepckeys_counts2displs(size, scounts, sdispls);

  VERBOSE_MOP(printf("%d: slsort_keys: 2. partitioning done\n", rank));


  /* alltoallv keys */
  VERBOSE_MOP(printf("%d: slsort_keys: 3. alltoallv\n", rank));

  TSTART(talltoall);
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  TSTOP(talltoall);
  for (rdispls[0] = 0, i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  TSTART(talltoallv);
  MPI_Alltoallv(keys, scounts, sdispls, pepckeys_sl_key_type_mpi, keys2, rcounts, rdispls, pepckeys_sl_key_type_mpi, comm);
  TSTOP(talltoallv);

  *nout = rdispls[size - 1] + rcounts[size - 1];

  for (i = 0; i < size; ++i)
  {
    fscounts[i] = scounts[i];
    frcounts[i] = rcounts[i];
    fsdispls[i] = sdispls[i];
    frdispls[i] = rdispls[i];
  }

  VERBOSE_MOP(printf("%d: slsort_keys: 3. alltoallv done\n", rank));


  /* post sort local (or mergep) */
  VERBOSE_MOP(printf("%d: slsort_keys: 4. merge keys\n", rank));

#ifdef SORT_INSTEAD_OF_MERGE

  TSTART(tinitpost);

  /* init post indexes */
  for (i = 0; i < *nout; ++i) irnkl2[i] = i;

  pepckeys_elem_set_size(&s0, *nout);
  pepckeys_elem_set_max_size(&s0, *nmax);
  pepckeys_elem_set_keys(&s0, keys2);
  pepckeys_elem_set_indices(&s0, irnkl2);
  pepckeys_elem_set_data(&s0, work);

  pepckeys_elem_set_size(&s1, *nout);
  pepckeys_elem_set_max_size(&s1, *nmax);
  pepckeys_elem_set_keys(&s1, keys2);
  pepckeys_elem_set_indices(&s1, irnkl2);
  pepckeys_elem_set_data(&s1, work);

  TSTOP(tinitpost);

  TSTART(tpostmerge);
  pepckeys_sort_radix(&s0, NULL, -1, -1, -1);
  TSTOP(tpostmerge);

#else

  TSTART(tinitpost);

  /* init post indexes */
  for (i = 0; i < *nout; ++i) irnkl[i] = i;

  memcpy(keys, keys2, *nout * sizeof(pepckeys_slkey_t));

  pepckeys_elem_set_size(&s0, *nout);
  pepckeys_elem_set_max_size(&s0, *nmax);
  pepckeys_elem_set_keys(&s0, keys);
  pepckeys_elem_set_indices(&s0, irnkl);
  pepckeys_elem_set_data(&s0, work);

  pepckeys_elem_set_size(&s1, *nout);
  pepckeys_elem_set_max_size(&s1, *nmax);
  pepckeys_elem_set_keys(&s1, keys2);
  pepckeys_elem_set_indices(&s1, irnkl2);
  pepckeys_elem_set_data(&s1, work);

  TSTOP(tinitpost);

  TSTART(tpostmerge);
  pepckeys_mergep_heap_idx(&s0, &s1, size, frdispls, frcounts);
  TSTOP(tpostmerge);

#endif

  for (i = 0; i < s0.size; ++i) irnkl[irnkl2[i]] = i + 1;

  VERBOSE_MOP(printf("%d: slsort_keys: 4. merge keys done\n", rank));


#ifdef VALIDATE
  o = pepckeys_mpi_elements_validate_order(&s1, 1, size, rank, comm);
 #ifdef VERBOSE
  printf("%d: slsort_keys: global order: %s - local order: %s\n", rank, (!o)?"success":"FAILED", (!l)?"success":"failed");
 #endif
#endif


  pepckeys_mpi_datatypes_release();

  TSTOP(ttotal);

  VERBOSE_MOP(printf("%d: nout: %" FINT_TYPE_FMT "\n", rank, *nout));

#ifdef TIMING
  if (rank == 0)
  {
    printf("%d: slsort_keys: %f\n", rank, ttotal);
    printf("%d: slsort_keys: initpre: %f\n", rank, tinitpre);
    printf("%d: slsort_keys: presort: %f\n", rank, tpresort);
    printf("%d: slsort_keys: partition: %f\n", rank, tpartition);
    printf("%d: slsort_keys: alltoall: %f\n", rank, talltoall);
    printf("%d: slsort_keys: alltoallv: %f\n", rank, talltoallv);
    printf("%d: slsort_keys: initpost: %f\n", rank, tinitpost);
    printf("%d: slsort_keys: postmerge: %f\n", rank, tpostmerge);
  }
#endif
}

//
//void slsort_keys_(pepckeys_sl_globalcount_t *nin,                            /* IN */
//                 pepckeys_sl_globalcount_t *nmax,                            /* IN */
//                 pepckeys_slkey_t *keys,                                     /* INOUT */
//                 pepckeys_sldata0_t *work,                                   /* INOUT */
//                 pepckeys_slindex_t *balance_weight,                         /* IN */
//                 pepckeys_sldata0_t *max_imbalance,                          /* IN */
//                 pepckeys_sl_globalcount_t *nout,                            /* OUT */
//                 pepckeys_slindex_t *indxl, pepckeys_slindex_t *irnkl,       /* OUT */
//                 pepckeys_slindex_t *fscounts, pepckeys_slindex_t *frcounts, /* OUT */
//                 pepckeys_slindex_t *fsdispls, pepckeys_slindex_t *frdispls, /* OUT */
//                 pepckeys_slkey_t *keys2,                                    /* SCRATCH */
//                 pepckeys_slindex_t *irnkl2,                                 /* SCRATCH */
//                 MPI_Fint *fsize, MPI_Fint *frank, MPI_Fint *fcomm)          /* IN */
//{
//  slsort_keys(nin, nmax, keys, work, balance_weight, max_imbalance, nout, indxl, irnkl, fscounts, frcounts, fsdispls, frdispls, keys2, irnkl2, fsize, frank, fcomm);
//}
//
