
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "sl_pepckeys.h"
#include "sl_pepcparts.h"

#include "fortran2c_types.h"


typedef FINT_TYPE_C finteger_t;
#define finteger_mpi  FINT_TYPE_MPI
#define finteger_fmt  FINT_TYPE_FMT

/*#define MAX_IMBALANCE  0.01*/

/*#define MPI_PARTITION_RADIX_2GROUPS*/
/*#define MPI_PARTITION_RADIX_NGROUPS  2*/

/*#define MPI_PARTITION_SAMPLE*/

/*#define MERGE_AND_UNPACK*/

#define PART_MINMAX

/*#define VERBOSE*/
/*#define VALIDATE*/
/*#define BORDER_STATS*/
/*#define RECEIVE_STATS*/
/*#define TIMING*/
/*#define TIMING_ROW*/


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


void slsort_parts_(finteger_t *, finteger_t *, pepcparts_slkey_t *, pepcparts_sldata0_t *, pepcparts_sldata1_t *, pepcparts_sldata2_t *, pepcparts_sldata3_t *, pepcparts_sldata4_t *, pepcparts_sldata5_t *,
                   pepcparts_sldata6_t *, pepcparts_sldata7_t *, pepcparts_sldata8_t *, pepcparts_sldata9_t *, pepcparts_sldata10_t *, pepcparts_sldata11_t *, pepcparts_sldata12_t *,
                   finteger_t *, double *, finteger_t *, finteger_t *, finteger_t *, finteger_t *, finteger_t *, finteger_t *,
                   void *, void *, pepcparts_slkey_t *, pepcparts_sldata8_t *, finteger_t *, finteger_t *, finteger_t *);

#pragma weak slsort_parts_ = slsort_parts
void slsort_parts(finteger_t *n,                                                             /* INOUT */
                  finteger_t *nmax,                                                          /* IN */
                  pepcparts_slkey_t *keys,                                                   /* INOUT */
                  pepcparts_sldata0_t *x,                                                    /* INOUT */
                  pepcparts_sldata1_t *y,                                                    /* INOUT */
                  pepcparts_sldata2_t *z,                                                    /* INOUT */
                  pepcparts_sldata3_t *ux,                                                   /* INOUT */
                  pepcparts_sldata4_t *uy,                                                   /* INOUT */
                  pepcparts_sldata5_t *uz,                                                   /* INOUT */
                  pepcparts_sldata6_t *q,                                                    /* INOUT */
                  pepcparts_sldata7_t *m,                                                    /* INOUT */
                  pepcparts_sldata8_t *work,                                                 /* INOUT */
                  pepcparts_sldata9_t  *ex,                                                  /* INOUT */
                  pepcparts_sldata10_t *ey,                                                  /* INOUT */
                  pepcparts_sldata11_t *ez,                                                  /* INOUT */
                  pepcparts_sldata12_t *pelabel,                                             /* INOUT */
                  finteger_t *balance_weight,                                                /* IN */
                  double *max_imbalance,                                                     /* IN */
                  finteger_t *indxl, finteger_t *irnkl,                                      /* OUT */
                  finteger_t *fscounts, finteger_t *frcounts,                                /* OUT */
                  finteger_t *fsdispls, finteger_t *frdispls,                                /* OUT */
                  void *parts0, void *parts1,                                                /* SCRATCH */
                  pepcparts_slkey_t *keys2, pepcparts_sldata8_t *work2, finteger_t *irnkl2,  /* SCRATCH */
                  finteger_t *fsize, finteger_t *frank)                                      /* IN */
{
  int size = *fsize;
  int rank = *frank;
  MPI_Comm comm = MPI_COMM_WORLD;

  typedef pepckeys_slint_t slint_t;
#define slint_fmt pepckeys_sl_int_type_fmt

  slint_t i;

  pepckeys_elements_t k0, k1;
  pepckeys_partcond_t pc;

  pepcparts_elements_t d0;
  pepcparts_packed_elements_t pd0;
  MPI_Datatype pdt;

  finteger_t nin = *n;

#ifdef MAX_IMBALANCE
  double imba = MAX_IMBALANCE;
#else
  double imba = *max_imbalance;
#endif

  int scounts[size], rcounts[size], sdispls[size], rdispls[size];

#define indxl2  irnkl2

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
  double ttotal;
  double tinitindxl, tcopy, tpresort, tpartition, tpack, talltoall,talltoallv, tmergeunpack, tmakeindices;
#ifndef MERGE_AND_UNPACK
  double tunpack;
#endif
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
  printf("%d: slsort_parts with %d processes\n", rank, size);
  printf("%d:  n: %" finteger_fmt "\n", rank, *n);
  printf("%d:  nmax: %" finteger_fmt "\n", rank, *nmax);
  printf("%d:  balance_weight: %" finteger_fmt "\n", rank, *balance_weight);
  printf("%d:  max_imbalance: %f\n", rank, *max_imbalance);

  printf("%d:  sizeof(integer) = %d\n", rank, (int) sizeof(finteger_t));
  printf("%d:  sizeof(integer*8) = %d\n", rank, (int) sizeof(FINT8_TYPE_C));
  printf("%d:  sizeof(pepckeys_key) = %d\n", rank, (int) sizeof(pepckeys_slkey_t));
  printf("%d:  sizeof(pepckeys_index) = %d\n", rank, (int) sizeof(pepckeys_slindex_t));
  printf("%d:  sizeof(pepcparts_key) = %d\n", rank, (int) sizeof(pepcparts_slkey_t));
  printf("%d:  sizeof(pepcparts_index) = %d\n", rank, (int) sizeof(pepcparts_slindex_t));
#endif

  pepckeys_sl_mpi_rank = rank;
  pepcparts_sl_mpi_rank = rank;

  pepckeys_mpi_datatypes_init();
  pepcparts_mpi_datatypes_init();

  /* init pre indexes */
  TSTART(tinitindxl);
  for (i = 0; i < *n; ++i) indxl[i] = i;
  TSTOP(tinitindxl);

  pepckeys_elem_set_size(&k0, *n);
  pepckeys_elem_set_max_size(&k0, *nmax);
  pepckeys_elem_set_keys(&k0, keys);
  pepckeys_elem_set_indices(&k0, indxl);
  pepckeys_elem_set_data(&k0, work);

  pepckeys_elem_set_size(&k1, *n);
  pepckeys_elem_set_max_size(&k1, *nmax);
  pepckeys_elem_set_keys(&k1, keys2);
  pepckeys_elem_set_indices(&k1, indxl2); /* just a dummy array, indxl2 is never used */
  pepckeys_elem_set_data(&k1, work2);

  TSTART(tcopy);
  pepckeys_elements_ncopy(&k0, &k1, *n);
  TSTOP(tcopy);

  /* pre sort keys (+ indices and work) */
  VERBOSE_MOP(printf("%d: slsort_parts: 1. sort keys\n", rank));

  TSTART(tpresort);
  pepckeys_sort_radix(&k0, NULL, SORT_RHIGH, SORT_RLOW, SORT_RWIDTH);
  TSTOP(tpresort);

/*  for (i = 0; i < k0.size; ++i) printf("%d  %d  %f\n", rank, (int) i, work[i]);*/

  VERBOSE_MOP(printf("%d: slsort_parts: 1. sort keys done\n", rank));


#ifdef VALIDATE
  l = pepckeys_elements_validate_order(&k0, 1);
#endif


  /* partitioning */
  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning\n", rank));

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
  }
#endif

/*  pepckeys_mseg_finalize_mode = SL_MSEG_FM_ALLORNOTHING;*/

  TSTART(tpartition);

#if defined(MPI_PARTITION_SAMPLE)

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: regular sampling\n", rank));

/*  pepckeys_mss_root = -1;*/

  pepckeys_mpi_partition_sample_regular(&k0, &pc, scounts, NULL, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_2GROUPS)

  nsubs = 2;

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: 2groups\n", rank));

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_2groups(&k0, &pc, sub_comms[1], NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#elif defined(MPI_PARTITION_RADIX_NGROUPS)

  nsubs = (MPI_PARTITION_RADIX_NGROUPS <= max_nsubs)?MPI_PARTITION_RADIX_NGROUPS:max_nsubs;

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: ngroups (%" slint_fmt ")\n", rank, nsubs));

  pepckeys_elem_set_indices(&k0, indxl2); /* ..._ngroups requires indices that can be modified (indxl2 is not used otherwise) */

  pepckeys_mpi_subgroups_create(nsubs, sub_comms, sub_sizes, sub_ranks, size, rank, comm);
  pepckeys_mpi_partition_exact_radix_ngroups(&k0, &pc, nsubs, sub_comms, NULL, PART_RHIGH, PART_RLOW, PART_RWIDTH, scounts, NULL, size, rank, comm);
  pepckeys_mpi_subgroups_delete(nsubs, sub_comms, size, rank, comm);

#else

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning: direct\n", rank));

  pepckeys_mpi_partition_exact_radix(&k0, &pc, PART_RHIGH, PART_RLOW, PART_RWIDTH, SL_SORTED_IN, scounts, NULL, size, rank, comm);

/*  if (rank == 0) printf("average finish round: %f\n", pepckeys_mseg_info_finish_rounds_avg);*/

#endif

  TSTOP(tpartition);

  pepckeys_counts2displs(size, scounts, sdispls);

  VERBOSE_MOP(printf("%d: slsort_parts: 2. partitioning done\n", rank));


  /* pack */
  VERBOSE_MOP(printf("%d: slsort_parts: 3. pack\n", rank));

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys2);
/*  pepcparts_elem_set_indices(&d0, indxl); */  /* indices are not required for packaging */
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work2, ex, ey, ez, pelabel);

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts0);

  TSTART(tpack);
  pepcparts_elements_pack_indexed(&d0, &pd0, (pepcparts_slindex_t *) indxl, NULL);
  TSTOP(tpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 3. pack done\n", rank));


  /* alltoallv */
  VERBOSE_MOP(printf("%d: slsort_parts: 4. alltoallv\n", rank));

  pepcparts_mpi_elements_packed_datatype_create(&pdt, 0);

  TSTART(talltoall);
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  TSTOP(talltoall);

  for (rdispls[0] = 0, i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];

#ifdef RECEIVE_STATS
  printf("%d: slsort_parts: total receive count: %d (nmax: %" finteger_fmt ")%s\n", rank, rdispls[size - 1] + rcounts[size - 1], *nmax, (rdispls[size - 1] + rcounts[size - 1] > *nmax)?" ERROR!!!":"");

  double w, sweights[size], rweights[size];
  slint_t j;
  for (j = 0; j < size; ++j)
  {
    sweights[j] = 0;
    for (i = sdispls[j]; i < sdispls[j]+scounts[j]; ++i) sweights[j] += work[i];
  }
  MPI_Alltoall(sweights, 1, MPI_DOUBLE, rweights, 1, MPI_DOUBLE, comm);
  w = 0;
  for (i = 0; i < size; ++i) w += rweights[i];

/*  printf("%d: slsort_parts: total receive weight: %f\n", rank, w);*/

  double cw[2], cw_min[2], cw_max[2], cw_sum[2];
  cw[0] = rdispls[size - 1] + rcounts[size - 1];
  cw[1] = w;
  MPI_Reduce(cw, cw_min, 2, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(cw, cw_max, 2, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(cw, cw_sum, 2, MPI_DOUBLE, MPI_SUM, 0, comm);
  if (rank == 0) printf("%d: slsort_parts: receive stats: %d  %d  %d  /  %f  %f  %f\n", rank, (int) cw_min[0], (int) cw_max[0], (int) (cw_sum[0] / size), cw_min[1], cw_max[1], cw_sum[1] / size);
#endif

  TSTART(talltoallv);
  MPI_Alltoallv(parts0, scounts, sdispls, pdt, parts1, rcounts, rdispls, pdt, comm);
  TSTOP(talltoallv);

  pepcparts_mpi_elements_packed_datatype_destroy(&pdt);

  *n = rdispls[size - 1] + rcounts[size - 1];

  for (i = 0; i < size; ++i)
  {
    fscounts[i] = scounts[i];
    frcounts[i] = rcounts[i];
    fsdispls[i] = sdispls[i];
    frdispls[i] = rdispls[i];
  }

  VERBOSE_MOP(printf("%d: slsort_parts: 4. alltoallv done\n", rank));


#ifdef MERGE_AND_UNPACK

  /* fused merge and unpack (indices created during merge) */
  VERBOSE_MOP(printf("%d: slsort_parts: 5. merge and unpack\n", rank));

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts1);

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys);
  pepcparts_elem_set_indices(&d0, irnkl2);
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);

  TSTART(tmergeunpack);
  pepcparts_mergep_heap_unpack_idx(&pd0, &d0, size, frdispls, frcounts);
  TSTOP(tmergeunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5. merge and unpack done\n", rank));

  TSTART(tmakeindices);
  for (i = 0; i < nin; ++i) ++indxl[i];
  for (i = 0; i < *n; ++i) irnkl[irnkl2[i]] = i + 1;
  TSTOP(tmakeindices);

#else

  /* fused merge and unpack of keys (indices created during merge) */
  VERBOSE_MOP(printf("%d: slsort_parts: 5a. merge and unpack indices\n", rank));

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts1);

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
/*  pepcparts_elem_set_keys(&d0, keys);*/
  pepcparts_elem_set_indices(&d0, irnkl2);
/*  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);*/ 

  TSTART(tmergeunpack);
  pepcparts_mergep_heap_unpack_idxonly(&pd0, &d0, size, frdispls, frcounts);
  TSTOP(tmergeunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5a. merge and unpack indices done\n", rank));

  /* unpack */
  VERBOSE_MOP(printf("%d: slsort_parts: 5b. unpack\n", rank));

  for (i = 0; i < *n; ++i) irnkl[irnkl2[i]] = i;

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys);
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);

  TSTART(tunpack);
  pepcparts_elements_unpack_indexed(&pd0, &d0, NULL, irnkl);
  TSTOP(tunpack);

  VERBOSE_MOP(printf("%d: slsort_parts: 5b. unpack data done\n", rank));

  TSTART(tmakeindices);
  for (i = 0; i < nin; ++i) ++indxl[i];
  for (i = 0; i < *n; ++i) ++irnkl[i];
  TSTOP(tmakeindices);

#endif


#ifdef VALIDATE
  o = pepcparts_mpi_elements_validate_order(&d0, 1, size, rank, comm);
 #ifdef VERBOSE
  printf("%d: slsort_parts: global order: %s - local order: %s\n", rank, (!o)?"success":"FAILED", (!l)?"success":"failed");
  printf("%d: %" pepcparts_sl_key_type_fmt " - %" pepcparts_sl_key_type_fmt "\n", rank, d0.keys[0], d0.keys[d0.size - 1]);
 #endif
#endif


  pepcparts_mpi_datatypes_release();

  TSTOP(ttotal);

  VERBOSE_MOP(printf("%d: out: n = %" FINT_TYPE_FMT "\n", rank, *n));

#ifdef BORDER_STATS
  pepcparts_sl_key_type_c lmm[2], gmm[2 * size];
  
  pepckeys_slint_t nb;
  double b, bsum, bmin, bmax;

  if (d0.size > 0)
  {
    lmm[0] = d0.keys[0];
    lmm[1] = d0.keys[d0.size - 1];

  } else lmm[0] = lmm[1] = 0;

  MPI_Gather(lmm, 2 * pepcparts_sl_key_size_mpi, pepcparts_sl_key_type_mpi, gmm, 2 * pepcparts_sl_key_size_mpi, pepcparts_sl_key_type_mpi, 0, comm);
  
  if (rank == 0)
  {
    nb = 0;
    bsum = 0.0;
    bmin = 100.0;
    bmax = 0.0;

/*    for (i = 0; i < size; ++i)
    {
      if (i > 0) printf("%" pepckeys_sl_int_type_fmt "  %llX -> %llX\n", i, gmm[i * 2 + 0], gmm[i * 2 - 1] ^ gmm[i * 2 + 0]);
      else printf("%" pepckeys_sl_int_type_fmt "  %llX\n", i, gmm[i * 2 + 0]);
      printf("%" pepckeys_sl_int_type_fmt "  %llX\n", i, gmm[i * 2 + 1]);
    }*/

/*    printf("%d: borders:\n", rank);*/
    for (i = 0; i < size - 1; ++i)
    {
      if (gmm[i * 2 + 1] > 0 && gmm[(i + 1) * 2] > 0 && (gmm[i * 2 + 1] ^ gmm[(i + 1) * 2]) > 0)
      {
        ++nb;
        b = log((double) (gmm[i * 2 + 1] ^ gmm[(i + 1) * 2])) / log(2.0);
/*        printf("%" pepckeys_sl_int_type_fmt "  %" pepcparts_sl_key_type_fmt "  %f\n", i, gmm[i * 2 + 1] ^ gmm[(i + 1) * 2], b);*/

      } else b = 0;
      
      bsum += b;
      if (b < bmin) bmin = b;
      if (b > bmax) bmax = b;
    }

    printf("%d: borders: avg: %f, min: %f, max: %f\n", rank, (nb)?(bsum / (double) nb):-1.0, bmin, bmax);
  }
#endif

#ifdef TIMING
  if (rank == 0)
  {
# ifndef TIMING_ROW
    printf("%d: slsort_parts: %f\n", rank, ttotal);
    printf("%d: slsort_parts: initindxl: %f\n", rank, tinitindxl);
    printf("%d: slsort_parts: copy: %f\n", rank, tcopy);
    printf("%d: slsort_parts: presort: %f\n", rank, tpresort);
    printf("%d: slsort_parts: partition: %f\n", rank, tpartition);
    printf("%d: slsort_parts: pack: %f\n", rank, tpack);
    printf("%d: slsort_parts: alltoall: %f\n", rank, talltoall);
    printf("%d: slsort_parts: alltoallv: %f\n", rank, talltoallv);
    printf("%d: slsort_parts: mergeunpack: %f\n", rank, tmergeunpack);
#  ifndef MERGE_AND_UNPACK
    printf("%d: slsort_parts: unpack: %f\n", rank, tunpack);
#  endif
    printf("%d: slsort_parts: makeindices: %f\n", rank, tmakeindices);
# else
    printf("%" slint_fmt "  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", ntotal, ttotal, tinitindxl, tcopy, tpresort, tpartition, tpack, talltoall, talltoallv, tmergeunpack, tunpack, tmakeindices);
# endif
  }
#endif
}
