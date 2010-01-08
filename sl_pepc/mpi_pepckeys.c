
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "sl_pepckeys.h"

#include "fortran2c_types.h"


typedef FINT_TYPE_C finteger_t;
#define finteger_mpi  FINT_TYPE_MPI
#define finteger_fmt  FINT_TYPE_FMT

/*#define MAX_IMBALANCE  0.01*/

/*#define MPI_PARTITION_RADIX_OLD*/

/*#define SORT_INSTEAD_OF_MERGE*/

/*#define VERBOSE
#define VALIDATE
#define TIMING*/


#ifdef TIMING
# define TSTART(tid)  tid = MPI_Wtime()
# define TSTOP(tid)   tid = MPI_Wtime() - tid
#else
# define TSTART(tid)  do { } while (0)
# define TSTOP(tid)   do { } while (0)
#endif

extern int pepckeys_sl_mpi_rank;


void slsort_keys_(finteger_t *, finteger_t *, pepckeys_slkey_t *, pepckeys_sldata0_t *, finteger_t *, double *, finteger_t *, finteger_t *, finteger_t *,
                  finteger_t *, finteger_t *, finteger_t *, finteger_t *, pepckeys_slkey_t *, finteger_t *, finteger_t *, finteger_t *);

#pragma weak slsort_keys_ = slsort_keys
void slsort_keys(finteger_t *nin,                                       /* IN */
                 finteger_t *nmax,                                      /* IN */
                 pepckeys_slkey_t *keys,                                /* INOUT */
                 pepckeys_sldata0_t *workload,                          /* INOUT */
                 finteger_t *balance_weight,                            /* IN */
                 double *max_imbalance,                                 /* IN */
                 finteger_t *nout,                                      /* OUT */
                 finteger_t *indxl, finteger_t *irnkl,                  /* OUT */
                 finteger_t *fscounts, finteger_t *frcounts,            /* OUT */
                 finteger_t *fsdispls, finteger_t *frdispls,            /* OUT */
                 pepckeys_slkey_t *keys2,                               /* SCRATCH */
                 finteger_t *irnkl2,                                    /* SCRATCH */
                 finteger_t *fsize, finteger_t *frank)                  /* IN */
{
  int size = *fsize;
  int rank = *frank;
  MPI_Comm comm = MPI_COMM_WORLD;

  pepckeys_elements_t s0, s1;
  pepckeys_slint_t i;
  pepckeys_partcond_t pc;

#ifdef MAX_IMBALANCE
  double imba = MAX_IMBALANCE;
#else
  double imba = *max_imbalance;
#endif

#define scounts  fscounts
#define rcounts  frcounts
#define sdispls  fsdispls
#define rdispls  frdispls

#ifdef VALIDATE
  pepckeys_slint_t o, l;
#endif

#ifdef TIMING
  double ttotal, tinitpre, tpresort, tpartition, talltoall, talltoallv, tinitpost, tpostmerge;
#endif

  TSTART(ttotal);

#ifdef VERBOSE
  printf("%d: slsort_keys with %d processes\n", rank, size);
  printf("%d:  nin: %" finteger_fmt "\n", rank, *nin);
  printf("%d:  nmax: %" finteger_fmt "\n", rank, *nmax);
  printf("%d:  balance_weight: %" finteger_fmt "\n", rank, *balance_weight);
  printf("%d:  max_imbalance: %f\n", rank, *max_imbalance);

  printf("%d:  sizeof(integer) = %d\n", rank, sizeof(finteger_t));
  printf("%d:  sizeof(key) = %d\n", rank, sizeof(pepckeys_slkey_t));
  printf("%d:  sizeof(index) = %d\n", rank, sizeof(pepckeys_slindex_t));
#endif

  pepckeys_sl_mpi_rank = rank;

  pepckeys_mpi_datatypes_init();

  TSTART(tinitpre);

  /* init pre indexes */
  for (i = 0; i < *nin; ++i) indxl[i] = i + 1;

  s0.size = *nin;
  s0.max_size = *nmax;
  s0.keys = keys;
  s0.indices = indxl;
  s0.data0 = workload;

  TSTOP(tinitpre);

  /* pre sort local */
#ifdef VERBOSE
  printf("%d: slsort_keys: 1. sort keys\n", rank);
#endif

  TSTART(tpresort);
  pepckeys_sort_radix(&s0, NULL, -1, -1, -1);
  TSTOP(tpresort);

#ifdef VERBOSE
  printf("%d: slsort_keys: 1. sort keys done\n", rank);
#endif

#ifdef VALIDATE
  l = pepckeys_elements_validate_order(&s0, 1);
#endif


  /* partitioning */
#ifdef VERBOSE
  printf("%d: slsort_keys: 2. partitioning\n", rank);
#endif

  /* no fixed min/max borders for this partition */
  pc.min_cpart = 0.0;
  pc.max_cpart = -1.0;
  pc.min_wpart = 0.0;
  pc.max_wpart = -1.0;

#ifdef MPI_PARTITION_RADIX_OLD
  pc.weighted = (*balance_weight != 0);
  pc.min_count = -(1.0 - imba);
  pc.max_count = -(1.0 + imba);
  pc.min_weight = -(1.0 - imba);
  pc.max_weight = -(1.0 + imba);

  TSTART(tpartition);
  pepckeys_mpi_partition_radix_old(&s0, &pc, -1, -1, 3, scounts, sdispls, size, rank, comm);
  TSTOP(tpartition);
#else
  if (*balance_weight)
  {
    pc.min_count = 0;
    pc.max_count = *nmax;
    pc.min_weight = -(1.0 - imba);
    pc.max_weight = -(1.0 + imba);

  } else
  {
    pc.min_count = -(1.0 - imba);
    pc.max_count = -(1.0 + imba);
    pc.min_weight = 0;
    pc.max_weight = -size; /* max = size x avg. = total */

    if (pc.max_count > *nmax) pc.max_count = *nmax;
  }

  TSTART(tpartition);
  pepckeys_mpi_partition_radix(&s0, &pc, -1, -1, 3, scounts, sdispls, size, rank, comm);
  TSTOP(tpartition);
#endif


#ifdef VERBOSE
  printf("%d: slsort_keys: 2. partitioning done\n", rank);
#endif
  
  /* alltoallv keys */
#ifdef VERBOSE
  printf("%d: slsort_keys: 3. alltoallv\n", rank);
#endif

  TSTART(talltoall);
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  TSTOP(talltoall);
  for (rdispls[0] = 0, i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  TSTART(talltoallv);
  MPI_Alltoallv(keys, scounts, sdispls, pepckeys_sl_key_type_mpi, keys2, rcounts, rdispls, pepckeys_sl_key_type_mpi, comm);
  TSTOP(talltoallv);

#ifdef VERBOSE
  printf("%d: slsort_keys: 3. alltoallv done\n", rank);
#endif

  *nout = rdispls[size - 1] + rcounts[size - 1];

  /* post sort local (or mergep) */
#ifdef VERBOSE
  printf("%d: slsort_keys: 4. merge keys\n", rank);
#endif

#ifdef SORT_INSTEAD_OF_MERGE

  TSTART(tinitpost);

  /* init post indexes */
  for (i = 0; i < *nout; ++i) irnkl2[i] = i;

  s0.size = *nout;
  s0.max_size = *nmax;
  s0.keys = keys2;
  s0.indices = irnkl2;
  s0.data0 = workload;

  s1.size = *nout;
  s1.max_size = *nmax;
  s1.keys = keys2;
  s1.indices = irnkl2;
  s1.data0 = workload;

  TSTOP(tinitpost);

  TSTART(tpostmerge);
  pepckeys_sort_radix(&s0, NULL, -1, -1, -1);
  TSTOP(tpostmerge);

#else

  TSTART(tinitpost);

  /* init post indexes */
  for (i = 0; i < *nout; ++i) irnkl[i] = i;

  memcpy(keys, keys2, *nout * sizeof(pepckeys_slkey_t));

  s0.size = *nout;
  s0.max_size = *nmax;
  s0.keys = keys;
  s0.indices = irnkl;
  s0.data0 = workload;

  s1.size = *nout;
  s1.max_size = *nmax;
  s1.keys = keys2;
  s1.indices = irnkl2;
  s1.data0 = workload;

  TSTOP(tinitpost);

  TSTART(tpostmerge);
  pepckeys_mergep_heap(&s0, &s1, size, rdispls, rcounts);
  TSTOP(tpostmerge);

#endif
  
  for (i = 0; i < s0.size; ++i) irnkl[irnkl2[i]] = i + 1;

#ifdef VERBOSE
  printf("%d: slsort_keys: 4. merge keys done\n", rank);
#endif

#ifdef VALIDATE
  o = pepckeys_mpi_elements_validate_order(&s1, 1, size, rank, comm);
 #ifdef VERBOSE
  printf("%d: slsort_keys: global order: %s - local order: %s\n", rank, (!o)?"success":"FAILED", (!l)?"success":"failed");
 #endif
#endif

  pepckeys_mpi_datatypes_release();

  TSTOP(ttotal);

#ifdef VERBOSE
  printf("%d: nout: %" FINT_TYPE_FMT "\n", rank, *nout);
#endif

#ifdef TIMING
  if (rank == 0)
  {
#ifdef MPI_PARTITION_RADIX_OLD
    printf("%d: slsort_keys: %f (old)\n", rank, ttotal);
#else
    printf("%d: slsort_keys: %f\n", rank, ttotal);
#endif
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
