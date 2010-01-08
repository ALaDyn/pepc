
#include <stdio.h>
#include <mpi.h>

#include "sl_pepckeys.h"
#include "sl_pepcparts.h"

#include "fortran2c_types.h"


typedef FINT_TYPE_C finteger_t;
#define finteger_mpi  FINT_TYPE_MPI
#define finteger_fmt  FINT_TYPE_FMT

/*#define MAX_IMBALANCE  0.01*/

/*#define MPI_PARTITION_RADIX_OLD*/

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

#define scounts  fscounts
#define rcounts  frcounts
#define sdispls  fsdispls
#define rdispls  frdispls

#define indxl2  irnkl2

#ifdef VALIDATE
  slint_t o, l;
#endif

#ifdef TIMING
  double ttotal;
  double tinitindxl, tcopy, tpresort, tpartition, tpack, talltoall,talltoallv, tmergeunpack, tmakeindices;
#endif

#ifdef VERBOSE
  printf("%d: slsort_parts with %d processes\n", rank, size);
  printf("%d:  n: %" finteger_fmt "\n", rank, *n);
  printf("%d:  nmax: %" finteger_fmt "\n", rank, *nmax);
  printf("%d:  balance_weight: %" finteger_fmt "\n", rank, *balance_weight);
  printf("%d:  max_imbalance: %f\n", rank, *max_imbalance);

  printf("%d:  sizeof(integer) = %d\n", rank, sizeof(finteger_t));
  printf("%d:  sizeof(integer*8) = %d\n", rank, sizeof(FINT8_TYPE_C));
  printf("%d:  sizeof(pepckeys_key) = %d\n", rank, sizeof(pepckeys_slkey_t));
  printf("%d:  sizeof(pepckeys_index) = %d\n", rank, sizeof(pepckeys_slindex_t));
  printf("%d:  sizeof(pepcparts_key) = %d\n", rank, sizeof(pepcparts_slkey_t));
  printf("%d:  sizeof(pepcparts_index) = %d\n", rank, sizeof(pepcparts_slindex_t));
#endif

  TSTART(ttotal);

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
#ifdef VERBOSE
  printf("%d: slsort_parts: 1. sort keys\n", rank);
#endif

  TSTART(tpresort);
  pepckeys_sort_radix(&k0, NULL, -1, -1, -1);
  TSTOP(tpresort);

#ifdef VERBOSE
  printf("%d: slsort_parts: 1. sort keys done\n", rank);
#endif


#ifdef VALIDATE
  l = pepckeys_elements_validate_order(&k0, 1);
#endif


  /* partitioning */
#ifdef VERBOSE
  printf("%d: slsort_parts: 2. partitioning\n", rank);
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
  pepckeys_mpi_partition_radix_old(&k0, &pc, -1, -1, 3, scounts, sdispls, size, rank, comm);
/*  pepckeys_mpi_partition_radix_old(&k0, &pc, 59, 15, 3, scounts, sdispls, size, rank, comm);*/
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
  pepckeys_mpi_partition_radix(&k0, &pc, -1, -1, 3, scounts, sdispls, size, rank, comm);
/*  pepckeys_mpi_partition_radix(&k0, &pc, 59, 15, 3, scounts, sdispls, size, rank, comm);*/
  TSTOP(tpartition);
#endif


#ifdef VERBOSE
  printf("%d: slsort_parts: 2. partitioning done\n", rank);
#endif


  /* pack */
#ifdef VERBOSE
  printf("%d: slsort_parts: 3. pack\n", rank);
#endif

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

#ifdef VERBOSE
  printf("%d: slsort_parts: 3. pack done\n", rank);
#endif


  /* alltoallv */
#ifdef VERBOSE
  printf("%d: slsort_parts: 4. alltoallv\n", rank);
#endif

  pepcparts_mpi_elements_packed_datatype_create(&pdt, 0);

  TSTART(talltoall);
  MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, comm);
  TSTOP(talltoall);
  for (rdispls[0] = 0, i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
  TSTART(talltoallv);
  MPI_Alltoallv(parts0, scounts, sdispls, pdt, parts1, rcounts, rdispls, pdt, comm);
  TSTOP(talltoallv);

  pepcparts_mpi_elements_packed_datatype_destroy(&pdt);

  *n = rdispls[size - 1] + rcounts[size - 1];
  
#ifdef VERBOSE
  printf("%d: slsort_parts: 4. alltoallv done\n", rank);
#endif


  /* fused merge and unpack (indices during merge) */
#ifdef VERBOSE
  printf("%d: slsort_parts: 5. merge and unpack\n", rank);
#endif

  pepcparts_pelem_set_size(&pd0, *n);
  pepcparts_pelem_set_max_size(&pd0, *nmax);
  pepcparts_pelem_set_elements(&pd0, parts1);

  pepcparts_elem_set_size(&d0, *n);
  pepcparts_elem_set_max_size(&d0, *nmax);
  pepcparts_elem_set_keys(&d0, keys);
  pepcparts_elem_set_indices(&d0, irnkl2);
  pepcparts_elem_set_data(&d0, x, y, z, ux, uy, uz, q, m, work, ex, ey, ez, pelabel);

  TSTART(tmergeunpack);
  pepcparts_mergep_heap_unpack(&pd0, &d0, size, rdispls, rcounts);
  TSTOP(tmergeunpack);

#ifdef VERBOSE
  printf("%d: slsort_parts: 5. merge and unpack done\n", rank);
#endif


  TSTART(tmakeindices);
  for (i = 0; i < nin; ++i) ++indxl[i];
  for (i = 0; i < *n; ++i) irnkl[irnkl2[i]] = i + 1;
  TSTOP(tmakeindices);


#ifdef VALIDATE
  o = pepcparts_mpi_elements_validate_order(&d0, 1, size, rank, comm);
 #ifdef VERBOSE
  printf("%d: slsort_parts: global order: %s - local order: %s\n", rank, (!o)?"success":"FAILED", (!l)?"success":"failed");
 #endif
#endif


  pepcparts_mpi_datatypes_release();

  TSTOP(ttotal);

#ifdef VERBOSE  
  printf("%d: out: n = %" FINT_TYPE_FMT "\n", rank, *n);
#endif

#ifdef TIMING
  if (rank == 0)
  {
    printf("%d: slsort_parts: %f\n", rank, ttotal);
    printf("%d: slsort_parts: initindxl: %f\n", rank, tinitindxl);
    printf("%d: slsort_parts: copy: %f\n", rank, tcopy);
    printf("%d: slsort_parts: presort: %f\n", rank, tpresort);
    printf("%d: slsort_parts: partition: %f\n", rank, tpartition);
    printf("%d: slsort_parts: pack: %f\n", rank, tpack);
    printf("%d: slsort_parts: alltoall: %f\n", rank, talltoall);
    printf("%d: slsort_parts: alltoallv: %f\n", rank, talltoallv);
    printf("%d: slsort_parts: mergeunpack: %f\n", rank, tmergeunpack);
    printf("%d: slsort_parts: makeindices: %f\n", rank, tmakeindices);
  }
#endif
}
