/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_select_exact_radix.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


#include "sl_common.h"


slint_t mpi_select_exact_radix(elements_t *s, slint_t nelements, slint_t nparts, partcond_t *pconds, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_radix */
{
  binning_t bm;

  binning_radix_create(&bm, rhigh, rlow, rwidth, sorted|SL_SORTED_IN);

#if 1
  splitter_t sp;
  sp.displs = sdispls;

  mpi_select_exact_generic_bulk(s, nelements, nparts, pconds, &bm, &sp, size, rank, comm);
#else
  mpi_select_exact_generic_bulk_old(s, nelements, nparts, pconds, &bm, sdispls, size, rank, comm);
#endif

  binning_radix_destroy(&bm);
  
  return 0;
}


slint_t mpi_select_exact_radix_grouped(elements_t *s, slint_t nelements, partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, slint_t rhigh, slint_t rlow, slint_t rwidth, slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_exact_radix_grouped */
{
  slint_t npconds = -1;
  partcond_t *pconds;


  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, NULL, &npconds, size, rank, comm);

  pconds = z_alloca(npconds, sizeof(partcond_t));
  
  mpi_gather_partconds_grouped(pcond, pcond_comm, group_comm, pconds, &npconds, size, rank, comm);

  mpi_select_exact_radix(s, nelements, npconds, pconds, rhigh, rlow, rwidth, sorted, sdispls, size, rank, comm);
  
  z_freea(pconds);
  
  return 0;
}
