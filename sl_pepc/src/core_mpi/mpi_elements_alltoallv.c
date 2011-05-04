/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_elements_alltoallv.c
 *  timestamp: 2010-12-19 13:24:45 +0100
 *  
 */


#include "sl_common.h"


slint_t mpi_elements_alltoallv_db(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_db */
{
#define xelem_call \
  MPI_Alltoallv(xelem_buf(sbuf), scounts, sdispls, xelem_mpi_datatype, xelem_buf(rbuf), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  return 0;
}


slint_t mpi_elements_alltoallv_ip(elements_t *sbuf, elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_alltoallv_ip */
{
  return 0;
}
