/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_linsplit.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


#include "sl_common.h"


#define DEBUG_OR_NOT  (sl_mpi_rank == -1)


slint_t mpi_linsplit(MPI_Comm comm_in, slkey_pure_t *keys_in, MPI_Comm *comms_out, slint_t *parity, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_linsplit */
{
  slint_t local_parity;


  comms_out[0] = comms_out[1] = MPI_COMM_NULL;

  Z_TRACE_IF(DEBUG_OR_NOT, "comm_in = %" zcomm_fmt "\n", zcomm_val(comm_in));

  if (comm_in == MPI_COMM_NULL) return 0;

  local_parity = key_pure_cmp_ne(keys_in[0], keys_in[1]);

  *parity = 0;
  MPI_Exscan(&local_parity, parity, 1, int_mpi_datatype, MPI_SUM, comm_in);

  *parity = (*parity % 2);

  if (local_parity)
  {
    MPI_Comm_split(comm_in, keys_in[*parity], rank, &comms_out[*parity]);
    MPI_Comm_split(comm_in, keys_in[!(*parity)], rank, &comms_out[!(*parity)]);

  } else
  {
    MPI_Comm_split(comm_in, (*parity == 0)?keys_in[0]:MPI_UNDEFINED, rank, &comms_out[0]);
    MPI_Comm_split(comm_in, (*parity == 1)?keys_in[1]:MPI_UNDEFINED, rank, &comms_out[1]);
  }

  return 0;
}


#ifdef key_integer

slint_t mpi_linsplit_radix(slkey_pure_t klow, slkey_pure_t khigh, MPI_Comm *comm0, MPI_Comm *comm1, int size, int rank, MPI_Comm comm)
{
  slint_t parity, parity_prefix;
  int comm0_size, comm1_size;
  
  parity = (khigh != klow);
  
  parity_prefix = 0;
  MPI_Exscan(&parity, &parity_prefix, 1, int_mpi_datatype, MPI_SUM, comm);

  MPI_Comm_split(comm, (parity_prefix % 2)?klow:khigh, rank, comm0);

  MPI_Comm_split(comm, (parity_prefix % 2)?khigh:klow, rank, comm1);

  MPI_Comm_size(*comm0, &comm0_size);
  MPI_Comm_size(*comm1, &comm1_size);

  if (klow == khigh)
  {
    if (comm0_size < comm1_size) MPI_Comm_free(comm0);
    else MPI_Comm_free(comm1);
  }

  return 0;
}

#endif /* key_integer */


slint_t mpi_linsplit2(MPI_Comm comm_in, slkey_pure_t *keys_in, MPI_Comm *comms_out, slint_t *parity, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_linsplit2 */
{
  slint_t local_parity;


  comms_out[0] = comms_out[1] = MPI_COMM_NULL;

  Z_TRACE_IF(DEBUG_OR_NOT, "comm_in = %" zcomm_fmt "\n", zcomm_val(comm_in));

  if (comm_in == MPI_COMM_NULL) return 0;

  local_parity = key_pure_cmp_ne(keys_in[0], keys_in[1]);

  parity[0] = 0;
  MPI_Exscan(&local_parity, &parity[0], 1, int_mpi_datatype, MPI_SUM, comm_in);

  parity[0] = (parity[0] % 2);

  MPI_Comm_split(comm_in, (parity[0] == 0 || local_parity)?keys_in[(parity[0])?1:0]:MPI_UNDEFINED, rank, &comms_out[0]);
  MPI_Comm_split(comm_in, (parity[0] == 1 || local_parity)?keys_in[(parity[0])?0:1]:MPI_UNDEFINED, rank, &comms_out[1]);

  return 0;
}
