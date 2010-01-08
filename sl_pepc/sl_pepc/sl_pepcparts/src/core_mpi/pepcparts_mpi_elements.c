/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_elements.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


slint pepcparts_mpi_elements_init_keys_from_file(elements_t *s, char *filename, slint from, slint to, slint const_bytes_per_line, slint root, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_elements_init_keys_from_file */
{
  slint i, m, n, p, q;

  elements_t e;
  
  int send_counts[size], send_dipls[size];

  n = to - from + 1;

  if (root < 0 || size == 1)
  {
    m = (n / size) + ((n % size) > rank);
    p = ((int) (n / size)) * rank + xmin(rank, n % size);

    q = pepcparts_elements_init_keys_from_file(s, 1, filename, p, p + m - 1, const_bytes_per_line);

  } else
  {
    if (rank == root)
    {
      pepcparts_elements_init_keys_from_file(&e, 0, filename, from, to, const_bytes_per_line);

      for (i = 0; i < size; i++)
      {
        send_counts[i] = (n / size) + ((n % size) > i);
        send_counts[i] *= key_size_mpi;

        send_dipls[i] = ((int) (n / size)) * i + xmin(i, n % size);
        send_dipls[i] *= key_size_mpi;
      }
    }

    m = (n / size) + ((n % size) > rank);

    pepcparts_elements_alloc(s, m, 1, 1);

    MPI_Scatterv(e.keys, send_counts, send_dipls, key_type_mpi, s->keys, m, key_type_mpi, root, comm);
    
    q = m;
    
    if (rank == root) pepcparts_elements_free(&e);
  }

  return q;
}


slint pepcparts_mpi_elements_validate_order(elements_t *s, slint n, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_elements_validate_order */
{
  slint local_result = 0, global_result;
  slkey_pure_t pure_keys[2];
  MPI_Status status;

  if (size < 0) MPI_Comm_size(comm, &size);
  if (rank < 0) MPI_Comm_rank(comm, &rank);

  if (size > 1)
  {
    /* send lowest key to the left neighbor */
    pure_keys[0] = key_purify(s[0].keys[0]);

    if (rank == 0) MPI_Recv(&pure_keys[1], 1, pepcparts_key_pure_mpi_datatype, rank + 1, 0, comm, &status);
    else if (rank + 1 == size) MPI_Send(&pure_keys[0], 1, pepcparts_key_pure_mpi_datatype, rank - 1, 0, comm);
    else MPI_Sendrecv(&pure_keys[0], 1, pepcparts_key_pure_mpi_datatype, rank - 1, 0, &pure_keys[1], 1, pepcparts_key_pure_mpi_datatype, rank + 1, 0, comm, &status);

    if (rank + 1 < size) local_result += (key_pure_cmp_gt(key_purify(s[n - 1].keys[s[n - 1].size - 1]), pure_keys[1]) != 0);

    /* send highest key to the right neighbor */
    pure_keys[0] = key_purify(s[n - 1].keys[s[n - 1].size - 1]);

    if (rank == 0) MPI_Send(&pure_keys[0], 1, pepcparts_key_pure_mpi_datatype, rank + 1, 0, comm);
    else if (rank + 1 == size) MPI_Recv(&pure_keys[1], 1, pepcparts_key_pure_mpi_datatype, rank - 1, 0, comm, &status);
    else MPI_Sendrecv(&pure_keys[0], 1, pepcparts_key_pure_mpi_datatype, rank + 1, 0, &pure_keys[1], 1, pepcparts_key_pure_mpi_datatype, rank - 1, 0, comm, &status);

    if (rank > 0) local_result += (key_pure_cmp_lt(key_purify(s[0].keys[0]), pure_keys[1]) != 0);

    /* reduce the local results of the validation between neighbor-processes */
    MPI_Allreduce(&local_result, &global_result, 1, sl_int_type_mpi, MPI_MAX, comm);
    /* exit if the validation fails */
    if (global_result > 0) return 1;
  }

  /* start the process-internal validation */
  local_result = (pepcparts_elements_validate_order(s, n) != 0);
  /* reduce the local results */
  MPI_Allreduce(&local_result, &global_result, 1, sl_int_type_mpi, MPI_MAX, comm);

  return global_result;
}


unsigned short pepcparts_mpi_cs16(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_cs16 */
{
  unsigned short crc_local, crc_global = 0;

  crc_local = pepcparts_cs_crc16(s, n, keys, data);

  MPI_Allreduce(&crc_local, &crc_global, 1, MPI_UNSIGNED_SHORT, MPI_SUM, comm);

  return crc_global;
}


unsigned int pepcparts_mpi_cs32(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_cs32 */
{
  unsigned int crc_local, crc_global = 0;

  crc_local = pepcparts_cs_crc32(s, n, keys, data);

  MPI_Allreduce(&crc_local, &crc_global, 1, MPI_UNSIGNED, MPI_SUM, comm);

  return crc_global;
}
