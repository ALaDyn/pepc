/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_rebalance.c
 *  timestamp: 2009-11-12 08:41:36 +0100
 *  
 */


#include "sl_common.h"


slint_t pepckeys_mpi_rebalance(elements_t *s0, elements_t *s1, slint_t stable, slint_t *dst_size, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_rebalance */
{
  const int tag = 1;
  int exit_code = 0;
  int i;
  slint_t s, t, u, v;
  slint_t *sizes, local_sizes[2];
  int *sendcounts, *senddispls, *recvcounts, *recvdispls;
  slint_t send_rank, send_size, recv_rank, recv_size;
  MPI_Status status;

  
  sizes = sl_alloc(size * 2, sizeof(slint_t));

  local_sizes[0] = s0->size;
  if (dst_size) local_sizes[1] = *dst_size; else local_sizes[1] = s0->size;

  MPI_Allgather(local_sizes, 2, sl_int_type_mpi, sizes, 2, sl_int_type_mpi, comm);

  s = t = 0;
  for (i = 0; i < size; ++i)
  {
    s += sizes[i * 2];
    if (sizes[i * 2 + 1] < 0) ++t; else s -= sizes[i * 2 + 1];
  }

  /* correct oversized 'sizes' */
  if (s < 0)
  {
    for (i = size - 1; i >= 0; --i)
    if (sizes[i * 2 + 1] > 0)
    {
      v = xmin(-s, sizes[i * 2 + 1]);
      sizes[i * 2 + 1] -= v;
      s += v;
    }
  }

  /* FIXME: APP_ASSERT(s >= 0); */

  /* correct autosize 'sizes' */
  if (t > 0)
  {
    u = s % t;
    for (i = 0; i < size; ++i)
    if (sizes[i * 2 + 1] < 0)
    {
      sizes[i * 2 + 1] = s / t;
      if (u > 0)
      {
        ++sizes[i * 2 + 1];
        --u;
      }
    }

  } else if (s > 0)
  {
    fprintf(stderr, "%d: pepckeys_mpi_rebalance: error: destination sizes too small (%" sl_int_type_fmt " > 0)\n", rank, s);
    exit_code = -1;
    goto free_and_exit;
  }

/*  printf("%d here: sizes = [", rank);
  for (i = 0; i < size; ++i) printf(" %" sl_int_type_fmt ":%" sl_int_type_fmt " ", sizes[i * 2], sizes[i * 2 + 1]);
  printf("]\n");*/

  sendcounts = sl_alloc(4 * size, sizeof(int));
  senddispls = sendcounts + 1 * size;
  recvcounts = sendcounts + 2 * size;
  recvdispls = sendcounts + 3 * size;

  memset(sendcounts, 0, size * sizeof(int));
  memset(recvcounts, 0, size * sizeof(int));

  send_rank = 0;
  send_size = -1;
  recv_rank = 0;
  recv_size = -1;

  while (send_rank < size && recv_rank < size)
  {
    if (send_size < 0)
    {
      send_size = sizes[send_rank * 2];
      if (!stable) send_size -= sizes[send_rank * 2 + 1];
      if (send_size < 0) send_size = 0;
    }

    if (recv_size < 0)
    {
      recv_size = sizes[recv_rank * 2 + 1];
      if (!stable) recv_size -= sizes[recv_rank * 2];
      if (recv_size < 0) recv_size = 0;
    }

    s = xmin(send_size, recv_size);

    if (send_rank == rank) sendcounts[recv_rank] += s;
    if (recv_rank == rank) recvcounts[send_rank] += s;

    send_size -= s;
    recv_size -= s;

    if (send_size <= 0)
    {
      ++send_rank;
      send_size = -1;
    }

    if (recv_size <= 0)
    {
      ++recv_rank;
      recv_size = -1;
    }
  }
    
  /* FIXME: APP_ASSERT(send_rank >= size && recv_rank >= size); */

  if (stable) senddispls[0] = recvdispls[0] = 0;
  else
  {
    senddispls[0] = sizes[rank * 2 + 1];
    recvdispls[0] = sizes[rank * 2];
  }
  for (i = 1; i < size; ++i)
  {
    senddispls[i] = senddispls[i - 1] + sendcounts[i - 1];
    recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
  }
  if (!stable)
  {
    sendcounts[rank] = recvcounts[rank] = xmin(sizes[rank * 2], sizes[rank * 2 + 1]);
    senddispls[rank] = recvdispls[rank] = 0;
  }

/*  printf("%d here: sendcounts = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", sendcounts[i]);
  printf("]\n");

  printf("%d here: senddispls = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", senddispls[i]);
  printf("]\n");

  printf("%d here: recvcounts = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", recvcounts[i]);
  printf("]\n");

  printf("%d here: recvdispls = [", rank);
  for (i = 0; i < size; ++i) printf(" %d ", recvdispls[i]);
  printf("]\n");*/

  
  if (s1) /* out-of-place */
  {
#ifdef LOCAL_NCOPY
    if (!stable)
    {
      elem_ncopy_at(s0, senddispls[rank], s1, recvdispls[rank], sendcounts[rank]);
      sendcounts[rank] = recvcounts[rank] = 0;
    }
#endif

#define xelem_call \
    MPI_Alltoallv(xelem_buf(s0), sendcounts, senddispls, xelem_mpi_datatype, xelem_buf(s1), recvcounts, recvdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  } else /* in-place */
  {
    /* unstable case: local elements don't need to move  */
    if (!stable) sendcounts[rank] = recvcounts[rank] = 0;

    /* send/recv graph is a dag! (next exchange with point-to-points is deadlock free!?!) */
    
    /* send right */
    for (i = rank + 1; i < size; ++i)
    if (sendcounts[i] > 0)
    {
#define xelem_call \
      MPI_Send(xelem_buf_at(s0, senddispls[i]), sendcounts[i], xelem_mpi_datatype, i, tag, comm);
#include "sl_xelem_call.h"
    }

    /* move right (local) */
    if (sendcounts[rank] > 0 && senddispls[rank] < recvdispls[rank]) elem_nmove_at(s0, senddispls[rank], s0, recvdispls[rank], sendcounts[rank]);
    
    /* receive left */
    for (i = rank - 1; i >= 0; --i)
    if (recvcounts[i] > 0)
    {
#define xelem_call \
      MPI_Recv(xelem_buf_at(s0, recvdispls[i]), recvcounts[i], xelem_mpi_datatype, i, tag, comm, &status);
#include "sl_xelem_call.h"
    }

    /* send left */
    for (i = rank - 1; i >= 0; --i)
    if (sendcounts[i] > 0)
    {
#define xelem_call \
      MPI_Send(xelem_buf_at(s0, senddispls[i]), sendcounts[i], xelem_mpi_datatype, i, tag, comm);
#include "sl_xelem_call.h"
    }

    /* move left (local) */
    if (sendcounts[rank] > 0 && senddispls[rank] > recvdispls[rank]) elem_nmove_at(s0, senddispls[rank], s0, recvdispls[rank], sendcounts[rank]);

    /* receive right */
    for (i = rank + 1; i < size; ++i)
    if (recvcounts[i] > 0)
    {
#define xelem_call \
      MPI_Recv(xelem_buf_at(s0, recvdispls[i]), recvcounts[i], xelem_mpi_datatype, i, tag, comm, &status);
#include "sl_xelem_call.h"
    }
  }

  if (dst_size) *dst_size = sizes[rank * 2 + 1];

  if (s1) s1->size = sizes[rank * 2 + 1];
  else s0->size = sizes[rank * 2 + 1];

  sl_free(sendcounts);

free_and_exit:
  sl_free(sizes);

  return exit_code;
}
