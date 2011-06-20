/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_rebalance.c
 *  
 */


#include "sl_common.h"


/*#define LOCAL_NCOPY*/
#define NONBLOCKING

#define MPI_REBALANCE_ALLTOALLV_TAG  1


slint_t mpi_rebalance(elements_t *s0, elements_t *s1, slint_t stable, slint_t *dst_size, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_rebalance */
{
  int exit_code = 0;
  int i;
  slint_t s, t, u, v;
  slint_t *sizes, local_sizes[2];
  int *sendcounts, *senddispls, *recvcounts, *recvdispls;
  slint_t send_rank, send_size, recv_rank, recv_size;

  
  sizes = z_alloc(size * 2, sizeof(slint_t));

  local_sizes[0] = s0->size;
  if (dst_size) local_sizes[1] = *dst_size; else local_sizes[1] = s0->size;

  MPI_Allgather(local_sizes, 2, int_mpi_datatype, sizes, 2, int_mpi_datatype, comm);

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
      v = z_min(-s, sizes[i * 2 + 1]);
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
    fprintf(stderr, "%d: mpi_rebalance: error: destination sizes too small (%" sl_int_type_fmt " > 0)\n", rank, s);
    exit_code = -1;
    goto free_and_exit;
  }

/*  printf("%d here: sizes = [", rank);
  for (i = 0; i < size; ++i) printf(" %" sl_int_type_fmt ":%" sl_int_type_fmt " ", sizes[i * 2], sizes[i * 2 + 1]);
  printf("]\n");*/

  sendcounts = z_alloc(4 * size, sizeof(int));
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

    s = z_min(send_size, recv_size);

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
    sendcounts[rank] = recvcounts[rank] = z_min(sizes[rank * 2], sizes[rank * 2 + 1]);
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

  } else /* in-place */
  {
    /* unstable case: local elements don't need to move  */
    if (!stable) sendcounts[rank] = recvcounts[rank] = 0;
  }
  
  mpi_rebalance_alltoallv(s0, sendcounts, senddispls, s1, recvcounts, recvdispls, size, rank, comm);


  if (dst_size) *dst_size = sizes[rank * 2 + 1];

  if (s1) s1->size = sizes[rank * 2 + 1];
  else s0->size = sizes[rank * 2 + 1];

  z_free(sendcounts);

free_and_exit:
  z_free(sizes);

  return exit_code;
}


slint_t mpi_rebalance_alltoallv(elements_t *sbuf, int *scounts, int *sdispls, elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_rebalance_alltoallv */
{
  slint_t i;

#ifdef NONBLOCKING
  slint_t j, nreqs;
  const slint_t request_per_element = 1 + data_n;
  MPI_Request reqs[(size - 1) * request_per_element];
  MPI_Status stats[(size - 1) * request_per_element];
#else
  MPI_Status status;
#endif
  
  if (sbuf != NULL) /* out-of-place */
  {
#ifdef LOCAL_NCOPY
    if (scounts[rank] > 0) elem_ncopy_at(sbuf, sdispls[rank], rbuf, rdispls[rank], scounts[rank]);
    scounts[rank] = rcounts[rank] = 0;
#endif

#define xelem_call \
    MPI_Alltoallv(xelem_buf(sbuf), scounts, sdispls, xelem_mpi_datatype, xelem_buf(rbuf), rcounts, rdispls, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

  } else /* in-place */
  {
    sbuf = rbuf;

#ifdef NONBLOCKING
    /* send */
    nreqs = 0;
    for (i = 1; i < size; ++i)
    {
      j = (rank + i) % size;
      if (scounts[j] > 0)
      {
#define xelem_call \
        MPI_Isend(xelem_buf_at(sbuf, sdispls[j]), scounts[j], xelem_mpi_datatype, j, MPI_REBALANCE_ALLTOALLV_TAG, comm, &reqs[nreqs]); ++nreqs;
#include "sl_xelem_call.h"
      }
    }
    MPI_Waitall(nreqs, reqs, stats);
    
    /* local move */
    if (scounts[rank] > 0) elem_nmove_at(sbuf, sdispls[rank], sbuf, rdispls[rank], scounts[rank]);
    
    /* receive */
    nreqs = 0;
    for (i = 1; i < size; ++i)
    {
      j = (rank - i + size) % size;
      if (rcounts[j] > 0)
      {
#define xelem_call \
        MPI_Irecv(xelem_buf_at(rbuf, rdispls[j]), rcounts[j], xelem_mpi_datatype, j, MPI_REBALANCE_ALLTOALLV_TAG, comm, &reqs[nreqs]); ++nreqs;
#include "sl_xelem_call.h"
      }
    }
    MPI_Waitall(nreqs, reqs, stats);

#else /* NONBLOCKING */

    /* send right */
    for (i = rank + 1; i < size; ++i)
    if (scounts[i] > 0)
    {
#define xelem_call \
      MPI_Send(xelem_buf_at(sbuf, sdispls[i]), scounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm);
#include "sl_xelem_call.h"
    }

    /* move right (local) */
    if (scounts[rank] > 0 && sdispls[rank] < rdispls[rank]) elem_nmove_at(sbuf, sdispls[rank], sbuf, rdispls[rank], scounts[rank]);
    
    /* receive left */
    for (i = rank - 1; i >= 0; --i)
    if (rcounts[i] > 0)
    {
#define xelem_call \
      MPI_Recv(xelem_buf_at(sbuf, rdispls[i]), rcounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm, &status);
#include "sl_xelem_call.h"
    }

    /* send left */
    for (i = rank - 1; i >= 0; --i)
    if (scounts[i] > 0)
    {
#define xelem_call \
      MPI_Send(xelem_buf_at(sbuf, sdispls[i]), scounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm);
#include "sl_xelem_call.h"
    }

    /* move left (local) */
    if (scounts[rank] > 0 && sdispls[rank] > rdispls[rank]) elem_nmove_at(sbuf, sdispls[rank], sbuf, rdispls[rank], scounts[rank]);

    /* receive right */
    for (i = rank + 1; i < size; ++i)
    if (rcounts[i] > 0)
    {
#define xelem_call \
      MPI_Recv(xelem_buf_at(sbuf, rdispls[i]), rcounts[i], xelem_mpi_datatype, i, MPI_REBALANCE_ALLTOALLV_TAG, comm, &status);
#include "sl_xelem_call.h"
    }

#endif /* NONBLOCKING */
  }

  return 0;
}


#undef LOCAL_NCOPY
#undef NONBLOCKING
#undef MPI_REBALANCE_ALLTOALLV_TAG
