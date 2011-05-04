/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_elements.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


/* sl_macro ME_TRACE_IF */

#include "sl_common.h"


void *me_sendrecv_replace_mem = NULL;  /* sl_global, sl_var me_sendrecv_replace_mem */
slint_t me_sendrecv_replace_memsize = 0;  /* sl_global, sl_var me_sendrecv_replace_memsize */
slint_t me_sendrecv_replace_mpi_maxsize = INT_MAX / 2;  /* sl_global, sl_var me_sendrecv_replace_mpi_maxsize */

#ifndef ME_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define ME_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define ME_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint mpi_elements_init_keys_from_file(elements_t *s, char *filename, slint from, slint to, slint const_bytes_per_line, slint root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_init_keys_from_file */
{
  slint i, m, n, p, q;

  elements_t e;
  
  int send_counts[size], send_dipls[size];

  n = to - from + 1;

  if (root < 0 || size == 1)
  {
    m = (n / size) + ((n % size) > rank);
    p = ((int) (n / size)) * rank + z_min(rank, n % size);

    q = elements_init_keys_from_file(s, 1, filename, p, p + m - 1, const_bytes_per_line);

  } else
  {
    if (rank == root)
    {
      elements_init_keys_from_file(&e, 0, filename, from, to, const_bytes_per_line);

      for (i = 0; i < size; i++)
      {
        send_counts[i] = (n / size) + ((n % size) > i);
        send_counts[i] *= key_size_mpi;

        send_dipls[i] = ((int) (n / size)) * i + z_min(i, n % size);
        send_dipls[i] *= key_size_mpi;
      }
    }

    m = (n / size) + ((n % size) > rank);

    elements_alloc(s, m, SLCM_ALL);

    MPI_Scatterv(e.keys, send_counts, send_dipls, key_type_mpi, s->keys, m, key_type_mpi, root, comm);
    
    q = m;
    
    if (rank == root) elements_free(&e);
  }

  return q;
}


slint mpi_elements_validate_order(elements_t *s, slint n, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_validate_order */
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

    if (rank == 0) MPI_Recv(&pure_keys[1], 1, pkey_mpi_datatype, rank + 1, 0, comm, &status);
    else if (rank + 1 == size) MPI_Send(&pure_keys[0], 1, pkey_mpi_datatype, rank - 1, 0, comm);
    else MPI_Sendrecv(&pure_keys[0], 1, pkey_mpi_datatype, rank - 1, 0, &pure_keys[1], 1, pkey_mpi_datatype, rank + 1, 0, comm, &status);

    if (rank + 1 < size) local_result += (key_pure_cmp_gt(key_purify(s[n - 1].keys[s[n - 1].size - 1]), pure_keys[1]) != 0);

    /* send highest key to the right neighbor */
    pure_keys[0] = key_purify(s[n - 1].keys[s[n - 1].size - 1]);

    if (rank == 0) MPI_Send(&pure_keys[0], 1, pkey_mpi_datatype, rank + 1, 0, comm);
    else if (rank + 1 == size) MPI_Recv(&pure_keys[1], 1, pkey_mpi_datatype, rank - 1, 0, comm, &status);
    else MPI_Sendrecv(&pure_keys[0], 1, pkey_mpi_datatype, rank + 1, 0, &pure_keys[1], 1, pkey_mpi_datatype, rank - 1, 0, comm, &status);

    if (rank > 0) local_result += (key_pure_cmp_lt(key_purify(s[0].keys[0]), pure_keys[1]) != 0);

    /* reduce the local results of the validation between neighbor-processes */
    MPI_Allreduce(&local_result, &global_result, 1, int_mpi_datatype, MPI_MAX, comm);
    /* exit if the validation fails */
    if (global_result > 0) return 1;
  }

  /* start the process-internal validation */
  local_result = (elements_validate_order(s, n) != 0);
  /* reduce the local results */
  MPI_Allreduce(&local_result, &global_result, 1, int_mpi_datatype, MPI_MAX, comm);

  return global_result;
}


slint_t mpi_linear_exchange_pure_keys(slkey_pure_t *in, slkey_pure_t *out, int size, int rank, MPI_Comm comm)  /* sl_proto, sl_func mpi_linear_exchange_pure_keys */
{
  MPI_Status status;

  out[0] = in[0];
  out[1] = in[1];

  /* exchange to the left */
  MPI_Sendrecv(&in[0], 1, pkey_mpi_datatype, (rank - 1 >= 0)?(rank - 1):MPI_PROC_NULL, 0,
               &out[1], 1, pkey_mpi_datatype, (rank + 1 < size)?(rank + 1):MPI_PROC_NULL, 0,
               comm, &status);

  /* exchange to the right */
  MPI_Sendrecv(&in[1], 1, pkey_mpi_datatype, (rank + 1 < size)?(rank + 1):MPI_PROC_NULL, 0,
               &out[0], 1, pkey_mpi_datatype, (rank - 1 >= 0)?(rank - 1):MPI_PROC_NULL, 0,
               comm, &status);

  if (rank - 1 < 0) out[0] = in[0];
  if (rank + 1 >= size) out[1] = in[1];

  return 0;
}


slint_t mpi_elements_check_order(elements_t *s, slint_t nelements, slint_t *orders, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_check_order */
{
  slint_t _orders[2], my_orders[2];
  slkey_pure_t my_keys[2], other_keys[2];

  if (orders == NULL) orders = _orders;

  /* check local order */
  my_orders[0] = (elements_validate_order(s, nelements) == 0);

  /* check global order (independent of local order) */
  elements_get_minmax_keys(s, nelements, my_keys);
  mpi_linear_exchange_pure_keys(my_keys, other_keys, size, rank, comm);
  my_orders[1] = (key_pure_cmp_le(other_keys[0], my_keys[0]) && key_pure_cmp_le(my_keys[1], other_keys[1]));

  MPI_Allreduce(my_orders, orders, 2, int_mpi_datatype, MPI_LAND, comm);

  return (orders[0] && orders[1]);
}


slint_t mpi_check_global_order(slkey_pure_t local_min, slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_check_global_order */
{
  slint_t local_order, global_order;
  slkey_pure_t my_keys[2], other_keys[2];


  my_keys[0] = local_min;
  my_keys[1] = local_max;
  
  mpi_linear_exchange_pure_keys(my_keys, other_keys, size, rank, comm);

  local_order = (key_pure_cmp_le(other_keys[0], my_keys[0]) && key_pure_cmp_le(my_keys[1], other_keys[1]));

  global_order = -1;
  
  if (root < 0) MPI_Allreduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_LAND, comm);
  else MPI_Reduce(&local_order, &global_order, 1, int_mpi_datatype, MPI_LAND, root, comm);
  
  return global_order;
}


slint_t mpi_elements_digest_sum(elements_t *s, slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_digest_sum */
{
  unsigned int lsum;

#ifdef ELEMENTS_DIGEST_SUM_ADD
# define MPI_SUM_OP  MPI_SUM
#else
# define MPI_SUM_OP  MPI_BXOR
#endif

  elements_digest_sum(s, nelements, components, &lsum);

  MPI_Allreduce(&lsum, sum, 1, MPI_UNSIGNED, MPI_SUM_OP, comm);

  return 0;
}


/* deprecated */
unsigned int mpi_cs32(elements_t *s, slint n, slint keys, slint data, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_cs32 */
{
  unsigned int crc32;

  mpi_elements_digest_sum(s, n, (keys?(SLCM_KEYS):0)|(data?(SLCM_DATA):0), &crc32, size, rank, comm);

  return crc32;
}


slint_t mpi_elements_digest_hash(elements_t *s, slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_digest_hash */
{
  slint_t hash_size = elements_digest_hash(NULL, 0, 0, NULL);

  if (!hash) return hash_size;

  elements_digest_hash(s, nelements, components, hash);
  
  MPI_Allreduce(MPI_IN_PLACE, hash, hash_size, MPI_BYTE, MPI_BXOR, comm);

  return 0;
}



slint_t mpi_elements_get_counts(elements_t *s, slint_t *clocal, slint_t *cglobal, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_get_counts */
{
  slint_t lc, gc;


  if (!clocal) clocal = &lc;
  if (!cglobal) cglobal = &gc;
 
  *clocal = s->size;

  if (root < 0) MPI_Allreduce(clocal, cglobal, 1, int_mpi_datatype, MPI_SUM, comm);
  else MPI_Reduce(clocal, cglobal, 1, int_mpi_datatype, MPI_SUM, root, comm);

  return *cglobal;
}


slweight_t mpi_elements_get_weights(elements_t *s, slweight_t *wlocal, slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_get_weights */
{
  slweight_t lw, gw;


  if (!wlocal) wlocal = &lw;
  if (!wglobal) wglobal = &gw;

#ifdef elem_weight
  *wlocal = elements_get_weight(s);

  if (root < 0) MPI_Allreduce(wlocal, wglobal, 1, weight_mpi_datatype, MPI_SUM, comm);
  else MPI_Reduce(wlocal, wglobal, 1, weight_mpi_datatype, MPI_SUM, root, comm);
#else
  *wlocal = *wglobal = 0.0;
#endif
  return *wglobal;
}


slint_t mpi_elements_get_counts_and_weights(elements_t *s, slint_t nelements, slint_t *counts, slweight_t *weights, int root, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_get_counts_and_weights */
{
  slint_t i;
  slweight_t lcw[2], gcw[2];


  lcw[0] = lcw[1] = gcw[0] = gcw[1] = 0.0;

  for (i = 0; i < nelements; ++i)
  {
    lcw[0] += s[i].size;
#ifdef elem_weight
    lcw[1] += elements_get_weight(&s[i]);
#endif
  }

  if (root < 0)
  {
    MPI_Allreduce(lcw, gcw, elem_weight_ifelse(2, 1), weight_mpi_datatype, MPI_SUM, comm);

  } else
  {
    MPI_Reduce(lcw, gcw, elem_weight_ifelse(2, 1), weight_mpi_datatype, MPI_SUM, root, comm);
  }
  
  counts[0] = lcw[0];
  counts[1] = gcw[0];
#ifdef elem_weight
  weights[0] = lcw[1];
  weights[1] = gcw[1];
#endif
    
  return 0;
}


slint_t mpi_elements_sendrecv_replace(elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm)  /* sl_proto, sl_func mpi_elements_sendrecv_replace */
{
  MPI_Status status;
  elements_t xs;
  int offset, current, maxc;

  slint_t replace, min_byte;
  slint_t smaxsize, rmaxsize;


  min_byte = key_byte;
#define xelem_call \
  if (min_byte < xelem_byte) min_byte = xelem_byte;
#include "sl_xelem_call.h"

  Z_TRACE_IF(ME_TRACE_IF, "me_sendrecv_replace_memsize: %" slint_fmt " vs. min_byte: %" slint_fmt "", me_sendrecv_replace_memsize, min_byte);

#if 1
  
  if (me_sendrecv_replace_mem == NULL || me_sendrecv_replace_memsize < min_byte)
  {
    replace = 1;
    rmaxsize = me_sendrecv_replace_mpi_maxsize;

  } else
  {
    replace = 0;
    rmaxsize = me_sendrecv_replace_memsize;
  }

  /* make sure total MPI message sizes fit into "int" ("workaround" (at least) necessary for JUROPA) */
  rmaxsize = z_min(rmaxsize, INT_MAX / 2);

  MPI_Sendrecv(&rmaxsize, 1, int_mpi_datatype, dest, sendtag, &smaxsize, 1, int_mpi_datatype, source, recvtag, comm, &status);

  Z_ASSERT(smaxsize != 0 && rmaxsize != 0);

#define xelem_call \
  maxc = z_min(smaxsize, rmaxsize) / xelem_byte; \
  Z_TRACE_IF(ME_TRACE_IF, xelem_name_str ": count: %d, maxc: %d", count, maxc); \
\
  offset = 0; \
  xelem_buf(&xs) = me_sendrecv_replace_mem; \
\
  while (offset < count) \
  { \
    current = z_min(count - offset, maxc); \
    Z_TRACE_IF(ME_TRACE_IF, "current: %d", current); \
    if (replace) \
    { \
      Z_TRACE_IF(ME_TRACE_IF, "sendrecv_replace: %d @ %d to %d / from %d", current, offset, dest, source); \
      MPI_Sendrecv_replace(xelem_buf_at(s, offset), current, xelem_mpi_datatype, dest, sendtag, source, recvtag, comm, &status); \
    } else \
    { \
      Z_TRACE_IF(ME_TRACE_IF, "sendrecv: %d @ %d to %d from %d", current, offset, dest, source); \
      MPI_Sendrecv(xelem_buf_at(s, offset), current, xelem_mpi_datatype, dest, sendtag, xelem_buf(&xs), current, xelem_mpi_datatype, source, recvtag, comm, &status); \
      Z_TRACE_IF(ME_TRACE_IF, "ncopy: %d from 0 to %d", current, offset); \
      xelem_ncopy_at(&xs, 0, s, offset, current); \
    } \
    offset += current; \
  }
#include "sl_xelem_call.h"

#else

  if (me_sendrecv_replace_mem == NULL || me_sendrecv_replace_memsize < min_byte)
  {
#define xelem_call \
    maxc = me_sendrecv_replace_mpi_maxsize / xelem_byte; \
    offset = 0; \
    while (offset < count) \
    { \
      current = z_min(maxc, count - offset); \
      MPI_Sendrecv_replace(xelem_buf_at(s, offset), current, xelem_mpi_datatype, dest, sendtag, source, recvtag, comm, &status); \
      offset += current; \
    }
#include "sl_xelem_call.h"
    return 0;
  }

#define xelem_call \
  maxc = (INT_MAX / 2) / xelem_byte; \
  xs.size = z_min(maxc, me_sendrecv_replace_memsize / xelem_byte); \
  xelem_buf(&xs) = me_sendrecv_replace_mem; \
  offset = 0; \
  while (offset < count) \
  { \
    current = z_min(xs.size, count - offset); \
    MPI_Sendrecv(xelem_buf_at(s, offset), current, xelem_mpi_datatype, dest, sendtag, xelem_buf(&xs), current, xelem_mpi_datatype, source, recvtag, comm, &status); \
    xelem_ncopy_at(&xs, 0, s, offset, current); \
    offset += current; \
  }
#include "sl_xelem_call.h"

#endif

  return 0;
}


unsigned int mpi_elements_crc32(elements_t *s, slint_t n, slint_t keys, slint_t data, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_elements_crc32 */
{
  unsigned int crc32_local, crc32_global;
  
  
  crc32_local = elements_crc32(s, 1, keys, data);

  MPI_Allreduce(&crc32_local, &crc32_global, 1, MPI_UNSIGNED, MPI_SUM, comm);

  return crc32_global;
}
