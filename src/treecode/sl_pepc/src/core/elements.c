/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/elements.c
 *  
 */


/* sl_macro E_TRACE_IF */

#include "sl_common.h"


#ifndef E_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define E_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define E_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t elements_alloc(elements_t *s, slint_t nelements, slcint_t components) /* sl_proto, sl_func elements_alloc */
{
  slint_t failed = 0;


  if (s == NULL) return -1;

  elem_null(s);

  if (nelements == 0) return 0;

  s->size = s->max_size = nelements;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = z_alloc(nelements, xelem_size_c * sizeof(xelem_type_c)); \
  if (xelem_buf(s) == NULL) failed = 1; \
}
#include "sl_xelem_call.h"

  if (failed)
  {
    elements_free(s);
    return -1;
  }
  
  return 0;
}


slint_t elements_free(elements_t *s) /* sl_proto, sl_func elements_free */
{
  if (s == NULL) return -1;

#define xelem_call          z_free(xelem_buf(s));
#include "sl_xelem_call.h"

  elem_null(s);

  return 0;
}


slint_t elements_alloc2(elements_t *s, slint_t nelements, slint_t keys, slint_t indices, slint_t data, slint_t weights) /* sl_proto, sl_func elements_alloc2 */
{
  if (s == NULL) return -1;

  elem_null(s);

  if (nelements == 0) return 0;

  s->size = s->max_size = nelements;
  if (keys) s->keys = z_alloc(nelements, sizeof(slkey_t));

#ifdef SL_INDEX
  if (indices) s->indices = z_alloc(nelements, sizeof(slindex_t));
#endif

  /* FIXME: allocate _only_ the weight data components */
  if (weights) data = 1;
  
  if (data)
  {
#define xelem_key_not
#define xelem_index_not
#define xelem_call          xelem_buf(s) = z_alloc(nelements, xelem_size_c * sizeof(xelem_type_c));
#include "sl_xelem_call.h"
  }

  if ((keys == 0) == (s->keys == NULL))
#ifdef SL_INDEX
  if ((indices == 0) == (s->indices == NULL))
#endif
#define xelem_key_not
#define xelem_index_not
#define xelem_call          if ((data == 0) == (xelem_buf(s) == NULL))
#include "sl_xelem_call.h"
    return 0;

  /* a required allocation failed, free all */
  elements_free(s);

  return -1;
}


slint_t elements_alloc_old(elements_t *s, slint_t nelements, slint_t keys, slint_t data) /* sl_proto, sl_func elements_alloc_old */
{
  return elements_alloc2(s, nelements, keys, 1, data, 1);
}


/* FIXME: alignment! */
slint_t elements_alloc_from_blocks(elements_t *s, slint_t nblocks, void **blocks, slint_t *blocksizes, slint_t alignment, slint_t nmax, slcint_t components) /* sl_proto, sl_func elements_alloc_from_blocks */
{
  slint_t csizes[2 + data_nmax], ctargets[2 + data_nmax];
  slint_t cbytes[nblocks], nelements[nblocks];
  slint_t max_nelements, max_i, cur_min;

  slint_t i, j, n;


  if (s == NULL) return -1;

  elem_null(s);

  for (i = 0; i < nblocks; ++i)
  {
    cbytes[i] = 0;

    Z_TRACE_IF(E_TRACE_IF, "block %" slint_fmt ": %" slint_fmt " @ %p", i, blocksizes[i], blocks[i]);
  }

  n = 1;
  i = 0;
#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  csizes[i] = xelem_byte; \
  ctargets[i] = 0; \
  cbytes[0] += csizes[i]; \
  n *= nblocks; \
  ++i; \
}
#include "sl_xelem_call.h"

  Z_TRACE_IF(E_TRACE_IF, "n: %" slint_fmt, n);

  if (cbytes[0] > 0) nelements[0] = (blocks[0])?(blocksizes[0] / cbytes[0]):0;
  else nelements[0] = -1;

  max_nelements = nelements[0];
  max_i = 0;

  Z_TRACE_IF(E_TRACE_IF, "nelements[0]: %" slint_fmt ", max nelements: %" slint_fmt, nelements[0], max_nelements);

  /* find best combination (scales with nblocks^ncomponents) */
  for (i = 1; i < n; ++i)
  {
    if (nmax >= 0 && max_nelements >= nmax) break;

    j = -1;
    do
    {
      ++j;

      cbytes[ctargets[j]] -= csizes[j];
      
      if (cbytes[ctargets[j]] > 0) nelements[ctargets[j]] = (blocks[ctargets[j]])?(blocksizes[ctargets[j]] / cbytes[ctargets[j]]):0;
      else nelements[ctargets[j]] = -1;

      ctargets[j] = (ctargets[j] + 1) % nblocks;

      cbytes[ctargets[j]] += csizes[j];

      if (cbytes[ctargets[j]] > 0) nelements[ctargets[j]] = (blocks[ctargets[j]])?(blocksizes[ctargets[j]] / cbytes[ctargets[j]]):0;
      else nelements[ctargets[j]] = -1;

    } while (ctargets[j] == 0);

    cur_min = -1;
    for (j = 0; j < nblocks; ++j) if (cur_min < 0) cur_min = nelements[j]; else if (nelements[j] >= 0) cur_min = z_min(cur_min, nelements[j]);

    if (cur_min > max_nelements)
    {
      max_nelements = cur_min;
      max_i = i;
    }

    Z_TRACE_IF(E_TRACE_IF, "%" slint_fmt ": cur: %" slint_fmt ", max: %" slint_fmt "(%" slint_fmt ")", i, cur_min, max_nelements, max_i);
  }

  Z_ASSERT(max_nelements >= 0);

  Z_TRACE_IF(E_TRACE_IF, "max i: %" slint_fmt ", max nelements: %" slint_fmt, max_i, max_nelements);

  if (nmax >= 0 && max_nelements > nmax) max_nelements = nmax;

  if (max_nelements == 0) return 0;

  s->size = s->max_size = max_nelements;

#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  xelem_buf(s) = (xelem_type_c *) blocks[max_i % nblocks]; \
  blocks[max_i % nblocks] = xelem_buf_at(s, max_nelements); \
  max_i /= nblocks; \
}
#include "sl_xelem_call.h"

  return 0;
}


slint_t elements_alloc_from_block2(elements_t *s, void *block, slint_t blocksize, slint_t alignment, slint_t nmax, slint_t keys, slint_t indices, slint_t data, slint_t weights) /* sl_proto, sl_func elements_alloc_from_block2 */
{
  return elements_alloc_from_blocks(s, 1, &block, &blocksize, alignment, nmax, ((keys)?SLCM_KEYS:0)|((indices)?SLCM_INDICES:0)|((data)?SLCM_DATA:0)|((weights)?SLCM_WEIGHTS:0));
}


/* FIXME: alignment! */
slint_t elements_alloc_from_block(elements_t *s, void *block, slint_t blocksize, slint_t alignment, slint_t nmax) /* sl_proto, sl_func elements_alloc_from_block */
{
  return elements_alloc_from_blocks(s, 1, &block, &blocksize, alignment, nmax, SLCM_ALL);
}


slint_t elements_copy(elements_t *s, elements_t *d) /* sl_proto, sl_func elements_copy */
{
  elem_copy(s, d);
  
  return 0;
}


slint_t elements_copy_at(elements_t *s, slint_t sat, elements_t *d, slint_t dat) /* sl_proto, sl_func elements_copy_at */
{
  elem_copy_at(s, sat, d, dat);
  
  return 0;
}


slint_t elements_ncopy(elements_t *s, elements_t *d, slint_t n) /* sl_proto, sl_func elements_ncopy */
{
  elem_ncopy(s, d, n);
  
  return 0;
}


slint_t elements_nmove(elements_t *s, elements_t *d, slint_t n) /* sl_proto, sl_func elements_nmove */
{
  elem_nmove(s, d, n);
  
  return 0;
}


slint_t elements_printf(elements_t *s, const char *prefix) /* sl_proto, sl_func elements_printf */
{
  if (s == NULL) return -1;

  printf("%s: [%" sl_int_type_fmt
#define xelem_call          ", %p"
#include "sl_xelem_call.h"
    "]\n", prefix, s->size
#define xelem_call          , xelem_buf(s)
#include "sl_xelem_call.h"
    );

  return 0;
}


slint_t elements_extract(elements_t *src, slint_t nelements, elements_t *dst0, elements_t *dst1) /* sl_proto, sl_func elements_extract */
{
  elements_t s;

  if (src == NULL) return -1;

  s = *src;

  if (dst0 != NULL)
  {
    elem_assign(&s, dst0);
    dst0->size = z_min(s.size, nelements);

    if (dst0->size <= 0) elem_null(dst0);
  }

  if (dst1 != NULL)
  {
    elem_assign_at(&s, nelements, dst1);
    dst1->size = z_max(s.size - nelements, 0);

    if (dst1->size <= 0) elem_null(dst1);
  }

  return 0;
}


slint_t elements_touch(elements_t *s) /* sl_proto, sl_func elements_touch */
{
  elements_t _s, end, t;

  elements_alloc(&t, 1, SLCM_ALL);

  elem_assign_at(s, s->size, &end);

  for (elem_assign(s, &_s); _s.keys < end.keys; elem_inc(&_s)) elem_copy(&_s, &t);

  elements_free(&t);

  return 0;
}


slint_t elements_digest_sum(elements_t *s, slint_t nelements, slcint_t components, unsigned int *sum) /* sl_proto, sl_func elements_digest_sum */
{
  slint_t i, j;
  unsigned int ssum = 0;

#ifdef ELEMENTS_DIGEST_SUM_ADD
# define SUM_OP(_a_,_b_)  _a_ + _b_
#else
# define SUM_OP(_a_,_b_)  _a_ ^ _b_
#endif

  *sum = 0;
  for (j = 0; j < nelements; j++)
  {
    for (i = 0; i < s[j].size; i++)
    {
#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  z_digest_sum_buffer((const void *) xelem_buf_at(&s[j], i), xelem_byte, (void *) &ssum); \
  *sum = SUM_OP(*sum, ssum); \
}
#include "sl_xelem_call.h"
    }
  }

#undef SUM_OP

  return 0;
}


/* deprecated */
unsigned int elements_crc32(elements_t *s, slint nelements, slint_t keys, slint_t data) /* sl_proto, sl_func elements_crc32 */
{
  unsigned int crc32;
  
  elements_digest_sum(s, nelements, (keys?(SLCM_KEYS):0)|(data?(SLCM_DATA):0), &crc32);

  return crc32;
}


slint_t elements_digest_hash(elements_t *s, slint_t nelements, slcint_t components, void *hash) /* sl_proto, sl_func elements_digest_hash */
{
  slint_t i, j;

  if (!hash) return z_digest_hash_read(NULL, NULL);

  void *hdl;

  z_digest_hash_open(&hdl);

  for (j = 0; j < nelements; j++)
  {
    for (i = 0; i < s[j].size; i++)
    {
#define xelem_call  if ((components & xelem_cm) || ((components & SLCM_WEIGHTS) && xelem_weight)) \
{ \
  z_digest_hash_write(hdl, (const void *) xelem_buf_at(&s[j], i), xelem_byte); \
}
#include "sl_xelem_call.h"
    }
  }

  z_digest_hash_read(hdl, hash);

  z_digest_hash_close(hdl);

  return 0;
}


slint_t elements_random_exchange(elements_t *s, slint_t rounds, elements_t *xs) /* sl_proto, sl_func elements_random_exchange */
{
  slint_t i, j, k = 0;
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL || xs->size < 1)
  {
    xs = &txs;
    elements_alloc(xs, 1, SLCM_ALL);
  }

  j = 0;
  elem_copy(s, xs);

  for (i = 0; i < rounds; i++)
  {
    k = z_rand() % s->size;
    elem_copy_at(s, k, s, j);
    j = k;
  }

  elem_copy_at(xs, 0, s, k);

  if (xs == &txs) elements_free(xs);

  return 0;
}


slint_t elements_keys_init_seed(unsigned long s) /* sl_proto, sl_func elements_keys_init_seed */
{
#ifdef key_val_srand
  key_val_srand(s);
#endif

  z_urandom_seed(s);
  z_nrandom_seed(s);

  return 0;
}

#define KEY_SET_VARIABLES() \
  slint_t ks_j, ks_m; \
  slkey_pure_t *ks_pk; \
  key_set_f ks_func; \
  key_set_data_t ks_data

#define KEY_SET_SET_INIT(_d_)       Z_NOP()
#define KEY_SET_SET(_k_, _d_)       Z_MOP(*(_k_) = *((slkey_pure_t *) (_d_));)

#define KEY_SET_SET_FUNC_INIT(_d_)  Z_MOP(ks_func = (key_set_f) ((void **) (_d_))[0]; ks_data = (key_set_data_t) ((void **) (_d_))[1];)
#define KEY_SET_SET_FUNC(_k_, _d_)  Z_MOP(ks_func(_k_, ks_data);)

#define KEY_SET_RAND_INIT(_d_)      Z_MOP( \
  if (_d_) { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand_minmax is missing!"); \
    ks_pk = (slkey_pure_t *) (_d_); \
  } else { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand is missing!"); \
  })
#define KEY_SET_RAND(_k_, _d_)      Z_MOP(if (_d_) *(_k_) = key_val_rand_minmax(ks_pk[0], ks_pk[1]); else *(_k_) = key_val_rand();)

#define KEY_SET_RAND_QUAD_INIT(_d_)  Z_MOP( \
  if (_d_) { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand_minmax is missing!"); \
    ks_pk = (slkey_pure_t *) (_d_); \
  } else { \
    if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand is missing!"); \
  })
#define KEY_SET_RAND_QUAD(_k_, _d_)  Z_MOP(if (_d_) *(_k_) = key_val_rand_minmax(ks_pk[0], ks_pk[1]); else *(_k_) = key_val_rand(); *(_k_) *= *(_k_);)

#ifdef key_integer
#define KEY_SET_RAND_AND_INIT(_d_)  Z_MOP( \
  ks_m = *((slint_t *) (_d_)); \
  if (!have_key_val_rand_minmax) Z_ERROR("key_val_rand is missing!");)
#define KEY_SET_RAND_AND(_k_, _d_)  Z_MOP(*(_k_) = key_val_rand(); for (ks_j = 0; ks_j < ks_m; ++ks_j) *(_k_) &= key_val_rand();)
#endif

#define KEY_SET_URAND_INIT(_d_)  Z_MOP(ks_pk = (slkey_pure_t *) (_d_);)
#define KEY_SET_URAND(_k_, _d_)  Z_MOP(*(_k_) = (slkey_pure_t) (ks_pk[0] + z_urandom() * (ks_pk[1] - ks_pk[0]) + 0.5);)

#define KEY_SET_NRAND_INIT(_d_)  Z_MOP(ks_pk = (slkey_pure_t *) (_d_);)
#define KEY_SET_NRAND(_k_, _d_)  Z_MOP(*(_k_) = (slkey_pure_t) (ks_pk[2] + z_nrandom() * ks_pk[3]); *(_k_) = z_minmax(ks_pk[0], *(_k_), ks_pk[1]);)


slint_t elements_keys_init(elements_t *s, keys_init_type_t t, keys_init_data_t d) /* sl_proto, sl_func elements_keys_init */
{
  slint_t i;

  KEY_SET_VARIABLES();

  if (s == NULL) return -1;

  switch (t)
  {
    case SL_EKIT_SET:
      KEY_SET_SET_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_SET(&s->keys[i], d);
      break;
    case SL_EKIT_SET_FUNC:
      KEY_SET_SET_FUNC_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_SET_FUNC(&s->keys[i], d);
      break;
    case SL_EKIT_RAND:
      KEY_SET_RAND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_RAND(&s->keys[i], d);
      break;
    case SL_EKIT_RAND_QUAD:
      KEY_SET_RAND_QUAD_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_RAND_QUAD(&s->keys[i], d);
      break;
#ifdef key_integer
    case SL_EKIT_RAND_AND:
      KEY_SET_RAND_AND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_RAND_AND(&s->keys[i], d);
      break;
#endif
    case SL_EKIT_URAND:
      KEY_SET_URAND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_URAND(&s->keys[i], d);
      break;
    case SL_EKIT_NRAND:
      KEY_SET_NRAND_INIT(d);
      for (i = 0; i < s->size; i++) KEY_SET_NRAND(&s->keys[i], d);
      break;
  }

  return 0;
}


slint_t elements_keys_init_randomized(elements_t *s, slint_t nkeys, keys_init_type_t t, keys_init_data_t d) /* sl_proto, sl_func elements_keys_init_randomized */
{
  slint_t i, j;

  KEY_SET_VARIABLES();

  if (s == NULL) return -1;

  i = 0;
  switch (t)
  {
    case SL_EKIT_SET:
      KEY_SET_SET_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_SET(&s->keys[i], d); }
      break;
    case SL_EKIT_SET_FUNC:
      KEY_SET_SET_FUNC_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_SET_FUNC(&s->keys[i], d); }
      break;
    case SL_EKIT_RAND:
      KEY_SET_RAND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_RAND(&s->keys[i], d); }
      break;
    case SL_EKIT_RAND_QUAD:
      KEY_SET_RAND_QUAD_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_RAND_QUAD(&s->keys[i], d); }
      break;
#ifdef key_integer
    case SL_EKIT_RAND_AND:
      KEY_SET_RAND_AND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_RAND_AND(&s->keys[i], d); }
      break;
#endif
    case SL_EKIT_URAND:
      KEY_SET_URAND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_URAND(&s->keys[i], d); }
      break;
    case SL_EKIT_NRAND:
      KEY_SET_NRAND_INIT(d);
      for (j = 0; j < nkeys; j++) { i = (i + z_rand()) % s->size; KEY_SET_NRAND(&s->keys[i], d); }
      break;
  }

  return 0;
}


#undef KEY_SET_VARIABLES
#undef KEY_SET_SET_INIT
#undef KEY_SET_SET
#undef KEY_SET_SET_FUNC_INIT
#undef KEY_SET_SET_FUNC
#undef KEY_SET_RAND_INIT
#undef KEY_SET_RAND
#undef KEY_SET_RAND_QUAD_INIT
#undef KEY_SET_RAND_QUAD
#undef KEY_SET_RAND_AND_INIT
#undef KEY_SET_RAND_AND
#undef KEY_SET_URAND_INIT
#undef KEY_SET_URAND
#undef KEY_SET_NRAND_INIT
#undef KEY_SET_NRAND


#define LINE_LENGTH  1024

slint_t elements_init_keys_from_file(elements_t *s, slint_t data, char *filename, slint_t from, slint_t to, slint_t const_bytes_per_line) /* sl_proto, sl_func elements_init_keys_from_file */
{
  FILE *inputfile;
  char buffer[LINE_LENGTH];
  slint_t i = 0, line = 0;
  slint_t bytes_per_line;
  slkey_pure_t inkey;

  elements_alloc(s, to - from + 1, SLCM_ALL|((data)?0:(~SLCM_DATA)));

/*  printf("opening '%s'\n", filename);*/

  inputfile = fopen(filename, "r");

/*  printf("inputfile = %p\n", inputfile);*/

  if (!inputfile) { return -1; }

#ifdef key_integer
  if (const_bytes_per_line)
  {
    fgets(buffer, LINE_LENGTH, inputfile);
    bytes_per_line = ftell(inputfile);
    rewind(inputfile);

    fseek(inputfile, from * bytes_per_line, SEEK_SET);

    line = from;
    
  } else while (line < from)
  {
    line++;
    if (!fgets(buffer, LINE_LENGTH, inputfile)) break;
  }

  while((i < s->size) && (line <= to))
  {
    if (!fgets(buffer, LINE_LENGTH, inputfile)) break;

#ifdef key_pure_type_fmt
    sscanf(buffer, "%" key_pure_type_fmt, &inkey);
    key_set_pure(&s->keys[i], inkey);
#endif

/*    printf("line: %d - input: '%s'", line, buffer);
    printf("i sscanf'd %d: %ld\n", r, s->keys[i]);*/

    ++line;

    ++i;
  }
#endif


  fclose(inputfile);

  return i;
}

#undef LINE_LENGTH


slint_t elements_save_keys_to_file(elements_t *s, char *filename) /* sl_proto, sl_func elements_save_keys_to_file */
{
  FILE *outputfile;
  slint_t i = 0;

  printf("opening '%s'\n", filename);

  outputfile = fopen(filename, "w");

  printf("inputfile = %p\n", outputfile);

  if (!outputfile) { return -1; }

  for (i = 0; i < s->size; ++i)
#ifdef key_pure_type_fmt
    fprintf(outputfile, "%" key_pure_type_fmt "\n", key_purify(s->keys[i]));
#else
    ;
#endif

  fclose(outputfile);

  return 0;
}


#define evo_body \
  slint_t i, j, k = 0, l = -1; \
  if (s == NULL) return -1; \
  for (j = 0; j < n; j++) \
  { \
    k++; \
    if ((l >= 0) && (s[j].size > 0)) \
    if (key_pure_cmp_lt(km(s[j].keys[0]), km(s[l].keys[s[l].size - 1]))) return j; \
    for (i = 1; i < s[j].size; i++, k++) \
    if (key_pure_cmp_lt(km(s[j].keys[i]), km(s[j].keys[i - 1]))) return k; \
    if (s[j].size > 0) l = j; \
  }

slint_t elements_validate_order(elements_t *s, slint_t n) /* sl_proto, sl_func elements_validate_order */
{

#define km(k) key_purify(k)

  evo_body

#undef km

  return 0;
}


slint_t elements_validate_order_bmask(elements_t *s, slint_t n, slkey_pure_t bmask) /* sl_proto, sl_func elements_validate_order_bmask */
{

#define  km(k) (key_purify(k) & bmask)

  evo_body

#undef km

  return 0;
}


slint_t elements_validate_order_weight(elements_t *s, slint_t n, slkey_pure_t weight) /* sl_proto, sl_func elements_validate_order_weight */
{

#define  km(k) (key_purify(k) / weight)

  evo_body

#undef km

  return 0;
}

#undef evo_body


#if defined(HAVE_GMP_H)
# include <gmp.h>
#endif


slint_t elements_keys_stats(elements_t *s, slkey_pure_t *stats) /* sl_proto, sl_func elements_keys_stats */
{
  slint_t i;
  slkey_pure_t k, kmin, kmax;


  if (s->size <= 0) return 0;

#ifdef HAVE_GMP_H
# ifdef key_integer
#  define stats_t                                    mpz_t
#  define stats_init(_x_)                            mpz_init(_x_)
#  define stats_free(_x_)                            mpz_clear(_x_)
#  define stats_set(_x_, _v_)                        Z_MOP(if (sizeof(_v_) <= sizeof(long)) \
                                                           { if (key_integer_unsigned) mpz_set_ui(_x_, (unsigned long) _v_); else mpz_set_si(_x_, (long) _v_); } \
                                                           else \
                                                           { if (key_integer_unsigned) z_gmp_mpz_set_ull(_x_, (unsigned long long) _v_); else z_gmp_mpz_set_sll(_x_, (long long) _v_); })
#  define stats_get(_x_, _type_)                     (_type_) ((sizeof(_type_) <= sizeof(long))? \
                                                               (key_integer_unsigned?mpz_get_ui(_x_):mpz_get_si(_x_)): \
                                                               (key_integer_unsigned?z_gmp_mpz_get_ull(_x_):z_gmp_mpz_get_sll(_x_)))
#  define stats_add(_x_, _v_)                        mpz_add(_x_, _x_, _v_);
#  define stats_addsqr(_x_, _v_, _t_)                Z_MOP(mpz_mul(_t_, _v_, _v_); mpz_add(_x_, _x_, _t_);)
#  define stats_avg(_sum_, _n_, _t_, _type_)         (_type_) ((mpz_tdiv_q_ui(_t_, _sum_, (unsigned long) _n_) < (_n_ / 2))? \
                                                               (stats_get(_t_, _type_)): \
                                                               (stats_get(_t_, _type_) + 1))
#  define stats_std(_sum_, _sqr_, _n_, _t_, _type_)  (_type_) (mpz_mul(_t_, _sum_, _sum_), \
                                                               mpz_tdiv_q_ui(_t_, _t_, (unsigned long) _n_), \
                                                               mpz_sub(_t_, _sqr_, _t_), \
                                                               mpz_tdiv_q_ui(_t_, _t_, (unsigned long) (_n_ - 1)), \
                                                               mpz_sqrt(_t_, _t_), stats_get(_t_, _type_))
#  define stats_print(_x_)                           gmp_printf(#_x_ ": %Zd\n", _x_)
# endif
#else
# define stats_t                                     long double
# define stats_init(_x_)                             Z_NOP()
# define stats_free(_x_)                             Z_NOP()
# define stats_set(_x_, _v_)                         _x_ = (stats_t) _v_
# define stats_get(_x_, _type_)                      (_type_) _x_
# define stats_add(_x_, _v_)                         _x_ += (stats_t) _v_
# define stats_addsqr(_x_, _v_, _t_)                 _x_ += (stats_t) _v_ * (stats_t) _v_
# define stats_avg(_sum_, _n_, _t_, _type_)          (_type_) (_sum_ / _n_ + 0.5)
# define stats_std(_sum_, _sqr_, _n_, _t_, _type_)   (_type_) sqrt((_sqr_ - ((_sum_ * _sum_) / (stats_t) _n_)) / (stats_t) (_n_ - 1))
# define stats_print(_x_)                            printf(#_x_ ": %Lf\n", _x_)
#endif

  stats_t x, xsum, xsqr, xt;

/*  elements_print_keys(s);*/

  kmin = kmax = key_get_pure(s->keys[0]);

  stats_init(x);
  stats_init(xsum);
  stats_init(xsqr);
  stats_init(xt);

  stats_set(x, 0);
  stats_set(xsum, 0);
  stats_set(xsqr, 0);
  stats_set(xt, 0);

  for (i = 0; i < s->size; i++)
  {
    k = key_get_pure(s->keys[i]);

/*    printf("key: %" sl_key_type_fmt "\n", k);*/

    if (k < kmin) kmin = k;
    if (k > kmax) kmax = k;

    stats_set(x, k);
    stats_add(xsum, x);
    stats_addsqr(xsqr, x, xt);

/*    stats_print(x);
    stats_print(xsum);
    stats_print(xsqr);*/
  }

  stats[SL_EKS_MIN] = kmin;
  stats[SL_EKS_MAX] = kmax;
  stats[SL_EKS_SUM] = stats_get(xsum, slkey_pure_t);
  stats[SL_EKS_AVG] = stats_avg(xsum, s->size, xt, slkey_pure_t);

/*  stats_print(ksum);
  stats_print(ksqr);*/

  if (s->size > 1) stats[SL_EKS_STD] = stats_std(xsum, xsqr, s->size, xt, slkey_pure_t);
  else stats[SL_EKS_STD] = 0;

  stats_free(x);
  stats_free(xsum);
  stats_free(xsqr);
  stats_free(xt);

  return 0;
}

#undef stats_t
#undef stats_init
#undef stats_free
#undef stats_set
#undef stats_get
#undef stats_add
#undef stats_addsqr
#undef stats_avg
#undef stats_std
#undef stats_print


slint_t elements_keys_stats_print(elements_t *s) /* sl_proto, sl_func elements_keys_stats_print */
{
  slkey_pure_t stats[SL_EKS_SIZE];


  elements_keys_stats(s, stats);

  printf("min: %" key_pure_type_fmt "\n", stats[SL_EKS_MIN]);
  printf("max: %" key_pure_type_fmt "\n", stats[SL_EKS_MAX]);
  printf("sum: %" key_pure_type_fmt "\n", stats[SL_EKS_SUM]);
  printf("avg: %" key_pure_type_fmt "\n", stats[SL_EKS_AVG]);
  printf("std: %" key_pure_type_fmt "\n", stats[SL_EKS_STD]);

  return 0;
}


slint_t elements_print_keys(elements_t *s) /* sl_proto, sl_func elements_print_keys */
{
  slint_t i;

  if (s == NULL) return -1;

  for (i = 0; i < s->size; i++)
  {
    printf(" [%3" sl_int_type_fmt "] @ %p = ", i, &s->keys[i]);
/*    key_printf(s->keys[i]);*/
#ifdef key_pure_type_fmt
    printf("%" key_pure_type_fmt "\n", key_purify(s->keys[i]));
#else
    printf("\n");
#endif
  }

  return 0;
}


slint_t elements_print_all(elements_t *s) /* sl_proto, sl_func elements_print_all */
{
  slint_t i;

  if (s == NULL) return -1;

  for (i = 0; i < s->size; i++)
  {
    printf(" [%3" slint_fmt "]", i);
#ifdef SL_INDEX
    printf(", idx @ %p = %" sl_index_type_fmt, &s->indices[i], s->indices[i]);
#endif
#ifdef key_pure_type_fmt
    printf(", key @ %p = %" key_pure_type_fmt, &s->keys[i], s->keys[i]);
#endif
#ifdef elem_weight
    printf(", weight = %" slweight_fmt, elem_weight(s, i));
#endif
    printf("\n");
  }

  return 0;
}


slweight_t elements_get_weight(elements_t *s) /* sl_proto, sl_func elements_get_weight */
{
#ifdef elem_weight
  slint_t i;
  slweight_t w = 0.0;

  for (i = 0; i < s->size; ++i) w += elem_weight(s, i);

  return w;
#else
  return 0.0;
#endif
}


slint_t elements_get_minmax_keys(elements_t *s, slint_t nelements, slkey_pure_t *minmaxkeys) /* sl_proto, sl_func elements_get_minmax_keys */
{
  slint_t i, j;

  if (s == NULL || nelements < 1) return -1;
  
  minmaxkeys[0] = minmaxkeys[1] = s[0].keys[0];
  
  for (j = 0; j < nelements; ++j)
  for (i = 0; i < s[j].size; ++i)
  {
    if (key_pure_cmp_lt(key_purify(s[j].keys[i]), minmaxkeys[0])) minmaxkeys[0] = key_purify(s[j].keys[i]);
    if (key_pure_cmp_gt(key_purify(s[j].keys[i]), minmaxkeys[1])) minmaxkeys[1] = key_purify(s[j].keys[i]);
  }

  return 0;
}
