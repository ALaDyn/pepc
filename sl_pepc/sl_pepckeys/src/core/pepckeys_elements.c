/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/elements.c
 *  timestamp: 2009-11-19 12:48:09 +0100
 *  
 */


#include "sl_common.h"


slint_t pepckeys_elements_alloc(elements_t *s, slint_t nelements, slint_t keys, slint_t data) /* pepckeys_sl_proto, sl_func pepckeys_elements_alloc */
{
  if (s == NULL) return -1;

  elem_null(s);

  if (nelements == 0) return 0;

  s->size = s->max_size = nelements;
  if (keys) s->keys = sl_alloc(nelements, sizeof(slkey_t));

#ifdef SL_INDEX
  s->indices = sl_alloc(nelements, sizeof(slindex_t));
#endif

  if (data)
  {
#define xelem_key_not
#define xelem_call          xelem_buf(s) = sl_alloc(nelements, xelem_size_c * sizeof(xelem_type_c));
#include "sl_xelem_call.h"
  }

  if ((keys == 0) == (s->keys == NULL))
#define xelem_key_not
#define xelem_call          if ((data == 0) == (xelem_buf(s) == NULL))
#include "sl_xelem_call.h"
    return 0;

  /* a required allocation failed, free all */
  pepckeys_elements_free(s);

  return -1;
}


slint_t pepckeys_elements_free(elements_t *s) /* pepckeys_sl_proto, sl_func pepckeys_elements_free */
{
  if (s == NULL) return -1;

#define xelem_call          sl_free(xelem_buf(s));
#include "sl_xelem_call.h"

  elem_null(s);

  return 0;
}


/* FIXME: alignment! */
slint_t pepckeys_elements_alloc_from_block(elements_t *s, void *block, slint_t blocksize, slint_t alignment) /* pepckeys_sl_proto, sl_func pepckeys_elements_alloc_from_block */
{
  slint_t maxn = blocksize / elem_byte;
  
  char *base = (char *) block;

  s->size = maxn;

#define xelem_call          xelem_buf(s) = (xelem_type_c *) base; \
                            base += (maxn * xelem_byte);
#include "sl_xelem_call.h"

/*  pepckeys_elements_printf(s);*/

  return 0;
}


slint_t pepckeys_elements_copy(elements_t *s, elements_t *d) /* pepckeys_sl_proto, sl_func pepckeys_elements_copy */
{
  elem_copy(s, d);
  
  return 0;
}


slint_t pepckeys_elements_copy_at(elements_t *s, slint_t sat, elements_t *d, slint_t dat) /* pepckeys_sl_proto, sl_func pepckeys_elements_copy_at */
{
  elem_copy_at(s, sat, d, dat);
  
  return 0;
}


slint_t pepckeys_elements_ncopy(elements_t *s, elements_t *d, slint_t n) /* pepckeys_sl_proto, sl_func pepckeys_elements_ncopy */
{
  elem_ncopy(s, d, n);
  
  return 0;
}


slint_t pepckeys_elements_nmove(elements_t *s, elements_t *d, slint_t n) /* pepckeys_sl_proto, sl_func pepckeys_elements_nmove */
{
  elem_nmove(s, d, n);
  
  return 0;
}


slint_t pepckeys_elements_printf(elements_t *s) /* pepckeys_sl_proto, sl_func pepckeys_elements_printf */
{
  if (s == NULL) return -1;

  printf("[%" sl_int_type_fmt
#define xelem_call          ", %p"
#include "sl_xelem_call.h"
    "]\n", s->size
#define xelem_call          , xelem_buf(s)
#include "sl_xelem_call.h"
    );

  return 0;
}


slint_t pepckeys_elements_extract(elements_t *src, slint_t nelements, elements_t *dst0, elements_t *dst1) /* pepckeys_sl_proto, sl_func pepckeys_elements_extract */
{
  elements_t s;

  if (src == NULL) return -1;

  s = *src;

  if (dst0 != NULL)
  {
    elem_assign(&s, dst0);
    dst0->size = xmin(s.size, nelements);

    if (dst0->size <= 0) elem_null(dst0);
  }

  if (dst1 != NULL)
  {
    elem_assign_at(&s, nelements, dst1);
    dst1->size = xmax(s.size - nelements, 0);

    if (dst1->size <= 0) elem_null(dst1);
  }

  return 0;
}


slint_t pepckeys_elements_touch(elements_t *s) /* pepckeys_sl_proto, sl_func pepckeys_elements_touch */
{
  elements_t _s, end, t;

  pepckeys_elements_alloc(&t, 1, 1, 1);

  elem_assign_at(s, s->size, &end);

  for (elem_assign(s, &_s); _s.keys < end.keys; elem_inc(&_s)) elem_copy(&_s, &t);

  pepckeys_elements_free(&t);

  return 0;
}


slint_t pepckeys_elements_random_exchange(elements_t *s, slint_t rounds, elements_t *xs) /* pepckeys_sl_proto, sl_func pepckeys_elements_random_exchange */
{
  slint_t i, j, k = 0;
  elements_t txs;

  if (s == NULL) return -1;

  if (xs == NULL)
  {
    xs = &txs;
    pepckeys_elements_alloc(xs, 1, 1, 1);
  }

  j = 0;
  elem_copy(s, xs);

  for (i = 0; i < rounds; i++)
  {
    k = sl_rand() % s->size;
    elem_copy_at(s, k, s, j);
    j = k;
  }

  elem_copy_at(xs, 0, s, k);

  if (xs == &txs) pepckeys_elements_free(xs);

  return 0;
}


slint_t pepckeys_elements_init_keys(elements_t *s, slint_t dtype, slint_t _min, slint_t _max) /* pepckeys_sl_proto, sl_func pepckeys_elements_init_keys */
{
  slint_t i;

  if (s == NULL) return -1;

#ifdef key_integer
  switch (dtype)
  {
    case 1:
      for (i = 0; i < s->size; i++) key_purify(s->keys[i]) = (key_pure_type_c) sl_rand_minmax(_min, _max);
      break;
    case 2:
/*      for (i = 0; i < s->size; i++) s->keys[i] = sl_key(_min + (i % (_max - _min + 1)));*/
      break;
    case 3:
/*      for (i = 0; i < s->size; i++) s->keys[i] = sl_key(_min + iround((double) i / (double) (s->size - 1) * (double) (_max - _min)));*/
      break;
   }
#endif

   return 0;
}


#define LINE_LENGTH  1024

slint_t pepckeys_elements_init_keys_from_file(elements_t *s, slint_t data, char *filename, slint_t from, slint_t to, slint_t const_bytes_per_line) /* pepckeys_sl_proto, sl_func pepckeys_elements_init_keys_from_file */
{
  FILE *inputfile;
  char buffer[LINE_LENGTH];
  slint_t i = 0, line = 0;
  slint_t bytes_per_line;

  pepckeys_elements_alloc(s, to - from + 1, 1, data);

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
    sscanf(buffer, "%" key_pure_type_fmt, &key_purify(s->keys[i]));
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


slint_t pepckeys_elements_save_keys_to_file(elements_t *s, char *filename) /* pepckeys_sl_proto, sl_func pepckeys_elements_save_keys_to_file */
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

slint_t pepckeys_elements_validate_order(elements_t *s, slint_t n) /* pepckeys_sl_proto, sl_func pepckeys_elements_validate_order */
{

#define km(k) key_purify(k)

  evo_body

#undef km

  return 0;
}


slint_t pepckeys_elements_validate_order_bmask(elements_t *s, slint_t n, slkey_pure_t bmask) /* pepckeys_sl_proto, sl_func pepckeys_elements_validate_order_bmask */
{

#define  km(k) (key_purify(k) & bmask)

  evo_body

#undef km

  return 0;
}


slint_t pepckeys_elements_validate_order_weight(elements_t *s, slint_t n, slkey_pure_t weight) /* pepckeys_sl_proto, sl_func pepckeys_elements_validate_order_weight */
{

#define  km(k) (key_purify(k) / weight)

  evo_body

#undef km

  return 0;
}

#undef evo_body


slint_t pepckeys_elements_print_keys(elements_t *s) /* pepckeys_sl_proto, sl_func pepckeys_elements_print_keys */
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
