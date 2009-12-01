/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/elements_packed.c
 *  timestamp: 2009-11-23 13:06:23 +0100
 *  
 */


#include "sl_common.h"


slint_t pepckeys_elements_alloc_packed(packed_elements_t *s, slint_t nelements) /* pepckeys_sl_proto, sl_func pepckeys_elements_alloc_packed */
{
  if (s == NULL) return -1;

  s->size = s->max_size = 0;
  s->elements = NULL;

  if (nelements == 0) return 0;

  s->size = nelements;
  
  s->elements = sl_alloc(nelements, sizeof(packed_element_t));

  if (s->elements != NULL) return 0;

  pepckeys_elements_free_packed(s);

  return -1;
}


slint_t pepckeys_elements_free_packed(packed_elements_t *s) /* pepckeys_sl_proto, sl_func pepckeys_elements_free_packed */
{
  if (s == NULL) return -1;

  sl_free(s->elements);

  s->size = s->max_size = 0;

  return 0;
}


slint_t pepckeys_elements_pack_indexed(elements_t *s, packed_elements_t *d, slindex_t *rindx, slindex_t *windx) /* pepckeys_sl_proto, sl_func pepckeys_elements_pack_indexed */
{
  slint_t i;

  if (s->size > d->size) return -1;

#define DO_INDEXED \
  for (i = 0; i < s->size; ++i) \
  { \
    key_copy_at(s->keys, READ_INDEX(i), &d->elements[WRITE_INDEX(i)].key, 0); \
    data_copy_at(s, READ_INDEX(i), &d->elements[WRITE_INDEX(i)], 0); \
  }

  if (rindx == NULL)
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }

  } else
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }
  }

#undef DO_INDEXED

  return 0;
}


slint_t pepckeys_elements_pack(elements_t *s, packed_elements_t *d) /* pepckeys_sl_proto, sl_func pepckeys_elements_pack */
{
  return pepckeys_elements_pack_indexed(s, d, NULL, NULL);
}


slint_t pepckeys_elements_unpack_indexed(packed_elements_t *s, elements_t *d, slindex_t *rindx, slindex_t *windx) /* pepckeys_sl_proto, sl_func pepckeys_elements_unpack_indexed */
{
  slint_t i;

  if (s->size > d->size) return -1;

#define DO_INDEXED \
  for (i = 0; i < s->size; ++i) \
  { \
    key_copy_at(&s->elements[READ_INDEX(i)].key, 0, d->keys, WRITE_INDEX(i)); \
    data_copy_at(&s->elements[READ_INDEX(i)], 0, d, WRITE_INDEX(i)); \
  }

  if (rindx == NULL)
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   _i_
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }

  } else
  {
    if (windx == NULL)
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  _i_
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX

    } else
    {
#define READ_INDEX(_i_)   rindx[_i_]
#define WRITE_INDEX(_i_)  windx[_i_]
       DO_INDEXED
#undef READ_INDEX
#undef WRITE_INDEX
    }
  }

#undef DO_INDEXED

  return 0;
}


slint_t pepckeys_elements_unpack(packed_elements_t *s, elements_t *d) /* pepckeys_sl_proto, sl_func pepckeys_elements_unpack */
{
  return pepckeys_elements_unpack_indexed(s, d, NULL, NULL);
}


slint_t pepckeys_elements_unpack_keys(packed_elements_t *s, slkey_t *k) /* pepckeys_sl_proto, sl_func pepckeys_elements_unpack_keys */
{
  slint_t i;

  for (i = 0; i < s->size; ++i) key_copy_at(&s->elements[i].key, 0, k, i);
}
