/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/elements_packed.c
 *  timestamp: 2011-02-10 21:20:49 +0100
 *  
 */


#include "sl_common.h"


slint_t elements_alloc_packed(packed_elements_t *s, slint_t nelements) /* sl_proto, sl_func elements_alloc_packed */
{
  if (s == NULL) return -1;

  s->size = s->max_size = 0;
  s->elements = NULL;

  if (nelements == 0) return 0;

  s->size = nelements;
  
  s->elements = z_alloc(nelements, pelem_byte);

  if (s->elements != NULL) return 0;

  elements_free_packed(s);

  return -1;
}


slint_t elements_free_packed(packed_elements_t *s) /* sl_proto, sl_func elements_free_packed */
{
  if (s == NULL) return -1;

  z_free(s->elements);

  s->size = s->max_size = 0;

  return 0;
}


slint_t elements_pack_indexed(elements_t *s, packed_elements_t *d, slindex_t *rindx, slindex_t *windx) /* sl_proto, sl_func elements_pack_indexed */
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


slint_t elements_pack(elements_t *s, packed_elements_t *d) /* sl_proto, sl_func elements_pack */
{
  return elements_pack_indexed(s, d, NULL, NULL);
}


slint_t elements_unpack_indexed(packed_elements_t *s, elements_t *d, slindex_t *rindx, slindex_t *windx) /* sl_proto, sl_func elements_unpack_indexed */
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


slint_t elements_unpack(packed_elements_t *s, elements_t *d) /* sl_proto, sl_func elements_unpack */
{
  return elements_unpack_indexed(s, d, NULL, NULL);
}


slint_t elements_unpack_keys(packed_elements_t *s, slkey_t *k) /* sl_proto, sl_func elements_unpack_keys */
{
  slint_t i;

  for (i = 0; i < s->size; ++i) key_copy_at(&s->elements[i].key, 0, k, i);
  
  return 0;
}
