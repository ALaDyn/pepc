/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepcparts_mergep_heap.c
 *  timestamp: 2010-01-05 10:54:06 +0100
 *  
 */


#include "sl_common.h"


/* - heap gleich im dst-Feld -> heap umdrehen (Min bewegt sich aufs Ende zu)?
*/

#define HEAPIFY_UP() \
  while (j > 0) \
  { \
    k = (j - 1) / 2; \
    if (key_pure_cmp_le(heap_keys[k], heap_keys[j])) break; \
    hkx = heap_keys[k]; heap_keys[k] = heap_keys[j]; heap_keys[j] = hkx; \
    hsx = heap_sources[k]; heap_sources[k] = heap_sources[j]; heap_sources[j] = hsx; \
    j = k; \
  } \

#define HEAPIFY_DOWN() \
  while ((k = 2 * j + 1) < heap_size) \
  { \
    if (k + 1 < heap_size) \
    if (key_pure_cmp_gt(heap_keys[k], heap_keys[k + 1])) ++k; \
    if (key_pure_cmp_gt(heap_keys[k], heap_keys[j])) break; \
    hkx = heap_keys[k]; heap_keys[k] = heap_keys[j]; heap_keys[j] = hkx; \
    hsx = heap_sources[k]; heap_sources[k] = heap_sources[j]; heap_sources[j] = hsx; \
    j = k; \
  }


slint pepcparts_mergep_heap(elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts) /* pepcparts_sl_proto, sl_func pepcparts_mergep_heap */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k, l;
  
  slindex_t local_displs[p], local_counts[p];


  memcpy(local_displs, displs, p * sizeof(slindex_t));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(slindex_t));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_get_pure(s->keys[local_displs[i]]);
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_UP();
    HEAPIFY_DOWN();
  }

  l = 0;

  while (heap_size > 0)
  {
    i = heap_sources[0];
    
    /* copy min element */
    elem_copy_at(s, local_displs[i], d, l);

    ++l;

    --local_counts[i];
    ++local_displs[i];

    if (local_counts[i] > 0) heap_keys[0] = key_get_pure(s->keys[local_displs[i]]);
    else
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }
    
    j = 0;
    HEAPIFY_DOWN();
  }

  return 0;
}


slint pepcparts_mergep_heap_unpack(packed_elements_t *s, elements_t *d, slint_t p, slindex_t *displs, slindex_t *counts) /* pepcparts_sl_proto, sl_func pepcparts_mergep_heap_unpack */
{
  slkey_pure_t heap_keys[p], hkx;
  slint_t heap_sources[p], hsx;
  slint_t heap_size = 0;

  slint_t i, j, k, l;
  
  slindex_t local_displs[p], local_counts[p];


  memcpy(local_displs, displs, p * sizeof(slindex_t));

  if (counts == NULL)
  {
    for (i = 0; i < p - 1; ++i) local_counts[i] = local_displs[i + 1] - local_displs[i];
    local_counts[p - 1] = s->size - local_displs[p - 1];

  } else memcpy(local_counts, counts, p * sizeof(slindex_t));

  for (i = 0; i < p; ++i)
  if (local_counts[i] > 0)
  {
    heap_keys[heap_size] = key_get_pure(s->elements[local_displs[i]].key);
    heap_sources[heap_size] = i;

    j = heap_size;
    ++heap_size;

    HEAPIFY_UP();
    HEAPIFY_DOWN();
  }

  l = 0;

  while (heap_size > 0)
  {
    i = heap_sources[0];
    
    /* copy min element */
    pelem_unpack_at(s, local_displs[i], d, l);

#ifdef SL_INDEX
    d->indices[l] = local_displs[i];
#endif

    ++l;

    --local_counts[i];
    ++local_displs[i];

    if (local_counts[i] > 0) heap_keys[0] = key_get_pure(s->elements[local_displs[i]].key);
    else
    {
      heap_keys[0] = heap_keys[heap_size - 1];
      heap_sources[0] = heap_sources[heap_size - 1];
      --heap_size;
    }
    
    j = 0;
    HEAPIFY_DOWN();
  }

  return 0;
}
