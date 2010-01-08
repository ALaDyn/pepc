/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepckeys_sort_counting.c
 *  timestamp: 2009-10-12 09:08:26 +0200
 *  
 */


#include "sl_common.h"


void pepckeys_make_counts(elements_t *s, slint_t ncounts, slint_t *counts)  /* sl_func pepckeys_make_counts */
{
  slint_t i;
  key_type_c *k;

  memset(counts, 0, ncounts * sizeof(slint_t));

  for (i = 0, k = s->keys; i < s->size; ++i, ++k) ++counts[key_purify(*k)];
}


void pepckeys_make_displs(slint_t ncounts, slint_t *counts, slint_t *displs)  /* sl_func pepckeys_make_displs */
{
  slint_t i;

  displs[0] = 0;
  for (i = 1; i < ncounts; ++i) displs[i] = displs[i - 1] + counts[i - 1];
}


void pepckeys_make_counts2displs(slint_t ncounts, slint_t *countsdispls)  /* sl_func pepckeys_make_counts2displs */
{
  slint_t i, s, t;

  s = 0;
  for (i = 0; i < ncounts; ++i)
  {
    t = countsdispls[i];
    countsdispls[i] = s;
    s += t;
  }
}


slint_t pepckeys_sort_counting_use_displs(elements_t *s, elements_t *d, slint_t ndispls, slint_t *displs)  /* pepckeys_sl_proto, sl_func pepckeys_sort_counting_use_displs */
{
  slint_t i;

  for (i = 0; i < s->size; ++i)
  {
    elem_copy_at(s, i, d, displs[key_purify(s->keys[i])]);
    ++displs[key_purify(s->keys[i])];
  }

  d->size = s->size;
  
  return 0;
}


slint_t pepckeys_sort_counting_use_counts(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts)  /* pepckeys_sl_proto, sl_func pepckeys_sort_counting_use_counts */
{
  slint_t r;
  slint_t *displs = NULL;

  if (counts == NULL)
  {
    displs = sl_alloc(ncounts, sizeof(slint_t));
    pepckeys_make_counts(s, ncounts, displs);
    pepckeys_make_counts2displs(ncounts, displs);

  } else
  {
    displs = counts;

    if (ncounts < 0)
    {
      ncounts *= -1;
      displs += ncounts;
      pepckeys_make_displs(ncounts, counts, displs);

    } else pepckeys_make_counts2displs(ncounts, displs);
  }
  
  r = pepckeys_sort_counting_use_displs(s, d, ncounts, displs);

  if (counts == NULL) sl_free(displs);

  return r;
}


slint_t pepckeys_sort_counting_get_counts(elements_t *s, elements_t *d, slint_t ncounts, slint_t *counts)  /* pepckeys_sl_proto, sl_func pepckeys_sort_counting_get_counts */
{
  slint_t r;
  slint_t *displs = NULL;
  
  if (counts == NULL)
  {
    displs = sl_alloc(ncounts, sizeof(slint_t));
    pepckeys_make_counts(s, ncounts, displs);
    pepckeys_make_counts2displs(ncounts, displs);

  } else
  {
    if (ncounts < 0)
    {
      ncounts *= -1;
      displs = counts + ncounts;

    } else displs = sl_alloc(ncounts, sizeof(slint_t));

    pepckeys_make_counts(s, ncounts, counts);
    pepckeys_make_displs(ncounts, counts, displs);
  }

  r = pepckeys_sort_counting_use_displs(s, d, ncounts, displs);
  
  if (counts == NULL || displs != counts + ncounts) sl_free(displs);
  
  return r;
}


slint_t pepckeys_sort_counting(elements_t *s, elements_t *d, slint_t ncounts)  /* pepckeys_sl_proto, sl_func pepckeys_sort_counting */
{
  return pepckeys_sort_counting_use_counts(s, d, ncounts, NULL);
}
