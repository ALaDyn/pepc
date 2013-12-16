/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/splitx.c
 *  
 */


#include "sl_common.h"


slint_t splitx_radix(elements_t *s, elements_t *sx, slint_t nclasses, slint_t shl, slint_t *counts) /* sl_proto, sl_func splitx_radix */
{
#define max_nclasses  nclasses

  slkey_pure_t bit_mask = nclasses - 1;

  slint_t i, j, k, c[max_nclasses];
  elements_t xi, end, parts[max_nclasses];

  elem_assign_at(s, s->size, &end);
  
  
  if (counts == NULL)
  {
    counts = c;
  
    for (i = 0; i < nclasses; i++) counts[i] = 0;

    for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++counts[key_radix_key2class(key_purify(*xi.keys), shl, bit_mask)];
  }

  elem_assign(s, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], counts[i - 1], &parts[i]);;

  elem_assign(s, &end);
  for (i = 0; i < nclasses; i++)
  {
    elem_add(&end, counts[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = key_radix_key2class(key_purify(*xi.keys), shl, bit_mask);

      while (j != i)
      {
        k = key_radix_key2class(key_purify(*parts[j].keys), shl, bit_mask);

        if (k != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = k;
      }

      elem_inc(&xi);
    }
  }

  return 0;
}

#undef max_nclasses


slint split2_lt_ge(elements_t *s, slkey_pure_t *k, elements_t *t) /* sl_proto, sl_func split2_lt_ge */
{
  elements_t low, high;

  elem_assign(s, &low);
  elem_assign_at(s, s->size - 1, &high);

  while (1)
  {
    while (low.keys < high.keys)
    if (key_pure_cmp_lt(key_purify(*low.keys), *k)) elem_inc(&low); else break;

    while (low.keys < high.keys)
    if (key_pure_cmp_ge(key_purify(*high.keys), *k)) elem_dec(&high); else break;

    if (low.keys >= high.keys) break;

    elem_copy(&low, t);
    elem_copy(&high, &low);
    elem_copy(t, &high);
    elem_inc(&low);
    elem_dec(&high);
  }

  return (low.keys - s->keys) + (key_pure_cmp_lt(key_purify(*low.keys), *k));
}


slint split2_le_gt(elements_t *s, slkey_pure_t *k, elements_t *t) /* sl_proto, sl_func split2_le_gt */
{
  elements_t low, high;

  elem_assign(s, &low);
  elem_assign_at(s, s->size - 1, &high);

  while (1)
  {
    while (low.keys < high.keys)
    if (key_pure_cmp_le(key_purify(*low.keys), *k)) elem_inc(&low); else break;

    while (low.keys < high.keys)
    if (key_pure_cmp_gt(key_purify(*high.keys), *k)) elem_dec(&high); else break;

    if (low.keys >= high.keys) break;

    elem_copy(&low, t);
    elem_copy(&high, &low);
    elem_copy(t, &high);
    elem_inc(&low);
    elem_dec(&high);
  }

  return (low.keys - s->keys) + (key_pure_cmp_le(key_purify(*low.keys), *k));
}


slint split3_lt_eq_gt(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle) /* sl_proto, sl_func split3_lt_eq_gt */
{
  elements_t low, high, middle;

/*  slint i;*/

  *nlt = *nle = 0;

  if (s->size <= 0) return 0;

  elem_assign(s, &low);
  elem_assign(s, &middle);
  elem_assign_at(s, s->size, &high);

  while (low.keys < high.keys)
  if (key_pure_cmp_eq(key_purify(*low.keys), *k))
  {
    elem_xchange(&middle, &low, t);
    elem_inc(&middle);
    elem_inc(&low);

  } else if (key_pure_cmp_lt(key_purify(*low.keys), *k))
  {
    elem_inc(&low);

  } else
  {
    elem_dec(&high);
    elem_xchange(&high, &low, t);
  }

  *nle = low.keys - s->keys;
  *nlt = *nle - (middle.keys - s->keys);

  return 0;
}


slint split3_lt_eq_gt_old(elements_t *s, slkey_pure_t *k, elements_t *t, slint *nlt, slint *nle) /* sl_proto, sl_func split3_lt_eq_gt_old */
{
  elements_t low, high, middle;

/*  slint i;*/

  *nlt = *nle = 0;

  if (s->size <= 0) return 0;

  elem_assign(s, &low);
  elem_assign(s, &middle);
  elem_assign_at(s, s->size - 1, &high);

  while (1)
  {
    while (low.keys < high.keys)
    if (key_pure_cmp_lt(key_purify(*low.keys), *k)) elem_inc(&low);
    else if (key_pure_cmp_eq(key_purify(*low.keys), *k))
    {
      elem_copy(&middle, t);
      elem_copy(&low, &middle);
      elem_copy(t, &low);
      elem_inc(&low);
      elem_inc(&middle);
    } else break;

    while (low.keys < high.keys)
    if (key_pure_cmp_gt(key_purify(*high.keys), *k)) elem_dec(&high); else break;

    if (low.keys >= high.keys) break;

    if (key_pure_cmp_eq(key_purify(*high.keys), *k))
    {
      elem_copy(&middle, t);
      elem_copy(&high, &middle);
      elem_inc(&middle);

    } else elem_copy(&high, t);

    if (low.keys >= middle.keys)
    {
      elem_copy(&low, &high);
      elem_copy(t, &low);

    } else elem_copy(t, &high);

    elem_inc(&low);
    elem_dec(&high);
  }

  *nle = (low.keys - s->keys) + (key_pure_cmp_le(key_purify(*low.keys), *k));
  *nlt = *nle - (middle.keys - s->keys) - (key_pure_cmp_eq(key_purify(*middle.keys), *k));

  return 0;
}


slint split2_b(elements_t *s, elements_t *sx, slkey_pure_t bmask) /* sl_proto, sl_func split2_b */
{
  elements_t xl, xh;

  elem_assign(s, &xl);
  elem_assign_at(s, s->size - 1, &xh);

  while (1)
  {
    while (xl.keys < xh.keys)
    if (key_purify(*xl.keys) & bmask) break; else elem_inc(&xl);

    while (xl.keys < xh.keys)
    if (key_purify(*xh.keys) & bmask) elem_dec(&xh); else break;

    if (xl.keys >= xh.keys) break;

    elem_xchange(&xl, &xh, sx);
    elem_inc(&xl);
    elem_dec(&xh);
  }

  return xl.keys - s->keys + ((key_purify(*xl.keys) & bmask) == 0);
}


slint splitk_k2c_af(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data) /* sl_proto, sl_func splitk_k2c_af */
{
  slint i;
  elements_t xi, end;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  elements_t *parts = NULL;
#else
  elements_t parts[k];
#endif

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  parts = z_alloca(k, sizeof(elements_t));
#endif

  /* compute the target of every class */
  elem_assign_at(s, c[0], &parts[0]);
  for (i = 1; i < k; i++) elem_assign_at(&parts[i - 1], c[i], &parts[i]);

  /* permute the keys home */
  elem_assign_at(s, s->size, &end);
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_add(&xi, c[i]))
  {
    while (1)
    {
      i = (k2c)(xi.keys, 0, k2c_data);

      elem_dec(&parts[i]);

      if (xi.keys >= parts[i].keys) break;

      elem_xchange(&parts[i], &xi, sx);
    }
  }

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  z_freea(parts);
#endif

  return 0;
}


slint splitk_k2c(elements_t *s, elements_t *sx, slint k, slint *c, k2c_func k2c, void *k2c_data) /* sl_proto, sl_func splitk_k2c */
{
  slint i, j, l;
  elements_t xi, end;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  elements_t *parts = NULL;
#else
  elements_t parts[k];
#endif

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  parts = z_alloca(k, sizeof(elements_t));
#endif

  /* compute the target of every class */
  elem_assign(s, &parts[0]);
  for (i = 1; i < k; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

  elem_assign(s, &end);
  for (i = 0; i < k; i++)
  {
    elem_add(&end, c[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = (k2c)(xi.keys, 0, k2c_data);

      while (j != i)
      {
        l = (k2c)(parts[j].keys, 0, k2c_data);

        if (l != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = l;
      }

      elem_inc(&xi);
    }
  }

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  z_freea(parts);
#endif

  return 0;
}


slint splitk_k2c_count(elements_t *s, slint k, slint *c, k2c_func k2c, void *k2c_data) /* sl_proto, sl_func splitk_k2c_count */
{
  slint i;

  elements_t xi, end;

  elem_assign(s, &xi);
  elem_assign_at(s, s->size, &end);

  for (i = 0; i < k; i++) c[i] = 0;

  while (xi.keys < end.keys)
  {
    i = (k2c)(xi.keys, 0, k2c_data);

    c[i]++;

    elem_inc(&xi);
  }

  return 0;
}
