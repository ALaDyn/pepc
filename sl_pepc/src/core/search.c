/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/search.c
 *  timestamp: 2011-02-24 21:37:35 +0100
 *  
 */

/* simple routines for searching in sorted elements */


#include "sl_common.h"


slint sl_search_sequential_lt(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_sequential_lt */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign(s->keys, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s != _e)
  if (key_cmp_lt(*_s, *k)) key_inc(_s); else break;

  return s->size - (_e - _s);
}


slint sl_search_sequential_le(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_sequential_le */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign(s->keys, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s != _e)
  if (key_cmp_le(*_s, *k)) key_inc(_s); else break;

  return s->size - (_e - _s);
}


slint sl_search_sequential_gt(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_sequential_gt */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - 1, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s != _e)
  if (key_cmp_gt(*_s, *k)) key_dec(_s); else break;

  return s->size - (_s - _e);
}


slint sl_search_sequential_ge(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_sequential_ge */
{
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - 1, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s != _e)
  if (key_cmp_ge(*_s, *k)) key_dec(_s); else break;

  return s->size - (_s - _e);
}


/* max number of elements less than k */
/* index i with s->keys[i-1] < k <= s->keys[i] */
/* max. comparisons: \lceil \log(2,n+1) \rceil */
slint sl_search_binary_lt(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_binary_lt */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_cmp_le(*k, s->keys[mi])) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


/* max number of elements less than or equal k */
/* index i with s->keys[i-1] <= k < s->keys[i] */
/* max. comparisons: \lceil \log(2,n+1) \rceil */
slint sl_search_binary_le(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_binary_le */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_cmp_lt(*k, s->keys[mi])) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint sl_search_binary_gt(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_binary_gt */
{
  return s->size - sl_search_binary_le(s, k);
}


slint sl_search_binary_ge(elements_t *s, slkey_t *k) /* sl_proto, sl_func sl_search_binary_ge */
{
  return s->size - sl_search_binary_lt(s, k);
}


slint sl_search_binary_lt2(elements_t *s, slkey_pure_t k) /* sl_proto, sl_func sl_search_binary_lt2 */
{
  slint le, ri, mi;

/*  printf("searching for lt %" sl_key_pure_type_fmt " in keys:\n", k);
  elements_print_keys(s);*/

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_le(k, key_purify(s->keys[mi]))) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint sl_search_binary_le2(elements_t *s, slkey_pure_t k) /* sl_proto, sl_func sl_search_binary_le2 */
{
  slint le, ri, mi;

/*  printf("searching for le %" sl_key_pure_type_fmt " in keys:\n", k);
  elements_print_keys(s);*/

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_lt(k, key_purify(s->keys[mi]))) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint sl_search_binary_gt2(elements_t *s, slkey_pure_t k) /* sl_proto, sl_func sl_search_binary_gt2 */
{
  return s->size - sl_search_binary_le2(s, k);
}


slint sl_search_binary_ge2(elements_t *s, slkey_pure_t k) /* sl_proto, sl_func sl_search_binary_ge2 */
{
  return s->size - sl_search_binary_lt2(s, k);
}


slint_t sl_search_binary_lt_bmask(elements_t *s, slkey_pure_t k, slkey_pure_t bmask) /* sl_proto, sl_func sl_search_binary_lt_bmask */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_le(k, key_purify(s->keys[mi]) & bmask)) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint_t sl_search_binary_le_bmask(elements_t *s, slkey_pure_t k, slkey_pure_t bmask) /* sl_proto, sl_func sl_search_binary_le_bmask */
{
  slint le, ri, mi;

  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_lt(k, key_purify(s->keys[mi]) & bmask)) ri = mi - 1;
    else le = mi + 1;
  }

  return le;
}


slint_t sl_search_binary_sign_switch(elements_t *s) /* sl_proto, sl_func sl_search_binary_sign_switch */
{
  slint le, ri, mi;

  if (s->size <= 0) return 0;
  
  if (key_pure_cmp_ge(key_purify(s->keys[0]), 0) && key_pure_cmp_ge(key_purify(s->keys[s->size - 1]), 0)) return s->size;

  if (key_pure_cmp_lt(key_purify(s->keys[0]), 0) && key_pure_cmp_lt(key_purify(s->keys[s->size - 1]), 0)) return 0;
  
  le = 0;
  ri = s->size - 1;

  while (le <= ri)
  {
    mi = (le + ri) / 2;
    if (key_pure_cmp_ge(key_purify(s->keys[mi]), 0)) le = mi + 1;
    else ri = mi - 1;
  }
  
  return le;
}


slint sl_search_hybrid_lt(elements_t *s, slkey_t *k, slint t) /* sl_proto, sl_func sl_search_hybrid_lt */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, t - 1, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s < _e)
  if (key_cmp_lt(*_s, *k)) key_add(_s, t); else break;

  n = (_s - s->keys) - (t - 1);

  elem_assign_at(s, n, &x);
  x.size = (z_min(_s, _e) - s->keys) - n;

  return n + sl_search_binary_lt(&x, k);
}


slint sl_search_hybrid_le(elements_t *s, slkey_t *k, slint t) /* sl_proto, sl_func sl_search_hybrid_le */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, t - 1, _s);
  key_assign_at(s->keys, s->size, _e);

  while (_s < _e)
  if (key_cmp_le(*_s, *k)) key_add(_s, t); else break;

  n = (_s - s->keys) - (t - 1);

  elem_assign_at(s, n, &x);
  x.size = (z_min(_s, _e) - s->keys) - n;

  return n + sl_search_binary_le(&x, k);
}


slint sl_search_hybrid_gt(elements_t *s, slkey_t *k, slint t) /* sl_proto, sl_func sl_search_hybrid_gt */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - t, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s > _e)
  if (key_cmp_gt(*_s, *k)) key_sub(_s, t); else break;

  n = (z_max(_e, _s) - s->keys) + 1;

  elem_assign_at(s, n, &x);
  x.size = ((_s - s->keys) + t) - n;

  return s->size - (n + sl_search_binary_le(&x, k));
}


slint sl_search_hybrid_ge(elements_t *s, slkey_t *k, slint t) /* sl_proto, sl_func sl_search_hybrid_ge */
{
  slint n;
  elements_t x;
  slkey_t *_s, *_e;

  if (s->size <= 0) return 0;

  key_assign_at(s->keys, s->size - t, _s);
  key_assign_at(s->keys, -1, _e);

  while (_s > _e)
  if (key_cmp_ge(*_s, *k)) key_sub(_s, t); else break;

  n = (z_max(_e, _s) - s->keys) + 1;

  elem_assign_at(s, n, &x);
  x.size = ((_s - s->keys) + t) - n;

  return s->size - (n + sl_search_binary_lt(&x, k));
}
