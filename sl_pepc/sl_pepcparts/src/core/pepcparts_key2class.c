/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/key2class.c
 *  timestamp: 2009-10-12 09:08:26 +0200
 *  
 */


#include "sl_common.h"


slint pepcparts_key2class_equal(slkey_t *k, slint i, void *ci) /* pepcparts_sl_proto, sl_func pepcparts_key2class_equal */
{
  slint *param = ci;

  return i % *param;
}


slint pepcparts_key2class_split(slkey_t *k, slint i, void *ci) /* pepcparts_sl_proto, sl_func pepcparts_key2class_split */
{
  slint *param = ci;

  return i / *param;
}


slint pepcparts_key2class_split_keys(slkey_t *k, slint i, void *ci) /* pepcparts_sl_proto, sl_func pepcparts_key2class_split_keys */
{
  slint *param = ci;

#ifdef key_integer
  return key_purify(*k) / *param;
#else
  return 0;
#endif
}


slint pepcparts_key2class_random(slkey_t *k, slint i, void *ci) /* pepcparts_sl_proto, sl_func pepcparts_key2class_random */
{
  slint *param = ci;

  /* FIXME: but how? */
/*  srand(key_2int(*k));*/

  return iround(sl_rand01() * (double) (*param - 1));
}


slint pepcparts_key2class_ci_nocounts(slkey_t *k, slint i, void *ci) /* pepcparts_sl_proto, sl_func pepcparts_key2class_ci_nocounts */
{
  slint c = 0;
  classification_info *_ci = ci;

  /* process all split-keys */
  while (c < _ci->nclasses - 1)
  {
    /* if the key is less than the split-key[c] */
    if (key_pure_cmp_lt(key_purify(*k), _ci->keys[c])) break;

    ++c;
  }

  return c;
}


slint pepcparts_key2class_ci(slkey_t *k, slint i, void *ci) /* pepcparts_sl_proto, sl_func pepcparts_key2class_ci */
{
  slint c = 0;
  classification_info *_ci = ci;

  /* process all split-keys */
  while (c < _ci->nclasses - 1)
  {
    /* if the key is less than the split-key[c] */
    if (key_pure_cmp_lt(key_purify(*k), _ci->keys[c])) break;
    /* else if the key is equal to the split-key[c] */
    else if ((key_pure_cmp_eq(_ci->keys[c], key_purify(*k))) && (_ci->counts[c] > 0))
    {
      --_ci->counts[c];
      break;
    }
    ++c;
  }

  return c;
}
