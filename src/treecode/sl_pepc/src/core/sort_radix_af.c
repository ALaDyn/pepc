/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_radix_af.c
 *  
 */

/*
 - a fixed size of the local arrays is very import for performance, therefor we set the size to a given maximum (at compiletime)
 - all other parameters, like highest-/lowest-bit or bit-width can be chosen da runtime
*/


#include "sl_common.h"


#ifdef key_integer


#define insertsort
/*#define insertsort_finalize*/
/*#define insertsort_finalize_adaptive*/


/* insert sort */
slint rs_rec_insertsort_af(elements_t *s, elements_t *sx, slint rhigh, slint rlow) /* sl_func rs_rec_insertsort_af */
{
  slint i, j;

  slkey_pure_t bmask = 0;

  if (rhigh - rlow + 1 <= key_radix_high) bmask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);

  bmask = (bmask - 1) << rlow;

  for (i = 1; i < s->size; i++)
  {
    if (key_pure_cmp_lt(key_purify(s->keys[i]) & bmask, key_purify(s->keys[i - 1]) & bmask))
    {
      j = i - 1;
      elem_copy_at(s, i, sx, 0);

      do
      {
        elem_copy_at(s, j, s, j + 1);
        if (--j < 0) break;

      } while (key_pure_cmp_lt(key_purify(*sx->keys) & bmask, key_purify(s->keys[j]) & bmask));

      elem_copy_at(sx, 0, s, j + 1);
    }
  }

  return 0;
}


slint rs_rec_af(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth, slint *finalize) /* sl_func rs_rec_af */
{
#define max_nclasses (z_powof2_typed(sort_radix_width_max, slkey_pure_t))

  slkey_pure_t bit_mask, nclasses;

  slint i, current_width, c[max_nclasses];
  elements_t xi, end, parts[max_nclasses];

  elem_assign_at(s, s->size, &end);

  current_width = z_min(rwidth, rhigh - rlow + 1);
  rhigh -= current_width - 1;

  nclasses = z_powof2_typed(current_width, slkey_pure_t);
  bit_mask = nclasses - 1;


  /* zero all counter */
  for (i = 0; i < nclasses; i++) c[i] = 0;

  /* count the number of elements in every class */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++c[key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];

  /* compute the target of every class */
  elem_assign_at(s, c[0], &parts[0]);
  parts[0].size = c[0];
  for (i = 1; i < nclasses; i++)
  {
    elem_assign_at(&parts[i - 1], c[i], &parts[i]);
    parts[i].size = c[i];
  }

  /* permute the keys home */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_add(&xi, c[i]))
  {
    while (1)
    {
      i = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

      elem_dec(&parts[i]);

      if (xi.keys >= parts[i].keys) break;

      elem_xchange(&parts[i], &xi, sx);
    }
  }

  --rhigh;

  if (rhigh >= rlow)
  {
    elem_assign(s, &xi);
    for (i = 0; i < nclasses; i++)
    {
      xi.size = c[i];

#ifdef insertsort
      if (xi.size > sort_radix_threshold_rec) rs_rec_af(&xi, sx, rhigh, rlow, rwidth, finalize);

 #ifdef insertsort_finalize

  #ifdef insertsort_finalize_adaptive
      else if (xi.size > 1) *finalize = 1;
  #endif /* insertsort_finalize_adaptive */

 #else /* insertsort_finalize */
      else if (xi.size > 1) rs_rec_insertsort_af(&xi, sx, rhigh, rlow);
 #endif /* insertsort_finalize */

#else /* insertsort */
      if (xi.size > 1) rs_rec(&xi, sx, rhigh, rlow, rwidth, finalize);
#endif /* insertsort */

      elem_add(&xi, c[i]);
    }
  }

  return 0;
}

#undef max_nclasses


slint sort_radix_af(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth) /* sl_proto, sl_func sort_radix_af */
{
  elements_t _sx;

  slint finalize = 1;

#ifdef insertsort_finalize_adaptive
  finalize = 0;
#endif /* insertsort_finalize_adaptive */

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth <= 0) rwidth = sort_radix_width_default;

  rs_rec_af(s, sx, rhigh, rlow, z_min(rwidth, sort_radix_width_max), &finalize);

#ifdef insertsort_finalize
  if (sort_radix_threshold_rec > 1 && finalize) rs_rec_insertsort_af(s, sx, rhigh, rlow);
#endif /* insertsort_finalize */

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_radix);

  return 0;
}


#undef insertsort
#undef insertsort_finalize
#undef insertsort_finalize_adaptive


#endif /* key_integer */
