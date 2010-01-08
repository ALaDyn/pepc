/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepckeys_sort_radix.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/*
 - a fixed size of the local arrays is very import for performance, therefor we set the size to a given maximum (at compiletime)
 - all other parameters, like highest-/lowest-bit or bit-width can be chosen at runtime
*/


#include "sl_common.h"


#ifdef key_integer


/* configure tuneable */
#ifdef SL_TUNEABLE

 /* sort_radix_threshold_rec */
 int pepckeys_tuneable_sort_radix_threshold_rec = sort_radix_threshold_rec; /* sl_var pepckeys_tuneable_sort_radix_threshold_rec */
 #undef sort_radix_threshold_rec
 #define sort_radix_threshold_rec pepckeys_tuneable_sort_radix_threshold_rec

#endif


#define radix_low                 0
#define radix_high                (sizeof(slkey_pure_t) * 8 - 1)

#define radix_key2class(k, x, y)  (((k) >> (x)) & y)


#define insertsort
/*#define insertsort_finalize*/
/*#define insertsort_finalize_adaptive*/


/* insert sort */
slint pepckeys_rs_rec_insertsort(elements_t *s, elements_t *sx, slint rhigh, slint rlow) /* sl_func pepckeys_rs_rec_insertsort */
{
  slint i, j;

  slkey_pure_t bmask = 0;

  if (rhigh - rlow + 1 <= radix_high) bmask = powof2_typed(rhigh - rlow + 1, slkey_pure_t);

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


slint pepckeys_rs_rec(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth, slint *finalize) /* sl_func pepckeys_rs_rec */
{
#define max_nclasses (powof2_typed(sort_radix_width_max, slkey_pure_t))

  slkey_pure_t bit_mask, nclasses;

  slint i, j, k, current_width, c[max_nclasses];
  elements_t xi, end, parts[max_nclasses];

  elem_assign_at(s, s->size, &end);

  current_width = xmin(rwidth, rhigh - rlow + 1);
  rhigh -= current_width - 1;

  nclasses = powof2_typed(current_width, slkey_pure_t);
  bit_mask = nclasses - 1;


  /* zero all counter */
  for (i = 0; i < nclasses; i++) c[i] = 0;

  /* count the number of elements in every class */
  for (elem_assign(s, &xi); xi.keys < end.keys; elem_inc(&xi)) ++c[radix_key2class(key_purify(*xi.keys), rhigh, bit_mask)];

  /* compute the target of every class */
  elem_assign(s, &parts[0]);
  for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

  /* split the elements */
  elem_assign(s, &end);
  for (i = 0; i < nclasses; i++)
  {
    elem_add(&end, c[i]);

    elem_assign(&parts[i], &xi);

    while (xi.keys < end.keys)
    {
      j = radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

      while (j != i)
      {
        k = radix_key2class(key_purify(*parts[j].keys), rhigh, bit_mask);

        if (k != j) elem_xchange(&xi, &parts[j], sx);

        elem_inc(&parts[j]);

        j = k;
      }

      elem_inc(&xi);
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
      if (xi.size > sort_radix_threshold_rec) pepckeys_rs_rec(&xi, sx, rhigh, rlow, rwidth, finalize);

 #ifdef insertsort_finalize

  #ifdef insertsort_finalize_adaptive
      else if (xi.size > 1) *finalize = 1;
  #endif /* insertsort_finalize_adaptive */

 #else /* insertsort_finalize */
      else if (xi.size > 1) pepckeys_rs_rec_insertsort(&xi, sx, rhigh, rlow);
 #endif /* insertsort_finalize */

#else /* insertsort */
      if (xi.size > 1) pepckeys_rs_rec(&xi, sx, rhigh, rlow, rwidth, finalize);
#endif /* insertsort */

      elem_add(&xi, c[i]);
    }
  }

  return 0;
}


slint pepckeys_sort_radix(elements_t *s, elements_t *sx, slint rhigh, slint rlow, slint rwidth) /* pepckeys_sl_proto, sl_func pepckeys_sort_radix */
{
  elements_t _sx;

  slint finalize = 1;

#ifdef insertsort_finalize_adaptive
  finalize = 0;
#endif /* insertsort_finalize_adaptive */

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix);

  if (sx == NULL)
  {
    sx = &_sx;
    pepckeys_elements_alloc(sx, 1, 1, 1);
  }

  if (rhigh < 0) rhigh = radix_high;
  if (rlow < 0) rlow = radix_low;
  if (rwidth <= 0) rwidth = sort_radix_width_default;

  pepckeys_rs_rec(s, sx, rhigh, rlow, xmin(rwidth, sort_radix_width_max), &finalize);

#ifdef insertsort_finalize
  if (sort_radix_threshold_rec > 1 && finalize) pepckeys_rs_rec_insertsort(s, sx, rhigh, rlow);
#endif /* insertsort_finalize */

  if (sx == &_sx) pepckeys_elements_free(sx);

  rti_tstop(rti_tid_sort_radix);

  return 0;
}


#endif /* key_integer */
