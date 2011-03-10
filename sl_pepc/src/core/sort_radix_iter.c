/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sort_radix_iter.c
 *  timestamp: 2011-02-08 21:32:53 +0100
 *  
 */


#include "sl_common.h"


#ifdef key_integer


/* configure tuneable */
#ifdef SL_TUNEABLE

 /* sort_radix_threshold_iter */
 int tuneable_sort_radix_threshold_iter = sort_radix_threshold_iter; /* sl_var tuneable_sort_radix_threshold_iter */
 #undef sort_radix_threshold_iter
 #define sort_radix_threshold_iter tuneable_sort_radix_threshold_iter

#endif


#define insertsort
/*#define insertsort_finalize*/


/* insert sort */
slint rs_iter_insertsort(elements_t *s, elements_t *sx, slint presorted, slint rhigh, slint rlow) /* sl_func rs_iter_insertsort */
{
  slint i, j;

  slkey_pure_t class_mask = 0, bmask = 0;


  if (presorted && (rhigh < key_radix_high)) class_mask = ~(z_powof2_typed(rhigh + 1, slkey_pure_t) - 1);

  if (rhigh - rlow + 1 <= key_radix_high) bmask = z_powof2_typed(rhigh - rlow + 1, slkey_pure_t);

  bmask = class_mask | (bmask - 1) << rlow;

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


slint rs_iter(elements_t *s, elements_t *sx, slint presorted, slint rhigh, slint rlow, slint rwidth) /* sl_func rs_iter */
{
#define max_nclasses (z_powof2_typed(sort_radix_width_max, slkey_pure_t))

  slkey_pure_t class_mask, bit_mask, current_class, nclasses;

  slint i, j, k, current_width, c[max_nclasses];
  elements_t xi, xj, end, parts[max_nclasses];
  
  slint rhigh_old = rhigh;


  if (presorted && (rhigh < key_radix_high)) class_mask = ~(z_powof2_typed(rhigh + 1, slkey_pure_t) - 1);
  else class_mask = 0;

  elem_assign_at(s, s->size, &end);

  while (rhigh >= rlow)
  {
    current_width = z_min(rwidth, rhigh - rlow + 1);
    rhigh -= current_width - 1;

    nclasses = z_powof2_typed(current_width, slkey_pure_t);
    bit_mask = nclasses - 1;

    elem_assign(s, &xi);

    while (xi.keys < end.keys)
    {
      elem_assign(&xi, &xj);

      current_class = class_mask & key_purify(*xi.keys);

      /* zero all counters */
      for (i = 0; i < nclasses; ++i) c[i] = 0;

      /* counting the number of sub-class elements in the current class */
      do
      {
        i = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

        ++c[i];

        elem_inc(&xi);

        /* stop if we are at the end */
        if (xi.keys >= end.keys) break;

      } while (current_class == (class_mask & key_purify(*xi.keys)));

      /* if the current class has enough elements, perform split */
      if (xi.keys - xj.keys > 
#ifdef insertsort
        sort_radix_threshold_iter
#else
        1
#endif
        )
      {
        /* compute the target of every sub-class */
        elem_assign(&xj, &parts[0]);
        for (i = 1; i < nclasses; i++) elem_assign_at(&parts[i - 1], c[i - 1], &parts[i]);;

        for (i = 0; i < nclasses; i++)
        {
          elem_add(&xj, c[i]);
          elem_assign(&parts[i], &xi);

          while (xi.keys < xj.keys)
          {
            j = key_radix_key2class(key_purify(*xi.keys), rhigh, bit_mask);

            while (j != i)
            {
              k = key_radix_key2class(key_purify(*parts[j].keys), rhigh, bit_mask);

              if (k != j) elem_xchange(&xi, &parts[j], sx);

              elem_inc(&parts[j]);

              j = k;
            }

            elem_inc(&xi);
          }
        }
      }
#ifdef insertsort
 #ifndef insertsort_finalize
        else
      {
        xj.size = xi.keys - xj.keys;

        rs_iter_insertsort(&xj, sx, presorted, rhigh_old, rlow);
      }
 #endif
#endif
    }

    class_mask |= (bit_mask << rhigh);

    --rhigh;
  }

  return 0;
}


slint sort_radix_iter(elements_t *s, elements_t *sx, slint presorted, slint rhigh, slint rlow, slint rwidth) /* sl_proto, sl_func sort_radix_iter */
{
  elements_t _sx;

  if (s == NULL) return -1;

  if (s->size < 2) return 0;

  rti_tstart(rti_tid_sort_radix_iter);

  if (sx == NULL || sx->size < 1)
  {
    sx = &_sx;
    elements_alloc(sx, 1, SLCM_ALL);
  }

  if (rhigh < 0) rhigh = key_radix_high;
  if (rlow < 0) rlow = key_radix_low;
  if (rwidth <= 0) rwidth = sort_radix_width_default;

  rs_iter(s, sx, presorted, rhigh, rlow, z_min(rwidth, sort_radix_width_max));

#ifdef insertsort
 #ifdef insertsort_finalize
  if (sort_radix_threshold_iter > 1) rs_iter_insertsort(s, sx, presorted, rhigh, rlow);
 #endif
#endif

  if (sx == &_sx) elements_free(sx);

  rti_tstop(rti_tid_sort_radix_iter);

  return 0;
}


#endif /* key_integer */
