/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sl_common.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


declare_ts_temp
declare_rti_env


slint pepckeys_ilog2c(slint x) /* pepckeys_sl_proto, sl_func pepckeys_ilog2c */
{
  slint l = 0;
  while (x /= 2) ++l;
  return l;
}


slint pepckeys_ilog2f(slint x) /* pepckeys_sl_proto, sl_func pepckeys_ilog2f */
{
  slint l = 0;
  if (x <= 1) return 0;
  --x;
  while (x /= 2) ++l;
  return l + 1;
}


slint pepckeys_print_bits(slint v) /* pepckeys_sl_proto, sl_func pepckeys_sl_proto pepckeys_print_bits */
{
  slint i;

  for (i = sizeof(v) * 8 - 1; i >= 0; i--) printf("%d", (v & (1L << i)) != 0);

  return 0;
}


slint pepckeys_pivot_random(elements_t *s) /* pepckeys_sl_proto, sl_func pepckeys_pivot_random */
{
  return sl_rand() % s->size;
}


/* calculates the intersection of [a0_start,a0_end] with [a1_start,a1_end]
   - storing the intersection in [ia_start,ia_end]
   - returning 'ia_end - ia_start + 1'
     x > 0: intersection consists of x elements
     x = 0: no intersection, areas touching
     x < 0: no intersection, x elements inbetween */
/*slint intersect_areas(slint a0_start, slint a0_end, slint a1_start, slint a1_end, slint *ia_start, slint *ia_end)
{
  slint temp_ia_start, temp_ia_end;

  if (ia_start == NULL) ia_start = &temp_ia_start;
  if (ia_end == NULL) ia_end = &temp_ia_end;

  *ia_start = xmax(a0_start, a1_start);
  *ia_end = xmin(a0_end, a1_end);

  return *ia_end - *ia_start + 1;
}*/

/* calculates the remaining area(s) when discarding [a1_start,a1_end] from [a0_start,a0_end]
   - storing the remaining area(s) in [a2_start,a2_end] and [a3_start,a3_end]
   - 'a2_end < a1_start && a1_end < a3_start' is always true
   - returning the added size of [a2_start,a2_end] and [a3_start,a3_end] */
/*slint subtract_areas(slint a0_start, slint a0_end, slint a1_start, slint a1_end, slint *a2_start, slint *a2_end, slint *a3_start, slint *a3_end)
{
  slint temp_a2_start, temp_a2_end, temp_a3_start, temp_a3_end;

  if (a2_start == NULL) a2_start = &temp_a2_start;
  if (a2_end == NULL) a2_end = &temp_a2_end;
  if (a3_start == NULL) a3_start = &temp_a3_start;
  if (a3_end == NULL) a3_end = &temp_a3_end;

  *a2_start = a0_start;
  *a2_end = xmin(a1_start - 1, a0_end);

  *a3_start = xmax(a1_end + 1, a0_start);
  *a3_end = a0_end;

  return xmax(0, *a2_end - *a2_start + 1) + xmax(0, *a3_end - *a3_start + 1);
}*/
