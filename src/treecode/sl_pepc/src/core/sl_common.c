/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sl_common.c
 *  
 */


#include "sl_common.h"


rti rti_env;  /* sl_global, sl_var rti_env */

int sl_mpi_rank_dummy = -2;  /* sl_global, sl_var sl_mpi_rank_dummy */


slint ilog2c(slint x) /* sl_proto, sl_func ilog2c */
{
  slint l = 0;
  while (x /= 2) ++l;
  return l;
}


slint ilog2f(slint x) /* sl_proto, sl_func ilog2f */
{
  slint l = 0;
  if (x <= 1) return 0;
  --x;
  while (x /= 2) ++l;
  return l + 1;
}


slint print_bits(slint v) /* sl_proto, sl_func print_bits */
{
  slint i;

  for (i = sizeof(v) * 8 - 1; i >= 0; i--) printf("%d", (v & (1L << i)) != 0);

  return 0;
}


slint pivot_random(elements_t *s) /* sl_proto, sl_func pivot_random */
{
  return z_rand() % s->size;
}


slint_t counts2displs(slint_t n, int *counts, int *displs) /* sl_proto, sl_func counts2displs */
{
  slint_t i;

  displs[0] = 0;
  for (i = 1; i < n; ++i) displs[i] = displs[i - 1] + counts[i - 1];
  
  return 0;
}


slint_t displs2counts(slint_t n, int *displs, int *counts, slint_t total_counts) /* sl_proto, sl_func displs2counts */
{
  slint_t i;

  if (counts)
  {
    for (i = 0; i < n - 1; ++i) counts[i] = displs[i + 1] - displs[i];
    counts[n - 1] = total_counts - displs[n - 1];

  } else
  {
    for (i = 0; i < n - 1; ++i) displs[i] = displs[i + 1] - displs[i];
    displs[n - 1] = total_counts - displs[n - 1];
  }
  
  return 0;
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

  *ia_start = z_max(a0_start, a1_start);
  *ia_end = z_min(a0_end, a1_end);

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
  *a2_end = z_min(a1_start - 1, a0_end);

  *a3_start = z_max(a1_end + 1, a0_start);
  *a3_end = a0_end;

  return z_max(0, *a2_end - *a2_start + 1) + z_max(0, *a3_end - *a3_start + 1);
}*/
