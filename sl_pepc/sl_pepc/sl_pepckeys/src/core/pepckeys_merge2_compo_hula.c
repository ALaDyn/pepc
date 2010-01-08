/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/pepckeys_merge2_compo_hula.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/*
   "Practical In-Place Merging"
   Bing-Chao Huang, Michael A. Langston
   Communications of the ACM, March 1988, Volume 31, Number 3

   Modified for practicability!
*/


#include "sl_common.h"


#define sort_again(s, sx)    pepckeys_sort_radix(s, sx, -1, -1, -1);

#define the_merge2_01(s0, s1, sx)       pepckeys_merge2_basic_sbin_01(s0, s1, sx)
#define the_merge2_0X_1(s0, s1, sx, t)  pepckeys_merge2_basic_straight_0X_1(s0, s1, sx, t)
#define the_merge2_X0_1(s0, s1, sx, t)  pepckeys_merge2_basic_straight_X0_1(s0, s1, sx, t)


slint pepckeys_merge2_compo_hula_(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_func pepckeys_merge2_compo_hula_ */
{
  slint i, j, k;

  slint n0, n1, t0, t1, e0, e1;

  elements_t wa, _src0, _src1, _sx, undersized;

  slint s, m, m0, m1, s0mod, u;

  elements_t _s0, _s1;

  elements_t sa, sb, sc, sd, sf, sg, x;

  elements_t first, second, last, current_head, current_tail, lowest_head, lowest_tail;

  s = (slint) floor(sqrt(s0->size + s1->size));

  if (s0->size < s || s1->size < s)
  {
    /* do an easier merge with bisection and rotates */
    the_merge2_01(s0, s1, sx);

/*    printf("easier merge2 done\n");*/

    return 0;
  }

  elem_assign(s0, &_s0);
  elem_assign(s1, &_s1);

  /* bring the size of s1 to a mulitple of s */
  pepckeys_merge2_simplify_s1(&_s0, &_s1, sx, s1->size % s);

/*  printf("_s0:\n"); pepckeys_elements_printf(&_s0); pepckeys_elements_print_keys(&_s0);
  printf("_s1:\n"); pepckeys_elements_printf(&_s1); pepckeys_elements_print_keys(&_s1);*/

  /* now: every list contains at least s elements and the size of s1 is a multiple of s */

  /** detect the s biggest elements **/
  elem_assign_at(&_s0, _s0.size - 1, &sa);
  elem_assign_at(&_s1, _s1.size - 1, &sb);

  m0 = m1 = 0;
  m = s;
  while (m-- > 0)
  if (key_cmp_ge(*sa.keys, *sb.keys))
  {
    elem_dec(&sa);
    ++m0;
  } else
  {
    elem_dec(&sb);
    ++m1;
  }

/*  printf("highest %d elements: %d from s0, %d from s1\n", s, m0, m1);*/

  /* now A & B are the workarea, C & D are the biggest elements */
  elem_inc(&sa); sa.size = m0;
  elem_inc(&sb); sb.size = m1;
  elem_assign_at(&_s0, _s0.size - s, &sc); sc.size = m1;
  elem_assign_at(&_s1, _s1.size - s, &sd); sd.size = m0;

/*  printf("sd:"); pepckeys_elements_printf(&sd);
  printf("sc:"); pepckeys_elements_printf(&sc);
  printf("sb:"); pepckeys_elements_printf(&sb);*/

  /* merge C and D forming block E at the end of s1 */
  the_merge2_0X_1(&sd, &sc, &sb, sx);

/*  printf("_s0:\n"); pepckeys_elements_printf(&_s0); pepckeys_elements_print_keys(&_s0);
  printf("_s1:\n"); pepckeys_elements_printf(&_s1); pepckeys_elements_print_keys(&_s1);*/

  /* after the merge, the workarea is in one piece at the end of s0 */
  elem_assign(&sc, &wa); wa.size = s;

/*  printf("wa:\n"); pepckeys_elements_printf(&wa); pepckeys_elements_print_keys(&wa);*/

  /* handle the undersized block at the beginning of s0 */
  s0mod = _s0.size % s;
  if (s0mod > 0)
  {
    elem_assign(&_s0, &sf); sf.size = s0mod; /* block F */
    elem_assign(&_s1, &sg); sg.size = s; /* block G */
    elem_assign_at(&wa, s - s0mod, &x); x.size = s0mod; /* ending part of the workarea */

/*    printf("sg:"); pepckeys_elements_printf(&sg);
    printf("sf:"); pepckeys_elements_printf(&sf);
    printf("x:"); pepckeys_elements_printf(&x);*/

    /* merge F & G forming H & I */
    the_merge2_X0_1(&sg, &sf, &x, sx);

/*    printf("_s0:\n"); pepckeys_elements_printf(&_s0); pepckeys_elements_print_keys(&_s0);
    printf("_s1:\n"); pepckeys_elements_printf(&_s1); pepckeys_elements_print_keys(&_s1);*/

    /* exchanging H (replaced F) and the part of the workarea at the front of s0 */
    pepckeys_elem_nxchange(&sf, &x, s0mod, sx);

/*    printf("_s0:\n"); pepckeys_elements_printf(&_s0); pepckeys_elements_print_keys(&_s0);
    printf("_s1:\n"); pepckeys_elements_printf(&_s1); pepckeys_elements_print_keys(&_s1);*/

    /* skip the finished undersized block at the front of s0 */
    elem_add(&_s0, s0mod); _s0.size -= s0mod;
  }

  /* bring the workarea to the front (behind an already finished undersized block) */
  pepckeys_elem_nxchange_ro0(&_s0, &wa, s, sx);
  elem_assign(&_s0, &wa); wa.size = s;

  /** now: all undersized blocks are handled, the sizes of _s0 and _s1 are multiples of s **/

/*  printf("merging elements of size %d @ %p and %d @ %p\n", _s0.size, _s0.keys, _s1.size, _s1.keys);
  printf("_s0:\n"); pepckeys_elements_printf(&_s0); pepckeys_elements_print_keys(&_s0);
  printf("_s1:\n"); pepckeys_elements_printf(&_s1); pepckeys_elements_print_keys(&_s1);*/


  /** sort the blocks **/

  /* the first block starts behind the workarea */
  elem_assign_at(&_s0, s, &first);
  elem_assign_at(&_s1, _s1.size - s, &last);

  /* while more the one block remains */
  while (first.keys < last.keys)
  {
    /* take the last block as lowest */
    elem_assign(&last, &lowest_head);
    elem_assign_at(&lowest_head, s - 1, &lowest_tail);

    /* start search at the first block */
    elem_assign(&first, &current_head);
    elem_assign_at(&current_head, s - 1, &current_tail);

/*    printf("starting, first @ %p\n", first.keys);
    printf("starting, lowest @ %p\n", lowest_head.keys);
    printf("starting, current @ [%p,%p]\n", current_head.keys, current_tail.keys);
    printf("starting, last @ %p\n", last.keys);*/

    while (current_tail.keys < last.keys)
    {
/*      printf("lowest: [%d,%d]\n", *lowest_head.keys, *lowest_tail.keys);
      printf("current: [%d,%d]\n", *current_head.keys, *current_tail.keys);*/

      if (key_cmp_lt(*current_tail.keys, *lowest_tail.keys) || (key_cmp_eq(*current_tail.keys, *lowest_tail.keys) && key_cmp_lt(*current_head.keys, *lowest_head.keys)))
      {
        elem_assign(&current_head, &lowest_head);
        elem_assign(&current_tail, &lowest_tail);

/*        printf("taken!\n");*/
      }

      elem_add(&current_head, s);
      elem_add(&current_tail, s);
    }

/*    printf("first %d @ %p\n", *first.keys, first.keys);
    printf("lowest %d @ %p\n", *lowest_head.keys, lowest_head.keys);*/

    /* bring the lowest block to the front (if not already there) */
    if (first.keys != lowest_head.keys) pepckeys_elem_nxchange(&first, &lowest_head, s, sx);

/*    printf("exchange done\n");*/

    /* continue sort on the remaining blocks */
    elem_add(&first, s);
  }

/*  printf("after blocksort\n");
  printf("_s0:\n"); pepckeys_elements_printf(&_s0); pepckeys_elements_print_keys(&_s0);
  printf("_s1:\n"); pepckeys_elements_printf(&_s1); pepckeys_elements_print_keys(&_s1);
  printf("wa:\n"); pepckeys_elements_printf(&wa); pepckeys_elements_print_keys(&wa);*/


  /** merge the blocks **/

  elem_assign_at(&_s0, s, &first); first.size = s;
  elem_assign_at(&_s1, _s1.size, &last);


  while (first.keys < last.keys)
  {
    elem_assign(&first, &lowest_head);
    elem_assign_at(&first, first.size - 1, &lowest_tail);

    elem_assign_at(&first, first.size, &second); second.size = s;

    while (second.keys < last.keys)
    if (key_cmp_le(*lowest_tail.keys, *second.keys))
    {
      first.size += s;
      elem_add(&lowest_tail, s);
      elem_add(&second, s);

    } else break;

    /* there is no more second block to merge with */
    if (second.keys >= last.keys) break;

/*    printf("merging blocks:\n");
    printf("first:\n"); pepckeys_elements_printf(&first); pepckeys_elements_print_keys(&first);
    printf("second:\n"); pepckeys_elements_printf(&second); pepckeys_elements_print_keys(&second);
    printf("wa:\n"); pepckeys_elements_printf(&wa); pepckeys_elements_print_keys(&wa);*/

    /* blockmerge */
    u = pepckeys_merge2_basic_straight_X0_1u(&first, &second, &wa, sx);

    elem_assign_at(&second, s - u, &first); first.size = u;
    elem_assign_at(&first, -s, &wa); wa.size = s;

/*    printf("after merging blocks:\n");
    printf("first:\n"); pepckeys_elements_printf(&first); pepckeys_elements_print_keys(&first);
    printf("wa:\n"); pepckeys_elements_printf(&wa); pepckeys_elements_print_keys(&wa);*/
  }

  /* bring the workarea to the end */
  pepckeys_elem_nxchange_ro0(&first, &wa, first.size, sx);
  elem_add(&wa, first.size);

/*  printf("starting final sort:\n");
  printf("wa:\n"); pepckeys_elements_printf(&wa); pepckeys_elements_print_keys(&wa);*/

  /* sort the workarea */
  sort_again(&wa, sx);

  return 0;
}


slint pepckeys_merge2_compo_hula(elements_t *s0, elements_t *s1, elements_t *xs) /* pepckeys_sl_proto, sl_func pepckeys_merge2_compo_hula */
{
  elements_t txs;

  if ((s0 == NULL) || (s1 == NULL)) return -1;

  /* if one location is empty, we are finished */
  if ((s0->size == 0) || (s1->size == 0)) return 0;

  /* calloc mode? */
  if (xs == NULL)
  {
    xs = &txs;

    /* this is really "in-place", need exact one additional element of eXtraspace */
    pepckeys_elements_alloc(xs, 1, 1, 1);
  }

  pepckeys_merge2_compo_hula_(s0, s1, xs);

  /* was in calloc mode? */
  if (xs == &txs) pepckeys_elements_free(xs);

  return 0;
}
