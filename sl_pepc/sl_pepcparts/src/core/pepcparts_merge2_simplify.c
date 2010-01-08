/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_simplify.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */


#include "sl_common.h"


#define the_merge2_01(s0, s1, sx)  pepcparts_merge2_basic_sseq_01(s0, s1, sx)


/*slint merge2_simplify_s0(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements)
{
  return 0;
}*/


slint pepcparts_merge2_simplify_s1(elements_t *s0, elements_t *s1, elements_t *sx, slint s1elements) /* pepcparts_sl_proto, sl_func pepcparts_merge2_simplify_s1 */
{
  slint m, m0, m1;

  elements_t _s0, _s1, x;


  s1elements = xmin(s1elements, s1->size);

/*  printf("simplifying %d elements from s1\n", s1elements);*/

  if (s1elements == 0) return 0;

  /** find the s1elements highest elements of s0 and s1 **/

  m = s1elements;
  m0 = m1 = 0;

  elem_assign_at(s0, s0->size - 1, &_s0);
  elem_assign_at(s1, s1->size - 1, &_s1);

  while (m-- > 0 && _s0.keys >= s0->keys)
  if (key_cmp_ge(*_s0.keys, *_s1.keys))
  {
    elem_dec(&_s0);
    ++m0;
  } else
  {
    elem_dec(&_s1);
  }

  elem_inc(&_s0);

  m1 = s1elements - m0;
  elem_assign_at(s1, s1->size - m1, &_s1);

/*  printf("highest %d elements, %d from s0, %d from s1\n", s1elements, m0, m1);*/

  /** bring the s1elements highest elements to the end of s1 **/

  elem_assign_at(s1, s1->size - s1elements, &x);

  pepcparts_elem_nxchange(&_s0, &x, m0, sx);

  /* merge the highest elements */
  x.size = m0;
  _s1.size = m1;

/*  printf("merging x(%d) @ %p & _s1(%d) @ %p\n", x.size, x.keys, _s1.size, _s1.keys);*/

  the_merge2_01(&x, &_s1, sx);

/*  pepcparts_elements_print_keys(s1);*/

  /* merge the highest elements of s0 */
  elem_assign(s0, &x); x.size -= m0;

  _s0.size = m0;

/*  printf("merging x(%d) @ %p & _s0(%d) @ %p\n", x.size, x.keys, _s0.size, _s0.keys);*/

  the_merge2_01(&x, &_s0, sx);

/*  pepcparts_elements_print_keys(s0);*/

  s1->size -= s1elements;

  return 0;
}
