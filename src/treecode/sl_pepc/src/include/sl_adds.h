/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_adds.h
 *  
 */


#ifndef __SL_ADDS_H__
#define __SL_ADDS_H__


/* sl_macro elem_set_size elem_set_max_size elem_set_keys elem_set_indices */
#define elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))

/* sl_macro pelem_set_size pelem_set_max_size pelem_set_elements */
#define pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))


#endif /* __SL_ADDS_H__ */
