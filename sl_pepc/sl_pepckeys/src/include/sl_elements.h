/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_elements.h
 *  timestamp: 2009-11-19 18:15:31 +0100
 *  
 */

/* SL_DATA_IGNORE -> prevents that the present data is processed */


#ifndef __SL_ELEMENTS_H__
#define __SL_ELEMENTS_H__


/* prevent the following functions from processing the DATAs (even though they may exist) */
#ifdef SL_DATA_IGNORE
 /* disable the single DATAs */
/* DATAX_TEMPLATE_BEGIN */
 #undef SL_DATA0
 #undef SL_DATA1
 #undef SL_DATA2
 #undef SL_DATA3
 #undef SL_DATA4
 #undef SL_DATA5
 #undef SL_DATA6
 #undef SL_DATA7
 #undef SL_DATA8
 #undef SL_DATA9
 #undef SL_DATA10
 #undef SL_DATA11
 #undef SL_DATA12
 #undef SL_DATA13
 #undef SL_DATA14
 #undef SL_DATA15
 #undef SL_DATA16
 #undef SL_DATA17
 #undef SL_DATA18
 #undef SL_DATA19
/* DATAX_TEMPLATE_END */
#endif /* SL_DATA_IGNORE */

#undef SL_DATA

#define ARRAY1_COPY(_s_, _d_)  ((_d_)[0] = (_s_)[0])
#define ARRAY2_COPY(_s_, _d_)  ARRAY1_COPY(_s_, _d_), ((_d_)[1] = (_s_)[1])
#define ARRAY3_COPY(_s_, _d_)  ARRAY2_COPY(_s_, _d_), ((_d_)[2] = (_s_)[2])
#define ARRAY4_COPY(_s_, _d_)  ARRAY3_COPY(_s_, _d_), ((_d_)[3] = (_s_)[3])
#define ARRAY5_COPY(_s_, _d_)  ARRAY4_COPY(_s_, _d_), ((_d_)[4] = (_s_)[4])
#define ARRAY6_COPY(_s_, _d_)  ARRAY5_COPY(_s_, _d_), ((_d_)[5] = (_s_)[5])
#define ARRAY7_COPY(_s_, _d_)  ARRAY6_COPY(_s_, _d_), ((_d_)[6] = (_s_)[6])
#define ARRAY8_COPY(_s_, _d_)  ARRAY7_COPY(_s_, _d_), ((_d_)[7] = (_s_)[7])
#define ARRAY9_COPY(_s_, _d_)  ARRAY8_COPY(_s_, _d_), ((_d_)[8] = (_s_)[8])

#include "sl_key.h"

#include "sl_index.h"

#include "sl_data_singles.h"

#include "sl_data.h"


#define elem_n                                          (key_n + index_n + data_n)
#define elem_byte                                       (key_byte + index_byte + data_byte)

#define elem_assign(_s_, _d_)                           (*(_d_) = *(_s_))
#define elem_assign_at(_s_, _sat_, _d_)                 ((_d_)->size = (_s_)->size - (_sat_), (_d_)->max_size = (_s_)->max_size - (_sat_), \
                                                         key_assign_at((_s_)->keys, _sat_, (_d_)->keys) \
                                                         cc_index_assign_at((_s_)->indices, _sat_, (_d_)->indices) \
                                                         cc_data_assign_at(_s_, _sat_, _d_))
#define elem_null(_e_)                                  ((_e_)->size = 0, (_e_)->max_size = 0, \
                                                         key_null((_e_)->keys) \
                                                         cc_index_null((_e_)->indices) \
                                                         cc_data_null(_e_))
#define elem_inc(_e_)                                   (key_inc((_e_)->keys) \
                                                         cc_index_inc((_e_)->indices) \
                                                         cc_data_inc(_e_))
#define elem_dec(_e_)                                   (key_dec((_e_)->keys) \
                                                         cc_index_dec((_e_)->indices) \
                                                         cc_data_dec(_e_))
#define elem_add(_e_, _n_)                              (key_add((_e_)->keys, _n_) \
                                                         cc_index_add((_e_)->indices, _n_) \
                                                         cc_data_add(_e_, _n_))
#define elem_sub(_e_, _n_)                              (key_sub((_e_)->keys, _n_) \
                                                         cc_index_sub((_e_)->indices, _n_) \
                                                         cc_data_sub(_e_, _n_))

#define elem_copy(_s_, _d_)                             (key_copy((_s_)->keys, (_d_)->keys) \
                                                         cc_index_copy((_s_)->indices, (_d_)->indices) \
                                                         cc_data_copy(_s_, _d_))
#define elem_ncopy(_s_, _d_, _n_)                       (key_ncopy((_s_)->keys, (_d_)->keys, _n_) \
                                                         cc_index_ncopy((_s_)->indices, (_d_)->indices, _n_) \
                                                         cc_data_ncopy(_s_, _d_, _n_))
#define elem_nmove(_s_, _d_, _n_)                       (key_nmove((_s_)->keys, (_d_)->keys, _n_) \
                                                         cc_index_nmove((_s_)->indices, (_d_)->indices, _n_) \
                                                         cc_data_nmove(_s_, _d_, _n_))

#define elem_copy_at(_s_, _sat_, _d_, _dat_)            (key_copy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_) \
                                                         cc_index_copy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_) \
                                                         cc_data_copy_at(_s_, _sat_, _d_, _dat_))
#define elem_ncopy_at(_s_, _sat_, _d_, _dat_, _n_)      (key_ncopy_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_) \
                                                         cc_index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_) \
                                                         cc_data_ncopy_at(_s_, _sat_, _d_, _dat_, _n_))
#define elem_nmove_at(_s_, _sat_, _d_, _dat_, _n_)      (key_nmove_at((_s_)->keys, _sat_, (_d_)->keys, _dat_, _n_) \
                                                         cc_index_ncopy_at((_s_)->indices, _sat_, (_d_)->indices, _dat_, _n_) \
                                                         cc_data_nmove_at(_s_, _sat_, _d_, _dat_, _n_))

#define elem_xchange(_e0_, _e1_, _t_)                   (key_xchange((_e0_)->keys, (_e1_)->keys, (_t_)->keys) \
                                                         cc_index_xchange((_e0_)->indices, (_e1_)->indices, (_t_)->indices) \
                                                         cc_data_xchange(_e0_, _e1_, _t_))
#define elem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  (key_xchange_at((_e0_)->keys, _at0_, (_e1_)->keys, _at1_, (_t_)->keys) \
                                                         cc_index_xchange_at((_e0_)->indices, _at0_, (_e1_)->indices, _at1_, (_t_)->indices) \
                                                         cc_data_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_))

/*#define elem_xchange(_e0_, _e1_, _t_)                 (elem_copy(_e0_, _t_), elem_copy(_e1_, _e0_), elem_copy(_t_, _e1_))
#define elem_xchange_at(_e0_, _at0_, _e1_, _at1_, _t_)  (elem_copy_at(_e0_, _at0_, _t_, 0), elem_copy_at(_e1_, _at1_, _e0_, _at0_), elem_copy_at(_t_, 0, _e1_, _at1_))*/

#define elem_weight(_e_, _at_)                          sl_elem_weight((_e_), (_at_))


#endif /* __SL_ELEMENTS_H__ */
