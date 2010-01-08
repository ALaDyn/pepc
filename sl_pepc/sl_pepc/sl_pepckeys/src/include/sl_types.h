/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_types.h
 *  timestamp: 2009-12-03 09:12:31 +0100
 *  
 */


#ifndef __SL_TYPES_H__
#define __SL_TYPES_H__


/* sl_type slint_t slint */
typedef sl_int_type_c slint_t, slint;  /* deprecated 'slint' */

/* sl_type slindex_t */
typedef sl_index_type_c slindex_t;

/* sl_type slkey_t */
typedef sl_key_type_c slkey_t;

/* sl_type slkey_pure_t */
typedef sl_key_pure_type_c slkey_pure_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type sldata0_t */
#ifdef sl_data0_type_c
typedef sl_data0_type_c sldata0_t;
#endif
/* sl_type sldata1_t */
#ifdef sl_data1_type_c
typedef sl_data1_type_c sldata1_t;
#endif
/* sl_type sldata2_t */
#ifdef sl_data2_type_c
typedef sl_data2_type_c sldata2_t;
#endif
/* sl_type sldata3_t */
#ifdef sl_data3_type_c
typedef sl_data3_type_c sldata3_t;
#endif
/* sl_type sldata4_t */
#ifdef sl_data4_type_c
typedef sl_data4_type_c sldata4_t;
#endif
/* sl_type sldata5_t */
#ifdef sl_data5_type_c
typedef sl_data5_type_c sldata5_t;
#endif
/* sl_type sldata6_t */
#ifdef sl_data6_type_c
typedef sl_data6_type_c sldata6_t;
#endif
/* sl_type sldata7_t */
#ifdef sl_data7_type_c
typedef sl_data7_type_c sldata7_t;
#endif
/* sl_type sldata8_t */
#ifdef sl_data8_type_c
typedef sl_data8_type_c sldata8_t;
#endif
/* sl_type sldata9_t */
#ifdef sl_data9_type_c
typedef sl_data9_type_c sldata9_t;
#endif
/* sl_type sldata10_t */
#ifdef sl_data10_type_c
typedef sl_data10_type_c sldata10_t;
#endif
/* sl_type sldata11_t */
#ifdef sl_data11_type_c
typedef sl_data11_type_c sldata11_t;
#endif
/* sl_type sldata12_t */
#ifdef sl_data12_type_c
typedef sl_data12_type_c sldata12_t;
#endif
/* sl_type sldata13_t */
#ifdef sl_data13_type_c
typedef sl_data13_type_c sldata13_t;
#endif
/* sl_type sldata14_t */
#ifdef sl_data14_type_c
typedef sl_data14_type_c sldata14_t;
#endif
/* sl_type sldata15_t */
#ifdef sl_data15_type_c
typedef sl_data15_type_c sldata15_t;
#endif
/* sl_type sldata16_t */
#ifdef sl_data16_type_c
typedef sl_data16_type_c sldata16_t;
#endif
/* sl_type sldata17_t */
#ifdef sl_data17_type_c
typedef sl_data17_type_c sldata17_t;
#endif
/* sl_type sldata18_t */
#ifdef sl_data18_type_c
typedef sl_data18_type_c sldata18_t;
#endif
/* sl_type sldata19_t */
#ifdef sl_data19_type_c
typedef sl_data19_type_c sldata19_t;
#endif
/* DATAX_TEMPLATE_END */

/* sl_type _elements_t elements_t */
typedef struct _elements_t
{
  slint_t size, max_size;
  slkey_t *keys;

#ifdef SL_INDEX
  slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef SL_DATA0
  sldata0_t *data0;
#endif
#ifdef SL_DATA1
  sldata1_t *data1;
#endif
#ifdef SL_DATA2
  sldata2_t *data2;
#endif
#ifdef SL_DATA3
  sldata3_t *data3;
#endif
#ifdef SL_DATA4
  sldata4_t *data4;
#endif
#ifdef SL_DATA5
  sldata5_t *data5;
#endif
#ifdef SL_DATA6
  sldata6_t *data6;
#endif
#ifdef SL_DATA7
  sldata7_t *data7;
#endif
#ifdef SL_DATA8
  sldata8_t *data8;
#endif
#ifdef SL_DATA9
  sldata9_t *data9;
#endif
#ifdef SL_DATA10
  sldata10_t *data10;
#endif
#ifdef SL_DATA11
  sldata11_t *data11;
#endif
#ifdef SL_DATA12
  sldata12_t *data12;
#endif
#ifdef SL_DATA13
  sldata13_t *data13;
#endif
#ifdef SL_DATA14
  sldata14_t *data14;
#endif
#ifdef SL_DATA15
  sldata15_t *data15;
#endif
#ifdef SL_DATA16
  sldata16_t *data16;
#endif
#ifdef SL_DATA17
  sldata17_t *data17;
#endif
#ifdef SL_DATA18
  sldata18_t *data18;
#endif
#ifdef SL_DATA19
  sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} elements_t;


/* sl_type _packed_element_t packed_element_t */
typedef struct _packed_element_t
{
  slkey_t key;

#ifdef SL_PACKED_INDEX
  slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef SL_DATA0
  sldata0_t data0[sl_data0_size_c];
#endif
#ifdef SL_DATA1
  sldata1_t data1[sl_data1_size_c];
#endif
#ifdef SL_DATA2
  sldata2_t data2[sl_data2_size_c];
#endif
#ifdef SL_DATA3
  sldata3_t data3[sl_data3_size_c];
#endif
#ifdef SL_DATA4
  sldata4_t data4[sl_data4_size_c];
#endif
#ifdef SL_DATA5
  sldata5_t data5[sl_data5_size_c];
#endif
#ifdef SL_DATA6
  sldata6_t data6[sl_data6_size_c];
#endif
#ifdef SL_DATA7
  sldata7_t data7[sl_data7_size_c];
#endif
#ifdef SL_DATA8
  sldata8_t data8[sl_data8_size_c];
#endif
#ifdef SL_DATA9
  sldata9_t data9[sl_data9_size_c];
#endif
#ifdef SL_DATA10
  sldata10_t data10[sl_data10_size_c];
#endif
#ifdef SL_DATA11
  sldata11_t data11[sl_data11_size_c];
#endif
#ifdef SL_DATA12
  sldata12_t data12[sl_data12_size_c];
#endif
#ifdef SL_DATA13
  sldata13_t data13[sl_data13_size_c];
#endif
#ifdef SL_DATA14
  sldata14_t data14[sl_data14_size_c];
#endif
#ifdef SL_DATA15
  sldata15_t data15[sl_data15_size_c];
#endif
#ifdef SL_DATA16
  sldata16_t data16[sl_data16_size_c];
#endif
#ifdef SL_DATA17
  sldata17_t data17[sl_data17_size_c];
#endif
#ifdef SL_DATA18
  sldata18_t data18[sl_data18_size_c];
#endif
#ifdef SL_DATA19
  sldata19_t data19[sl_data19_size_c];
#endif
/* DATAX_TEMPLATE_END */

} packed_element_t;


/* sl_type _packed_elements_t packed_elements_t */
typedef struct _packed_elements_t
{
  slint_t size, max_size;
  
  packed_element_t *elements;
  
} packed_elements_t;


/* sl_type _classification_info_t classification_info_t classification_info */
typedef struct _classification_info_t
{
  slint_t nclasses;
  slkey_pure_t *keys;
  slint_t *counts;
  slint_t *masks;

  /* */
  slint_t *all_local_sizes;
  slint_t *local_lt_eq_counts;
  slint_t *all_local_lt_eq_counts;

} classification_info_t, classification_info;  /* deprecated 'classification_info' */


/* key2class, sl_type key2class_f */
typedef slint_t (*key2class_f)(slkey_t *, slint, void *);

/* pivot-element, sl_type pivot_f */
typedef slint_t (*pivot_f)(elements_t *);

/* sorting-network, sl_type sortnet_f sortnet_data_t */
typedef void *sortnet_data_t;
typedef slint_t (*sortnet_f)(slint_t size, slint_t rank, slint_t stage, sortnet_data_t snd, slint_t *up);

/* merge2, sl_type merge2x_f merge2X_f */
typedef slint_t (*merge2x_f)(elements_t *s0, elements_t *s1, elements_t *sx);
typedef slint_t (*merge2X_f)(elements_t *s0, elements_t *s1, elements_t *sx, elements_t *t);


/* deprecated, sl_type k2c_func pivot_func sn_func m2x_func m2X_func */
typedef key2class_f k2c_func;
typedef pivot_f pivot_func;
typedef sortnet_f sn_func;
typedef merge2x_f m2x_func;
typedef merge2X_f m2X_func;


/* partition conditions, sl_type _partcond_t partcond_t */
typedef struct _partcond_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} partcond_t;

#endif /* __SL_TYPES_H__ */
