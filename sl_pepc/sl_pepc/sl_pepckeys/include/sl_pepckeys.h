
#ifndef __SL_PEPCKEYS_H__
#define __SL_PEPCKEYS_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#include "sl_rti.h"


#include "fortran2c_types.h"


/* enable MPI */
#define SL_USE_MPI

/* enable runtime_informations */
#define SL_USE_RTI
#define SL_USE_RTI_TIM


/* standard (SL) integer data type */
#define pepckeys_sl_int_type_c          long
#define pepckeys_sl_int_type_mpi        MPI_LONG
#define pepckeys_sl_int_size_mpi        1
#define pepckeys_sl_int_type_fmt        "ld"


/* index data type */
#define pepckeys_sl_index_type_c        FINT_TYPE_C
#define pepckeys_sl_index_type_mpi      FINT_TYPE_MPI
#define pepckeys_sl_index_size_mpi      1
#define pepckeys_sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define pepckeys_SL_INDEX


/* keys */
#define pepckeys_sl_key_type_c          FINT8_TYPE_C
#define pepckeys_sl_key_type_mpi        FINT8_TYPE_MPI
#define pepckeys_sl_key_size_mpi        1
#define pepckeys_sl_key_type_fmt        FINT8_TYPE_FMT
#define pepckeys_sl_key_integer

/* data0: work loads */
#define pepckeys_SL_DATA0
#define pepckeys_sl_data0_type_c        FREAL8_TYPE_C
#define pepckeys_sl_data0_size_c        1
#define pepckeys_sl_data0_type_mpi      FREAL8_TYPE_MPI
#define pepckeys_sl_data0_size_mpi      1

#define pepckeys_sl_elem_weight(e, at)  ((e)->data0[at])
/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_config_intern.h
 *  timestamp: 2009-12-03 09:11:41 +0100
 *  
 */




/* override SL_USE_MPI from sl_config.h */
#ifdef SL_USE_MPI_IGNORE
# undef SL_USE_MPI
#endif

#ifdef SL_USE_MPI_FORCE
# ifndef SL_USE_MPI
#  define SL_USE_MPI
# endif
#endif


/* override inlining */
#ifdef NO_INLINE
# ifndef inline
#  define inline
# endif
#endif


#ifndef pepckeys_SL_INDEX
# undef SL_PACKED_INDEX
#endif


/* if no special, given, primary and heavy used integer-type ... */
#ifndef pepckeys_sl_int_type_c
  /* ... use a default one */
# define pepckeys_sl_int_type_c               long      /* sl_macro */
# undef pepckeys_sl_int_type_mpi
# define pepckeys_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef pepckeys_sl_int_size_mpi
# define pepckeys_sl_int_size_mpi             1         /* sl_macro */
# undef pepckeys_sl_int_type_fmt
# define pepckeys_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepckeys_sl_int_type_mpi) || !defined(pepckeys_sl_int_size_mpi)
#   error "pepckeys_sl_int_type_mpi and/or pepckeys_sl_int_size_mpi missing"
#  endif
# endif
# ifndef pepckeys_sl_int_type_fmt
#  error "pepckeys_sl_int_type_fmt macro is missing, using d as default"
#  define pepckeys_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef pepckeys_sl_index_type_c
 /* ... use the primary integer type */
# define pepckeys_sl_index_type_c             pepckeys_sl_int_type_c    /* sl_macro */
# undef pepckeys_sl_index_type_mpi
# define pepckeys_sl_index_type_mpi           pepckeys_sl_int_type_mpi  /* sl_macro */
# undef pepckeys_sl_index_size_mpi
# define pepckeys_sl_index_size_mpi           pepckeys_sl_int_size_mpi  /* sl_macro */
# undef pepckeys_sl_index_type_fmt
# define pepckeys_sl_index_type_fmt           pepckeys_sl_int_type_fmt  /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepckeys_sl_index_type_mpi) || !defined(pepckeys_sl_index_size_mpi)
#   error "pepckeys_sl_index_type_mpi and/or pepckeys_sl_index_size_mpi missing"
#  endif
# endif
# ifndef pepckeys_sl_index_type_fmt
#  error "pepckeys_sl_index_type_fmt macro is missing, using d as default"
#  define pepckeys_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef pepckeys_sl_key_pure_type_c
# define pepckeys_sl_key_pure_type_c          pepckeys_sl_key_type_c  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_type_mpi
# define pepckeys_sl_key_pure_type_mpi        pepckeys_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_size_mpi
# define pepckeys_sl_key_pure_size_mpi        pepckeys_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_type_fmt
# ifdef pepckeys_sl_key_type_fmt
#  define pepckeys_sl_key_pure_type_fmt       pepckeys_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef pepckeys_sl_key_purify
 #define pepckeys_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_get_pure
 #define pepckeys_sl_key_get_pure(k)          (pepckeys_sl_key_purify(k))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_set_pure
 #define pepckeys_sl_key_set_pure(k, p)       (pepckeys_sl_key_purify(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef pepckeys_sl_key_pure_cmp_eq
 #define pepckeys_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_ne
 #define pepckeys_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_lt
 #define pepckeys_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_le
 #define pepckeys_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_gt
 #define pepckeys_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_pure_cmp_ge
 #define pepckeys_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef pepckeys_sl_key_cmp_eq
 #define pepckeys_sl_key_cmp_eq(k0, k1)       (pepckeys_sl_key_pure_cmp_eq(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_ne
 #define pepckeys_sl_key_cmp_ne(k0, k1)       (pepckeys_sl_key_pure_cmp_ne(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_lt
 #define pepckeys_sl_key_cmp_lt(k0, k1)       (pepckeys_sl_key_pure_cmp_lt(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_le
 #define pepckeys_sl_key_cmp_le(k0, k1)       (pepckeys_sl_key_pure_cmp_le(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_gt
 #define pepckeys_sl_key_cmp_gt(k0, k1)       (pepckeys_sl_key_pure_cmp_gt(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_cmp_ge
 #define pepckeys_sl_key_cmp_ge(k0, k1)       (pepckeys_sl_key_pure_cmp_ge(pepckeys_sl_key_purify(k0), pepckeys_sl_key_purify(k1)))  /* sl_macro */
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
/* sl_macro pepckeys_SL_DATA0_IGNORE */
#ifdef pepckeys_SL_DATA0_IGNORE
# undef pepckeys_SL_DATA0
#endif
/* sl_macro pepckeys_SL_DATA1_IGNORE */
#ifdef pepckeys_SL_DATA1_IGNORE
# undef pepckeys_SL_DATA1
#endif
/* sl_macro pepckeys_SL_DATA2_IGNORE */
#ifdef pepckeys_SL_DATA2_IGNORE
# undef pepckeys_SL_DATA2
#endif
/* sl_macro pepckeys_SL_DATA3_IGNORE */
#ifdef pepckeys_SL_DATA3_IGNORE
# undef pepckeys_SL_DATA3
#endif
/* sl_macro pepckeys_SL_DATA4_IGNORE */
#ifdef pepckeys_SL_DATA4_IGNORE
# undef pepckeys_SL_DATA4
#endif
/* sl_macro pepckeys_SL_DATA5_IGNORE */
#ifdef pepckeys_SL_DATA5_IGNORE
# undef pepckeys_SL_DATA5
#endif
/* sl_macro pepckeys_SL_DATA6_IGNORE */
#ifdef pepckeys_SL_DATA6_IGNORE
# undef pepckeys_SL_DATA6
#endif
/* sl_macro pepckeys_SL_DATA7_IGNORE */
#ifdef pepckeys_SL_DATA7_IGNORE
# undef pepckeys_SL_DATA7
#endif
/* sl_macro pepckeys_SL_DATA8_IGNORE */
#ifdef pepckeys_SL_DATA8_IGNORE
# undef pepckeys_SL_DATA8
#endif
/* sl_macro pepckeys_SL_DATA9_IGNORE */
#ifdef pepckeys_SL_DATA9_IGNORE
# undef pepckeys_SL_DATA9
#endif
/* sl_macro pepckeys_SL_DATA10_IGNORE */
#ifdef pepckeys_SL_DATA10_IGNORE
# undef pepckeys_SL_DATA10
#endif
/* sl_macro pepckeys_SL_DATA11_IGNORE */
#ifdef pepckeys_SL_DATA11_IGNORE
# undef pepckeys_SL_DATA11
#endif
/* sl_macro pepckeys_SL_DATA12_IGNORE */
#ifdef pepckeys_SL_DATA12_IGNORE
# undef pepckeys_SL_DATA12
#endif
/* sl_macro pepckeys_SL_DATA13_IGNORE */
#ifdef pepckeys_SL_DATA13_IGNORE
# undef pepckeys_SL_DATA13
#endif
/* sl_macro pepckeys_SL_DATA14_IGNORE */
#ifdef pepckeys_SL_DATA14_IGNORE
# undef pepckeys_SL_DATA14
#endif
/* sl_macro pepckeys_SL_DATA15_IGNORE */
#ifdef pepckeys_SL_DATA15_IGNORE
# undef pepckeys_SL_DATA15
#endif
/* sl_macro pepckeys_SL_DATA16_IGNORE */
#ifdef pepckeys_SL_DATA16_IGNORE
# undef pepckeys_SL_DATA16
#endif
/* sl_macro pepckeys_SL_DATA17_IGNORE */
#ifdef pepckeys_SL_DATA17_IGNORE
# undef pepckeys_SL_DATA17
#endif
/* sl_macro pepckeys_SL_DATA18_IGNORE */
#ifdef pepckeys_SL_DATA18_IGNORE
# undef pepckeys_SL_DATA18
#endif
/* sl_macro pepckeys_SL_DATA19_IGNORE */
#ifdef pepckeys_SL_DATA19_IGNORE
# undef pepckeys_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* default element weights */
#ifndef pepckeys_sl_elem_weight
# define pepckeys_sl_elem_weight(e, at)       1  /* sl_macro */
#endif


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_types.h
 *  timestamp: 2009-12-03 09:12:31 +0100
 *  
 */




/* sl_type pepckeys_slint_t pepckeys_slint */
typedef pepckeys_sl_int_type_c pepckeys_slint_t, pepckeys_slint;  /* deprecated 'pepckeys_slint' */

/* sl_type pepckeys_slindex_t */
typedef pepckeys_sl_index_type_c pepckeys_slindex_t;

/* sl_type pepckeys_slkey_t */
typedef pepckeys_sl_key_type_c pepckeys_slkey_t;

/* sl_type pepckeys_slkey_pure_t */
typedef pepckeys_sl_key_pure_type_c pepckeys_slkey_pure_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type pepckeys_sldata0_t */
#ifdef pepckeys_sl_data0_type_c
typedef pepckeys_sl_data0_type_c pepckeys_sldata0_t;
#endif
/* sl_type pepckeys_sldata1_t */
#ifdef pepckeys_sl_data1_type_c
typedef pepckeys_sl_data1_type_c pepckeys_sldata1_t;
#endif
/* sl_type pepckeys_sldata2_t */
#ifdef pepckeys_sl_data2_type_c
typedef pepckeys_sl_data2_type_c pepckeys_sldata2_t;
#endif
/* sl_type pepckeys_sldata3_t */
#ifdef pepckeys_sl_data3_type_c
typedef pepckeys_sl_data3_type_c pepckeys_sldata3_t;
#endif
/* sl_type pepckeys_sldata4_t */
#ifdef pepckeys_sl_data4_type_c
typedef pepckeys_sl_data4_type_c pepckeys_sldata4_t;
#endif
/* sl_type pepckeys_sldata5_t */
#ifdef pepckeys_sl_data5_type_c
typedef pepckeys_sl_data5_type_c pepckeys_sldata5_t;
#endif
/* sl_type pepckeys_sldata6_t */
#ifdef pepckeys_sl_data6_type_c
typedef pepckeys_sl_data6_type_c pepckeys_sldata6_t;
#endif
/* sl_type pepckeys_sldata7_t */
#ifdef pepckeys_sl_data7_type_c
typedef pepckeys_sl_data7_type_c pepckeys_sldata7_t;
#endif
/* sl_type pepckeys_sldata8_t */
#ifdef pepckeys_sl_data8_type_c
typedef pepckeys_sl_data8_type_c pepckeys_sldata8_t;
#endif
/* sl_type pepckeys_sldata9_t */
#ifdef pepckeys_sl_data9_type_c
typedef pepckeys_sl_data9_type_c pepckeys_sldata9_t;
#endif
/* sl_type pepckeys_sldata10_t */
#ifdef pepckeys_sl_data10_type_c
typedef pepckeys_sl_data10_type_c pepckeys_sldata10_t;
#endif
/* sl_type pepckeys_sldata11_t */
#ifdef pepckeys_sl_data11_type_c
typedef pepckeys_sl_data11_type_c pepckeys_sldata11_t;
#endif
/* sl_type pepckeys_sldata12_t */
#ifdef pepckeys_sl_data12_type_c
typedef pepckeys_sl_data12_type_c pepckeys_sldata12_t;
#endif
/* sl_type pepckeys_sldata13_t */
#ifdef pepckeys_sl_data13_type_c
typedef pepckeys_sl_data13_type_c pepckeys_sldata13_t;
#endif
/* sl_type pepckeys_sldata14_t */
#ifdef pepckeys_sl_data14_type_c
typedef pepckeys_sl_data14_type_c pepckeys_sldata14_t;
#endif
/* sl_type pepckeys_sldata15_t */
#ifdef pepckeys_sl_data15_type_c
typedef pepckeys_sl_data15_type_c pepckeys_sldata15_t;
#endif
/* sl_type pepckeys_sldata16_t */
#ifdef pepckeys_sl_data16_type_c
typedef pepckeys_sl_data16_type_c pepckeys_sldata16_t;
#endif
/* sl_type pepckeys_sldata17_t */
#ifdef pepckeys_sl_data17_type_c
typedef pepckeys_sl_data17_type_c pepckeys_sldata17_t;
#endif
/* sl_type pepckeys_sldata18_t */
#ifdef pepckeys_sl_data18_type_c
typedef pepckeys_sl_data18_type_c pepckeys_sldata18_t;
#endif
/* sl_type pepckeys_sldata19_t */
#ifdef pepckeys_sl_data19_type_c
typedef pepckeys_sl_data19_type_c pepckeys_sldata19_t;
#endif
/* DATAX_TEMPLATE_END */

/* sl_type pepckeys__elements_t pepckeys_elements_t */
typedef struct pepckeys__elements_t
{
  pepckeys_slint_t size, max_size;
  pepckeys_slkey_t *keys;

#ifdef pepckeys_SL_INDEX
  pepckeys_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
  pepckeys_sldata0_t *data0;
#endif
#ifdef pepckeys_SL_DATA1
  pepckeys_sldata1_t *data1;
#endif
#ifdef pepckeys_SL_DATA2
  pepckeys_sldata2_t *data2;
#endif
#ifdef pepckeys_SL_DATA3
  pepckeys_sldata3_t *data3;
#endif
#ifdef pepckeys_SL_DATA4
  pepckeys_sldata4_t *data4;
#endif
#ifdef pepckeys_SL_DATA5
  pepckeys_sldata5_t *data5;
#endif
#ifdef pepckeys_SL_DATA6
  pepckeys_sldata6_t *data6;
#endif
#ifdef pepckeys_SL_DATA7
  pepckeys_sldata7_t *data7;
#endif
#ifdef pepckeys_SL_DATA8
  pepckeys_sldata8_t *data8;
#endif
#ifdef pepckeys_SL_DATA9
  pepckeys_sldata9_t *data9;
#endif
#ifdef pepckeys_SL_DATA10
  pepckeys_sldata10_t *data10;
#endif
#ifdef pepckeys_SL_DATA11
  pepckeys_sldata11_t *data11;
#endif
#ifdef pepckeys_SL_DATA12
  pepckeys_sldata12_t *data12;
#endif
#ifdef pepckeys_SL_DATA13
  pepckeys_sldata13_t *data13;
#endif
#ifdef pepckeys_SL_DATA14
  pepckeys_sldata14_t *data14;
#endif
#ifdef pepckeys_SL_DATA15
  pepckeys_sldata15_t *data15;
#endif
#ifdef pepckeys_SL_DATA16
  pepckeys_sldata16_t *data16;
#endif
#ifdef pepckeys_SL_DATA17
  pepckeys_sldata17_t *data17;
#endif
#ifdef pepckeys_SL_DATA18
  pepckeys_sldata18_t *data18;
#endif
#ifdef pepckeys_SL_DATA19
  pepckeys_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} pepckeys_elements_t;


/* sl_type pepckeys__packed_element_t pepckeys_packed_element_t */
typedef struct pepckeys__packed_element_t
{
  pepckeys_slkey_t key;

#ifdef SL_PACKED_INDEX
  pepckeys_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
  pepckeys_sldata0_t data0[pepckeys_sl_data0_size_c];
#endif
#ifdef pepckeys_SL_DATA1
  pepckeys_sldata1_t data1[pepckeys_sl_data1_size_c];
#endif
#ifdef pepckeys_SL_DATA2
  pepckeys_sldata2_t data2[pepckeys_sl_data2_size_c];
#endif
#ifdef pepckeys_SL_DATA3
  pepckeys_sldata3_t data3[pepckeys_sl_data3_size_c];
#endif
#ifdef pepckeys_SL_DATA4
  pepckeys_sldata4_t data4[pepckeys_sl_data4_size_c];
#endif
#ifdef pepckeys_SL_DATA5
  pepckeys_sldata5_t data5[pepckeys_sl_data5_size_c];
#endif
#ifdef pepckeys_SL_DATA6
  pepckeys_sldata6_t data6[pepckeys_sl_data6_size_c];
#endif
#ifdef pepckeys_SL_DATA7
  pepckeys_sldata7_t data7[pepckeys_sl_data7_size_c];
#endif
#ifdef pepckeys_SL_DATA8
  pepckeys_sldata8_t data8[pepckeys_sl_data8_size_c];
#endif
#ifdef pepckeys_SL_DATA9
  pepckeys_sldata9_t data9[pepckeys_sl_data9_size_c];
#endif
#ifdef pepckeys_SL_DATA10
  pepckeys_sldata10_t data10[pepckeys_sl_data10_size_c];
#endif
#ifdef pepckeys_SL_DATA11
  pepckeys_sldata11_t data11[pepckeys_sl_data11_size_c];
#endif
#ifdef pepckeys_SL_DATA12
  pepckeys_sldata12_t data12[pepckeys_sl_data12_size_c];
#endif
#ifdef pepckeys_SL_DATA13
  pepckeys_sldata13_t data13[pepckeys_sl_data13_size_c];
#endif
#ifdef pepckeys_SL_DATA14
  pepckeys_sldata14_t data14[pepckeys_sl_data14_size_c];
#endif
#ifdef pepckeys_SL_DATA15
  pepckeys_sldata15_t data15[pepckeys_sl_data15_size_c];
#endif
#ifdef pepckeys_SL_DATA16
  pepckeys_sldata16_t data16[pepckeys_sl_data16_size_c];
#endif
#ifdef pepckeys_SL_DATA17
  pepckeys_sldata17_t data17[pepckeys_sl_data17_size_c];
#endif
#ifdef pepckeys_SL_DATA18
  pepckeys_sldata18_t data18[pepckeys_sl_data18_size_c];
#endif
#ifdef pepckeys_SL_DATA19
  pepckeys_sldata19_t data19[pepckeys_sl_data19_size_c];
#endif
/* DATAX_TEMPLATE_END */

} pepckeys_packed_element_t;


/* sl_type pepckeys__packed_elements_t pepckeys_packed_elements_t */
typedef struct pepckeys__packed_elements_t
{
  pepckeys_slint_t size, max_size;
  
  pepckeys_packed_element_t *elements;
  
} pepckeys_packed_elements_t;


/* sl_type pepckeys__classification_info_t pepckeys_classification_info_t pepckeys_classification_info */
typedef struct pepckeys__classification_info_t
{
  pepckeys_slint_t nclasses;
  pepckeys_slkey_pure_t *keys;
  pepckeys_slint_t *counts;
  pepckeys_slint_t *masks;

  /* */
  pepckeys_slint_t *all_local_sizes;
  pepckeys_slint_t *local_lt_eq_counts;
  pepckeys_slint_t *all_local_lt_eq_counts;

} pepckeys_classification_info_t, pepckeys_classification_info;  /* deprecated 'pepckeys_classification_info' */


/* key2class, sl_type pepckeys_key2class_f */
typedef pepckeys_slint_t (*pepckeys_key2class_f)(pepckeys_slkey_t *, pepckeys_slint, void *);

/* pivot-element, sl_type pepckeys_pivot_f */
typedef pepckeys_slint_t (*pepckeys_pivot_f)(pepckeys_elements_t *);

/* sorting-network, sl_type pepckeys_sortnet_f pepckeys_sortnet_data_t */
typedef void *pepckeys_sortnet_data_t;
typedef pepckeys_slint_t (*pepckeys_sortnet_f)(pepckeys_slint_t size, pepckeys_slint_t rank, pepckeys_slint_t stage, pepckeys_sortnet_data_t snd, pepckeys_slint_t *up);

/* merge2, sl_type pepckeys_merge2x_f pepckeys_merge2X_f */
typedef pepckeys_slint_t (*pepckeys_merge2x_f)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
typedef pepckeys_slint_t (*pepckeys_merge2X_f)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);


/* deprecated, sl_type pepckeys_k2c_func pepckeys_pivot_func pepckeys_sn_func pepckeys_m2x_func pepckeys_m2X_func */
typedef pepckeys_key2class_f pepckeys_k2c_func;
typedef pepckeys_pivot_f pepckeys_pivot_func;
typedef pepckeys_sortnet_f pepckeys_sn_func;
typedef pepckeys_merge2x_f pepckeys_m2x_func;
typedef pepckeys_merge2X_f pepckeys_m2X_func;


/* partition conditions, sl_type pepckeys__partcond_t pepckeys_partcond_t */
typedef struct pepckeys__partcond_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} pepckeys_partcond_t;

/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_adds.h
 *  timestamp: 2009-11-23 13:06:23 +0100
 *  
 */




/* sl_macro pepckeys_elem_set_size pepckeys_elem_set_max_size pepckeys_elem_set_keys pepckeys_elem_set_indices */
#define pepckeys_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepckeys_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepckeys_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define pepckeys_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))

/* sl_macro pepckeys_pelem_set_size pepckeys_pelem_set_max_size pepckeys_pelem_set_elements */
#define pepckeys_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepckeys_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepckeys_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))


#ifndef SL_FRACRES
# define SL_FRACRES       1000000.0
# define SL_FRAC2INT(p)   ((pepckeys_slint_t) ((p) * -SL_FRACRES))
# define SL_INT2FRAC(i)   (((double) (i)) / -SL_FRACRES)
# define SL_ISFRAC(i)     ((i) < 0)
#endif


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos.h
 *  timestamp: 2010-01-05 17:56:41 +0100
 *  
 */



/* src/core/checksum_crc.c */
unsigned short pepckeys_cs_crc16(pepckeys_elements_t *s, pepckeys_slint n, pepckeys_slint keys, pepckeys_slint data);
unsigned int pepckeys_cs_crc32(pepckeys_elements_t *s, pepckeys_slint n, pepckeys_slint keys, pepckeys_slint data);

/* src/core/elements.c */
pepckeys_slint_t pepckeys_elements_alloc(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t keys, pepckeys_slint_t data);
pepckeys_slint_t pepckeys_elements_free(pepckeys_elements_t *s);
pepckeys_slint_t pepckeys_elements_alloc_from_block(pepckeys_elements_t *s, void *block, pepckeys_slint_t blocksize, pepckeys_slint_t alignment);
pepckeys_slint_t pepckeys_elements_copy(pepckeys_elements_t *s, pepckeys_elements_t *d);
pepckeys_slint_t pepckeys_elements_copy_at(pepckeys_elements_t *s, pepckeys_slint_t sat, pepckeys_elements_t *d, pepckeys_slint_t dat);
pepckeys_slint_t pepckeys_elements_ncopy(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t n);
pepckeys_slint_t pepckeys_elements_nmove(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t n);
pepckeys_slint_t pepckeys_elements_printf(pepckeys_elements_t *s);
pepckeys_slint_t pepckeys_elements_extract(pepckeys_elements_t *src, pepckeys_slint_t nelements, pepckeys_elements_t *dst0, pepckeys_elements_t *dst1);
pepckeys_slint_t pepckeys_elements_touch(pepckeys_elements_t *s);
pepckeys_slint_t pepckeys_elements_random_exchange(pepckeys_elements_t *s, pepckeys_slint_t rounds, pepckeys_elements_t *xs);
pepckeys_slint_t pepckeys_elements_init_keys(pepckeys_elements_t *s, pepckeys_slint_t dtype, pepckeys_slint_t _min, pepckeys_slint_t _max);
pepckeys_slint_t pepckeys_elements_init_keys_from_file(pepckeys_elements_t *s, pepckeys_slint_t data, char *filename, pepckeys_slint_t from, pepckeys_slint_t to, pepckeys_slint_t const_bytes_per_line);
pepckeys_slint_t pepckeys_elements_save_keys_to_file(pepckeys_elements_t *s, char *filename);
pepckeys_slint_t pepckeys_elements_validate_order(pepckeys_elements_t *s, pepckeys_slint_t n);
pepckeys_slint_t pepckeys_elements_validate_order_bmask(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t pepckeys_elements_validate_order_weight(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slkey_pure_t weight);
pepckeys_slint_t pepckeys_elements_print_keys(pepckeys_elements_t *s);

/* src/core/elements_packed.c */
pepckeys_slint_t pepckeys_elements_alloc_packed(pepckeys_packed_elements_t *s, pepckeys_slint_t nelements);
pepckeys_slint_t pepckeys_elements_free_packed(pepckeys_packed_elements_t *s);
pepckeys_slint_t pepckeys_elements_pack_indexed(pepckeys_elements_t *s, pepckeys_packed_elements_t *d, pepckeys_slindex_t *rindx, pepckeys_slindex_t *windx);
pepckeys_slint_t pepckeys_elements_pack(pepckeys_elements_t *s, pepckeys_packed_elements_t *d);
pepckeys_slint_t pepckeys_elements_unpack_indexed(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slindex_t *rindx, pepckeys_slindex_t *windx);
pepckeys_slint_t pepckeys_elements_unpack(pepckeys_packed_elements_t *s, pepckeys_elements_t *d);
pepckeys_slint_t pepckeys_elements_unpack_keys(pepckeys_packed_elements_t *s, pepckeys_slkey_t *k);

/* src/core/key2class.c */
pepckeys_slint pepckeys_key2class_equal(pepckeys_slkey_t *k, pepckeys_slint i, void *ci);
pepckeys_slint pepckeys_key2class_split(pepckeys_slkey_t *k, pepckeys_slint i, void *ci);
pepckeys_slint pepckeys_key2class_split_keys(pepckeys_slkey_t *k, pepckeys_slint i, void *ci);
pepckeys_slint pepckeys_key2class_random(pepckeys_slkey_t *k, pepckeys_slint i, void *ci);
pepckeys_slint pepckeys_key2class_ci_nocounts(pepckeys_slkey_t *k, pepckeys_slint i, void *ci);
pepckeys_slint pepckeys_key2class_ci(pepckeys_slkey_t *k, pepckeys_slint i, void *ci);

/* src/core/merge2_basic.c */
pepckeys_slint pepckeys_merge2_basic_01_x(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_m2x_func _x0_1, pepckeys_m2x_func _0x_1);
pepckeys_slint pepckeys_merge2_basic_01_X(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t, pepckeys_m2X_func _X0_1, pepckeys_m2X_func _0X_1);

/* src/core/merge2_basic_auto.c */
pepckeys_slint pepckeys_merge2_basic_auto_01_x(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);

/* src/core/merge2_basic_sbin.c */
pepckeys_slint pepckeys_merge2_basic_sbin_x0_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_sbin_0x_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_sbin_01_x(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_sbin_01(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);

/* src/core/merge2_basic_shyb.c */
pepckeys_slint pepckeys_merge2_basic_shyb_x0_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_shyb_0x_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_shyb_01_x(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_shyb_01(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);

/* src/core/merge2_basic_sseq.c */
pepckeys_slint pepckeys_merge2_basic_sseq_x0_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_sseq_0x_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_sseq_01_x(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_sseq_01(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);

/* src/core/merge2_basic_straight.c */
pepckeys_slint pepckeys_merge2_basic_straight_x0_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_straight_0x_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_straight_01_x(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_straight_x_0_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint pepckeys_merge2_basic_straight_X0_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint pepckeys_merge2_basic_straight_0X_1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint pepckeys_merge2_basic_straight_01_X(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint pepckeys_merge2_basic_straight_X0_1u(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);

/* src/core/pepckeys_merge2_compo_hula.c */
pepckeys_slint pepckeys_merge2_compo_hula(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs);

/* src/core/pepckeys_merge2_compo_tridgell.c */
pepckeys_slint pepckeys_merge2_compo_tridgell(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);

/* src/core/pepckeys_merge2_memory_adaptive.c */
pepckeys_slint pepckeys_merge2_memory_adaptive(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);

/* src/core/merge2_simplify.c */
pepckeys_slint pepckeys_merge2_simplify_s1(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_slint s1elements);

/* src/core/pepckeys_mergep_heap.c */
pepckeys_slint pepckeys_mergep_heap(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);
pepckeys_slint pepckeys_mergep_heap_unpack(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);

/* src/core/search.c */
pepckeys_slint pepckeys_sl_search_sequential_lt(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_sequential_le(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_sequential_gt(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_sequential_ge(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_binary_lt(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_binary_le(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_binary_gt(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_binary_ge(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint pepckeys_sl_search_hybrid_lt(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);
pepckeys_slint pepckeys_sl_search_hybrid_le(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);
pepckeys_slint pepckeys_sl_search_hybrid_gt(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);
pepckeys_slint pepckeys_sl_search_hybrid_ge(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);

/* src/core/sl_common.c */
pepckeys_slint pepckeys_ilog2c(pepckeys_slint x);
pepckeys_slint pepckeys_ilog2f(pepckeys_slint x);
pepckeys_slint pepckeys_print_bits(pepckeys_slint v);
pepckeys_slint pepckeys_pivot_random(pepckeys_elements_t *s);

/* src/core/sl_elem.c */
pepckeys_slint_t pepckeys_elem_set_data(pepckeys_elements_t *e, ...);
pepckeys_slint_t pepckeys_elem_reverse(pepckeys_elements_t *e, pepckeys_elements_t *t);
pepckeys_slint_t pepckeys_elem_nxchange_at(pepckeys_elements_t *e0, pepckeys_slint_t at0, pepckeys_elements_t *e1, pepckeys_slint_t at1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t pepckeys_elem_nxchange(pepckeys_elements_t *e0, pepckeys_elements_t *e1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t pepckeys_elem_nxchange_ro0(pepckeys_elements_t *e0, pepckeys_elements_t *e1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t pepckeys_elem_rotate(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t pepckeys_elem_rotate_ro0(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t pepckeys_elem_rotate_ro1(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);

/* src/core/pepckeys_sort_counting.c */
pepckeys_slint_t pepckeys_sort_counting_use_displs(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ndispls, pepckeys_slint_t *displs);
pepckeys_slint_t pepckeys_sort_counting_use_counts(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts, pepckeys_slint_t *counts);
pepckeys_slint_t pepckeys_sort_counting_get_counts(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts, pepckeys_slint_t *counts);
pepckeys_slint_t pepckeys_sort_counting(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts);

/* src/core/pepckeys_sort_heap.c */
pepckeys_slint pepckeys_sort_heap(pepckeys_elements_t *s, pepckeys_elements_t *xs);

/* src/core/pepckeys_sort_insert.c */
pepckeys_slint pepckeys_sort_insert(pepckeys_elements_t *s, pepckeys_elements_t *sx);

/* src/core/sort_permute.c */
pepckeys_slint pepckeys_sort_permute_forward(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint *perm, pepckeys_slint offset, pepckeys_slint mask_bit);
pepckeys_slint pepckeys_sort_permute_backward(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint *perm, pepckeys_slint offset, pepckeys_slint mask_bit);

/* src/core/pepckeys_sort_quick.c */
pepckeys_slint pepckeys_sort_quick(pepckeys_elements_t *s, pepckeys_elements_t *xs);

/* src/core/pepckeys_sort_radix.c */
pepckeys_slint pepckeys_sort_radix(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint rhigh, pepckeys_slint rlow, pepckeys_slint rwidth);

/* src/core/pepckeys_sort_radix_1bit.c */
pepckeys_slint pepckeys_sort_radix_1bit(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint rhigh, pepckeys_slint rlow);

/* src/core/pepckeys_sort_radix_af.c */
pepckeys_slint pepckeys_sort_radix_af(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint rhigh, pepckeys_slint rlow, pepckeys_slint rwidth);

/* src/core/pepckeys_sort_radix_iter.c */
pepckeys_slint pepckeys_sort_radix_iter(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint presorted, pepckeys_slint rhigh, pepckeys_slint rlow, pepckeys_slint rwidth);

/* src/core/sortnet.c */
pepckeys_slint pepckeys_sn_hypercube_lh(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint pepckeys_sn_hypercube_hl(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint pepckeys_sn_odd_even_trans(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint pepckeys_sn_batcher(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint pepckeys_sn_bitonic(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint pepckeys_sn_connected(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);

/* src/core/splitx.c */
pepckeys_slint pepckeys_split2_lt_ge(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t);
pepckeys_slint pepckeys_split2_le_gt(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t);
pepckeys_slint pepckeys_split3_lt_eq_gt(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t, pepckeys_slint *nlt, pepckeys_slint *nle);
pepckeys_slint pepckeys_split3_lt_eq_gt_old(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t, pepckeys_slint *nlt, pepckeys_slint *nle);
pepckeys_slint pepckeys_split2_b(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slkey_pure_t bmask);
pepckeys_slint pepckeys_splitk_k2c_af(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);
pepckeys_slint pepckeys_splitk_k2c(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);
pepckeys_slint pepckeys_splitk_k2c_count(pepckeys_elements_t *s, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI

/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  timestamp: 2010-01-05 17:56:41 +0100
 *  
 */



/* src/core_mpi/mpi_common.c */
pepckeys_slint_t pepckeys_mpi_datatypes_init();
pepckeys_slint_t pepckeys_mpi_datatypes_release();

/* src/core_mpi/mpi_elements.c */
pepckeys_slint pepckeys_mpi_elements_init_keys_from_file(pepckeys_elements_t *s, char *filename, pepckeys_slint from, pepckeys_slint to, pepckeys_slint const_bytes_per_line, pepckeys_slint root, int size, int rank, MPI_Comm comm);
pepckeys_slint pepckeys_mpi_elements_validate_order(pepckeys_elements_t *s, pepckeys_slint n, int size, int rank, MPI_Comm comm);
unsigned short pepckeys_mpi_cs16(pepckeys_elements_t *s, pepckeys_slint n, pepckeys_slint keys, pepckeys_slint data, int size, int rank, MPI_Comm comm);
unsigned int pepckeys_mpi_cs32(pepckeys_elements_t *s, pepckeys_slint n, pepckeys_slint keys, pepckeys_slint data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
pepckeys_slint_t pepckeys_mpi_elements_packed_datatype_create(MPI_Datatype *pdt, pepckeys_slint_t structured);
pepckeys_slint_t pepckeys_mpi_elements_packed_datatype_destroy(MPI_Datatype *pdt);

/* src/core_mpi/pepckeys_mpi_find_exact.c */
pepckeys_slint_t pepckeys_mpi_find_exact_equal(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *ex_start, pepckeys_slint_t *ex_size, int size, int rank, MPI_Comm comm);
pepckeys_slint_t pepckeys_mpi_find_exact(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *dst_size, pepckeys_slint_t *ex_start, pepckeys_slint_t *ex_sizes, pepckeys_slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_find_exact_old.c */
pepckeys_slint_t pepckeys_mpi_find_exact_old(pepckeys_elements_t *s, pepckeys_slint_t counterpart, pepckeys_slint_t high, pepckeys_elements_t *xs, pepckeys_slint_t *start, pepckeys_slint_t *end, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_merge2.c */
pepckeys_slint_t pepckeys_mpi_merge2(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *dst_size, pepckeys_merge2x_f m2, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_merge2_old.c */
pepckeys_slint pepckeys_mpi_merge2_old(pepckeys_elements_t *s, pepckeys_slint counterpart, pepckeys_slint high, pepckeys_m2x_func m2, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_mergek.c */
pepckeys_slint_t pepckeys_mpi_mergek_equal(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t pepckeys_mpi_mergek(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_mergek_old.c */
pepckeys_slint pepckeys_mpi_mergek_old(pepckeys_elements_t *s, pepckeys_sn_func sn, void *snp, pepckeys_m2x_func m2, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_partition_joink.c */
pepckeys_slint pepckeys_mpi_partition_joink(pepckeys_elements_t *s, pepckeys_slint *sizes, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_partition_radix.c */
pepckeys_slint_t pepckeys_mpi_partition_radix(pepckeys_elements_t *s, pepckeys_partcond_t *pc, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_partition_radix_old.c */
pepckeys_slint_t pepckeys_mpi_partition_radix_old(pepckeys_elements_t *s, pepckeys_partcond_t *pc, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_rebalance.c */
pepckeys_slint_t pepckeys_mpi_rebalance(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_slint_t stable, pepckeys_slint_t *dst_size, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sample_ci.c */
pepckeys_slint pepckeys_mpi_sample_ci_init(pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);
pepckeys_slint pepckeys_mpi_sample_ci_duplicate(pepckeys_classification_info *ci_src, pepckeys_classification_info *ci_dup, int size, int rank, MPI_Comm comm);
pepckeys_slint pepckeys_mpi_sample_ci_free(pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);
pepckeys_slint pepckeys_mpi_sample_ci_print(pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sample_complete.c */
pepckeys_slint pepckeys_mpi_sample_complete(pepckeys_elements_t *s, pepckeys_slint threshold, pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sample_permutation.c */
pepckeys_slint pepckeys_mpi_sample_permutation(pepckeys_elements_t *s, pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sample_precise_counts.c */
pepckeys_slint pepckeys_mpi_sample_precise_counts(pepckeys_elements_t *s, pepckeys_slint threshold, pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sample_select_qs.c */
pepckeys_slint pepckeys_mpi_sample_select_qs(pepckeys_elements_t *s, pepckeys_slint threshold, pepckeys_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_select_qs.c */
pepckeys_slint pepckeys_mpi_select_qs(pepckeys_elements_t *s, pepckeys_slint n, pepckeys_slint *iths, pepckeys_pivot_func pi, pepckeys_slint threshold, pepckeys_slkey_pure_t *e, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sort_merge.c */
pepckeys_slint_t pepckeys_mpi_sort_merge(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_splitk.c */
pepckeys_slint pepckeys_mpi_splitk(pepckeys_elements_t *s, pepckeys_k2c_func k2c, void *ci, pepckeys_elements_t *sx, pepckeys_elements_t *sa, pepckeys_slint *nne, pepckeys_slint *nue, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_splitk_dummy.c */
pepckeys_slint pepckeys_mpi_splitk_dummy(pepckeys_elements_t *s, pepckeys_k2c_func k2c, void *ci, pepckeys_elements_t *sx, pepckeys_slint *send_stats, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#endif /* __SL_PEPCKEYS_H__ */
