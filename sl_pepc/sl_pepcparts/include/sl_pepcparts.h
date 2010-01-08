
#ifndef __SL_PEPCPARTS_H__
#define __SL_PEPCPARTS_H__

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
#define pepcparts_sl_int_type_c          long
#define pepcparts_sl_int_type_mpi        MPI_LONG
#define pepcparts_sl_int_size_mpi        1
#define pepcparts_sl_int_type_fmt        "ld"


/* index data type */
#define pepcparts_sl_index_type_c        FINT_TYPE_C
#define pepcparts_sl_index_type_mpi      FINT_TYPE_MPI
#define pepcparts_sl_index_size_mpi      1
#define pepcparts_sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define pepcparts_SL_INDEX


/* keys */
#define pepcparts_sl_key_type_c          FINT8_TYPE_C
#define pepcparts_sl_key_type_mpi        FINT8_TYPE_MPI
#define pepcparts_sl_key_size_mpi        1
#define pepcparts_sl_key_type_fmt        FINT8_TYPE_FMT
#define pepcparts_sl_key_integer

/* data0: x */
#define pepcparts_SL_DATA0
#define pepcparts_sl_data0_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data0_size_c        1
#define pepcparts_sl_data0_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data0_size_mpi      1

/* data1: y */
#define pepcparts_SL_DATA1
#define pepcparts_sl_data1_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data1_size_c        1
#define pepcparts_sl_data1_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data1_size_mpi      1

/* data2: z */
#define pepcparts_SL_DATA2
#define pepcparts_sl_data2_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data2_size_c        1
#define pepcparts_sl_data2_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data2_size_mpi      1

/* data3: ux */
#define pepcparts_SL_DATA3
#define pepcparts_sl_data3_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data3_size_c        1
#define pepcparts_sl_data3_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data3_size_mpi      1

/* data4: uy */
#define pepcparts_SL_DATA4
#define pepcparts_sl_data4_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data4_size_c        1
#define pepcparts_sl_data4_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data4_size_mpi      1

/* data5: uz */
#define pepcparts_SL_DATA5
#define pepcparts_sl_data5_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data5_size_c        1
#define pepcparts_sl_data5_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data5_size_mpi      1

/* data6: q */
#define pepcparts_SL_DATA6
#define pepcparts_sl_data6_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data6_size_c        1
#define pepcparts_sl_data6_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data6_size_mpi      1

/* data7: m */
#define pepcparts_SL_DATA7
#define pepcparts_sl_data7_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data7_size_c        1
#define pepcparts_sl_data7_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data7_size_mpi      1

/* data8: work */
#define pepcparts_SL_DATA8
#define pepcparts_sl_data8_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data8_size_c        1
#define pepcparts_sl_data8_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data8_size_mpi      1

/* data9: ex */
#define pepcparts_SL_DATA9
#define pepcparts_sl_data9_type_c        FREAL8_TYPE_C
#define pepcparts_sl_data9_size_c        1
#define pepcparts_sl_data9_type_mpi      FREAL8_TYPE_MPI
#define pepcparts_sl_data9_size_mpi      1

/* data10: ey */
#define pepcparts_SL_DATA10
#define pepcparts_sl_data10_type_c       FREAL8_TYPE_C
#define pepcparts_sl_data10_size_c       1
#define pepcparts_sl_data10_type_mpi     FREAL8_TYPE_MPI
#define pepcparts_sl_data10_size_mpi     1

/* data11: ez */
#define pepcparts_SL_DATA11
#define pepcparts_sl_data11_type_c       FREAL8_TYPE_C
#define pepcparts_sl_data11_size_c       1
#define pepcparts_sl_data11_type_mpi     FREAL8_TYPE_MPI
#define pepcparts_sl_data11_size_mpi     1

/* data12: pelabel */
#define pepcparts_SL_DATA12
#define pepcparts_sl_data12_type_c       FINT_TYPE_C
#define pepcparts_sl_data12_size_c       1
#define pepcparts_sl_data12_type_mpi     FINT_TYPE_MPI
#define pepcparts_sl_data12_size_mpi     1


/* weighted elements */
#define pepcparts_sl_elem_weight(e, at)  ((e)->data8[at])
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


#ifndef pepcparts_SL_INDEX
# undef SL_PACKED_INDEX
#endif


/* if no special, given, primary and heavy used integer-type ... */
#ifndef pepcparts_sl_int_type_c
  /* ... use a default one */
# define pepcparts_sl_int_type_c               long      /* sl_macro */
# undef pepcparts_sl_int_type_mpi
# define pepcparts_sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef pepcparts_sl_int_size_mpi
# define pepcparts_sl_int_size_mpi             1         /* sl_macro */
# undef pepcparts_sl_int_type_fmt
# define pepcparts_sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepcparts_sl_int_type_mpi) || !defined(pepcparts_sl_int_size_mpi)
#   error "pepcparts_sl_int_type_mpi and/or pepcparts_sl_int_size_mpi missing"
#  endif
# endif
# ifndef pepcparts_sl_int_type_fmt
#  error "pepcparts_sl_int_type_fmt macro is missing, using d as default"
#  define pepcparts_sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef pepcparts_sl_index_type_c
 /* ... use the primary integer type */
# define pepcparts_sl_index_type_c             pepcparts_sl_int_type_c    /* sl_macro */
# undef pepcparts_sl_index_type_mpi
# define pepcparts_sl_index_type_mpi           pepcparts_sl_int_type_mpi  /* sl_macro */
# undef pepcparts_sl_index_size_mpi
# define pepcparts_sl_index_size_mpi           pepcparts_sl_int_size_mpi  /* sl_macro */
# undef pepcparts_sl_index_type_fmt
# define pepcparts_sl_index_type_fmt           pepcparts_sl_int_type_fmt  /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepcparts_sl_index_type_mpi) || !defined(pepcparts_sl_index_size_mpi)
#   error "pepcparts_sl_index_type_mpi and/or pepcparts_sl_index_size_mpi missing"
#  endif
# endif
# ifndef pepcparts_sl_index_type_fmt
#  error "pepcparts_sl_index_type_fmt macro is missing, using d as default"
#  define pepcparts_sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef pepcparts_sl_key_pure_type_c
# define pepcparts_sl_key_pure_type_c          pepcparts_sl_key_type_c  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_type_mpi
# define pepcparts_sl_key_pure_type_mpi        pepcparts_sl_key_type_mpi  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_size_mpi
# define pepcparts_sl_key_pure_size_mpi        pepcparts_sl_key_size_mpi  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_type_fmt
# ifdef pepcparts_sl_key_type_fmt
#  define pepcparts_sl_key_pure_type_fmt       pepcparts_sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef pepcparts_sl_key_purify
 #define pepcparts_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_get_pure
 #define pepcparts_sl_key_get_pure(k)          (pepcparts_sl_key_purify(k))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_set_pure
 #define pepcparts_sl_key_set_pure(k, p)       (pepcparts_sl_key_purify(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef pepcparts_sl_key_pure_cmp_eq
 #define pepcparts_sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_ne
 #define pepcparts_sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_lt
 #define pepcparts_sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_le
 #define pepcparts_sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_gt
 #define pepcparts_sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_pure_cmp_ge
 #define pepcparts_sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef pepcparts_sl_key_cmp_eq
 #define pepcparts_sl_key_cmp_eq(k0, k1)       (pepcparts_sl_key_pure_cmp_eq(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_ne
 #define pepcparts_sl_key_cmp_ne(k0, k1)       (pepcparts_sl_key_pure_cmp_ne(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_lt
 #define pepcparts_sl_key_cmp_lt(k0, k1)       (pepcparts_sl_key_pure_cmp_lt(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_le
 #define pepcparts_sl_key_cmp_le(k0, k1)       (pepcparts_sl_key_pure_cmp_le(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_gt
 #define pepcparts_sl_key_cmp_gt(k0, k1)       (pepcparts_sl_key_pure_cmp_gt(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_cmp_ge
 #define pepcparts_sl_key_cmp_ge(k0, k1)       (pepcparts_sl_key_pure_cmp_ge(pepcparts_sl_key_purify(k0), pepcparts_sl_key_purify(k1)))  /* sl_macro */
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
/* sl_macro pepcparts_SL_DATA0_IGNORE */
#ifdef pepcparts_SL_DATA0_IGNORE
# undef pepcparts_SL_DATA0
#endif
/* sl_macro pepcparts_SL_DATA1_IGNORE */
#ifdef pepcparts_SL_DATA1_IGNORE
# undef pepcparts_SL_DATA1
#endif
/* sl_macro pepcparts_SL_DATA2_IGNORE */
#ifdef pepcparts_SL_DATA2_IGNORE
# undef pepcparts_SL_DATA2
#endif
/* sl_macro pepcparts_SL_DATA3_IGNORE */
#ifdef pepcparts_SL_DATA3_IGNORE
# undef pepcparts_SL_DATA3
#endif
/* sl_macro pepcparts_SL_DATA4_IGNORE */
#ifdef pepcparts_SL_DATA4_IGNORE
# undef pepcparts_SL_DATA4
#endif
/* sl_macro pepcparts_SL_DATA5_IGNORE */
#ifdef pepcparts_SL_DATA5_IGNORE
# undef pepcparts_SL_DATA5
#endif
/* sl_macro pepcparts_SL_DATA6_IGNORE */
#ifdef pepcparts_SL_DATA6_IGNORE
# undef pepcparts_SL_DATA6
#endif
/* sl_macro pepcparts_SL_DATA7_IGNORE */
#ifdef pepcparts_SL_DATA7_IGNORE
# undef pepcparts_SL_DATA7
#endif
/* sl_macro pepcparts_SL_DATA8_IGNORE */
#ifdef pepcparts_SL_DATA8_IGNORE
# undef pepcparts_SL_DATA8
#endif
/* sl_macro pepcparts_SL_DATA9_IGNORE */
#ifdef pepcparts_SL_DATA9_IGNORE
# undef pepcparts_SL_DATA9
#endif
/* sl_macro pepcparts_SL_DATA10_IGNORE */
#ifdef pepcparts_SL_DATA10_IGNORE
# undef pepcparts_SL_DATA10
#endif
/* sl_macro pepcparts_SL_DATA11_IGNORE */
#ifdef pepcparts_SL_DATA11_IGNORE
# undef pepcparts_SL_DATA11
#endif
/* sl_macro pepcparts_SL_DATA12_IGNORE */
#ifdef pepcparts_SL_DATA12_IGNORE
# undef pepcparts_SL_DATA12
#endif
/* sl_macro pepcparts_SL_DATA13_IGNORE */
#ifdef pepcparts_SL_DATA13_IGNORE
# undef pepcparts_SL_DATA13
#endif
/* sl_macro pepcparts_SL_DATA14_IGNORE */
#ifdef pepcparts_SL_DATA14_IGNORE
# undef pepcparts_SL_DATA14
#endif
/* sl_macro pepcparts_SL_DATA15_IGNORE */
#ifdef pepcparts_SL_DATA15_IGNORE
# undef pepcparts_SL_DATA15
#endif
/* sl_macro pepcparts_SL_DATA16_IGNORE */
#ifdef pepcparts_SL_DATA16_IGNORE
# undef pepcparts_SL_DATA16
#endif
/* sl_macro pepcparts_SL_DATA17_IGNORE */
#ifdef pepcparts_SL_DATA17_IGNORE
# undef pepcparts_SL_DATA17
#endif
/* sl_macro pepcparts_SL_DATA18_IGNORE */
#ifdef pepcparts_SL_DATA18_IGNORE
# undef pepcparts_SL_DATA18
#endif
/* sl_macro pepcparts_SL_DATA19_IGNORE */
#ifdef pepcparts_SL_DATA19_IGNORE
# undef pepcparts_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* default element weights */
#ifndef pepcparts_sl_elem_weight
# define pepcparts_sl_elem_weight(e, at)       1  /* sl_macro */
#endif


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_types.h
 *  timestamp: 2009-12-03 09:12:31 +0100
 *  
 */




/* sl_type pepcparts_slint_t pepcparts_slint */
typedef pepcparts_sl_int_type_c pepcparts_slint_t, pepcparts_slint;  /* deprecated 'pepcparts_slint' */

/* sl_type pepcparts_slindex_t */
typedef pepcparts_sl_index_type_c pepcparts_slindex_t;

/* sl_type pepcparts_slkey_t */
typedef pepcparts_sl_key_type_c pepcparts_slkey_t;

/* sl_type pepcparts_slkey_pure_t */
typedef pepcparts_sl_key_pure_type_c pepcparts_slkey_pure_t;

/* DATAX_TEMPLATE_BEGIN */
/* sl_type pepcparts_sldata0_t */
#ifdef pepcparts_sl_data0_type_c
typedef pepcparts_sl_data0_type_c pepcparts_sldata0_t;
#endif
/* sl_type pepcparts_sldata1_t */
#ifdef pepcparts_sl_data1_type_c
typedef pepcparts_sl_data1_type_c pepcparts_sldata1_t;
#endif
/* sl_type pepcparts_sldata2_t */
#ifdef pepcparts_sl_data2_type_c
typedef pepcparts_sl_data2_type_c pepcparts_sldata2_t;
#endif
/* sl_type pepcparts_sldata3_t */
#ifdef pepcparts_sl_data3_type_c
typedef pepcparts_sl_data3_type_c pepcparts_sldata3_t;
#endif
/* sl_type pepcparts_sldata4_t */
#ifdef pepcparts_sl_data4_type_c
typedef pepcparts_sl_data4_type_c pepcparts_sldata4_t;
#endif
/* sl_type pepcparts_sldata5_t */
#ifdef pepcparts_sl_data5_type_c
typedef pepcparts_sl_data5_type_c pepcparts_sldata5_t;
#endif
/* sl_type pepcparts_sldata6_t */
#ifdef pepcparts_sl_data6_type_c
typedef pepcparts_sl_data6_type_c pepcparts_sldata6_t;
#endif
/* sl_type pepcparts_sldata7_t */
#ifdef pepcparts_sl_data7_type_c
typedef pepcparts_sl_data7_type_c pepcparts_sldata7_t;
#endif
/* sl_type pepcparts_sldata8_t */
#ifdef pepcparts_sl_data8_type_c
typedef pepcparts_sl_data8_type_c pepcparts_sldata8_t;
#endif
/* sl_type pepcparts_sldata9_t */
#ifdef pepcparts_sl_data9_type_c
typedef pepcparts_sl_data9_type_c pepcparts_sldata9_t;
#endif
/* sl_type pepcparts_sldata10_t */
#ifdef pepcparts_sl_data10_type_c
typedef pepcparts_sl_data10_type_c pepcparts_sldata10_t;
#endif
/* sl_type pepcparts_sldata11_t */
#ifdef pepcparts_sl_data11_type_c
typedef pepcparts_sl_data11_type_c pepcparts_sldata11_t;
#endif
/* sl_type pepcparts_sldata12_t */
#ifdef pepcparts_sl_data12_type_c
typedef pepcparts_sl_data12_type_c pepcparts_sldata12_t;
#endif
/* sl_type pepcparts_sldata13_t */
#ifdef pepcparts_sl_data13_type_c
typedef pepcparts_sl_data13_type_c pepcparts_sldata13_t;
#endif
/* sl_type pepcparts_sldata14_t */
#ifdef pepcparts_sl_data14_type_c
typedef pepcparts_sl_data14_type_c pepcparts_sldata14_t;
#endif
/* sl_type pepcparts_sldata15_t */
#ifdef pepcparts_sl_data15_type_c
typedef pepcparts_sl_data15_type_c pepcparts_sldata15_t;
#endif
/* sl_type pepcparts_sldata16_t */
#ifdef pepcparts_sl_data16_type_c
typedef pepcparts_sl_data16_type_c pepcparts_sldata16_t;
#endif
/* sl_type pepcparts_sldata17_t */
#ifdef pepcparts_sl_data17_type_c
typedef pepcparts_sl_data17_type_c pepcparts_sldata17_t;
#endif
/* sl_type pepcparts_sldata18_t */
#ifdef pepcparts_sl_data18_type_c
typedef pepcparts_sl_data18_type_c pepcparts_sldata18_t;
#endif
/* sl_type pepcparts_sldata19_t */
#ifdef pepcparts_sl_data19_type_c
typedef pepcparts_sl_data19_type_c pepcparts_sldata19_t;
#endif
/* DATAX_TEMPLATE_END */

/* sl_type pepcparts__elements_t pepcparts_elements_t */
typedef struct pepcparts__elements_t
{
  pepcparts_slint_t size, max_size;
  pepcparts_slkey_t *keys;

#ifdef pepcparts_SL_INDEX
  pepcparts_slindex_t *indices;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
  pepcparts_sldata0_t *data0;
#endif
#ifdef pepcparts_SL_DATA1
  pepcparts_sldata1_t *data1;
#endif
#ifdef pepcparts_SL_DATA2
  pepcparts_sldata2_t *data2;
#endif
#ifdef pepcparts_SL_DATA3
  pepcparts_sldata3_t *data3;
#endif
#ifdef pepcparts_SL_DATA4
  pepcparts_sldata4_t *data4;
#endif
#ifdef pepcparts_SL_DATA5
  pepcparts_sldata5_t *data5;
#endif
#ifdef pepcparts_SL_DATA6
  pepcparts_sldata6_t *data6;
#endif
#ifdef pepcparts_SL_DATA7
  pepcparts_sldata7_t *data7;
#endif
#ifdef pepcparts_SL_DATA8
  pepcparts_sldata8_t *data8;
#endif
#ifdef pepcparts_SL_DATA9
  pepcparts_sldata9_t *data9;
#endif
#ifdef pepcparts_SL_DATA10
  pepcparts_sldata10_t *data10;
#endif
#ifdef pepcparts_SL_DATA11
  pepcparts_sldata11_t *data11;
#endif
#ifdef pepcparts_SL_DATA12
  pepcparts_sldata12_t *data12;
#endif
#ifdef pepcparts_SL_DATA13
  pepcparts_sldata13_t *data13;
#endif
#ifdef pepcparts_SL_DATA14
  pepcparts_sldata14_t *data14;
#endif
#ifdef pepcparts_SL_DATA15
  pepcparts_sldata15_t *data15;
#endif
#ifdef pepcparts_SL_DATA16
  pepcparts_sldata16_t *data16;
#endif
#ifdef pepcparts_SL_DATA17
  pepcparts_sldata17_t *data17;
#endif
#ifdef pepcparts_SL_DATA18
  pepcparts_sldata18_t *data18;
#endif
#ifdef pepcparts_SL_DATA19
  pepcparts_sldata19_t *data19;
#endif
/* DATAX_TEMPLATE_END */

} pepcparts_elements_t;


/* sl_type pepcparts__packed_element_t pepcparts_packed_element_t */
typedef struct pepcparts__packed_element_t
{
  pepcparts_slkey_t key;

#ifdef SL_PACKED_INDEX
  pepcparts_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
  pepcparts_sldata0_t data0[pepcparts_sl_data0_size_c];
#endif
#ifdef pepcparts_SL_DATA1
  pepcparts_sldata1_t data1[pepcparts_sl_data1_size_c];
#endif
#ifdef pepcparts_SL_DATA2
  pepcparts_sldata2_t data2[pepcparts_sl_data2_size_c];
#endif
#ifdef pepcparts_SL_DATA3
  pepcparts_sldata3_t data3[pepcparts_sl_data3_size_c];
#endif
#ifdef pepcparts_SL_DATA4
  pepcparts_sldata4_t data4[pepcparts_sl_data4_size_c];
#endif
#ifdef pepcparts_SL_DATA5
  pepcparts_sldata5_t data5[pepcparts_sl_data5_size_c];
#endif
#ifdef pepcparts_SL_DATA6
  pepcparts_sldata6_t data6[pepcparts_sl_data6_size_c];
#endif
#ifdef pepcparts_SL_DATA7
  pepcparts_sldata7_t data7[pepcparts_sl_data7_size_c];
#endif
#ifdef pepcparts_SL_DATA8
  pepcparts_sldata8_t data8[pepcparts_sl_data8_size_c];
#endif
#ifdef pepcparts_SL_DATA9
  pepcparts_sldata9_t data9[pepcparts_sl_data9_size_c];
#endif
#ifdef pepcparts_SL_DATA10
  pepcparts_sldata10_t data10[pepcparts_sl_data10_size_c];
#endif
#ifdef pepcparts_SL_DATA11
  pepcparts_sldata11_t data11[pepcparts_sl_data11_size_c];
#endif
#ifdef pepcparts_SL_DATA12
  pepcparts_sldata12_t data12[pepcparts_sl_data12_size_c];
#endif
#ifdef pepcparts_SL_DATA13
  pepcparts_sldata13_t data13[pepcparts_sl_data13_size_c];
#endif
#ifdef pepcparts_SL_DATA14
  pepcparts_sldata14_t data14[pepcparts_sl_data14_size_c];
#endif
#ifdef pepcparts_SL_DATA15
  pepcparts_sldata15_t data15[pepcparts_sl_data15_size_c];
#endif
#ifdef pepcparts_SL_DATA16
  pepcparts_sldata16_t data16[pepcparts_sl_data16_size_c];
#endif
#ifdef pepcparts_SL_DATA17
  pepcparts_sldata17_t data17[pepcparts_sl_data17_size_c];
#endif
#ifdef pepcparts_SL_DATA18
  pepcparts_sldata18_t data18[pepcparts_sl_data18_size_c];
#endif
#ifdef pepcparts_SL_DATA19
  pepcparts_sldata19_t data19[pepcparts_sl_data19_size_c];
#endif
/* DATAX_TEMPLATE_END */

} pepcparts_packed_element_t;


/* sl_type pepcparts__packed_elements_t pepcparts_packed_elements_t */
typedef struct pepcparts__packed_elements_t
{
  pepcparts_slint_t size, max_size;
  
  pepcparts_packed_element_t *elements;
  
} pepcparts_packed_elements_t;


/* sl_type pepcparts__classification_info_t pepcparts_classification_info_t pepcparts_classification_info */
typedef struct pepcparts__classification_info_t
{
  pepcparts_slint_t nclasses;
  pepcparts_slkey_pure_t *keys;
  pepcparts_slint_t *counts;
  pepcparts_slint_t *masks;

  /* */
  pepcparts_slint_t *all_local_sizes;
  pepcparts_slint_t *local_lt_eq_counts;
  pepcparts_slint_t *all_local_lt_eq_counts;

} pepcparts_classification_info_t, pepcparts_classification_info;  /* deprecated 'pepcparts_classification_info' */


/* key2class, sl_type pepcparts_key2class_f */
typedef pepcparts_slint_t (*pepcparts_key2class_f)(pepcparts_slkey_t *, pepcparts_slint, void *);

/* pivot-element, sl_type pepcparts_pivot_f */
typedef pepcparts_slint_t (*pepcparts_pivot_f)(pepcparts_elements_t *);

/* sorting-network, sl_type pepcparts_sortnet_f pepcparts_sortnet_data_t */
typedef void *pepcparts_sortnet_data_t;
typedef pepcparts_slint_t (*pepcparts_sortnet_f)(pepcparts_slint_t size, pepcparts_slint_t rank, pepcparts_slint_t stage, pepcparts_sortnet_data_t snd, pepcparts_slint_t *up);

/* merge2, sl_type pepcparts_merge2x_f pepcparts_merge2X_f */
typedef pepcparts_slint_t (*pepcparts_merge2x_f)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
typedef pepcparts_slint_t (*pepcparts_merge2X_f)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);


/* deprecated, sl_type pepcparts_k2c_func pepcparts_pivot_func pepcparts_sn_func pepcparts_m2x_func pepcparts_m2X_func */
typedef pepcparts_key2class_f pepcparts_k2c_func;
typedef pepcparts_pivot_f pepcparts_pivot_func;
typedef pepcparts_sortnet_f pepcparts_sn_func;
typedef pepcparts_merge2x_f pepcparts_m2x_func;
typedef pepcparts_merge2X_f pepcparts_m2X_func;


/* partition conditions, sl_type pepcparts__partcond_t pepcparts_partcond_t */
typedef struct pepcparts__partcond_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} pepcparts_partcond_t;

/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_adds.h
 *  timestamp: 2009-11-23 13:06:23 +0100
 *  
 */




/* sl_macro pepcparts_elem_set_size pepcparts_elem_set_max_size pepcparts_elem_set_keys pepcparts_elem_set_indices */
#define pepcparts_elem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepcparts_elem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepcparts_elem_set_keys(_e_, _k_)      ((_e_)->keys = (_k_))
#define pepcparts_elem_set_indices(_e_, _i_)   ((_e_)->indices = (_i_))

/* sl_macro pepcparts_pelem_set_size pepcparts_pelem_set_max_size pepcparts_pelem_set_elements */
#define pepcparts_pelem_set_size(_e_, _s_)      ((_e_)->size = (_s_))
#define pepcparts_pelem_set_max_size(_e_, _s_)  ((_e_)->max_size = (_s_))
#define pepcparts_pelem_set_elements(_e_, _l_)  ((_e_)->elements = (_l_))


#ifndef SL_FRACRES
# define SL_FRACRES       1000000.0
# define SL_FRAC2INT(p)   ((pepcparts_slint_t) ((p) * -SL_FRACRES))
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
unsigned short pepcparts_cs_crc16(pepcparts_elements_t *s, pepcparts_slint n, pepcparts_slint keys, pepcparts_slint data);
unsigned int pepcparts_cs_crc32(pepcparts_elements_t *s, pepcparts_slint n, pepcparts_slint keys, pepcparts_slint data);

/* src/core/elements.c */
pepcparts_slint_t pepcparts_elements_alloc(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t keys, pepcparts_slint_t data);
pepcparts_slint_t pepcparts_elements_free(pepcparts_elements_t *s);
pepcparts_slint_t pepcparts_elements_alloc_from_block(pepcparts_elements_t *s, void *block, pepcparts_slint_t blocksize, pepcparts_slint_t alignment);
pepcparts_slint_t pepcparts_elements_copy(pepcparts_elements_t *s, pepcparts_elements_t *d);
pepcparts_slint_t pepcparts_elements_copy_at(pepcparts_elements_t *s, pepcparts_slint_t sat, pepcparts_elements_t *d, pepcparts_slint_t dat);
pepcparts_slint_t pepcparts_elements_ncopy(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t n);
pepcparts_slint_t pepcparts_elements_nmove(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t n);
pepcparts_slint_t pepcparts_elements_printf(pepcparts_elements_t *s);
pepcparts_slint_t pepcparts_elements_extract(pepcparts_elements_t *src, pepcparts_slint_t nelements, pepcparts_elements_t *dst0, pepcparts_elements_t *dst1);
pepcparts_slint_t pepcparts_elements_touch(pepcparts_elements_t *s);
pepcparts_slint_t pepcparts_elements_random_exchange(pepcparts_elements_t *s, pepcparts_slint_t rounds, pepcparts_elements_t *xs);
pepcparts_slint_t pepcparts_elements_init_keys(pepcparts_elements_t *s, pepcparts_slint_t dtype, pepcparts_slint_t _min, pepcparts_slint_t _max);
pepcparts_slint_t pepcparts_elements_init_keys_from_file(pepcparts_elements_t *s, pepcparts_slint_t data, char *filename, pepcparts_slint_t from, pepcparts_slint_t to, pepcparts_slint_t const_bytes_per_line);
pepcparts_slint_t pepcparts_elements_save_keys_to_file(pepcparts_elements_t *s, char *filename);
pepcparts_slint_t pepcparts_elements_validate_order(pepcparts_elements_t *s, pepcparts_slint_t n);
pepcparts_slint_t pepcparts_elements_validate_order_bmask(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t pepcparts_elements_validate_order_weight(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slkey_pure_t weight);
pepcparts_slint_t pepcparts_elements_print_keys(pepcparts_elements_t *s);

/* src/core/elements_packed.c */
pepcparts_slint_t pepcparts_elements_alloc_packed(pepcparts_packed_elements_t *s, pepcparts_slint_t nelements);
pepcparts_slint_t pepcparts_elements_free_packed(pepcparts_packed_elements_t *s);
pepcparts_slint_t pepcparts_elements_pack_indexed(pepcparts_elements_t *s, pepcparts_packed_elements_t *d, pepcparts_slindex_t *rindx, pepcparts_slindex_t *windx);
pepcparts_slint_t pepcparts_elements_pack(pepcparts_elements_t *s, pepcparts_packed_elements_t *d);
pepcparts_slint_t pepcparts_elements_unpack_indexed(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slindex_t *rindx, pepcparts_slindex_t *windx);
pepcparts_slint_t pepcparts_elements_unpack(pepcparts_packed_elements_t *s, pepcparts_elements_t *d);
pepcparts_slint_t pepcparts_elements_unpack_keys(pepcparts_packed_elements_t *s, pepcparts_slkey_t *k);

/* src/core/key2class.c */
pepcparts_slint pepcparts_key2class_equal(pepcparts_slkey_t *k, pepcparts_slint i, void *ci);
pepcparts_slint pepcparts_key2class_split(pepcparts_slkey_t *k, pepcparts_slint i, void *ci);
pepcparts_slint pepcparts_key2class_split_keys(pepcparts_slkey_t *k, pepcparts_slint i, void *ci);
pepcparts_slint pepcparts_key2class_random(pepcparts_slkey_t *k, pepcparts_slint i, void *ci);
pepcparts_slint pepcparts_key2class_ci_nocounts(pepcparts_slkey_t *k, pepcparts_slint i, void *ci);
pepcparts_slint pepcparts_key2class_ci(pepcparts_slkey_t *k, pepcparts_slint i, void *ci);

/* src/core/merge2_basic.c */
pepcparts_slint pepcparts_merge2_basic_01_x(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_m2x_func _x0_1, pepcparts_m2x_func _0x_1);
pepcparts_slint pepcparts_merge2_basic_01_X(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t, pepcparts_m2X_func _X0_1, pepcparts_m2X_func _0X_1);

/* src/core/merge2_basic_auto.c */
pepcparts_slint pepcparts_merge2_basic_auto_01_x(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);

/* src/core/merge2_basic_sbin.c */
pepcparts_slint pepcparts_merge2_basic_sbin_x0_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_sbin_0x_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_sbin_01_x(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_sbin_01(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);

/* src/core/merge2_basic_shyb.c */
pepcparts_slint pepcparts_merge2_basic_shyb_x0_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_shyb_0x_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_shyb_01_x(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_shyb_01(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);

/* src/core/merge2_basic_sseq.c */
pepcparts_slint pepcparts_merge2_basic_sseq_x0_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_sseq_0x_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_sseq_01_x(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_sseq_01(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);

/* src/core/merge2_basic_straight.c */
pepcparts_slint pepcparts_merge2_basic_straight_x0_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_straight_0x_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_straight_01_x(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_straight_x_0_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint pepcparts_merge2_basic_straight_X0_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint pepcparts_merge2_basic_straight_0X_1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint pepcparts_merge2_basic_straight_01_X(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint pepcparts_merge2_basic_straight_X0_1u(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);

/* src/core/pepcparts_merge2_compo_hula.c */
pepcparts_slint pepcparts_merge2_compo_hula(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs);

/* src/core/pepcparts_merge2_compo_tridgell.c */
pepcparts_slint pepcparts_merge2_compo_tridgell(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);

/* src/core/pepcparts_merge2_memory_adaptive.c */
pepcparts_slint pepcparts_merge2_memory_adaptive(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);

/* src/core/merge2_simplify.c */
pepcparts_slint pepcparts_merge2_simplify_s1(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_slint s1elements);

/* src/core/pepcparts_mergep_heap.c */
pepcparts_slint pepcparts_mergep_heap(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);
pepcparts_slint pepcparts_mergep_heap_unpack(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);

/* src/core/search.c */
pepcparts_slint pepcparts_sl_search_sequential_lt(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_sequential_le(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_sequential_gt(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_sequential_ge(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_binary_lt(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_binary_le(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_binary_gt(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_binary_ge(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint pepcparts_sl_search_hybrid_lt(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);
pepcparts_slint pepcparts_sl_search_hybrid_le(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);
pepcparts_slint pepcparts_sl_search_hybrid_gt(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);
pepcparts_slint pepcparts_sl_search_hybrid_ge(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);

/* src/core/sl_common.c */
pepcparts_slint pepcparts_ilog2c(pepcparts_slint x);
pepcparts_slint pepcparts_ilog2f(pepcparts_slint x);
pepcparts_slint pepcparts_print_bits(pepcparts_slint v);
pepcparts_slint pepcparts_pivot_random(pepcparts_elements_t *s);

/* src/core/sl_elem.c */
pepcparts_slint_t pepcparts_elem_set_data(pepcparts_elements_t *e, ...);
pepcparts_slint_t pepcparts_elem_reverse(pepcparts_elements_t *e, pepcparts_elements_t *t);
pepcparts_slint_t pepcparts_elem_nxchange_at(pepcparts_elements_t *e0, pepcparts_slint_t at0, pepcparts_elements_t *e1, pepcparts_slint_t at1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t pepcparts_elem_nxchange(pepcparts_elements_t *e0, pepcparts_elements_t *e1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t pepcparts_elem_nxchange_ro0(pepcparts_elements_t *e0, pepcparts_elements_t *e1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t pepcparts_elem_rotate(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t pepcparts_elem_rotate_ro0(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t pepcparts_elem_rotate_ro1(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);

/* src/core/pepcparts_sort_counting.c */
pepcparts_slint_t pepcparts_sort_counting_use_displs(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ndispls, pepcparts_slint_t *displs);
pepcparts_slint_t pepcparts_sort_counting_use_counts(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts, pepcparts_slint_t *counts);
pepcparts_slint_t pepcparts_sort_counting_get_counts(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts, pepcparts_slint_t *counts);
pepcparts_slint_t pepcparts_sort_counting(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts);

/* src/core/pepcparts_sort_heap.c */
pepcparts_slint pepcparts_sort_heap(pepcparts_elements_t *s, pepcparts_elements_t *xs);

/* src/core/pepcparts_sort_insert.c */
pepcparts_slint pepcparts_sort_insert(pepcparts_elements_t *s, pepcparts_elements_t *sx);

/* src/core/sort_permute.c */
pepcparts_slint pepcparts_sort_permute_forward(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint *perm, pepcparts_slint offset, pepcparts_slint mask_bit);
pepcparts_slint pepcparts_sort_permute_backward(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint *perm, pepcparts_slint offset, pepcparts_slint mask_bit);

/* src/core/pepcparts_sort_quick.c */
pepcparts_slint pepcparts_sort_quick(pepcparts_elements_t *s, pepcparts_elements_t *xs);

/* src/core/pepcparts_sort_radix.c */
pepcparts_slint pepcparts_sort_radix(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint rhigh, pepcparts_slint rlow, pepcparts_slint rwidth);

/* src/core/pepcparts_sort_radix_1bit.c */
pepcparts_slint pepcparts_sort_radix_1bit(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint rhigh, pepcparts_slint rlow);

/* src/core/pepcparts_sort_radix_af.c */
pepcparts_slint pepcparts_sort_radix_af(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint rhigh, pepcparts_slint rlow, pepcparts_slint rwidth);

/* src/core/pepcparts_sort_radix_iter.c */
pepcparts_slint pepcparts_sort_radix_iter(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint presorted, pepcparts_slint rhigh, pepcparts_slint rlow, pepcparts_slint rwidth);

/* src/core/sortnet.c */
pepcparts_slint pepcparts_sn_hypercube_lh(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint pepcparts_sn_hypercube_hl(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint pepcparts_sn_odd_even_trans(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint pepcparts_sn_batcher(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint pepcparts_sn_bitonic(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint pepcparts_sn_connected(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);

/* src/core/splitx.c */
pepcparts_slint pepcparts_split2_lt_ge(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t);
pepcparts_slint pepcparts_split2_le_gt(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t);
pepcparts_slint pepcparts_split3_lt_eq_gt(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t, pepcparts_slint *nlt, pepcparts_slint *nle);
pepcparts_slint pepcparts_split3_lt_eq_gt_old(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t, pepcparts_slint *nlt, pepcparts_slint *nle);
pepcparts_slint pepcparts_split2_b(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slkey_pure_t bmask);
pepcparts_slint pepcparts_splitk_k2c_af(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);
pepcparts_slint pepcparts_splitk_k2c(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);
pepcparts_slint pepcparts_splitk_k2c_count(pepcparts_elements_t *s, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI

/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  timestamp: 2010-01-05 17:56:41 +0100
 *  
 */



/* src/core_mpi/mpi_common.c */
pepcparts_slint_t pepcparts_mpi_datatypes_init();
pepcparts_slint_t pepcparts_mpi_datatypes_release();

/* src/core_mpi/mpi_elements.c */
pepcparts_slint pepcparts_mpi_elements_init_keys_from_file(pepcparts_elements_t *s, char *filename, pepcparts_slint from, pepcparts_slint to, pepcparts_slint const_bytes_per_line, pepcparts_slint root, int size, int rank, MPI_Comm comm);
pepcparts_slint pepcparts_mpi_elements_validate_order(pepcparts_elements_t *s, pepcparts_slint n, int size, int rank, MPI_Comm comm);
unsigned short pepcparts_mpi_cs16(pepcparts_elements_t *s, pepcparts_slint n, pepcparts_slint keys, pepcparts_slint data, int size, int rank, MPI_Comm comm);
unsigned int pepcparts_mpi_cs32(pepcparts_elements_t *s, pepcparts_slint n, pepcparts_slint keys, pepcparts_slint data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
pepcparts_slint_t pepcparts_mpi_elements_packed_datatype_create(MPI_Datatype *pdt, pepcparts_slint_t structured);
pepcparts_slint_t pepcparts_mpi_elements_packed_datatype_destroy(MPI_Datatype *pdt);

/* src/core_mpi/pepcparts_mpi_find_exact.c */
pepcparts_slint_t pepcparts_mpi_find_exact_equal(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *ex_start, pepcparts_slint_t *ex_size, int size, int rank, MPI_Comm comm);
pepcparts_slint_t pepcparts_mpi_find_exact(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *dst_size, pepcparts_slint_t *ex_start, pepcparts_slint_t *ex_sizes, pepcparts_slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_find_exact_old.c */
pepcparts_slint_t pepcparts_mpi_find_exact_old(pepcparts_elements_t *s, pepcparts_slint_t counterpart, pepcparts_slint_t high, pepcparts_elements_t *xs, pepcparts_slint_t *start, pepcparts_slint_t *end, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2.c */
pepcparts_slint_t pepcparts_mpi_merge2(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *dst_size, pepcparts_merge2x_f m2, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2_old.c */
pepcparts_slint pepcparts_mpi_merge2_old(pepcparts_elements_t *s, pepcparts_slint counterpart, pepcparts_slint high, pepcparts_m2x_func m2, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek.c */
pepcparts_slint_t pepcparts_mpi_mergek_equal(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t pepcparts_mpi_mergek(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek_old.c */
pepcparts_slint pepcparts_mpi_mergek_old(pepcparts_elements_t *s, pepcparts_sn_func sn, void *snp, pepcparts_m2x_func m2, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_joink.c */
pepcparts_slint pepcparts_mpi_partition_joink(pepcparts_elements_t *s, pepcparts_slint *sizes, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_radix.c */
pepcparts_slint_t pepcparts_mpi_partition_radix(pepcparts_elements_t *s, pepcparts_partcond_t *pc, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_radix_old.c */
pepcparts_slint_t pepcparts_mpi_partition_radix_old(pepcparts_elements_t *s, pepcparts_partcond_t *pc, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int *scounts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_rebalance.c */
pepcparts_slint_t pepcparts_mpi_rebalance(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_slint_t stable, pepcparts_slint_t *dst_size, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sample_ci.c */
pepcparts_slint pepcparts_mpi_sample_ci_init(pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);
pepcparts_slint pepcparts_mpi_sample_ci_duplicate(pepcparts_classification_info *ci_src, pepcparts_classification_info *ci_dup, int size, int rank, MPI_Comm comm);
pepcparts_slint pepcparts_mpi_sample_ci_free(pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);
pepcparts_slint pepcparts_mpi_sample_ci_print(pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_complete.c */
pepcparts_slint pepcparts_mpi_sample_complete(pepcparts_elements_t *s, pepcparts_slint threshold, pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_permutation.c */
pepcparts_slint pepcparts_mpi_sample_permutation(pepcparts_elements_t *s, pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_precise_counts.c */
pepcparts_slint pepcparts_mpi_sample_precise_counts(pepcparts_elements_t *s, pepcparts_slint threshold, pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sample_select_qs.c */
pepcparts_slint pepcparts_mpi_sample_select_qs(pepcparts_elements_t *s, pepcparts_slint threshold, pepcparts_classification_info *ci, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_select_qs.c */
pepcparts_slint pepcparts_mpi_select_qs(pepcparts_elements_t *s, pepcparts_slint n, pepcparts_slint *iths, pepcparts_pivot_func pi, pepcparts_slint threshold, pepcparts_slkey_pure_t *e, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sort_merge.c */
pepcparts_slint_t pepcparts_mpi_sort_merge(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_splitk.c */
pepcparts_slint pepcparts_mpi_splitk(pepcparts_elements_t *s, pepcparts_k2c_func k2c, void *ci, pepcparts_elements_t *sx, pepcparts_elements_t *sa, pepcparts_slint *nne, pepcparts_slint *nue, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_splitk_dummy.c */
pepcparts_slint pepcparts_mpi_splitk_dummy(pepcparts_elements_t *s, pepcparts_k2c_func k2c, void *ci, pepcparts_elements_t *sx, pepcparts_slint *send_stats, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#endif /* __SL_PEPCPARTS_H__ */
