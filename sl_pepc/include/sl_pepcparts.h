
#ifndef __SL_PEPCPARTS_H__
#define __SL_PEPCPARTS_H__

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif /* SL_USE_MPI */

#include "sl_rti.h"

#define SL_PROTO(_f_)  _f_

#include "fortran2c_types.h"


/* enable runtime_informations */
/*#define SL_USE_RTI
#define SL_USE_RTI_TIM*/


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

#define pepcparts_sl_data8_weight
/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_config_intern.h
 *  timestamp: 2011-02-14 10:22:21 +0100
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
# undef pepcparts_SL_PACKED_INDEX
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


/* if no special datatype for (intern) weight ... */
#ifndef pepcparts_sl_weight_type_c
 /* ... use the double */
# define pepcparts_sl_weight_type_c             double      /* sl_macro */
# undef pepcparts_sl_weight_type_mpi
# define pepcparts_sl_weight_type_mpi           MPI_DOUBLE  /* sl_macro */
# undef pepcparts_sl_weight_size_mpi
# define pepcparts_sl_weight_size_mpi           1           /* sl_macro */
# undef pepcparts_sl_weight_type_fmt
# define pepcparts_sl_weight_type_fmt           "f"         /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepcparts_sl_weight_type_mpi) || !defined(pepcparts_sl_weight_size_mpi)
#   error "pepcparts_sl_weight_type_mpi and/or pepcparts_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef pepcparts_sl_weight_type_fmt
#  error "pepcparts_sl_weight_type_fmt macro is missing, using f as default"
#  define pepcparts_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef pepcparts_sl_index_type_c
 /* ... use the primary integer type */
# define pepcparts_sl_index_type_c             pepcparts_sl_int_type_c
# undef pepcparts_sl_index_type_mpi
# define pepcparts_sl_index_type_mpi           pepcparts_sl_int_type_mpi
# undef pepcparts_sl_index_size_mpi
# define pepcparts_sl_index_size_mpi           pepcparts_sl_int_size_mpi
# undef pepcparts_sl_index_type_fmt
# define pepcparts_sl_index_type_fmt           pepcparts_sl_int_type_fmt
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
 /* key val -> key val */
 #define pepcparts_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define pepcparts_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef pepcparts_sl_key_set_pure
 /* key component pointer and key val */
 #define pepcparts_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
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


/* default random key */
#ifdef pepcparts_sl_key_integer
# if !defined(pepcparts_sl_key_val_srand) || !defined(pepcparts_sl_key_val_rand) || !defined(pepcparts_sl_key_val_rand_minmax)
#  undef pepcparts_sl_key_val_srand
#  undef pepcparts_sl_key_val_rand
#  undef pepcparts_sl_key_val_rand_minmax
#  define pepcparts_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define pepcparts_sl_key_val_rand()                     ((pepcparts_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define pepcparts_sl_key_val_rand_minmax(_min_, _max_)  ((pepcparts_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0_IGNORE
# undef pepcparts_SL_DATA0
#endif
#ifdef pepcparts_SL_DATA1_IGNORE
# undef pepcparts_SL_DATA1
#endif
#ifdef pepcparts_SL_DATA2_IGNORE
# undef pepcparts_SL_DATA2
#endif
#ifdef pepcparts_SL_DATA3_IGNORE
# undef pepcparts_SL_DATA3
#endif
#ifdef pepcparts_SL_DATA4_IGNORE
# undef pepcparts_SL_DATA4
#endif
#ifdef pepcparts_SL_DATA5_IGNORE
# undef pepcparts_SL_DATA5
#endif
#ifdef pepcparts_SL_DATA6_IGNORE
# undef pepcparts_SL_DATA6
#endif
#ifdef pepcparts_SL_DATA7_IGNORE
# undef pepcparts_SL_DATA7
#endif
#ifdef pepcparts_SL_DATA8_IGNORE
# undef pepcparts_SL_DATA8
#endif
#ifdef pepcparts_SL_DATA9_IGNORE
# undef pepcparts_SL_DATA9
#endif
#ifdef pepcparts_SL_DATA10_IGNORE
# undef pepcparts_SL_DATA10
#endif
#ifdef pepcparts_SL_DATA11_IGNORE
# undef pepcparts_SL_DATA11
#endif
#ifdef pepcparts_SL_DATA12_IGNORE
# undef pepcparts_SL_DATA12
#endif
#ifdef pepcparts_SL_DATA13_IGNORE
# undef pepcparts_SL_DATA13
#endif
#ifdef pepcparts_SL_DATA14_IGNORE
# undef pepcparts_SL_DATA14
#endif
#ifdef pepcparts_SL_DATA15_IGNORE
# undef pepcparts_SL_DATA15
#endif
#ifdef pepcparts_SL_DATA16_IGNORE
# undef pepcparts_SL_DATA16
#endif
#ifdef pepcparts_SL_DATA17_IGNORE
# undef pepcparts_SL_DATA17
#endif
#ifdef pepcparts_SL_DATA18_IGNORE
# undef pepcparts_SL_DATA18
#endif
#ifdef pepcparts_SL_DATA19_IGNORE
# undef pepcparts_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro pepcparts_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef pepcparts_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef pepcparts_sl_data0_weight
# undef pepcparts_sl_data1_weight
# undef pepcparts_sl_data2_weight
# undef pepcparts_sl_data3_weight
# undef pepcparts_sl_data4_weight
# undef pepcparts_sl_data5_weight
# undef pepcparts_sl_data6_weight
# undef pepcparts_sl_data7_weight
# undef pepcparts_sl_data8_weight
# undef pepcparts_sl_data9_weight
# undef pepcparts_sl_data10_weight
# undef pepcparts_sl_data11_weight
# undef pepcparts_sl_data12_weight
# undef pepcparts_sl_data13_weight
# undef pepcparts_sl_data14_weight
# undef pepcparts_sl_data15_weight
# undef pepcparts_sl_data16_weight
# undef pepcparts_sl_data17_weight
# undef pepcparts_sl_data18_weight
# undef pepcparts_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable pepcparts_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(pepcparts_sl_data0_weight) && !defined(pepcparts_SL_DATA0)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data1_weight) && !defined(pepcparts_SL_DATA1)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data2_weight) && !defined(pepcparts_SL_DATA2)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data3_weight) && !defined(pepcparts_SL_DATA3)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data4_weight) && !defined(pepcparts_SL_DATA4)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data5_weight) && !defined(pepcparts_SL_DATA5)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data6_weight) && !defined(pepcparts_SL_DATA6)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data7_weight) && !defined(pepcparts_SL_DATA7)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data8_weight) && !defined(pepcparts_SL_DATA8)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data9_weight) && !defined(pepcparts_SL_DATA9)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data10_weight) && !defined(pepcparts_SL_DATA10)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data11_weight) && !defined(pepcparts_SL_DATA11)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data12_weight) && !defined(pepcparts_SL_DATA12)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data13_weight) && !defined(pepcparts_SL_DATA13)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data14_weight) && !defined(pepcparts_SL_DATA14)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data15_weight) && !defined(pepcparts_SL_DATA15)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data16_weight) && !defined(pepcparts_SL_DATA16)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data17_weight) && !defined(pepcparts_SL_DATA17)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data18_weight) && !defined(pepcparts_SL_DATA18)
# undef pepcparts_sl_elem_weight
#endif
#if defined(pepcparts_sl_data19_weight) && !defined(pepcparts_SL_DATA19)
# undef pepcparts_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro pepcparts_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data0_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA1
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data1_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA2
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data2_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA3
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data3_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA4
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data4_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA5
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data5_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA6
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data6_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA7
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data7_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA8
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data8_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA9
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data9_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA10
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data10_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA11
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data11_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA12
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data12_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA13
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data13_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA14
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data14_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA15
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data15_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA16
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data16_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA17
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data17_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA18
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data18_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepcparts_SL_DATA19
# ifdef pepcparts_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepcparts_sl_data19_flex
#   define pepcparts_FLECKS_GUARD
#  endif
# endif
#endif
/* DATAX_TEMPLATE_END */


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_types.h
 *  timestamp: 2011-03-03 13:24:38 +0100
 *  
 */




/* sl_type pepcparts_slint_t pepcparts_slint */
typedef pepcparts_sl_int_type_c pepcparts_slint_t, pepcparts_slint;  /* deprecated 'pepcparts_slint' */

#define pepcparts_slint_fmt   pepcparts_sl_int_type_fmt    /* sl_macro */

/* sl_type pepcparts_slindex_t */
typedef pepcparts_sl_index_type_c pepcparts_slindex_t;

#define pepcparts_sindex_fmt  pepcparts_sl_index_type_fmt  /* sl_macro */

/* sl_type pepcparts_slkey_t */
typedef pepcparts_sl_key_type_c pepcparts_slkey_t;

/* sl_type pepcparts_slkey_pure_t pepcparts_slpkey_t */
typedef pepcparts_sl_key_pure_type_c pepcparts_slkey_pure_t, pepcparts_slpkey_t;

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

/* sl_type pepcparts_slweight_t */
typedef pepcparts_sl_weight_type_c pepcparts_slweight_t;

#define pepcparts_slweight_fmt  pepcparts_sl_weight_type_fmt  /* sl_macro */

/* sl_type pepcparts__slpwkey_t pepcparts_slpwkey_t */
typedef struct pepcparts__slpwkey_t
{
  pepcparts_slpkey_t pkey;
  pepcparts_slweight_t weight;

} pepcparts_slpwkey_t;


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

#ifdef pepcparts_SL_PACKED_INDEX
  pepcparts_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepcparts_SL_DATA0
# ifdef pepcparts_sl_data0_flex
  pepcparts_sldata0_t data0[];
# else
  pepcparts_sldata0_t data0[pepcparts_sl_data0_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA1
# ifdef pepcparts_sl_data1_flex
  pepcparts_sldata1_t data1[];
# else
  pepcparts_sldata1_t data1[pepcparts_sl_data1_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA2
# ifdef pepcparts_sl_data2_flex
  pepcparts_sldata2_t data2[];
# else
  pepcparts_sldata2_t data2[pepcparts_sl_data2_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA3
# ifdef pepcparts_sl_data3_flex
  pepcparts_sldata3_t data3[];
# else
  pepcparts_sldata3_t data3[pepcparts_sl_data3_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA4
# ifdef pepcparts_sl_data4_flex
  pepcparts_sldata4_t data4[];
# else
  pepcparts_sldata4_t data4[pepcparts_sl_data4_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA5
# ifdef pepcparts_sl_data5_flex
  pepcparts_sldata5_t data5[];
# else
  pepcparts_sldata5_t data5[pepcparts_sl_data5_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA6
# ifdef pepcparts_sl_data6_flex
  pepcparts_sldata6_t data6[];
# else
  pepcparts_sldata6_t data6[pepcparts_sl_data6_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA7
# ifdef pepcparts_sl_data7_flex
  pepcparts_sldata7_t data7[];
# else
  pepcparts_sldata7_t data7[pepcparts_sl_data7_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA8
# ifdef pepcparts_sl_data8_flex
  pepcparts_sldata8_t data8[];
# else
  pepcparts_sldata8_t data8[pepcparts_sl_data8_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA9
# ifdef pepcparts_sl_data9_flex
  pepcparts_sldata9_t data9[];
# else
  pepcparts_sldata9_t data9[pepcparts_sl_data9_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA10
# ifdef pepcparts_sl_data10_flex
  pepcparts_sldata10_t data10[];
# else
  pepcparts_sldata10_t data10[pepcparts_sl_data10_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA11
# ifdef pepcparts_sl_data11_flex
  pepcparts_sldata11_t data11[];
# else
  pepcparts_sldata11_t data11[pepcparts_sl_data11_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA12
# ifdef pepcparts_sl_data12_flex
  pepcparts_sldata12_t data12[];
# else
  pepcparts_sldata12_t data12[pepcparts_sl_data12_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA13
# ifdef pepcparts_sl_data13_flex
  pepcparts_sldata13_t data13[];
# else
  pepcparts_sldata13_t data13[pepcparts_sl_data13_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA14
# ifdef pepcparts_sl_data14_flex
  pepcparts_sldata14_t data14[];
# else
  pepcparts_sldata14_t data14[pepcparts_sl_data14_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA15
# ifdef pepcparts_sl_data15_flex
  pepcparts_sldata15_t data15[];
# else
  pepcparts_sldata15_t data15[pepcparts_sl_data15_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA16
# ifdef pepcparts_sl_data16_flex
  pepcparts_sldata16_t data16[];
# else
  pepcparts_sldata16_t data16[pepcparts_sl_data16_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA17
# ifdef pepcparts_sl_data17_flex
  pepcparts_sldata17_t data17[];
# else
  pepcparts_sldata17_t data17[pepcparts_sl_data17_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA18
# ifdef pepcparts_sl_data18_flex
  pepcparts_sldata18_t data18[];
# else
  pepcparts_sldata18_t data18[pepcparts_sl_data18_size_c];
# endif
#endif
#ifdef pepcparts_SL_DATA19
# ifdef pepcparts_sl_data19_flex
  pepcparts_sldata19_t data19[];
# else
  pepcparts_sldata19_t data19[pepcparts_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} pepcparts_packed_element_t;


/* sl_type pepcparts__packed_elements_t pepcparts_packed_elements_t */
typedef struct pepcparts__packed_elements_t
{
  pepcparts_slint_t size, max_size;
  
  pepcparts_packed_element_t *elements;
  
} pepcparts_packed_elements_t;


#ifndef SLCINT_T
#define SLCINT_T
typedef long long int slcint_t;
#define slcint_fmt  "ll"
/*#define slcint_sfx  LL*/
#endif


#define SLCM_KEYS     (((slcint_t) 1) << 0)
#define SLCM_INDICES  (((slcint_t) 1) << 1)
#define SLCM_WEIGHTS  (((slcint_t) 1) << 2)

/* DATAX_TEMPLATE_BEGIN */
#define SLCM_DATA0    (((slcint_t) 1) << (3+0))
#define SLCM_DATA1    (((slcint_t) 1) << (3+1))
#define SLCM_DATA2    (((slcint_t) 1) << (3+2))
#define SLCM_DATA3    (((slcint_t) 1) << (3+3))
#define SLCM_DATA4    (((slcint_t) 1) << (3+4))
#define SLCM_DATA5    (((slcint_t) 1) << (3+5))
#define SLCM_DATA6    (((slcint_t) 1) << (3+6))
#define SLCM_DATA7    (((slcint_t) 1) << (3+7))
#define SLCM_DATA8    (((slcint_t) 1) << (3+8))
#define SLCM_DATA9    (((slcint_t) 1) << (3+9))
#define SLCM_DATA10    (((slcint_t) 1) << (3+10))
#define SLCM_DATA11    (((slcint_t) 1) << (3+11))
#define SLCM_DATA12    (((slcint_t) 1) << (3+12))
#define SLCM_DATA13    (((slcint_t) 1) << (3+13))
#define SLCM_DATA14    (((slcint_t) 1) << (3+14))
#define SLCM_DATA15    (((slcint_t) 1) << (3+15))
#define SLCM_DATA16    (((slcint_t) 1) << (3+16))
#define SLCM_DATA17    (((slcint_t) 1) << (3+17))
#define SLCM_DATA18    (((slcint_t) 1) << (3+18))
#define SLCM_DATA19    (((slcint_t) 1) << (3+19))
/* DATAX_TEMPLATE_END */

#define SLCM_DATA     (((slcint_t) 0) \
  |SLCM_DATA0 \
  |SLCM_DATA1 \
  |SLCM_DATA2 \
  |SLCM_DATA3 \
  |SLCM_DATA4 \
  |SLCM_DATA5 \
  |SLCM_DATA6 \
  |SLCM_DATA7 \
  |SLCM_DATA8 \
  |SLCM_DATA9 \
  |SLCM_DATA10 \
  |SLCM_DATA11 \
  |SLCM_DATA12 \
  |SLCM_DATA13 \
  |SLCM_DATA14 \
  |SLCM_DATA15 \
  |SLCM_DATA16 \
  |SLCM_DATA17 \
  |SLCM_DATA18 \
  |SLCM_DATA19 \
  )

#define SLCM_ALL      (~((slcint_t) 0))


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


/* sl_type pepcparts__mergek_t pepcparts_mergek_t */
typedef struct pepcparts__mergek_t
{
  pepcparts_sortnet_f sn;
  pepcparts_sortnet_data_t snd;

  pepcparts_merge2x_f m2x;
  pepcparts_elements_t *sx;

} pepcparts_mergek_t;


/* sl_type pepcparts_keys_init_type_t pepcparts_keys_init_data_t */
typedef pepcparts_slint_t pepcparts_keys_init_type_t;
typedef void *pepcparts_keys_init_data_t;

/* sl_type pepcparts_key_set_data_t pepcparts_key_set_f */
typedef void *pepcparts_key_set_data_t;
typedef void (*pepcparts_key_set_f)(pepcparts_slkey_pure_t *k, pepcparts_key_set_data_t d);


#undef SL_EKIT_SET
#define SL_EKIT_SET         1
#undef SL_EKIT_SET_FUNC
#define SL_EKIT_SET_FUNC    2
#undef SL_EKIT_RAND
#define SL_EKIT_RAND        3
#undef SL_EKIT_RAND_QUAD
#define SL_EKIT_RAND_QUAD   4
#undef SL_EKIT_RAND_AND
#define SL_EKIT_RAND_AND    5
#undef SL_EKIT_NRAND
#define SL_EKIT_NRAND       6


#ifndef SL_EIK_OFFSET
# define SL_EIK_OFFSET     65536LL
#endif

#ifndef SL_EIK_SET
# define SL_EIK_SET        SL_EIK_OFFSET*1
#endif

#ifndef SL_EIK_RAND
# define SL_EIK_RAND       SL_EIK_OFFSET*2
#endif

#ifndef SL_EIK_RAND_QUAD
# define SL_EIK_RAND_QUAD  SL_EIK_OFFSET*3
#endif

#ifndef SL_EIK_RAND_AND
# define SL_EIK_RAND_AND   SL_EIK_OFFSET*4
#endif

#ifndef SL_EIK_RAND_NORM
# define SL_EIK_RAND_NORM  SL_EIK_OFFSET*5
#endif


/* pepcparts_elements_keys_stats */
#ifndef SL_EKS_MIN
# define SL_EKS_MIN   0
#endif

#ifndef SL_EKS_MAX
# define SL_EKS_MAX   1
#endif

#ifndef SL_EKS_SUM
# define SL_EKS_SUM   2
#endif

#ifndef SL_EKS_AVG
# define SL_EKS_AVG   3
#endif

#ifndef SL_EKS_STD
# define SL_EKS_STD   4
#endif

#ifndef SL_EKS_SIZE
# define SL_EKS_SIZE  5
#endif


#ifndef SL_SORTED_IN
# define SL_SORTED_IN   0x1LL
#endif

#ifndef SL_SORTED_OUT
# define SL_SORTED_OUT  0x2LL
#endif


#ifndef SL_MSEG_FM_EXACT
# define SL_MSEG_FM_EXACT         0
#endif
#ifndef SL_MSEG_FM_ALLORNOTHING
# define SL_MSEG_FM_ALLORNOTHING  1
#endif
#ifndef SL_MSEG_FM_MIDDLE
# define SL_MSEG_FM_MIDDLE        2
#endif


/* partition conditions, sl_type pepcparts__partcond2_t pepcparts_partcond2_t */
typedef struct pepcparts__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} pepcparts_partcond2_t;


#ifndef SLPC_COUNTS_MM
# define SLPC_COUNTS_MM   0x1
#endif
#ifndef SLPC_COUNTS_LH
# define SLPC_COUNTS_LH   0x2
#endif
#ifndef SLPC_WEIGHTS_MM
# define SLPC_WEIGHTS_MM  0x4
#endif
#ifndef SLPC_WEIGHTS_LH
# define SLPC_WEIGHTS_LH  0x8
#endif

/* partition conditions, sl_type pepcparts__partcond_t pepcparts_partcond_t pepcparts_partcond_p */
typedef struct pepcparts__partcond_t
{
  pepcparts_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} pepcparts_partcond_t, *pepcparts_partcond_p;


/* internal partition conditions, sl_type pepcparts__partcond_intern_t pepcparts_partcond_intern_t pepcparts_partcond_intern_p */
typedef struct pepcparts__partcond_intern_t
{
  pepcparts_slint_t pcm;
  pepcparts_slint_t count_min, count_max;
  pepcparts_slint_t count_low, count_high;
#ifdef elem_weight
  double weight_min, weight_max;
  double weight_low, weight_high;
#endif

} pepcparts_partcond_intern_t, *pepcparts_partcond_intern_p;


/* sl_type pepcparts__parttype_t pepcparts_parttype_t pepcparts_parttype_p */
typedef struct pepcparts__parttype_t
{

} pepcparts_parttype_t, *pepcparts_parttype_p;


/* generic binning method */

/* sl_type pepcparts__bin_t pepcparts_bin_t */
typedef struct pepcparts__bin_t
{
  pepcparts_elements_t s;

#ifdef elem_weight
  double weight;
#endif

} pepcparts_bin_t;


/* sl_type pepcparts__splitter_t pepcparts_splitter_t */
typedef struct pepcparts__splitter_t
{
  pepcparts_slint_t n;

  int *displs;
  pepcparts_slkey_pure_t *s;
  pepcparts_slint_t *sn;

} pepcparts_splitter_t;


struct pepcparts__binning_t;

/* sl_type pepcparts_binning_pre_f pepcparts_binning_exec_f pepcparts_binning_refine_f pepcparts_binning_hit_f pepcparts_binning_finalize_f pepcparts_binning_post_f */
typedef pepcparts_slint_t (*pepcparts_binning_pre_f)(struct pepcparts__binning_t *bm);
typedef pepcparts_slint_t (*pepcparts_binning_exec_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, double *counts, double *weights);
typedef pepcparts_slint_t (*pepcparts_binning_refine_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, double *counts, double *weights, pepcparts_splitter_t *sp, pepcparts_slint_t s, pepcparts_bin_t *new_bin);
typedef pepcparts_slint_t (*pepcparts_binning_hit_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, double *counts, pepcparts_splitter_t *sp, pepcparts_slint_t s);
typedef pepcparts_slint_t (*pepcparts_binning_finalize_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, double dcw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, double *lcw, pepcparts_splitter_t *sp, pepcparts_slint_t s);
typedef pepcparts_slint_t (*pepcparts_binning_post_f)(struct pepcparts__binning_t *bm);

#ifdef SL_DEPRECATED
/* sl_type pepcparts_binning_exec_pre_old_f pepcparts_binning_exec_post_old_f pepcparts_binning_refinable_old_f pepcparts_binning_refine_old_f */
typedef pepcparts_slint_t (*pepcparts_binning_exec_pre_old_f)(struct pepcparts__binning_t *bm);
typedef pepcparts_slint_t (*pepcparts_binning_exec_post_old_f)(struct pepcparts__binning_t *bm);
typedef pepcparts_slint_t (*pepcparts_binning_refinable_old_f)(struct pepcparts__binning_t *bm);
typedef pepcparts_slint_t (*pepcparts_binning_refine_old_f)(struct pepcparts__binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, double *counts, pepcparts_bin_t *new_bin);
#endif


/* sl_type pepcparts__binning_data_t pepcparts_binning_data_t */
typedef union pepcparts__binning_data_t
{
  struct
  {
    pepcparts_slint_t rhigh, rlow, rwidth;
    pepcparts_slint_t rcurrent;
    pepcparts_slkey_pure_t bit_mask;

    pepcparts_elements_t sx;

  } radix;

} pepcparts_binning_data_t;


/* sl_type pepcparts__binning_t pepcparts_binning_t */
typedef struct pepcparts__binning_t
{
  pepcparts_slint_t nbins, max_nbins;
  
  pepcparts_binning_pre_f pre;
  pepcparts_binning_exec_f exec;
  pepcparts_binning_refine_f refine;
  pepcparts_binning_hit_f hit;
  pepcparts_binning_finalize_f finalize;
  pepcparts_binning_post_f post;

  pepcparts_slint_t sorted;

#ifdef elem_weight
  pepcparts_slint_t doweights;
#endif

#ifdef SL_DEPRECATED
  pepcparts_binning_exec_pre_old_f exec_pre_old;
  pepcparts_binning_exec_post_old_f exec_post_old;
  pepcparts_binning_refinable_old_f refinable_old;
  pepcparts_binning_refine_old_f refine_old;
#endif

  pepcparts_binning_data_t bd;

} pepcparts_binning_t;


/* sl_type pepcparts__local_bins_t pepcparts_local_bins_t */
typedef struct pepcparts__local_bins_t
{
  pepcparts_binning_t *bm;

  pepcparts_slint_t nbins, max_nbins;
  pepcparts_slint_t nelements;

#ifdef elem_weight
  pepcparts_slint_t doweights, weight_factor;
#endif

  pepcparts_slint_t nbinnings, max_nbinnings;

  pepcparts_slint_t nbins_new, last_new_b, last_new_k;
  pepcparts_bin_t *bins, *bins_new;
  pepcparts_bin_t *bins0, *bins1;

  double *counts, *bin_counts;
#ifdef elem_weight
  double *weights, *bin_weights;
#endif

  pepcparts_slint_t *bcws;
  double *cws, *bin_cws;

  pepcparts_slint_t last_exec_b;

} pepcparts_local_bins_t;


/* sl_type pepcparts__global_bins_t pepcparts_global_bins_t */
typedef struct pepcparts__global_bins_t
{
  pepcparts_binning_t *bm;
  
  pepcparts_local_bins_t lb;
  
  double *counts;
#ifdef elem_weight
  double *weights;
#endif

  pepcparts_slint_t *bcws;
  double *cws;

} pepcparts_global_bins_t;


/* sl_macro pepcparts_WEIGHT_FACTOR */
#ifdef elem_weight
# define pepcparts_WEIGHT_FACTOR  2
#else
# define pepcparts_WEIGHT_FACTOR  1
#endif


/* sl_macro pepcparts_get1d pepcparts_get2d pepcparts_get3d pepcparts_get4d */
#define pepcparts_get1d(x0)                           (x0)
#define pepcparts_get2d(x1, d0, x0)                  ((x0) + (d0) *  (x1))
#define pepcparts_get3d(x2, d1, x1, d0, x0)          ((x0) + (d0) * ((x1) + (d1) *  (x2)))
#define pepcparts_get4d(x3, d2, x2, d1, x1, d0, x0)  ((x0) + (d0) * ((x1) + (d1) * ((x2) + (d2) * (x3))))


/* sl_macro pepcparts_lb_bin_count pepcparts_lb_bin_weight */
#define pepcparts_lb_bin_count(_lb_, _b_, _j_)    ((_lb_)->bins[(_b_) * (_lb_)->nelements + _j_].s.size)
#ifdef elem_weight
# define pepcparts_lb_bin_weight(_lb_, _b_, _j_)  ((_lb_)->bins[(_b_) * (_lb_)->nelements + _j_].weight)
#else
# define pepcparts_lb_bin_weight(_lb_, _b_, _j_)  0
#endif

/* sl_macro pepcparts_lb_bin_counts pepcparts_lb_bin_weights */
#ifdef elem_weight
# define pepcparts_lb_bin_counts(_lb_, _b_, _j_, _k_)   ((_lb_)->bin_cws + pepcparts_get4d((_lb_)->bcws[_b_], (_lb_)->weight_factor, 0, (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
# define pepcparts_lb_bin_weights(_lb_, _b_, _j_, _k_)  ((_lb_)->bin_cws + pepcparts_get4d((_lb_)->bcws[_b_], (_lb_)->weight_factor, 1, (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
#else
# define pepcparts_lb_bin_counts(_lb_, _b_, _j_, _k_)   ((_lb_)->bin_cws + pepcparts_get4d((_lb_)->bcws[_b_], 1, 0, (_lb_)->nelements, _j_, (_lb_)->bm->nbins, _k_))
# define pepcparts_lb_bin_weights(_lb_, _b_, _j_, _k_)  NULL
#endif

/* sl_macro pepcparts_lb_counts pepcparts_lb_weights */
#ifdef elem_weight
# define pepcparts_lb_counts(_lb_, _b_, _k_)   ((_lb_)->cws + pepcparts_get3d((_lb_)->bcws[_b_], (_lb_)->weight_factor, 0, (_lb_)->bm->nbins, (_k_)))
# define pepcparts_lb_weights(_lb_, _b_, _k_)  ((_lb_)->cws + pepcparts_get3d((_lb_)->bcws[_b_], (_lb_)->weight_factor, 1, (_lb_)->bm->nbins, (_k_)))
#else
# define pepcparts_lb_counts(_lb_, _b_, _k_)   ((_lb_)->cws + pepcparts_get3d((_lb_)->bcws[_b_], 1, 0, (_lb_)->bm->nbins, (_k_)))
# define pepcparts_lb_weights(_lb_, _b_, _k_)  NULL
#endif

/* sl_macro pepcparts_gb_counts pepcparts_gb_weights */
#ifdef elem_weight
# define pepcparts_gb_counts(_gb_, _b_, _k_)   ((_gb_)->cws + pepcparts_get3d((_gb_)->bcws[_b_], (_gb_)->lb.weight_factor, 0, (_gb_)->bm->nbins, (_k_)))
# define pepcparts_gb_weights(_gb_, _b_, _k_)  ((_gb_)->cws + pepcparts_get3d((_gb_)->bcws[_b_], (_gb_)->lb.weight_factor, 1, (_gb_)->bm->nbins, (_k_)))
#else
# define pepcparts_gb_counts(_gb_, _b_, _k_)   ((_gb_)->cws + pepcparts_get3d((_gb_)->bcws[_b_], 1, 0, (_gb_)->bm->nbins, (_k_)))
# define pepcparts_gb_weights(_gb_, _b_, _k_)  NULL
#endif


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_adds.h
 *  timestamp: 2011-02-08 18:10:11 +0100
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


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_globals.h
 *  timestamp: 2011-03-11 09:08:28 +0100
 *  
 */




/* src/include/sl_rti_intern.h */
extern rti pepcparts_rti_env;

/* src/core/pepcparts_sort_radix.c */
extern pepcparts_slint_t pepcparts_sa_ip_threshold;

/* src/core_mpi/mpi_common.c */
#ifdef SL_USE_MPI
extern MPI_Datatype pepcparts_int_mpi_datatype;
extern MPI_Datatype pepcparts_key_mpi_datatype;
extern MPI_Datatype pepcparts_pkey_mpi_datatype;
extern MPI_Datatype pepcparts_pwkey_mpi_datatype;
extern MPI_Datatype pepcparts_index_mpi_datatype;
extern MPI_Datatype pepcparts_weight_mpi_datatype;
extern MPI_Datatype pepcparts_data_mpi_datatype[];
#endif
#ifdef SL_USE_MPI
extern int pepcparts_sl_mpi_rank;
#endif

/* src/core_mpi/mpi_elements.c */
extern void *pepcparts_me_sendrecv_replace_mem;
extern pepcparts_slint_t pepcparts_me_sendrecv_replace_memsize;
extern pepcparts_slint_t pepcparts_me_sendrecv_replace_mpi_maxsize;

/* src/core_mpi/pepcparts_mpi_select_exact_generic.c */
extern int pepcparts_mseg_root;
extern double pepcparts_mseg_border_update_count_reduction;
extern double pepcparts_mseg_border_update_weight_reduction;
extern pepcparts_slint_t pepcparts_mseg_border_update_full;
extern pepcparts_slint_t pepcparts_mseg_info_rounds;
extern pepcparts_slint_t *pepcparts_mseg_info_finish_rounds;
extern double pepcparts_mseg_info_finish_rounds_avg;
extern pepcparts_slint_t pepcparts_mseg_binnings;
extern pepcparts_slint_t pepcparts_mseg_finalize_mode;

/* src/core_mpi/mpi_select_sample.c */
extern int pepcparts_mss_root;

/* src/core_mpi/pepcparts_mpi_sort_partition.c */
extern double pepcparts_msp_t[];
extern pepcparts_slint_t pepcparts_msp_sync;


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos.h
 *  timestamp: 2011-03-11 09:08:28 +0100
 *  
 */




/* src/core/binning.c */
pepcparts_slint_t SL_PROTO(pepcparts_binning_create)(pepcparts_local_bins_t *lb, pepcparts_slint_t max_nbins, pepcparts_slint_t max_nbinnings, pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t doweights, pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_destroy)(pepcparts_local_bins_t *lb);
pepcparts_slint_t SL_PROTO(pepcparts_binning_pre)(pepcparts_local_bins_t *lb);
pepcparts_slint_t SL_PROTO(pepcparts_binning_exec_reset)(pepcparts_local_bins_t *lb);
pepcparts_slint_t SL_PROTO(pepcparts_binning_exec)(pepcparts_local_bins_t *lb, pepcparts_slint_t b);
pepcparts_slint_t SL_PROTO(pepcparts_binning_refine)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_hit)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_finalize)(pepcparts_local_bins_t *lb, pepcparts_slint_t b, pepcparts_slweight_t dcw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slweight_t *lcw, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_post)(pepcparts_local_bins_t *lb);

/* src/core/binning_radix.c */
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_create)(pepcparts_binning_t *bm, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_destroy)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_pre)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_exec)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slweight_t *counts, pepcparts_slweight_t *weights);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_refine)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slweight_t *counts, pepcparts_slweight_t *weights, pepcparts_splitter_t *sp, pepcparts_slint_t s, pepcparts_bin_t *new_bin);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_hit)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slweight_t *counts, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_finalize)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slweight_t dcw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slweight_t *lcw, pepcparts_splitter_t *sp, pepcparts_slint_t s);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_post)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_exec_pre_old)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_exec_post_old)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_refinable_old)(pepcparts_binning_t *bm);
pepcparts_slint_t SL_PROTO(pepcparts_binning_radix_refine_old)(pepcparts_binning_t *bm, pepcparts_bin_t *bin, pepcparts_slint_t k, pepcparts_slweight_t *counts, pepcparts_bin_t *new_bin);

/* src/core/elements.c */
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_free)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc2)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t keys, pepcparts_slint_t indices, pepcparts_slint_t data, pepcparts_slint_t weights);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_old)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t keys, pepcparts_slint_t data);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_from_blocks)(pepcparts_elements_t *s, pepcparts_slint_t nblocks, void **blocks, pepcparts_slint_t *blocksizes, pepcparts_slint_t alignment, pepcparts_slint_t nmax, slcint_t components);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_from_block2)(pepcparts_elements_t *s, void *block, pepcparts_slint_t blocksize, pepcparts_slint_t alignment, pepcparts_slint_t nmax, pepcparts_slint_t keys, pepcparts_slint_t indices, pepcparts_slint_t data, pepcparts_slint_t weights);
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_from_block)(pepcparts_elements_t *s, void *block, pepcparts_slint_t blocksize, pepcparts_slint_t alignment, pepcparts_slint_t nmax);
pepcparts_slint_t SL_PROTO(pepcparts_elements_copy)(pepcparts_elements_t *s, pepcparts_elements_t *d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_copy_at)(pepcparts_elements_t *s, pepcparts_slint_t sat, pepcparts_elements_t *d, pepcparts_slint_t dat);
pepcparts_slint_t SL_PROTO(pepcparts_elements_ncopy)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_elements_nmove)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_elements_printf)(pepcparts_elements_t *s, const char *prefix);
pepcparts_slint_t SL_PROTO(pepcparts_elements_extract)(pepcparts_elements_t *src, pepcparts_slint_t nelements, pepcparts_elements_t *dst0, pepcparts_elements_t *dst1);
pepcparts_slint_t SL_PROTO(pepcparts_elements_touch)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_digest_sum)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(pepcparts_elements_crc32)(pepcparts_elements_t *s, pepcparts_slint nelements, pepcparts_slint_t keys, pepcparts_slint_t data);
pepcparts_slint_t SL_PROTO(pepcparts_elements_digest_hash)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, void *hash);
pepcparts_slint_t SL_PROTO(pepcparts_elements_random_exchange)(pepcparts_elements_t *s, pepcparts_slint_t rounds, pepcparts_elements_t *xs);
pepcparts_slint_t SL_PROTO(pepcparts_elements_init_keys2)(pepcparts_elements_t *s, pepcparts_slint_t dtype, pepcparts_slkey_pure_t key_min, pepcparts_slkey_pure_t key_max);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_init)(pepcparts_elements_t *s, pepcparts_keys_init_type_t t, pepcparts_keys_init_data_t d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_init_keys_from_file)(pepcparts_elements_t *s, pepcparts_slint_t data, char *filename, pepcparts_slint_t from, pepcparts_slint_t to, pepcparts_slint_t const_bytes_per_line);
pepcparts_slint_t SL_PROTO(pepcparts_elements_save_keys_to_file)(pepcparts_elements_t *s, char *filename);
pepcparts_slint_t SL_PROTO(pepcparts_elements_validate_order)(pepcparts_elements_t *s, pepcparts_slint_t n);
pepcparts_slint_t SL_PROTO(pepcparts_elements_validate_order_bmask)(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_elements_validate_order_weight)(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slkey_pure_t weight);
pepcparts_slint_t SL_PROTO(pepcparts_elements_keys_stats)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *stats);
pepcparts_slint_t SL_PROTO(pepcparts_elements_print_keys)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_print_all)(pepcparts_elements_t *s);
pepcparts_slweight_t SL_PROTO(pepcparts_elements_get_weight)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_get_minmax_keys)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slkey_pure_t *minmaxkeys);

/* src/core/elements_packed.c */
pepcparts_slint_t SL_PROTO(pepcparts_elements_alloc_packed)(pepcparts_packed_elements_t *s, pepcparts_slint_t nelements);
pepcparts_slint_t SL_PROTO(pepcparts_elements_free_packed)(pepcparts_packed_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_elements_pack_indexed)(pepcparts_elements_t *s, pepcparts_packed_elements_t *d, pepcparts_slindex_t *rindx, pepcparts_slindex_t *windx);
pepcparts_slint_t SL_PROTO(pepcparts_elements_pack)(pepcparts_elements_t *s, pepcparts_packed_elements_t *d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack_indexed)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slindex_t *rindx, pepcparts_slindex_t *windx);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d);
pepcparts_slint_t SL_PROTO(pepcparts_elements_unpack_keys)(pepcparts_packed_elements_t *s, pepcparts_slkey_t *k);

/* src/core/merge2_common.c */
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_auto_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_m2x_func _x0_1, pepcparts_m2x_func _0x_1);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_01_X)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t, pepcparts_m2X_func _X0_1, pepcparts_m2X_func _0X_1);
pepcparts_slint SL_PROTO(pepcparts_merge2_simplify_s1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_slint s1elements);
pepcparts_slint SL_PROTO(pepcparts_merge2_memory_adaptive)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);

/* src/core/merge2_hula.c */
pepcparts_slint SL_PROTO(pepcparts_merge2_compo_hula)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs);

/* src/core/merge2_search.c */
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sseq_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sseq_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sseq_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sseq_01)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sbin_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sbin_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sbin_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_sbin_01)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_shyb_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_shyb_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_shyb_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_shyb_01)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *t);

/* src/core/merge2_straight.c */
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_x0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_0x_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_01_x)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_x_0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_X0_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_0X_1)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_01_X)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_merge2_basic_straight_X0_1u)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx, pepcparts_elements_t *t);

/* src/core/merge2_tridgell.c */
pepcparts_slint SL_PROTO(pepcparts_merge2_compo_tridgell)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *sx);

/* src/core/mergep_2way.c */
pepcparts_slint_t SL_PROTO(pepcparts_mergep_2way_ip_int)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t p, int *displs, pepcparts_merge2x_f m2x);

/* src/core/mergep_heap.c */
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_int)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, int *displs, int *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_int_idx)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, int *displs, int *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_idx)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_unpack_idx)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_mergep_heap_unpack_idxonly)(pepcparts_packed_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t p, pepcparts_slindex_t *displs, pepcparts_slindex_t *counts);

/* src/core/search.c */
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_lt)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_le)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_gt)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_sequential_ge)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_lt)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_le)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_gt)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_ge)(pepcparts_elements_t *s, pepcparts_slkey_t *k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_lt2)(pepcparts_elements_t *s, pepcparts_slkey_pure_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_le2)(pepcparts_elements_t *s, pepcparts_slkey_pure_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_gt2)(pepcparts_elements_t *s, pepcparts_slkey_pure_t k);
pepcparts_slint SL_PROTO(pepcparts_sl_search_binary_ge2)(pepcparts_elements_t *s, pepcparts_slkey_pure_t k);
pepcparts_slint_t SL_PROTO(pepcparts_sl_search_binary_lt_bmask)(pepcparts_elements_t *s, pepcparts_slkey_pure_t k, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_sl_search_binary_le_bmask)(pepcparts_elements_t *s, pepcparts_slkey_pure_t k, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_sl_search_binary_sign_switch)(pepcparts_elements_t *s);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_lt)(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_le)(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_gt)(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);
pepcparts_slint SL_PROTO(pepcparts_sl_search_hybrid_ge)(pepcparts_elements_t *s, pepcparts_slkey_t *k, pepcparts_slint t);

/* src/core/sl_common.c */
pepcparts_slint SL_PROTO(pepcparts_ilog2c)(pepcparts_slint x);
pepcparts_slint SL_PROTO(pepcparts_ilog2f)(pepcparts_slint x);
long long SL_PROTO(pepcparts_sl_random64)();
pepcparts_slint SL_PROTO(pepcparts_print_bits)(pepcparts_slint v);
pepcparts_slint SL_PROTO(pepcparts_pivot_random)(pepcparts_elements_t *s);
pepcparts_slint_t SL_PROTO(pepcparts_counts2displs)(pepcparts_slint_t n, int *counts, int *displs);
pepcparts_slint_t SL_PROTO(pepcparts_displs2counts)(pepcparts_slint_t n, int *displs, int *counts, pepcparts_slint_t total_counts);

/* src/core/sl_elem.c */
pepcparts_slint_t SL_PROTO(pepcparts_elem_set_data)(pepcparts_elements_t *e, ...);
pepcparts_slint_t SL_PROTO(pepcparts_elem_reverse)(pepcparts_elements_t *e, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_nxchange_at)(pepcparts_elements_t *e0, pepcparts_slint_t at0, pepcparts_elements_t *e1, pepcparts_slint_t at1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_nxchange)(pepcparts_elements_t *e0, pepcparts_elements_t *e1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_nxchange_ro0)(pepcparts_elements_t *e0, pepcparts_elements_t *e1, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_rotate)(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_rotate_ro0)(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);
pepcparts_slint_t SL_PROTO(pepcparts_elem_rotate_ro1)(pepcparts_elements_t *e, pepcparts_slint_t m, pepcparts_slint_t n, pepcparts_elements_t *t);

/* src/core/pepcparts_sort_counting.c */
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting_use_displs)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ndispls, pepcparts_slint_t *displs);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting_use_counts)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts, pepcparts_slint_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting_get_counts)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts, pepcparts_slint_t *counts);
pepcparts_slint_t SL_PROTO(pepcparts_sort_counting)(pepcparts_elements_t *s, pepcparts_elements_t *d, pepcparts_slint_t ncounts);

/* src/core/pepcparts_sort_heap.c */
pepcparts_slint SL_PROTO(pepcparts_sort_heap)(pepcparts_elements_t *s, pepcparts_elements_t *xs);

/* src/core/pepcparts_sort_insert.c */
pepcparts_slint_t SL_PROTO(pepcparts_sort_insert_bmask_kernel)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slkey_pure_t bmask);
pepcparts_slint_t SL_PROTO(pepcparts_sort_insert)(pepcparts_elements_t *s, pepcparts_elements_t *sx);

/* src/core/sort_permute.c */
pepcparts_slint SL_PROTO(pepcparts_sort_permute_forward)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint *perm, pepcparts_slint offset, pepcparts_slint mask_bit);
pepcparts_slint SL_PROTO(pepcparts_sort_permute_backward)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint *perm, pepcparts_slint offset, pepcparts_slint mask_bit);

/* src/core/pepcparts_sort_quick.c */
pepcparts_slint SL_PROTO(pepcparts_sort_quick)(pepcparts_elements_t *s, pepcparts_elements_t *xs);

/* src/core/pepcparts_sort_radix.c */
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_ip)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_db)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_ma)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth);

/* src/core/pepcparts_sort_radix_1bit.c */
pepcparts_slint_t SL_PROTO(pepcparts_sort_radix_1bit_kernel)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow);
pepcparts_slint SL_PROTO(pepcparts_sort_radix_1bit)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow);

/* src/core/pepcparts_sort_radix_af.c */
pepcparts_slint SL_PROTO(pepcparts_sort_radix_af)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint rhigh, pepcparts_slint rlow, pepcparts_slint rwidth);

/* src/core/pepcparts_sort_radix_iter.c */
pepcparts_slint SL_PROTO(pepcparts_sort_radix_iter)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint presorted, pepcparts_slint rhigh, pepcparts_slint rlow, pepcparts_slint rwidth);

/* src/core/sortnet.c */
pepcparts_slint SL_PROTO(pepcparts_sn_hypercube_lh)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_hypercube_hl)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_odd_even_trans)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_odd)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_even)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_batcher)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_bitonic)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);
pepcparts_slint SL_PROTO(pepcparts_sn_connected)(pepcparts_slint size, pepcparts_slint rank, pepcparts_slint stage, void *snp, pepcparts_slint *up);

/* src/core/splitter.c */
pepcparts_slint_t SL_PROTO(pepcparts_splitter_reset)(pepcparts_splitter_t *sp);

/* src/core/splitx.c */
pepcparts_slint_t SL_PROTO(pepcparts_splitx_radix)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint_t nclasses, pepcparts_slint_t shl, pepcparts_slint_t *counts);
pepcparts_slint SL_PROTO(pepcparts_split2_lt_ge)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_split2_le_gt)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t);
pepcparts_slint SL_PROTO(pepcparts_split3_lt_eq_gt)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t, pepcparts_slint *nlt, pepcparts_slint *nle);
pepcparts_slint SL_PROTO(pepcparts_split3_lt_eq_gt_old)(pepcparts_elements_t *s, pepcparts_slkey_pure_t *k, pepcparts_elements_t *t, pepcparts_slint *nlt, pepcparts_slint *nle);
pepcparts_slint SL_PROTO(pepcparts_split2_b)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slkey_pure_t bmask);
pepcparts_slint SL_PROTO(pepcparts_splitk_k2c_af)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);
pepcparts_slint SL_PROTO(pepcparts_splitk_k2c)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);
pepcparts_slint SL_PROTO(pepcparts_splitk_k2c_count)(pepcparts_elements_t *s, pepcparts_slint k, pepcparts_slint *c, pepcparts_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI

/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  timestamp: 2011-03-11 09:08:28 +0100
 *  
 */




/* src/core_mpi/mpi_binning.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_create)(pepcparts_global_bins_t *gb, pepcparts_slint_t max_nbins, pepcparts_slint_t max_nbinnings, pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t doweights, pepcparts_binning_t *bm, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_destroy)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_pre)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_exec_reset)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_exec_local)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_exec_global)(pepcparts_global_bins_t *gb, pepcparts_slint_t root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_refine)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_hit)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slint_t k, pepcparts_splitter_t *sp, pepcparts_slint_t s, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_finalize)(pepcparts_global_bins_t *gb, pepcparts_slint_t b, pepcparts_slweight_t dcw, pepcparts_slint_t lc_min, pepcparts_slint_t lc_max, pepcparts_slweight_t *lcw, pepcparts_splitter_t *sp, pepcparts_slint_t s, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_binning_post)(pepcparts_global_bins_t *gb, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_common.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_datatypes_init)();
pepcparts_slint_t SL_PROTO(pepcparts_mpi_datatypes_release)();
pepcparts_slint_t SL_PROTO(pepcparts_mpi_get_grid_properties)(pepcparts_slint_t ndims, pepcparts_slint_t *dims, pepcparts_slint_t *pos, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_get_grid)(pepcparts_slint_t ndims, pepcparts_slint_t *dims, pepcparts_slint_t *pos, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_subgroups_create)(pepcparts_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_subgroups_delete)(pepcparts_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements.c */
pepcparts_slint SL_PROTO(pepcparts_mpi_elements_init_keys_from_file)(pepcparts_elements_t *s, char *filename, pepcparts_slint from, pepcparts_slint to, pepcparts_slint const_bytes_per_line, pepcparts_slint root, int size, int rank, MPI_Comm comm);
pepcparts_slint SL_PROTO(pepcparts_mpi_elements_validate_order)(pepcparts_elements_t *s, pepcparts_slint n, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_linear_exchange_pure_keys)(pepcparts_slkey_pure_t *in, pepcparts_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_check_order)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t *orders, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_check_global_order)(pepcparts_slkey_pure_t local_min, pepcparts_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_digest_sum)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(pepcparts_mpi_cs32)(pepcparts_elements_t *s, pepcparts_slint n, pepcparts_slint keys, pepcparts_slint data, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_digest_hash)(pepcparts_elements_t *s, pepcparts_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_get_counts)(pepcparts_elements_t *s, pepcparts_slint_t *clocal, pepcparts_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
pepcparts_slweight_t SL_PROTO(pepcparts_mpi_elements_get_weights)(pepcparts_elements_t *s, pepcparts_slweight_t *wlocal, pepcparts_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_get_counts_and_weights)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t *counts, pepcparts_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_sendrecv_replace)(pepcparts_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(pepcparts_mpi_elements_crc32)(pepcparts_elements_t *s, pepcparts_slint_t n, pepcparts_slint_t keys, pepcparts_slint_t data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_alltoallv.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_db)(pepcparts_elements_t *sbuf, int *scounts, int *sdispls, pepcparts_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_alltoallv_ip)(pepcparts_elements_t *sbuf, pepcparts_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, pepcparts_slint_t structured);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);

/* src/core_mpi/pepcparts_mpi_find_exact.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_find_exact_equal)(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *ex_start, pepcparts_slint_t *ex_size, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_find_exact)(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *dst_size, pepcparts_slint_t *ex_start, pepcparts_slint_t *ex_sizes, pepcparts_slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_linsplit.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_linsplit)(MPI_Comm comm_in, pepcparts_slkey_pure_t *keys_in, MPI_Comm *comms_out, pepcparts_slint_t *parity, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_linsplit2)(MPI_Comm comm_in, pepcparts_slkey_pure_t *keys_in, MPI_Comm *comms_out, pepcparts_slint_t *parity, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_merge2.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_merge2)(pepcparts_elements_t *s, pepcparts_slint_t other_rank, pepcparts_slint_t high_rank, pepcparts_slint_t *dst_size, pepcparts_merge2x_f m2, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_mergek.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek_equal)(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek)(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_mergek_equal2)(pepcparts_elements_t *s, pepcparts_sortnet_f sn, pepcparts_sortnet_data_t snd, pepcparts_merge2x_f m2x, pepcparts_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);

/* src/core_mpi/pepcparts_mpi_partition_exact_generic.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_generic)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_generic2)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_partition_exact_radix.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_radix)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_exact_radix_grouped.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_radix_ngroups)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, pepcparts_slint_t ngroups, MPI_Comm *group_comms, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_exact_radix_2groups)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, MPI_Comm group_comm, pepcparts_elements_t *sx, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_sample.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_partition_sample_regular)(pepcparts_elements_t *s, pepcparts_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_rebalance.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_rebalance)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_slint_t stable, pepcparts_slint_t *dst_size, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_rebalance_alltoallv)(pepcparts_elements_t *sbuf, int *scounts, int *sdispls, pepcparts_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_common.c */
pepcparts_slint_t SL_PROTO(pepcparts_init_partconds)(pepcparts_slint_t npconds, pepcparts_partcond_t *pconds, pepcparts_slint_t nparts, pepcparts_slint_t total_count, pepcparts_slweight_t total_weight);
pepcparts_slint_t SL_PROTO(pepcparts_init_partconds_intern)(pepcparts_slint_t npconds, pepcparts_partcond_intern_t *pci, pepcparts_partcond_t *pc, pepcparts_slint_t nparts, pepcparts_slint_t total_count, pepcparts_slweight_t total_weight);
pepcparts_slint_t SL_PROTO(pepcparts_merge_partconds)(pepcparts_partcond_t *pconds_in, pepcparts_slint_t npconds_in, pepcparts_partcond_t *pcond_out);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_gather_partconds_grouped)(pepcparts_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, pepcparts_partcond_t *pconds_out, pepcparts_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_gather_partconds)(pepcparts_partcond_t *pcond_in, pepcparts_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_allgather_partconds)(pepcparts_partcond_t *pcond_in, pepcparts_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_bcast_partconds)(pepcparts_slint_t npconds, pepcparts_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_post_check_partconds)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_post_check_partconds_intern)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_stats)(pepcparts_elements_t *s, pepcparts_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_select_exact_generic.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_generic_bulk)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_binning_t *bm, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_generic_grouped)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepcparts_binning_t *bm, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_generic)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_binning_t *bm, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_select_exact_radix.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_radix)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_exact_radix_grouped)(pepcparts_elements_t *s, pepcparts_slint_t nelements, pepcparts_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, pepcparts_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_sample.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_select_sample_regular)(pepcparts_elements_t *s, pepcparts_slint_t nparts, pepcparts_partcond_t *pconds, pepcparts_slint_t nsamples, pepcparts_splitter_t *sp, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sort_merge.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_merge)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_merge2)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t merge_type, pepcparts_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepcparts_mpi_sort_partition.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t part_type, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_radix)(pepcparts_elements_t *s0, pepcparts_elements_t *s1, pepcparts_elements_t *xs, pepcparts_slint_t part_type, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_exact_radix)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_partcond_t *pcond, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_exact_radix_ngroups)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_partcond_t *pcond, pepcparts_slint_t ngroups, MPI_Comm *group_comms, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_sort_partition_exact_radix_2groups)(pepcparts_elements_t *s, pepcparts_elements_t *sx, pepcparts_partcond_t *pcond, MPI_Comm group_comm, pepcparts_slint_t rhigh, pepcparts_slint_t rlow, pepcparts_slint_t rwidth, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_xcounts2ycounts.c */
pepcparts_slint_t SL_PROTO(pepcparts_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, pepcparts_slint_t ytotal, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_xcounts2ycounts_grouped)(int *xcounts, pepcparts_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
pepcparts_slint_t SL_PROTO(pepcparts_mpi_subxdispls2ycounts)(pepcparts_slint_t nsubs, int *sub_xdispls, pepcparts_slint_t *sub_sources, pepcparts_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_PEPCPARTS_H__ */
