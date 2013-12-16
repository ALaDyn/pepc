
#ifndef __SL_PEPCKEYS_H__
#define __SL_PEPCKEYS_H__

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
#define pepckeys_sl_int_type_c          long long
#define pepckeys_sl_int_type_mpi        MPI_LONG_LONG
#define pepckeys_sl_int_size_mpi        1
#define pepckeys_sl_int_type_fmt        "lld"


/* index data type, corresponds to kind_default */
#define pepckeys_sl_index_type_c        FINT_TYPE_C
#define pepckeys_sl_index_type_mpi      FINT_TYPE_MPI
#define pepckeys_sl_index_size_mpi      1
#define pepckeys_sl_index_type_fmt      FINT_TYPE_FMT

/* use indices */
#define pepckeys_SL_INDEX

/* global key number data type, corresponds to kind_particle */
#define pepckeys_sl_globalcount_type_c        FINT8_TYPE_C
#define pepckeys_sl_globalcount_type_mpi      FINT8_TYPE_MPI
#define pepckeys_sl_globalcount_size_mpi      1
#define pepckeys_sl_globalcount_type_fmt      FINT8_TYPE_FMT

/* keys, corresponds to kind_key */
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

/* weighted elements */
#define pepckeys_sl_elem_weight(e, at)  ((e)->data0[at])

#define pepckeys_sl_data0_weight
/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_config_intern.h
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


#ifndef pepckeys_SL_INDEX
# undef pepckeys_SL_PACKED_INDEX
#endif


/* if no special datatype for (sl default) integer ... */
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


/* if no special datatype for (intern) weight ... */
#ifndef pepckeys_sl_weight_type_c
 /* ... use (sl default) integer */
# define pepckeys_sl_weight_type_c             pepckeys_sl_int_type_c    /* sl_macro */
# undef pepckeys_sl_weight_type_mpi
# define pepckeys_sl_weight_type_mpi           pepckeys_sl_int_type_mpi  /* sl_macro */
# undef pepckeys_sl_weight_size_mpi
# define pepckeys_sl_weight_size_mpi           pepckeys_sl_int_size_mpi  /* sl_macro */
# undef pepckeys_sl_weight_type_fmt
# define pepckeys_sl_weight_type_fmt           pepckeys_sl_int_type_fmt  /* sl_macro */
# undef pepckeys_sl_weight_intequiv
# define pepckeys_sl_weight_intequiv                            /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(pepckeys_sl_weight_type_mpi) || !defined(pepckeys_sl_weight_size_mpi)
#   error "pepckeys_sl_weight_type_mpi and/or pepckeys_sl_weight_size_mpi missing"
#  endif
# endif
# ifndef pepckeys_sl_weight_type_fmt
#  error "pepckeys_sl_weight_type_fmt macro is missing, using f as default"
#  define pepckeys_sl_weight_type_fmt  "f"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef pepckeys_sl_index_type_c
 /* ... use the primary integer type */
# define pepckeys_sl_index_type_c             pepckeys_sl_int_type_c
# undef pepckeys_sl_index_type_mpi
# define pepckeys_sl_index_type_mpi           pepckeys_sl_int_type_mpi
# undef pepckeys_sl_index_size_mpi
# define pepckeys_sl_index_size_mpi           pepckeys_sl_int_size_mpi
# undef pepckeys_sl_index_type_fmt
# define pepckeys_sl_index_type_fmt           pepckeys_sl_int_type_fmt
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
 /* key val -> key val */
 #define pepckeys_sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_get_pure
 /* key component pointer -> key val pointer */
 #define pepckeys_sl_key_get_pure(k)          (k)  /* sl_macro */
#endif
#ifndef pepckeys_sl_key_set_pure
 /* key component pointer and key val */
 #define pepckeys_sl_key_set_pure(k, p)       (*(k) = p)  /* sl_macro */
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


/* default random key */
#ifdef pepckeys_sl_key_integer
# if !defined(pepckeys_sl_key_val_srand) || !defined(pepckeys_sl_key_val_rand) || !defined(pepckeys_sl_key_val_rand_minmax)
#  undef pepckeys_sl_key_val_srand
#  undef pepckeys_sl_key_val_rand
#  undef pepckeys_sl_key_val_rand_minmax
#  define pepckeys_sl_key_val_srand(_s_)                 z_srand(_s_)                                        /* sl_macro */
#  define pepckeys_sl_key_val_rand()                     ((pepckeys_sl_key_pure_type_c) z_rand())                     /* sl_macro */
#  define pepckeys_sl_key_val_rand_minmax(_min_, _max_)  ((pepckeys_sl_key_pure_type_c) z_rand_minmax(_min_, _max_))  /* sl_macro */
# endif
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0_IGNORE
# undef pepckeys_SL_DATA0
#endif
#ifdef pepckeys_SL_DATA1_IGNORE
# undef pepckeys_SL_DATA1
#endif
#ifdef pepckeys_SL_DATA2_IGNORE
# undef pepckeys_SL_DATA2
#endif
#ifdef pepckeys_SL_DATA3_IGNORE
# undef pepckeys_SL_DATA3
#endif
#ifdef pepckeys_SL_DATA4_IGNORE
# undef pepckeys_SL_DATA4
#endif
#ifdef pepckeys_SL_DATA5_IGNORE
# undef pepckeys_SL_DATA5
#endif
#ifdef pepckeys_SL_DATA6_IGNORE
# undef pepckeys_SL_DATA6
#endif
#ifdef pepckeys_SL_DATA7_IGNORE
# undef pepckeys_SL_DATA7
#endif
#ifdef pepckeys_SL_DATA8_IGNORE
# undef pepckeys_SL_DATA8
#endif
#ifdef pepckeys_SL_DATA9_IGNORE
# undef pepckeys_SL_DATA9
#endif
#ifdef pepckeys_SL_DATA10_IGNORE
# undef pepckeys_SL_DATA10
#endif
#ifdef pepckeys_SL_DATA11_IGNORE
# undef pepckeys_SL_DATA11
#endif
#ifdef pepckeys_SL_DATA12_IGNORE
# undef pepckeys_SL_DATA12
#endif
#ifdef pepckeys_SL_DATA13_IGNORE
# undef pepckeys_SL_DATA13
#endif
#ifdef pepckeys_SL_DATA14_IGNORE
# undef pepckeys_SL_DATA14
#endif
#ifdef pepckeys_SL_DATA15_IGNORE
# undef pepckeys_SL_DATA15
#endif
#ifdef pepckeys_SL_DATA16_IGNORE
# undef pepckeys_SL_DATA16
#endif
#ifdef pepckeys_SL_DATA17_IGNORE
# undef pepckeys_SL_DATA17
#endif
#ifdef pepckeys_SL_DATA18_IGNORE
# undef pepckeys_SL_DATA18
#endif
#ifdef pepckeys_SL_DATA19_IGNORE
# undef pepckeys_SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* sl_macro pepckeys_sl_elem_weight */


/* disable sl_dataX_weight if there is not weight */
#ifndef pepckeys_sl_elem_weight
/* DATAX_TEMPLATE_BEGIN */
# undef pepckeys_sl_data0_weight
# undef pepckeys_sl_data1_weight
# undef pepckeys_sl_data2_weight
# undef pepckeys_sl_data3_weight
# undef pepckeys_sl_data4_weight
# undef pepckeys_sl_data5_weight
# undef pepckeys_sl_data6_weight
# undef pepckeys_sl_data7_weight
# undef pepckeys_sl_data8_weight
# undef pepckeys_sl_data9_weight
# undef pepckeys_sl_data10_weight
# undef pepckeys_sl_data11_weight
# undef pepckeys_sl_data12_weight
# undef pepckeys_sl_data13_weight
# undef pepckeys_sl_data14_weight
# undef pepckeys_sl_data15_weight
# undef pepckeys_sl_data16_weight
# undef pepckeys_sl_data17_weight
# undef pepckeys_sl_data18_weight
# undef pepckeys_sl_data19_weight
/* DATAX_TEMPLATE_END */
#endif


/* disable pepckeys_sl_elem_weight if the weight component is missing */
/* DATAX_TEMPLATE_BEGIN */
#if defined(pepckeys_sl_data0_weight) && !defined(pepckeys_SL_DATA0)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data1_weight) && !defined(pepckeys_SL_DATA1)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data2_weight) && !defined(pepckeys_SL_DATA2)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data3_weight) && !defined(pepckeys_SL_DATA3)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data4_weight) && !defined(pepckeys_SL_DATA4)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data5_weight) && !defined(pepckeys_SL_DATA5)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data6_weight) && !defined(pepckeys_SL_DATA6)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data7_weight) && !defined(pepckeys_SL_DATA7)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data8_weight) && !defined(pepckeys_SL_DATA8)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data9_weight) && !defined(pepckeys_SL_DATA9)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data10_weight) && !defined(pepckeys_SL_DATA10)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data11_weight) && !defined(pepckeys_SL_DATA11)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data12_weight) && !defined(pepckeys_SL_DATA12)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data13_weight) && !defined(pepckeys_SL_DATA13)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data14_weight) && !defined(pepckeys_SL_DATA14)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data15_weight) && !defined(pepckeys_SL_DATA15)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data16_weight) && !defined(pepckeys_SL_DATA16)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data17_weight) && !defined(pepckeys_SL_DATA17)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data18_weight) && !defined(pepckeys_SL_DATA18)
# undef pepckeys_sl_elem_weight
#endif
#if defined(pepckeys_sl_data19_weight) && !defined(pepckeys_SL_DATA19)
# undef pepckeys_sl_elem_weight
#endif
/* DATAX_TEMPLATE_END */


/* verify that the flex component is the last (FIXME: only if packed is on?) */
/* sl_macro pepckeys_FLECKS_GUARD */
/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data0_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA1
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data1_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA2
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data2_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA3
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data3_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA4
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data4_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA5
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data5_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA6
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data6_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA7
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data7_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA8
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data8_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA9
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data9_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA10
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data10_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA11
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data11_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA12
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data12_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA13
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data13_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA14
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data14_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA15
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data15_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA16
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data16_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA17
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data17_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA18
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data18_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
#ifdef pepckeys_SL_DATA19
# ifdef pepckeys_FLECKS_GUARD
#  error "flexible data component is not the last data component!"
# else
#  ifdef pepckeys_sl_data19_flex
#   define pepckeys_FLECKS_GUARD
#  endif
# endif
#endif
/* DATAX_TEMPLATE_END */


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_types.h
 *  
 */




/* sl_type pepckeys_slint_t pepckeys_slint */
typedef pepckeys_sl_int_type_c pepckeys_slint_t, pepckeys_slint;  /* deprecated 'pepckeys_slint' */

#define pepckeys_slint_fmt   pepckeys_sl_int_type_fmt    /* sl_macro */

/* sl_type pepckeys_slindex_t */
typedef pepckeys_sl_index_type_c pepckeys_slindex_t;

#define pepckeys_sindex_fmt  pepckeys_sl_index_type_fmt  /* sl_macro */

/* sl_type pepckeys_slkey_t */
typedef pepckeys_sl_key_type_c pepckeys_slkey_t;

/* sl_type pepckeys_sl_globalcount */
typedef pepckeys_sl_globalcount_type_c pepckeys_sl_globalcount_t;

/* sl_type pepckeys_slkey_pure_t pepckeys_slpkey_t */
typedef pepckeys_sl_key_pure_type_c pepckeys_slkey_pure_t, pepckeys_slpkey_t;

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

/* sl_type pepckeys_slweight_t */
typedef pepckeys_sl_weight_type_c pepckeys_slweight_t;

#define pepckeys_slweight_fmt  pepckeys_sl_weight_type_fmt  /* sl_macro */

#if defined(pepckeys_sl_elem_weight) && defined(pepckeys_sl_weight_intequiv)
typedef pepckeys_sl_weight_type_c pepckeys_slcount_t;       /* sl_type pepckeys_slcount_t */
# define pepckeys_slcount_fmt  pepckeys_sl_weight_type_fmt  /* sl_macro */
#else
typedef pepckeys_sl_int_type_c pepckeys_slcount_t;
# define pepckeys_slcount_fmt  pepckeys_sl_int_type_fmt
#endif


/* sl_type pepckeys__slpwkey_t pepckeys_slpwkey_t */
typedef struct pepckeys__slpwkey_t
{
  pepckeys_slpkey_t pkey;
  pepckeys_slweight_t weight;

} pepckeys_slpwkey_t;


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

#ifdef pepckeys_SL_PACKED_INDEX
  pepckeys_slindex_t index;
#endif

/* DATAX_TEMPLATE_BEGIN */
#ifdef pepckeys_SL_DATA0
# ifdef pepckeys_sl_data0_flex
  pepckeys_sldata0_t data0[];
# else
  pepckeys_sldata0_t data0[pepckeys_sl_data0_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA1
# ifdef pepckeys_sl_data1_flex
  pepckeys_sldata1_t data1[];
# else
  pepckeys_sldata1_t data1[pepckeys_sl_data1_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA2
# ifdef pepckeys_sl_data2_flex
  pepckeys_sldata2_t data2[];
# else
  pepckeys_sldata2_t data2[pepckeys_sl_data2_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA3
# ifdef pepckeys_sl_data3_flex
  pepckeys_sldata3_t data3[];
# else
  pepckeys_sldata3_t data3[pepckeys_sl_data3_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA4
# ifdef pepckeys_sl_data4_flex
  pepckeys_sldata4_t data4[];
# else
  pepckeys_sldata4_t data4[pepckeys_sl_data4_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA5
# ifdef pepckeys_sl_data5_flex
  pepckeys_sldata5_t data5[];
# else
  pepckeys_sldata5_t data5[pepckeys_sl_data5_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA6
# ifdef pepckeys_sl_data6_flex
  pepckeys_sldata6_t data6[];
# else
  pepckeys_sldata6_t data6[pepckeys_sl_data6_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA7
# ifdef pepckeys_sl_data7_flex
  pepckeys_sldata7_t data7[];
# else
  pepckeys_sldata7_t data7[pepckeys_sl_data7_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA8
# ifdef pepckeys_sl_data8_flex
  pepckeys_sldata8_t data8[];
# else
  pepckeys_sldata8_t data8[pepckeys_sl_data8_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA9
# ifdef pepckeys_sl_data9_flex
  pepckeys_sldata9_t data9[];
# else
  pepckeys_sldata9_t data9[pepckeys_sl_data9_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA10
# ifdef pepckeys_sl_data10_flex
  pepckeys_sldata10_t data10[];
# else
  pepckeys_sldata10_t data10[pepckeys_sl_data10_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA11
# ifdef pepckeys_sl_data11_flex
  pepckeys_sldata11_t data11[];
# else
  pepckeys_sldata11_t data11[pepckeys_sl_data11_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA12
# ifdef pepckeys_sl_data12_flex
  pepckeys_sldata12_t data12[];
# else
  pepckeys_sldata12_t data12[pepckeys_sl_data12_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA13
# ifdef pepckeys_sl_data13_flex
  pepckeys_sldata13_t data13[];
# else
  pepckeys_sldata13_t data13[pepckeys_sl_data13_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA14
# ifdef pepckeys_sl_data14_flex
  pepckeys_sldata14_t data14[];
# else
  pepckeys_sldata14_t data14[pepckeys_sl_data14_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA15
# ifdef pepckeys_sl_data15_flex
  pepckeys_sldata15_t data15[];
# else
  pepckeys_sldata15_t data15[pepckeys_sl_data15_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA16
# ifdef pepckeys_sl_data16_flex
  pepckeys_sldata16_t data16[];
# else
  pepckeys_sldata16_t data16[pepckeys_sl_data16_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA17
# ifdef pepckeys_sl_data17_flex
  pepckeys_sldata17_t data17[];
# else
  pepckeys_sldata17_t data17[pepckeys_sl_data17_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA18
# ifdef pepckeys_sl_data18_flex
  pepckeys_sldata18_t data18[];
# else
  pepckeys_sldata18_t data18[pepckeys_sl_data18_size_c];
# endif
#endif
#ifdef pepckeys_SL_DATA19
# ifdef pepckeys_sl_data19_flex
  pepckeys_sldata19_t data19[];
# else
  pepckeys_sldata19_t data19[pepckeys_sl_data19_size_c];
# endif
#endif
/* DATAX_TEMPLATE_END */

} pepckeys_packed_element_t;


/* sl_type pepckeys__packed_elements_t pepckeys_packed_elements_t */
typedef struct pepckeys__packed_elements_t
{
  pepckeys_slint_t size, max_size;
  
  pepckeys_packed_element_t *elements;
  
} pepckeys_packed_elements_t;


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

/* sl_type pepckeys_tproc_f pepckeys__tproc_data_t pepckeys_tproc_data_t */
typedef struct pepckeys__tproc_data_t
{
  int *int_tprocs;
  int int_mask;

} pepckeys_tproc_data_t;

typedef int (*pepckeys_tproc_f)(pepckeys_elements_t *s0, pepckeys_slint_t x, void *data);

#ifndef SL_TPROC_INT
# define SL_TPROC_INT       ((void *) 1)
#endif
#ifndef SL_TPROC_INT_MASK
# define SL_TPROC_INT_MASK  ((void *) 2)
#endif


/* deprecated, sl_type pepckeys_k2c_func pepckeys_pivot_func pepckeys_sn_func pepckeys_m2x_func pepckeys_m2X_func */
typedef pepckeys_key2class_f pepckeys_k2c_func;
typedef pepckeys_pivot_f pepckeys_pivot_func;
typedef pepckeys_sortnet_f pepckeys_sn_func;
typedef pepckeys_merge2x_f pepckeys_m2x_func;
typedef pepckeys_merge2X_f pepckeys_m2X_func;


/* sl_type pepckeys__mergek_t pepckeys_mergek_t */
typedef struct pepckeys__mergek_t
{
  pepckeys_sortnet_f sn;
  pepckeys_sortnet_data_t snd;

  pepckeys_merge2x_f m2x;
  pepckeys_elements_t *sx;

} pepckeys_mergek_t;


/* sl_type pepckeys_keys_init_type_t pepckeys_keys_init_data_t */
typedef pepckeys_slint_t pepckeys_keys_init_type_t;
typedef void *pepckeys_keys_init_data_t;

/* sl_type pepckeys_key_set_data_t pepckeys_key_set_f */
typedef void *pepckeys_key_set_data_t;
typedef void (*pepckeys_key_set_f)(pepckeys_slkey_pure_t *k, pepckeys_key_set_data_t d);


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
#undef SL_EKIT_URAND
#define SL_EKIT_URAND       6
#undef SL_EKIT_NRAND
#define SL_EKIT_NRAND       7


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


/* pepckeys_elements_keys_stats */
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


/* partition conditions, sl_type pepckeys__partcond2_t pepckeys_partcond2_t */
typedef struct pepckeys__partcond2_t
{
  int weighted;
  double min_count, max_count;
  double min_weight, max_weight;
  double min_cpart, max_cpart;
  double min_wpart, max_wpart;

} pepckeys_partcond2_t;


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

/* partition conditions, sl_type pepckeys__partcond_t pepckeys_partcond_t pepckeys_partcond_p */
typedef struct pepckeys__partcond_t
{
  pepckeys_slint_t pcm;
  double count_min, count_max;
  double count_low, count_high;
  double weight_min, weight_max;
  double weight_low, weight_high;

} pepckeys_partcond_t, *pepckeys_partcond_p;


/* internal partition conditions, sl_type pepckeys__partcond_intern_t pepckeys_partcond_intern_t pepckeys_partcond_intern_p */
typedef struct pepckeys__partcond_intern_t
{
  pepckeys_slint_t pcm;
  pepckeys_slint_t count_min, count_max;
  pepckeys_slint_t count_low, count_high;
#ifdef elem_weight
  pepckeys_slweight_t weight_min, weight_max;
  pepckeys_slweight_t weight_low, weight_high;
#endif

} pepckeys_partcond_intern_t, *pepckeys_partcond_intern_p;


/* sl_type pepckeys__parttype_t pepckeys_parttype_t pepckeys_parttype_p */
typedef struct pepckeys__parttype_t
{
  pepckeys_slint_t type;

} pepckeys_parttype_t, *pepckeys_parttype_p;


/* generic binning method */

/* sl_type pepckeys__bin_t pepckeys_bin_t */
typedef struct pepckeys__bin_t
{
  pepckeys_elements_t s;

#ifdef elem_weight
  pepckeys_slweight_t weight;
#endif

} pepckeys_bin_t;


/* sl_type pepckeys__splitter_t pepckeys_splitter_t */
typedef struct pepckeys__splitter_t
{
  pepckeys_slint_t n;

  int *displs;
  pepckeys_slkey_pure_t *s;
  pepckeys_slint_t *sn;

} pepckeys_splitter_t;


struct pepckeys__binning_t;

/* sl_type pepckeys_binning_pre_f pepckeys_binning_exec_f pepckeys_binning_refine_f pepckeys_binning_hit_f pepckeys_binning_finalize_f pepckeys_binning_post_f */
typedef pepckeys_slint_t (*pepckeys_binning_pre_f)(struct pepckeys__binning_t *bm);
typedef pepckeys_slint_t (*pepckeys_binning_exec_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights);
typedef pepckeys_slint_t (*pepckeys_binning_refine_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights, pepckeys_splitter_t *sp, pepckeys_slint_t s, pepckeys_bin_t *new_bin);
typedef pepckeys_slint_t (*pepckeys_binning_hit_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_splitter_t *sp, pepckeys_slint_t s);
typedef pepckeys_slint_t (*pepckeys_binning_finalize_f)(struct pepckeys__binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s);
typedef pepckeys_slint_t (*pepckeys_binning_post_f)(struct pepckeys__binning_t *bm);


/* sl_type pepckeys__binning_data_t pepckeys_binning_data_t */
typedef union pepckeys__binning_data_t
{
  struct
  {
    pepckeys_slint_t rhigh, rlow, rwidth;
    pepckeys_slint_t rcurrent;
    pepckeys_slkey_pure_t bit_mask;

    pepckeys_elements_t sx;

  } radix;

} pepckeys_binning_data_t;


/* sl_type pepckeys__binning_t pepckeys_binning_t */
typedef struct pepckeys__binning_t
{
  pepckeys_slint_t nbins, max_nbins;
  
  pepckeys_binning_pre_f pre;
  pepckeys_binning_exec_f exec;
  pepckeys_binning_refine_f refine;
  pepckeys_binning_hit_f hit;
  pepckeys_binning_finalize_f finalize;
  pepckeys_binning_post_f post;

  pepckeys_slint_t sorted;

  pepckeys_slint_t docounts;
#ifdef elem_weight
  pepckeys_slint_t doweights;
#endif

  pepckeys_binning_data_t bd;

} pepckeys_binning_t;


/* sl_type pepckeys__local_bins_t pepckeys_local_bins_t */
typedef struct pepckeys__local_bins_t
{
  pepckeys_binning_t *bm;

  pepckeys_slint_t nbins, max_nbins;
  pepckeys_slint_t nelements;

  pepckeys_slint_t docounts;
#ifdef elem_weight
  pepckeys_slint_t doweights;
#endif

  pepckeys_slint_t nbinnings, max_nbinnings;

  pepckeys_slint_t nbins_new, last_new_b, last_new_k;
  pepckeys_bin_t *bins, *bins_new;
  pepckeys_bin_t *bins0, *bins1;

  pepckeys_slint_t *bcws;

#if defined(elem_weight) && defined(pepckeys_sl_weight_intequiv)
  pepckeys_slint_t cw_factor, w_index, bin_cw_factor;
  pepckeys_slweight_t *cws, *bin_cws;
#else
  pepckeys_slint_t *cs, *bin_cs;
# ifdef elem_weight
  pepckeys_slweight_t *ws, *bin_ws;
# endif
#endif

  pepckeys_slint_t last_exec_b;

} pepckeys_local_bins_t;


/* sl_type pepckeys__global_bins_t pepckeys_global_bins_t */
typedef struct pepckeys__global_bins_t
{
  pepckeys_binning_t *bm;
  
  pepckeys_local_bins_t lb;

  pepckeys_slint_t *bcws;

#if defined(elem_weight) && defined(pepckeys_sl_weight_intequiv)
  pepckeys_slweight_t *cws;
#else
  pepckeys_slint_t *cs;
# ifdef elem_weight
  pepckeys_slweight_t *ws;
# endif
#endif

} pepckeys_global_bins_t;


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_adds.h
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


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_globals.h
 *  
 */




/* src/core/sl_common.c */
extern rti pepckeys_rti_env;
extern int pepckeys_sl_mpi_rank_dummy;

/* src/core/pepckeys_sort_radix.c */
extern pepckeys_slint_t pepckeys_sr_ip_threshold;
extern pepckeys_slint_t pepckeys_sr_db_threshold;
extern pepckeys_slint_t pepckeys_sr_ma_threshold;

/* src/core_mpi/mpi_common.c */
#ifdef SL_USE_MPI
extern MPI_Datatype pepckeys_int_mpi_datatype;
extern MPI_Datatype pepckeys_key_mpi_datatype;
extern MPI_Datatype pepckeys_pkey_mpi_datatype;
extern MPI_Datatype pepckeys_pwkey_mpi_datatype;
extern MPI_Datatype pepckeys_index_mpi_datatype;
extern MPI_Datatype pepckeys_weight_mpi_datatype;
extern MPI_Datatype pepckeys_data_mpi_datatype[];
#endif
#ifdef SL_USE_MPI
extern int pepckeys_sl_mpi_rank;
#endif

/* src/core_mpi/mpi_elements.c */
extern void *pepckeys_me_sendrecv_replace_mem;
extern pepckeys_slint_t pepckeys_me_sendrecv_replace_memsize;
extern pepckeys_slint_t pepckeys_me_sendrecv_replace_mpi_maxsize;

/* src/core_mpi/pepckeys_mpi_elements_alltoall_specific.c */
extern double pepckeys_meas_t[];

/* src/core_mpi/pepckeys_mpi_select_exact_generic.c */
extern int pepckeys_mseg_root;
extern double pepckeys_mseg_border_update_count_reduction;
extern double pepckeys_mseg_border_update_weight_reduction;
extern pepckeys_slint_t pepckeys_mseg_forward_only;
extern pepckeys_slint_t pepckeys_mseg_info_rounds;
extern pepckeys_slint_t *pepckeys_mseg_info_finish_rounds;
extern double pepckeys_mseg_info_finish_rounds_avg;
extern pepckeys_slint_t pepckeys_mseg_binnings;
extern pepckeys_slint_t pepckeys_mseg_finalize_mode;

/* src/core_mpi/mpi_select_sample.c */
extern int pepckeys_mss_root;

/* src/core_mpi/pepckeys_mpi_sort_merge.c */
extern double pepckeys_msm_t[];
extern pepckeys_slint_t pepckeys_msm_sync;

/* src/core_mpi/pepckeys_mpi_sort_partition.c */
extern double pepckeys_msp_t[];
extern pepckeys_slint_t pepckeys_msp_sync;
extern pepckeys_partcond_t *pepckeys_msp_r_pc;

/* src/core_mpi/mpi_sort_special.c */
extern double pepckeys_mss_i_t[];
extern double pepckeys_mss_p_t[];
extern double pepckeys_mss_b_t[];
extern pepckeys_slint_t pepckeys_mss_sync;
extern pepckeys_slint_t pepckeys_mss_i_sync;
extern pepckeys_slint_t pepckeys_mss_p_sync;
extern pepckeys_slint_t pepckeys_mss_b_sync;


/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos.h
 *  
 */




/* src/core/binning.c */
pepckeys_slint_t SL_PROTO(pepckeys_binning_create)(pepckeys_local_bins_t *lb, pepckeys_slint_t max_nbins, pepckeys_slint_t max_nbinnings, pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t docounts, pepckeys_slint_t doweights, pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_binning_destroy)(pepckeys_local_bins_t *lb);
pepckeys_slint_t SL_PROTO(pepckeys_binning_pre)(pepckeys_local_bins_t *lb);
pepckeys_slint_t SL_PROTO(pepckeys_binning_exec_reset)(pepckeys_local_bins_t *lb);
pepckeys_slint_t SL_PROTO(pepckeys_binning_exec)(pepckeys_local_bins_t *lb, pepckeys_slint_t b);
pepckeys_slint_t SL_PROTO(pepckeys_binning_refine)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_hit)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_finalize)(pepckeys_local_bins_t *lb, pepckeys_slint_t b, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_post)(pepckeys_local_bins_t *lb);

/* src/core/binning_radix.c */
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_create)(pepckeys_binning_t *bm, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_destroy)(pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_pre)(pepckeys_binning_t *bm);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_exec)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_refine)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_slweight_t *weights, pepckeys_splitter_t *sp, pepckeys_slint_t s, pepckeys_bin_t *new_bin);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_hit)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t k, pepckeys_slcount_t *counts, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_finalize)(pepckeys_binning_t *bm, pepckeys_bin_t *bin, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s);
pepckeys_slint_t SL_PROTO(pepckeys_binning_radix_post)(pepckeys_binning_t *bm);

/* src/core/elements.c */
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_free)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc2)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t keys, pepckeys_slint_t indices, pepckeys_slint_t data, pepckeys_slint_t weights);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_old)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t keys, pepckeys_slint_t data);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_from_blocks)(pepckeys_elements_t *s, pepckeys_slint_t nblocks, void **blocks, pepckeys_slint_t *blocksizes, pepckeys_slint_t alignment, pepckeys_slint_t nmax, slcint_t components);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_from_block2)(pepckeys_elements_t *s, void *block, pepckeys_slint_t blocksize, pepckeys_slint_t alignment, pepckeys_slint_t nmax, pepckeys_slint_t keys, pepckeys_slint_t indices, pepckeys_slint_t data, pepckeys_slint_t weights);
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_from_block)(pepckeys_elements_t *s, void *block, pepckeys_slint_t blocksize, pepckeys_slint_t alignment, pepckeys_slint_t nmax);
pepckeys_slint_t SL_PROTO(pepckeys_elements_copy)(pepckeys_elements_t *s, pepckeys_elements_t *d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_copy_at)(pepckeys_elements_t *s, pepckeys_slint_t sat, pepckeys_elements_t *d, pepckeys_slint_t dat);
pepckeys_slint_t SL_PROTO(pepckeys_elements_ncopy)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_elements_nmove)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_elements_printf)(pepckeys_elements_t *s, const char *prefix);
pepckeys_slint_t SL_PROTO(pepckeys_elements_extract)(pepckeys_elements_t *src, pepckeys_slint_t nelements, pepckeys_elements_t *dst0, pepckeys_elements_t *dst1);
pepckeys_slint_t SL_PROTO(pepckeys_elements_touch)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_digest_sum)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, unsigned int *sum);
unsigned int SL_PROTO(pepckeys_elements_crc32)(pepckeys_elements_t *s, pepckeys_slint nelements, pepckeys_slint_t keys, pepckeys_slint_t data);
pepckeys_slint_t SL_PROTO(pepckeys_elements_digest_hash)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, void *hash);
pepckeys_slint_t SL_PROTO(pepckeys_elements_random_exchange)(pepckeys_elements_t *s, pepckeys_slint_t rounds, pepckeys_elements_t *xs);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init_seed)(unsigned long s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init)(pepckeys_elements_t *s, pepckeys_keys_init_type_t t, pepckeys_keys_init_data_t d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_init_randomized)(pepckeys_elements_t *s, pepckeys_slint_t nkeys, pepckeys_keys_init_type_t t, pepckeys_keys_init_data_t d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_init_keys_from_file)(pepckeys_elements_t *s, pepckeys_slint_t data, char *filename, pepckeys_slint_t from, pepckeys_slint_t to, pepckeys_slint_t const_bytes_per_line);
pepckeys_slint_t SL_PROTO(pepckeys_elements_save_keys_to_file)(pepckeys_elements_t *s, char *filename);
pepckeys_slint_t SL_PROTO(pepckeys_elements_validate_order)(pepckeys_elements_t *s, pepckeys_slint_t n);
pepckeys_slint_t SL_PROTO(pepckeys_elements_validate_order_bmask)(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_elements_validate_order_weight)(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slkey_pure_t weight);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_stats)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *stats);
pepckeys_slint_t SL_PROTO(pepckeys_elements_keys_stats_print)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_print_keys)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_print_all)(pepckeys_elements_t *s);
pepckeys_slweight_t SL_PROTO(pepckeys_elements_get_weight)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_get_minmax_keys)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slkey_pure_t *minmaxkeys);

/* src/core/elements_packed.c */
pepckeys_slint_t SL_PROTO(pepckeys_elements_alloc_packed)(pepckeys_packed_elements_t *s, pepckeys_slint_t nelements);
pepckeys_slint_t SL_PROTO(pepckeys_elements_free_packed)(pepckeys_packed_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_elements_pack_indexed)(pepckeys_elements_t *s, pepckeys_packed_elements_t *d, pepckeys_slindex_t *rindx, pepckeys_slindex_t *windx);
pepckeys_slint_t SL_PROTO(pepckeys_elements_pack)(pepckeys_elements_t *s, pepckeys_packed_elements_t *d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack_indexed)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slindex_t *rindx, pepckeys_slindex_t *windx);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d);
pepckeys_slint_t SL_PROTO(pepckeys_elements_unpack_keys)(pepckeys_packed_elements_t *s, pepckeys_slkey_t *k);

/* src/core/merge2_common.c */
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_auto_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_m2x_func _x0_1, pepckeys_m2x_func _0x_1);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_01_X)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t, pepckeys_m2X_func _X0_1, pepckeys_m2X_func _0X_1);
pepckeys_slint SL_PROTO(pepckeys_merge2_simplify_s1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_slint s1elements);
pepckeys_slint SL_PROTO(pepckeys_merge2_memory_adaptive)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);

/* src/core/merge2_hula.c */
pepckeys_slint SL_PROTO(pepckeys_merge2_compo_hula)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs);

/* src/core/merge2_search.c */
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sseq_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sseq_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sseq_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sseq_01)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sbin_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sbin_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sbin_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_sbin_01)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_shyb_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_shyb_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_shyb_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_shyb_01)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *t);

/* src/core/merge2_straight.c */
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_x0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_0x_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_01_x)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_x_0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_X0_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_0X_1)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_01_X)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_merge2_basic_straight_X0_1u)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx, pepckeys_elements_t *t);

/* src/core/merge2_tridgell.c */
pepckeys_slint SL_PROTO(pepckeys_merge2_compo_tridgell)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *sx);

/* src/core/mergep_2way.c */
pepckeys_slint_t SL_PROTO(pepckeys_mergep_2way_ip_int)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t p, int *displs, pepckeys_merge2x_f m2x);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_2way_ip_int_rec)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t p, int *displs, pepckeys_merge2x_f m2x);

/* src/core/mergep_heap.c */
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_int)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, int *displs, int *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_int_idx)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, int *displs, int *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_idx)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_unpack_idx)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_mergep_heap_unpack_idxonly)(pepckeys_packed_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t p, pepckeys_slindex_t *displs, pepckeys_slindex_t *counts);

/* src/core/search.c */
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_lt)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_le)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_gt)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_sequential_ge)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_lt)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_le)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_gt)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_ge)(pepckeys_elements_t *s, pepckeys_slkey_t *k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_lt2)(pepckeys_elements_t *s, pepckeys_slkey_pure_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_le2)(pepckeys_elements_t *s, pepckeys_slkey_pure_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_gt2)(pepckeys_elements_t *s, pepckeys_slkey_pure_t k);
pepckeys_slint SL_PROTO(pepckeys_sl_search_binary_ge2)(pepckeys_elements_t *s, pepckeys_slkey_pure_t k);
pepckeys_slint_t SL_PROTO(pepckeys_sl_search_binary_lt_bmask)(pepckeys_elements_t *s, pepckeys_slkey_pure_t k, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_sl_search_binary_le_bmask)(pepckeys_elements_t *s, pepckeys_slkey_pure_t k, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_sl_search_binary_sign_switch)(pepckeys_elements_t *s);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_lt)(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_le)(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_gt)(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);
pepckeys_slint SL_PROTO(pepckeys_sl_search_hybrid_ge)(pepckeys_elements_t *s, pepckeys_slkey_t *k, pepckeys_slint t);

/* src/core/sl_common.c */
pepckeys_slint SL_PROTO(pepckeys_ilog2c)(pepckeys_slint x);
pepckeys_slint SL_PROTO(pepckeys_ilog2f)(pepckeys_slint x);
pepckeys_slint SL_PROTO(pepckeys_print_bits)(pepckeys_slint v);
pepckeys_slint SL_PROTO(pepckeys_pivot_random)(pepckeys_elements_t *s);
pepckeys_slint_t SL_PROTO(pepckeys_counts2displs)(pepckeys_slint_t n, int *counts, int *displs);
pepckeys_slint_t SL_PROTO(pepckeys_displs2counts)(pepckeys_slint_t n, int *displs, int *counts, pepckeys_slint_t total_counts);

/* src/core/sl_elem.c */
pepckeys_slint_t SL_PROTO(pepckeys_elem_set_data)(pepckeys_elements_t *e, ...);
pepckeys_slint_t SL_PROTO(pepckeys_elem_reverse)(pepckeys_elements_t *e, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_nxchange_at)(pepckeys_elements_t *e0, pepckeys_slint_t at0, pepckeys_elements_t *e1, pepckeys_slint_t at1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_nxchange)(pepckeys_elements_t *e0, pepckeys_elements_t *e1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_nxchange_ro0)(pepckeys_elements_t *e0, pepckeys_elements_t *e1, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_rotate)(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_rotate_ro0)(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);
pepckeys_slint_t SL_PROTO(pepckeys_elem_rotate_ro1)(pepckeys_elements_t *e, pepckeys_slint_t m, pepckeys_slint_t n, pepckeys_elements_t *t);

/* src/core/pepckeys_sort_counting.c */
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting_use_displs)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ndispls, pepckeys_slint_t *displs);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting_use_counts)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts, pepckeys_slint_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting_get_counts)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts, pepckeys_slint_t *counts);
pepckeys_slint_t SL_PROTO(pepckeys_sort_counting)(pepckeys_elements_t *s, pepckeys_elements_t *d, pepckeys_slint_t ncounts);

/* src/core/pepckeys_sort_heap.c */
pepckeys_slint SL_PROTO(pepckeys_sort_heap)(pepckeys_elements_t *s, pepckeys_elements_t *xs);

/* src/core/pepckeys_sort_insert.c */
pepckeys_slint_t SL_PROTO(pepckeys_sort_insert_bmask_kernel)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slkey_pure_t bmask);
pepckeys_slint_t SL_PROTO(pepckeys_sort_insert)(pepckeys_elements_t *s, pepckeys_elements_t *sx);

/* src/core/sort_permute.c */
pepckeys_slint_t SL_PROTO(pepckeys_sort_permute_forward)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t *perm, pepckeys_slint_t offset, pepckeys_slint_t mask_bit);
pepckeys_slint_t SL_PROTO(pepckeys_sort_permute_backward)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t *perm, pepckeys_slint_t offset, pepckeys_slint_t mask_bit);

/* src/core/pepckeys_sort_quick.c */
pepckeys_slint SL_PROTO(pepckeys_sort_quick)(pepckeys_elements_t *s, pepckeys_elements_t *xs);

/* src/core/pepckeys_sort_radix.c */
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_ip)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_db)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_ma)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth);

/* src/core/pepckeys_sort_radix_1bit.c */
pepckeys_slint_t SL_PROTO(pepckeys_sort_radix_1bit_kernel)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow);
pepckeys_slint SL_PROTO(pepckeys_sort_radix_1bit)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow);

/* src/core/pepckeys_sort_radix_af.c */
pepckeys_slint SL_PROTO(pepckeys_sort_radix_af)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint rhigh, pepckeys_slint rlow, pepckeys_slint rwidth);

/* src/core/pepckeys_sort_radix_iter.c */
pepckeys_slint SL_PROTO(pepckeys_sort_radix_iter)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint presorted, pepckeys_slint rhigh, pepckeys_slint rlow, pepckeys_slint rwidth);

/* src/core/sortnet.c */
pepckeys_slint SL_PROTO(pepckeys_sn_hypercube_lh)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_hypercube_hl)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_odd_even_trans)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_odd)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_even)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_batcher)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_bitonic)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);
pepckeys_slint SL_PROTO(pepckeys_sn_connected)(pepckeys_slint size, pepckeys_slint rank, pepckeys_slint stage, void *snp, pepckeys_slint *up);

/* src/core/splitter.c */
pepckeys_slint_t SL_PROTO(pepckeys_splitter_reset)(pepckeys_splitter_t *sp);

/* src/core/splitx.c */
pepckeys_slint_t SL_PROTO(pepckeys_splitx_radix)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint_t nclasses, pepckeys_slint_t shl, pepckeys_slint_t *counts);
pepckeys_slint SL_PROTO(pepckeys_split2_lt_ge)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_split2_le_gt)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t);
pepckeys_slint SL_PROTO(pepckeys_split3_lt_eq_gt)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t, pepckeys_slint *nlt, pepckeys_slint *nle);
pepckeys_slint SL_PROTO(pepckeys_split3_lt_eq_gt_old)(pepckeys_elements_t *s, pepckeys_slkey_pure_t *k, pepckeys_elements_t *t, pepckeys_slint *nlt, pepckeys_slint *nle);
pepckeys_slint SL_PROTO(pepckeys_split2_b)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slkey_pure_t bmask);
pepckeys_slint SL_PROTO(pepckeys_splitk_k2c_af)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);
pepckeys_slint SL_PROTO(pepckeys_splitk_k2c)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);
pepckeys_slint SL_PROTO(pepckeys_splitk_k2c_count)(pepckeys_elements_t *s, pepckeys_slint k, pepckeys_slint *c, pepckeys_k2c_func k2c, void *k2c_data);


#ifdef SL_USE_MPI

/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_protos_mpi.h
 *  
 */




/* src/core_mpi/mpi_binning.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_create)(pepckeys_global_bins_t *gb, pepckeys_slint_t max_nbins, pepckeys_slint_t max_nbinnings, pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t docounts, pepckeys_slint_t doweights, pepckeys_binning_t *bm, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_destroy)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_pre)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_exec_reset)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_exec_local)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_exec_global)(pepckeys_global_bins_t *gb, pepckeys_slint_t root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_refine)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_hit)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t k, pepckeys_splitter_t *sp, pepckeys_slint_t s, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_finalize)(pepckeys_global_bins_t *gb, pepckeys_slint_t b, pepckeys_slint_t dc, pepckeys_slweight_t dw, pepckeys_slint_t lc_min, pepckeys_slint_t lc_max, pepckeys_slcount_t *lcs, pepckeys_slweight_t *lws, pepckeys_splitter_t *sp, pepckeys_slint_t s, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_binning_post)(pepckeys_global_bins_t *gb, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_common.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_datatypes_init)();
pepckeys_slint_t SL_PROTO(pepckeys_mpi_datatypes_release)();
pepckeys_slint_t SL_PROTO(pepckeys_mpi_get_grid_properties)(pepckeys_slint_t ndims, pepckeys_slint_t *dims, pepckeys_slint_t *pos, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_get_grid)(pepckeys_slint_t ndims, pepckeys_slint_t *dims, pepckeys_slint_t *pos, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_subgroups_create)(pepckeys_slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_subgroups_delete)(pepckeys_slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm);
int SL_PROTO(pepckeys_sl_MPI_Allreduce)(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank);

/* src/core_mpi/mpi_elements.c */
pepckeys_slint SL_PROTO(pepckeys_mpi_elements_init_keys_from_file)(pepckeys_elements_t *s, char *filename, pepckeys_slint from, pepckeys_slint to, pepckeys_slint const_bytes_per_line, pepckeys_slint root, int size, int rank, MPI_Comm comm);
pepckeys_slint SL_PROTO(pepckeys_mpi_elements_validate_order)(pepckeys_elements_t *s, pepckeys_slint n, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_linear_exchange_pure_keys)(pepckeys_slkey_pure_t *in, pepckeys_slkey_pure_t *out, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_check_order)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t *orders, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_check_global_order)(pepckeys_slkey_pure_t local_min, pepckeys_slkey_pure_t local_max, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_digest_sum)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, unsigned int *sum, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(pepckeys_mpi_cs32)(pepckeys_elements_t *s, pepckeys_slint n, pepckeys_slint keys, pepckeys_slint data, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_digest_hash)(pepckeys_elements_t *s, pepckeys_slint_t nelements, slcint_t components, void *hash, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_get_counts)(pepckeys_elements_t *s, pepckeys_slint_t *clocal, pepckeys_slint_t *cglobal, int root, int size, int rank, MPI_Comm comm);
pepckeys_slweight_t SL_PROTO(pepckeys_mpi_elements_get_weights)(pepckeys_elements_t *s, pepckeys_slweight_t *wlocal, pepckeys_slweight_t *wglobal, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_get_counts_and_weights)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t *counts, pepckeys_slweight_t *weights, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_sendrecv_replace)(pepckeys_elements_t *s, int count, int dest, int sendtag, int source, int recvtag, int size, int rank, MPI_Comm comm);
unsigned int SL_PROTO(pepckeys_mpi_elements_crc32)(pepckeys_elements_t *s, pepckeys_slint_t n, pepckeys_slint_t keys, pepckeys_slint_t data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_elements_alltoall_specific.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoall_specific)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_tproc_f tproc, void *data, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_alltoallv.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_db)(pepckeys_elements_t *sbuf, int *scounts, int *sdispls, pepckeys_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_alltoallv_ip)(pepckeys_elements_t *sbuf, pepckeys_elements_t *sx, int *scounts, int *sdispls, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_elements_packed.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_packed_datatype_create)(MPI_Datatype *pdt, pepckeys_slint_t structured);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_elements_packed_datatype_destroy)(MPI_Datatype *pdt);

/* src/core_mpi/pepckeys_mpi_find_exact.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_find_exact_equal)(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *ex_start, pepckeys_slint_t *ex_size, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_find_exact)(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *dst_size, pepckeys_slint_t *ex_start, pepckeys_slint_t *ex_sizes, pepckeys_slint_t *nx_move, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_linsplit.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_linsplit)(MPI_Comm comm_in, pepckeys_slkey_pure_t *keys_in, MPI_Comm *comms_out, pepckeys_slint_t *parity, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_linsplit_radix)(pepckeys_slkey_pure_t klow, pepckeys_slkey_pure_t khigh, MPI_Comm *comm0, MPI_Comm *comm1, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_linsplit2)(MPI_Comm comm_in, pepckeys_slkey_pure_t *keys_in, MPI_Comm *comms_out, pepckeys_slint_t *parity, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_merge2.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_merge2)(pepckeys_elements_t *s, pepckeys_slint_t other_rank, pepckeys_slint_t high_rank, pepckeys_slint_t *dst_size, pepckeys_merge2x_f m2, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_mergek.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_equal)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_sorted)(pepckeys_elements_t *s, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_mergek_equal2)(pepckeys_elements_t *s, pepckeys_sortnet_f sn, pepckeys_sortnet_data_t snd, pepckeys_merge2x_f m2x, pepckeys_elements_t *xs, int *sizes, int *ranks, MPI_Comm *comms);

/* src/core_mpi/pepckeys_mpi_partition_exact_generic.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_generic)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, pepckeys_binning_t *bm, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_partition_exact_radix.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_radix)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_exact_radix_grouped.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_radix_ngroups)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, pepckeys_slint_t ngroups, MPI_Comm *group_comms, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_exact_radix_2groups)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, MPI_Comm group_comm, pepckeys_elements_t *sx, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_partition_sample.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_partition_sample_regular)(pepckeys_elements_t *s, pepckeys_partcond_t *pcond, int *scounts, int *rcounts, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_rebalance.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_rebalance)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_slint_t stable, pepckeys_slint_t *dst_size, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_rebalance_alltoallv)(pepckeys_elements_t *sbuf, int *scounts, int *sdispls, pepckeys_elements_t *rbuf, int *rcounts, int *rdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_common.c */
pepckeys_slint_t SL_PROTO(pepckeys_init_partconds)(pepckeys_slint_t npconds, pepckeys_partcond_t *pconds, pepckeys_slint_t nparts, pepckeys_slint_t total_count, pepckeys_slweight_t total_weight);
pepckeys_slint_t SL_PROTO(pepckeys_init_partconds_intern)(pepckeys_slint_t npconds, pepckeys_partcond_intern_t *pci, pepckeys_partcond_t *pc, pepckeys_slint_t nparts, pepckeys_slint_t total_count, pepckeys_slweight_t total_weight);
pepckeys_slint_t SL_PROTO(pepckeys_merge_partconds)(pepckeys_partcond_t *pconds_in, pepckeys_slint_t npconds_in, pepckeys_partcond_t *pcond_out);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_gather_partconds_grouped)(pepckeys_partcond_t *pcond_in, MPI_Comm pcond_in_comm, MPI_Comm pconds_out_comm, pepckeys_partcond_t *pconds_out, pepckeys_slint_t *npconds_out, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_gather_partconds)(pepckeys_partcond_t *pcond_in, pepckeys_partcond_t *pconds_out, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_allgather_partconds)(pepckeys_partcond_t *pcond_in, pepckeys_partcond_t *pconds_out, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_bcast_partconds)(pepckeys_slint_t npconds, pepckeys_partcond_t *pconds, int root, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_post_check_partconds)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_post_check_partconds_intern)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_intern_t *pci, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_stats)(pepckeys_elements_t *s, pepckeys_slint_t nparts, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_select_exact_generic.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_generic_bulk)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_binning_t *bm, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_generic_grouped)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepckeys_binning_t *bm, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_generic)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_binning_t *bm, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_select_exact_radix.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_radix)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_exact_radix_grouped)(pepckeys_elements_t *s, pepckeys_slint_t nelements, pepckeys_partcond_t *pcond, MPI_Comm pcond_comm, MPI_Comm group_comm, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, pepckeys_slint_t sorted, int *sdispls, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_select_sample.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_select_sample_regular)(pepckeys_elements_t *s, pepckeys_slint_t nparts, pepckeys_partcond_t *pconds, pepckeys_slint_t nsamples, pepckeys_splitter_t *sp, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sort_merge.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_merge)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_merge2)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t merge_type, pepckeys_slint_t sort_type, double *times, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_merge_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t merge_type, pepckeys_slint_t sort_type, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);

/* src/core_mpi/pepckeys_mpi_sort_partition.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t part_type, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t part_type, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_exact_radix)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_partcond_t *pcond, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_exact_radix_ngroups)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_partcond_t *pcond, pepckeys_slint_t ngroups, MPI_Comm *group_comms, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_partition_exact_radix_2groups)(pepckeys_elements_t *s, pepckeys_elements_t *sx, pepckeys_partcond_t *pcond, MPI_Comm group_comm, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_sort_special.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_insert_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slpkey_t *mmkeys, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_presorted_radix)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slint_t merge_type, pepckeys_slint_t rhigh, pepckeys_slint_t rlow, pepckeys_slint_t rwidth, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_sort_back)(pepckeys_elements_t *s0, pepckeys_elements_t *s1, pepckeys_elements_t *xs, pepckeys_slpkey_t *lh, pepckeys_slint_t ntotal, int size, int rank, MPI_Comm comm);

/* src/core_mpi/mpi_xcounts2ycounts.c */
pepckeys_slint_t SL_PROTO(pepckeys_mpi_xcounts2ycounts_all2all)(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_xcounts2ycounts_sparse)(int *xcounts, int *ycounts, pepckeys_slint_t ytotal, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_xcounts2ycounts_grouped)(int *xcounts, pepckeys_slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm);
pepckeys_slint_t SL_PROTO(pepckeys_mpi_subxdispls2ycounts)(pepckeys_slint_t nsubs, int *sub_xdispls, pepckeys_slint_t *sub_sources, pepckeys_slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm);


#endif /* SL_USE_MPI */


#undef SL_PROTO
#endif /* __SL_PEPCKEYS_H__ */
