/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_config_intern.h
 *  timestamp: 2009-12-03 09:11:41 +0100
 *  
 */


#ifndef __SL_CONFIG_INTERN_H__
#define __SL_CONFIG_INTERN_H__


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


#ifndef SL_INDEX
# undef SL_PACKED_INDEX
#endif


/* if no special, given, primary and heavy used integer-type ... */
#ifndef sl_int_type_c
  /* ... use a default one */
# define sl_int_type_c               long      /* sl_macro */
# undef sl_int_type_mpi
# define sl_int_type_mpi             MPI_LONG  /* sl_macro */
# undef sl_int_size_mpi
# define sl_int_size_mpi             1         /* sl_macro */
# undef sl_int_type_fmt
# define sl_int_type_fmt             "ld"      /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(sl_int_type_mpi) || !defined(sl_int_size_mpi)
#   error "sl_int_type_mpi and/or sl_int_size_mpi missing"
#  endif
# endif
# ifndef sl_int_type_fmt
#  error "sl_int_type_fmt macro is missing, using d as default"
#  define sl_int_type_fmt  "d"
# endif
#endif


/* if no special datatype for indexes ... */
#ifndef sl_index_type_c
 /* ... use the primary integer type */
# define sl_index_type_c             sl_int_type_c    /* sl_macro */
# undef sl_index_type_mpi
# define sl_index_type_mpi           sl_int_type_mpi  /* sl_macro */
# undef sl_index_size_mpi
# define sl_index_size_mpi           sl_int_size_mpi  /* sl_macro */
# undef sl_index_type_fmt
# define sl_index_type_fmt           sl_int_type_fmt  /* sl_macro */
#else
  /* ... use the given one and check whether an mpi-type is present and required */
# ifdef SL_USE_MPI
#  if !defined(sl_index_type_mpi) || !defined(sl_index_size_mpi)
#   error "sl_index_type_mpi and/or sl_index_size_mpi missing"
#  endif
# endif
# ifndef sl_index_type_fmt
#  error "sl_index_type_fmt macro is missing, using d as default"
#  define sl_index_type_fmt  "d"
# endif
#endif


/* default pure keys */
#ifndef sl_key_pure_type_c
# define sl_key_pure_type_c          sl_key_type_c  /* sl_macro */
#endif
#ifndef sl_key_pure_type_mpi
# define sl_key_pure_type_mpi        sl_key_type_mpi  /* sl_macro */
#endif
#ifndef sl_key_pure_size_mpi
# define sl_key_pure_size_mpi        sl_key_size_mpi  /* sl_macro */
#endif
#ifndef sl_key_pure_type_fmt
# ifdef sl_key_type_fmt
#  define sl_key_pure_type_fmt       sl_key_type_fmt  /* sl_macro */
# endif
#endif

#ifndef sl_key_purify
 #define sl_key_purify(k)            (k)  /* sl_macro */
#endif
#ifndef sl_key_get_pure
 #define sl_key_get_pure(k)          (sl_key_purify(k))  /* sl_macro */
#endif
#ifndef sl_key_set_pure
 #define sl_key_set_pure(k, p)       (sl_key_purify(k) = p)  /* sl_macro */
#endif


/* default pure key comparisons */
#ifndef sl_key_pure_cmp_eq
 #define sl_key_pure_cmp_eq(k0, k1)  ((k0) == (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_ne
 #define sl_key_pure_cmp_ne(k0, k1)  ((k0) != (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_lt
 #define sl_key_pure_cmp_lt(k0, k1)  ((k0) < (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_le
 #define sl_key_pure_cmp_le(k0, k1)  ((k0) <= (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_gt
 #define sl_key_pure_cmp_gt(k0, k1)  ((k0) > (k1))  /* sl_macro */
#endif
#ifndef sl_key_pure_cmp_ge
 #define sl_key_pure_cmp_ge(k0, k1)  ((k0) >= (k1))  /* sl_macro */
#endif


/* default key comparisons */
#ifndef sl_key_cmp_eq
 #define sl_key_cmp_eq(k0, k1)       (sl_key_pure_cmp_eq(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_ne
 #define sl_key_cmp_ne(k0, k1)       (sl_key_pure_cmp_ne(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_lt
 #define sl_key_cmp_lt(k0, k1)       (sl_key_pure_cmp_lt(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_le
 #define sl_key_cmp_le(k0, k1)       (sl_key_pure_cmp_le(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_gt
 #define sl_key_cmp_gt(k0, k1)       (sl_key_pure_cmp_gt(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif
#ifndef sl_key_cmp_ge
 #define sl_key_cmp_ge(k0, k1)       (sl_key_pure_cmp_ge(sl_key_purify(k0), sl_key_purify(k1)))  /* sl_macro */
#endif


/* disable data components on request */
/* DATAX_TEMPLATE_BEGIN */
/* sl_macro SL_DATA0_IGNORE */
#ifdef SL_DATA0_IGNORE
# undef SL_DATA0
#endif
/* sl_macro SL_DATA1_IGNORE */
#ifdef SL_DATA1_IGNORE
# undef SL_DATA1
#endif
/* sl_macro SL_DATA2_IGNORE */
#ifdef SL_DATA2_IGNORE
# undef SL_DATA2
#endif
/* sl_macro SL_DATA3_IGNORE */
#ifdef SL_DATA3_IGNORE
# undef SL_DATA3
#endif
/* sl_macro SL_DATA4_IGNORE */
#ifdef SL_DATA4_IGNORE
# undef SL_DATA4
#endif
/* sl_macro SL_DATA5_IGNORE */
#ifdef SL_DATA5_IGNORE
# undef SL_DATA5
#endif
/* sl_macro SL_DATA6_IGNORE */
#ifdef SL_DATA6_IGNORE
# undef SL_DATA6
#endif
/* sl_macro SL_DATA7_IGNORE */
#ifdef SL_DATA7_IGNORE
# undef SL_DATA7
#endif
/* sl_macro SL_DATA8_IGNORE */
#ifdef SL_DATA8_IGNORE
# undef SL_DATA8
#endif
/* sl_macro SL_DATA9_IGNORE */
#ifdef SL_DATA9_IGNORE
# undef SL_DATA9
#endif
/* sl_macro SL_DATA10_IGNORE */
#ifdef SL_DATA10_IGNORE
# undef SL_DATA10
#endif
/* sl_macro SL_DATA11_IGNORE */
#ifdef SL_DATA11_IGNORE
# undef SL_DATA11
#endif
/* sl_macro SL_DATA12_IGNORE */
#ifdef SL_DATA12_IGNORE
# undef SL_DATA12
#endif
/* sl_macro SL_DATA13_IGNORE */
#ifdef SL_DATA13_IGNORE
# undef SL_DATA13
#endif
/* sl_macro SL_DATA14_IGNORE */
#ifdef SL_DATA14_IGNORE
# undef SL_DATA14
#endif
/* sl_macro SL_DATA15_IGNORE */
#ifdef SL_DATA15_IGNORE
# undef SL_DATA15
#endif
/* sl_macro SL_DATA16_IGNORE */
#ifdef SL_DATA16_IGNORE
# undef SL_DATA16
#endif
/* sl_macro SL_DATA17_IGNORE */
#ifdef SL_DATA17_IGNORE
# undef SL_DATA17
#endif
/* sl_macro SL_DATA18_IGNORE */
#ifdef SL_DATA18_IGNORE
# undef SL_DATA18
#endif
/* sl_macro SL_DATA19_IGNORE */
#ifdef SL_DATA19_IGNORE
# undef SL_DATA19
#endif
/* DATAX_TEMPLATE_END */


/* default element weights */
#ifndef sl_elem_weight
# define sl_elem_weight(e, at)       1  /* sl_macro */
#endif


#endif /* __SL_CONFIG_INTERN_H__ */
