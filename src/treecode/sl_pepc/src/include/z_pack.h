/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/z_pack.h
 *  
 */


#ifndef __Z_PACK_H__
#define __Z_PACK_H__


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include "z_pack_conf.h"

#ifdef Z_PACK_RENAME
# include "z_pack_rename.h"
#endif


#define Z_MOP(_mop_)  do { _mop_ } while (0)
#define Z_NOP()       Z_MOP()


#ifdef Z_PACK_MPI
# ifdef Z_MPI_COMM_PTR
#  define zcomm_fmt       "p"
#  define zcomm_val(_c_)  (_c_)
# else
#  define zcomm_fmt       "s"
#  define zcomm_val(_c_)  (((_c_) != MPI_COMM_NULL)?"valid":"null")
# endif
#endif


#ifdef Z_PACK_NUMERIC

#define z_max(_a_,_b_)            (((_a_)>(_b_))?(_a_):(_b_))
#define z_min(_a_,_b_)            (((_a_)<(_b_))?(_a_):(_b_))
#define z_max3(_a_,_b_,_c_)       z_max(_a_,z_max(_b_,_c_))
#define z_min3(_a_,_b_,_c_)       z_min(_a_,z_min(_b_,_c_))
#define z_minmax(_a_, _b_, _c_)   (((_b_)<(_a_))?(_a_):(((_b_)>(_c_))?(_c_):(_b_)))
#define z_abs(_a_)                (((_a_) >= 0)?(_a_):-(_a_))

#include <math.h>

#if __STDC_VERSION__ >= 199901L
# define z_round(_a_)             round(_a_)
#else
# define z_round(_a_)             (((_a_) >= 0)?floor((_a_) + 0.5):ceil((_a_) - 0.5))
#endif

#define z_powof2_typed(_a_, _t_)  (((_t_) 1) << (_a_))
#define z_powof2(_a_)             z_powof2_typed(_a_, z_int_t)

#define z_ispowof2(_a_)           ((_a_ & (_a_ - 1)) == 0)

#define z_get1d(_x0_)                                       (_x0_)
#define z_get2d(_x1_, _d0_, _x0_)                          ((_x0_) + (_d0_) *  (_x1_))
#define z_get3d(_x2_, _d1_, _x1_, _d0_, _x0_)              ((_x0_) + (_d0_) * ((_x1_) + (_d1_) *  (_x2_)))
#define z_get4d(_x3_, _d2_, _x2_, _d1_, _x1_, _d0_, _x0_)  ((_x0_) + (_d0_) * ((_x1_) + (_d1_) * ((_x2_) + (_d2_) * (_x3_))))

#endif /* Z_PACK_NUMERIC */


#ifdef Z_PACK_DEBUG

#include <stdio.h>

extern FILE *z_notice_fstream, *z_error_fstream, *z_debug_fstream;

#define Z_NOTICE_FSTREAM  (z_notice_fstream?z_notice_fstream:stdout)
#define Z_ERROR_FSTREAM   (z_error_fstream?z_error_fstream:stderr)
#ifdef Z_DEBUG_FSTREAM_STDERR
# define Z_DEBUG_FSTREAM  (z_debug_fstream?z_debug_fstream:stderr)
#else
# define Z_DEBUG_FSTREAM  (z_debug_fstream?z_debug_fstream:stdout)
#endif

#if !defined(Z_DEBUG_MESG_STR) || !defined(Z_DEBUG_MESG_ARG)
# undef Z_DEBUG_MESG_STR
# undef Z_DEBUG_MESG_ARG
# ifdef Z_PACK_MPI
#  define Z_DEBUG_MESG_STR  "%d: "
#  define Z_DEBUG_MESG_ARG  Z_PACK_MPI_RANK
# else
#  define Z_DEBUG_MESG_STR  "%s"
#  define Z_DEBUG_MESG_ARG  ""
# endif
#endif

#if !defined(Z_DEBUG_CODE_STR) || !defined(Z_DEBUG_CODE_ARG)
# undef Z_DEBUG_CODE_STR
# undef Z_DEBUG_CODE_ARG
# define Z_DEBUG_CODE_STR  "%s:%i:%s: "
# define Z_DEBUG_CODE_ARG(_fi_, _li_, _fu_)  _fi_, _li_, _fu_
#endif

#define Z_FPRINTF(_stream_, _format_, ...)                          fprintf(_stream_, _format_ "%s", __VA_ARGS__)
#define Z_FPRINTF_MESG(_stream_, _format_, ...)                     fprintf(_stream_, Z_DEBUG_MESG_STR _format_ "%s", Z_DEBUG_MESG_ARG, __VA_ARGS__)
#define Z_FPRINTF_CODE(_stream_, _format_, ...)                     fprintf(_stream_, Z_DEBUG_MESG_STR Z_DEBUG_CODE_STR _format_ "%s", Z_DEBUG_MESG_ARG, Z_DEBUG_CODE_ARG(__FILE__, __LINE__, __func__), __VA_ARGS__)
#define Z_FPRINTF_CODE_(_fi_, _li_, _fu_, _stream_, _format_, ...)  fprintf(_stream_, Z_DEBUG_MESG_STR Z_DEBUG_CODE_STR _format_ "%s", Z_DEBUG_MESG_ARG, Z_DEBUG_CODE_ARG(_fi_, _li_, _fu_), __VA_ARGS__)

#define Z_NOTICE(...)                              Z_FPRINTF_MESG(Z_NOTICE_FSTREAM, __VA_ARGS__, "\n")
#define Z_NOTICE_IF(_if_, ...)                     Z_MOP(if (_if_) Z_NOTICE(__VA_ARGS__);)
#define Z_NOTICE_(_fi_, _li_, _fu_, ...)           Z_FPRINTF_MESG_(_fi_, _li_, _fu_, Z_NOTICE_FSTREAM, __VA_ARGS__, "\n")
#define Z_NOTICE_IF_(_fi_, _li_, _fu_, _if_, ...)  Z_MOP(if (_if_) Z_NOTICE_(_fi_, _li_, _fu_, __VA_ARGS__);)

#define Z_ERROR(...)                              Z_FPRINTF_MESG(Z_ERROR_FSTREAM, __VA_ARGS__, "\n")
#define Z_ERROR_IF(_if_, ...)                     Z_MOP(if (_if_) Z_ERROR(__VA_ARGS__);)
#define Z_ERROR_(_fi_, _li_, _fu_, ...)           Z_FPRINTF_MESG_(_fi_, _li_, _fu_, Z_ERROR_FSTREAM, __VA_ARGS__, "\n")
#define Z_ERROR_IF_(_fi_, _li_, _fu_, _if_, ...)  Z_MOP(if (_if_) Z_ERROR_(_fi_, _li_, _fu_, __VA_ARGS__);)

#ifdef Z_DEBUG_LEVEL
# define Z_DEBUG_INTRO(_level_, ...)                        Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF_CODE(Z_DEBUG_FSTREAM, __VA_ARGS__, "");)
# define Z_DEBUG_CORE(_level_, ...)                         Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF(Z_DEBUG_FSTREAM, __VA_ARGS__, "");)
# define Z_DEBUG_OUTRO(_level_, ...)                        Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF(Z_DEBUG_FSTREAM, __VA_ARGS__, "\n");)
# define Z_DEBUG(_level_, ...)                              (((_level_) <= (Z_DEBUG_LEVEL))?(Z_FPRINTF_CODE(Z_DEBUG_FSTREAM, __VA_ARGS__, "\n")):0)
# define Z_DEBUG_IF(_if_, _level_, ...)                     (((_if_) && ((_level_) <= (Z_DEBUG_LEVEL)))?(Z_FPRINTF_CODE(Z_DEBUG_FSTREAM, __VA_ARGS__, "\n")):0)
# define Z_DEBUG_(_fi_, _li_, _fu_, _level_, ...)           Z_MOP(if ((_level_) <= (Z_DEBUG_LEVEL)) Z_FPRINTF_CODE_(_fi_, _li_, _fu_, Z_DEBUG_FSTREAM, __VA_ARGS__, "\n");)
# define Z_DEBUG_IF_(_fi_, _li_, _fu_, _if_, _level_, ...)  Z_MOP(if ((_if_) && ((_level_) <= (Z_DEBUG_LEVEL))) Z_FPRINTF_CODE_(_fi_, _li_, _fu_, Z_DEBUG_FSTREAM, __VA_ARGS__, "\n");)
#else
# define Z_DEBUG_INTRO(...)                                 Z_NOP()
# define Z_DEBUG_CORE(...)                                  Z_NOP()
# define Z_DEBUG_OUTRO(...)                                 Z_NOP()
# define Z_DEBUG(...)                                       Z_NOP()
# define Z_DEBUG_IF(...)                                    Z_NOP()
# define Z_DEBUG_(...)                                      Z_NOP()
# define Z_DEBUG_IF_(...)                                   Z_NOP()
#endif

#define Z_TRACE_LEVEL  3

#define Z_TRACE(...)                              Z_DEBUG(Z_TRACE_LEVEL, __VA_ARGS__)
#define Z_TRACE_IF(_if_, ...)                     Z_DEBUG_IF(_if_, Z_TRACE_LEVEL, __VA_ARGS__)
#define Z_TRACE_(_fi_, _li_, _fu_, ...)           Z_DEBUG_(_fi_, _li_, _fu_, Z_TRACE_LEVEL, __VA_ARGS__)
#define Z_TRACE_IF_(_fi_, _li_, _fu_, _if_, ...)  Z_DEBUG_IF_(_fi_, _li_, _fu_, _if_, Z_TRACE_LEVEL, __VA_ARGS__)

#define Z_TRACE_ARRAY(_i_, _n_, _ef_, _e_, ...) \
  Z_MOP(Z_DEBUG_INTRO(Z_TRACE_LEVEL, __VA_ARGS__); \
        for (_i_ = 0; _i_ < (_n_); ++_i_) Z_DEBUG_CORE(Z_TRACE_LEVEL, _ef_, _e_); \
        Z_DEBUG_OUTRO(Z_TRACE_LEVEL, "");)
#define Z_TRACE_ARRAY_IF(_if_, ...)  Z_MOP(if (_if_) Z_TRACE_ARRAY(__VA_ARGS__);)

#define Z_ASSERT_LEVEL  0

#define Z_ASSERT(_x_)           Z_MOP(if (_x_); else Z_DEBUG(Z_ASSERT_LEVEL, "ASSERT: '" #_x_ "' failed."); )
#define Z_ASSERT_IF(_if_, _x_)  Z_MOP(if (_if_) Z_ASSERT(_x_); )

#endif /* Z_PACK_DEBUG */


#ifdef Z_PACK_ALLOC

#include <stdlib.h>

#ifndef z_alloc_hook
# define z_alloc_hook_func(_n_, _s_, _p_, _fi_, _li_, _fu_)  (_p_)
#else
inline static void *z_alloc_hook_func(z_int_t n, z_int_t s, void *p, const char *fi, int li, const char *fu)
{
  z_alloc_hook(n, s, p, fi, li, fu);
  return p;
}
#endif
#ifndef z_free_hook
# define z_free_hook(_p_)
#endif

#ifdef Z_ALLOC_DEBUG
# define z_alloc(_n_, _s_)  z_alloc_hook_func((z_int_t) _n_, (z_int_t) _s_, calloc((_n_), (_s_)), __FILE__, __LINE__, __func__)
# define z_free(_p_)        Z_MOP(z_free_hook(_p_); free(_p_); _p_ = NULL;)
#else
# define z_alloc(_n_, _s_)  z_alloc_hook_func((z_int_t) _n_, (z_int_t) _s_, malloc((_n_) * (_s_)), __FILE__, __LINE__, __func__)
# define z_free(_p_)        Z_MOP(z_free_hook(_p_); free(_p_);)
#endif

#ifndef z_alloca_hook
# define z_alloca_hook_func(_n_, _s_, _p_, _fi_, _li_, _fu_)  (_p_)
#else
inline static void *z_alloca_hook_func(z_int_t n, z_int_t s, void *p, const char *fi, int li, const char *fu)
{
  z_alloca_hook(n, s, p, fi, li, fu);
  return p;
}
#endif
#ifndef z_freea_hook
# define z_freea_hook(_p_)
#endif

#include <alloca.h>

#define z_alloca(_n_, _s_)  z_alloca_hook_func((z_int_t) _n_, (z_int_t) _s_, alloca((_n_) * (_s_)), __FILE__, __LINE__, __func__)
#ifdef Z_ALLOC_DEBUG
# define z_freea(_p_)       Z_MOP(z_freea_hook(_p_); _p_ = NULL;)
#else
# define z_freea(_p_)       Z_MOP(z_freea_hook(_p_);)
#endif

#endif /* Z_PACK_ALLOC */


#ifdef Z_PACK_TIME

#ifdef Z_PACK_MPI
# include <mpi.h>
typedef double z_time_t;
# define z_time_save(_t_)            (_t_ = MPI_Wtime())
# define z_time_diff_s(_f_, _t_)     ((_t_) - (_f_))
# define z_time_get_s()              (MPI_Wtime())
inline static double z_time_wtime() { return MPI_Wtime(); }
#else
# include <sys/time.h>
typedef struct timeval z_time_t;
# define z_time_save(_t_)            (gettimeofday(&(_t_), NULL))
# define z_time_diff_s(_f_, _t_)     ((double) (((_t_).tv_sec - (_f_).tv_sec) + ((_t_).tv_usec - (_f_).tv_usec) / 1000000.0))
# define z_time_get_s()              (z_time_wtime())
double z_time_wtime();
#endif

#endif /* Z_PACK_TIME */


#ifdef Z_PACK_RANDOM

#ifndef z_rand
# undef z_srand
# undef Z_RAND_MIN
# undef Z_RAND_MAX
# if defined(HAVE_RANDOM)
#  include <stdlib.h>
#  define Z_RAND_MIN    0
#  define Z_RAND_MAX    RAND_MAX
#  define z_srand(_s_)  srandom(_s_)
#  define z_rand()      random()
# elif defined(HAVE_RAND)
#  include <stdlib.h>
#  define Z_RAND_MIN    0
#  define Z_RAND_MAX    RAND_MAX
#  define z_srand(_s_)  srand(_s_)
#  define z_rand()      rand()
# else
void z_srandom(unsigned long seed);
unsigned long z_random();
#  define Z_RAND_MIN    0
#  define Z_RAND_MAX    0xFFFFFFFF
#  define z_srand(_s_)  z_srandom(_s_)
#  define z_rand()      z_random()
# endif
#endif

#ifndef z_rand_minmax
# define z_rand_minmax(_min_, _max_)  (_min_ + ((double) (_max_ - _min_) * (z_rand() - Z_RAND_MIN) / (Z_RAND_MAX - Z_RAND_MIN)))
#endif

void z_srandom64(unsigned long seed);
long long z_random64();
long long z_random64_minmax(long long min, long long max);
unsigned long long z_random64u();
unsigned long long z_random64u_minmax(unsigned long long min, unsigned long long max);

void z_nrandom_seed(unsigned long s);
double z_nrandom();
void z_urandom_seed(unsigned long s);
double z_urandom();

#endif /* Z_PACK_RANDOM */


#ifdef Z_PACK_DIGEST

z_int_t z_digest_sum_buffer(const void *buffer, z_int_t length, void *sum);
#ifdef HAVE_GCRYPT_H
extern int z_digest_hash_gcrypt_algo;
#endif
void z_digest_hash_open(void **hdl);
void z_digest_hash_close(void *hdl);
void z_digest_hash_write(void *hdl, const void *buffer, z_int_t length);
z_int_t z_digest_hash_read(void *hdl, void *hash);

#endif /* Z_PACK_DIGEST */


#ifdef Z_PACK_CRC32

extern const z_int_t z_crc32_table_size;
extern const z_crc32_t z_crc32_table[];

void z_crc32_make_table(z_crc32_t *tbl);
void z_crc32_print_table(z_crc32_t *tbl);

z_crc32_t z_crc32_buffer_update(z_crc32_t crc, const void *buffer, z_int_t length);
z_crc32_t z_crc32_buffer(const void *buffer, z_int_t length);

#endif /* Z_CRC32 */


#if defined(Z_PACK_GMP) && defined(HAVE_GMP_H)

#ifdef HAVE_GMP_H
# include <gmp.h>
#endif

void z_gmp_mpz_set_ull(mpz_t z, unsigned long long v);
void z_gmp_mpz_set_sll(mpz_t z, long long v);
unsigned long long z_gmp_mpz_get_ull(mpz_t z);
long long z_gmp_mpz_get_sll(mpz_t z);

#endif /* Z_PACK_GMP && HAVE_GMP_H */


#endif /* __Z_PACK_H__ */
