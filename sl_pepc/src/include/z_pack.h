/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/z_pack.h
 *  timestamp: 2011-03-06 21:59:31 +0100
 *  
 */


#ifndef __Z_PACK_H__
#define __Z_PACK_H__


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif


#include "z_pack_conf.h"
#include "z_pack_rename.h"


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

#define z_round(_a_)              round(_a_)

#define z_powof2_typed(_a_, _t_)  (((_t_) 1) << (_a_))
#define z_powof2(_a_)             z_powof2_typed(_a_, z_int_t)

#endif /* Z_PACK_NUMERIC */


#ifdef Z_PACK_ALLOC

#ifndef cc_z_alloc_pre_hook
# define cc_z_alloc_pre_hook(_n_, _s_, _file_, _line_, _func_)
#endif
#ifndef z_alloc_post_hook
# define z_alloc_post_hook(_n_, _s_, _p_, _file_, _line_, _func_)  (_p_)
#endif
#ifndef z_free_hook
# define z_free_hook(_p_)
#endif

#ifdef Z_ALLOC_DEBUG
# define z_alloc(_n_, _s_)  (cc_z_alloc_pre_hook(_n_, _s_, __FILE__, __LINE__, __func__) z_alloc_post_hook(_n_, _s_, calloc((_n_), (_s_)), __FILE__, __LINE__, __func__))
# define z_free(_p_)        Z_MOP(z_free_hook(_p_); free(_p_); _p_ = NULL;)
#else
# define z_alloc(_n_, _s_)  (cc_z_alloc_pre_hook(_n_, _s_, __FILE__, __LINE__, __func__) z_alloc_post_hook(_n_, _s_, malloc((_n_) * (_s_)), __FILE__, __LINE__, __func__))
# define z_free(_p_)        Z_MOP(z_free_hook(_p_); free(_p_);)
#endif

#ifndef cc_z_alloca_pre_hook
# define cc_z_alloca_pre_hook(_n_, _s_, _file_, _line_, _func_)
#endif
#ifndef z_alloca_post_hook
# define z_alloca_post_hook(_n_, _s_, _p_, _file_, _line_, _func_)  (_p_)
#endif
#ifndef z_freea_hook
# define z_freea_hook(_p_)
#endif

#include <alloca.h>

#define z_alloca(_n_, _s_)  (cc_z_alloca_pre_hook(_n_, _s_, __FILE__, __LINE__, __func__) z_alloca_post_hook(_n_, _s_, alloca((_n_) * (_s_)), __FILE__, __LINE__, __func__))
#ifdef Z_ALLOC_DEBUG
# define z_freea(_p_)       Z_MOP(z_freea_hook(_p_); _p_ = NULL;)
#else
# define z_freea(_p_)       Z_MOP(z_freea_hook(_p_);)
#endif

#endif /* Z_PACK_ALLOC */


#ifdef Z_PACK_DEBUG

extern FILE *z_notice_fstream, *z_error_fstream, *z_debug_fstream;

#define Z_NOTICE_FSTREAM  (z_notice_fstream?z_notice_fstream:stdout)
#define Z_ERROR_FSTREAM   (z_error_fstream?z_error_fstream:stderr)
#ifdef Z_DEBUG_FSTREAM_STDERR
# define Z_DEBUG_FSTREAM  (z_debug_fstream?z_debug_fstream:stderr)
#else
# define Z_DEBUG_FSTREAM  (z_debug_fstream?z_debug_fstream:stdout)
#endif

#ifdef Z_PACK_MPI
# if !defined(Z_DEBUG_MPI_STR) || !defined(Z_DEBUG_MPI_PARAM)
#  undef Z_DEBUG_MPI_STR
#  define Z_DEBUG_MPI_STR    "%d: "
#  undef Z_DEBUG_MPI_PARAM
#  define Z_DEBUG_MPI_PARAM  Z_PACK_MPI_RANK
# endif
#else
# undef Z_DEBUG_MPI_STR
# define Z_DEBUG_MPI_STR     "%s"
# undef Z_DEBUG_MPI_PARAM
# define Z_DEBUG_MPI_PARAM   ""
#endif

#define Z_NOTICE(_format_, _args_...)           fprintf(Z_NOTICE_FSTREAM, Z_DEBUG_MPI_STR _format_ "\n", Z_DEBUG_MPI_PARAM, ##_args_)
#define Z_NOTICE_IF(_if_, _format_, _args_...)  Z_MOP(if (_if_) Z_NOTICE(_format_, ##_args_); )

#define Z_ERROR(_format_, _args_...)            fprintf(Z_ERROR_FSTREAM, Z_DEBUG_MPI_STR "%s:%i:%s: " _format_ "\n", Z_DEBUG_MPI_PARAM, __FILE__, __LINE__, __func__, ##_args_)
#define Z_ERROR_IF(_if_, _format_, _args_...)   Z_MOP(if (_if_) Z_ERROR(_format_, ##_args_); )

#ifdef Z_DEBUG_LEVEL
# define Z_DEBUG_INTRO(_level_) \
    Z_MOP(if (_level_ <= Z_DEBUG_LEVEL) fprintf(Z_DEBUG_FSTREAM, Z_DEBUG_MPI_STR "%s:%i:%s: ", Z_DEBUG_MPI_PARAM, __FILE__, __LINE__, __func__); )
# define Z_DEBUG_CORE(_level_, _format_, _args_...) \
    Z_MOP(if (_level_ <= Z_DEBUG_LEVEL) fprintf(Z_DEBUG_FSTREAM, _format_, ##_args_); )
# define Z_DEBUG_OUTRO(_level_) \
    Z_MOP(if (_level_ <= Z_DEBUG_LEVEL) fprintf(Z_DEBUG_FSTREAM, "\n"); )
# define Z_DEBUG(_level_, _format_, _args_... ) \
    ((_level_ <= Z_DEBUG_LEVEL)?(fprintf(Z_DEBUG_FSTREAM, Z_DEBUG_MPI_STR "%s:%i:%s: " _format_ "\n", Z_DEBUG_MPI_PARAM, __FILE__, __LINE__, __func__, ##_args_)):0)
# define Z_DEBUG_IF(_if_, _level_, _format_, _args_... ) \
    (((_if_) && (_level_ <= Z_DEBUG_LEVEL))?(fprintf(Z_DEBUG_FSTREAM, Z_DEBUG_MPI_STR "%s:%i:%s: " _format_ "\n", Z_DEBUG_MPI_PARAM, __FILE__, __LINE__, __func__, ##_args_)):0)
#else
# define Z_DEBUG_INTRO(_x_...)                  Z_NOP()
# define Z_DEBUG_CORE(_x_...)                   Z_NOP()
# define Z_DEBUG_OUTRO(_x_...)                  Z_NOP()
# define Z_DEBUG(_x_...)                        0
# define Z_DEBUG_IF(_x_...)                     0
#endif

#define Z_ASSERT(_x_)                           Z_MOP(if (_x_); else Z_DEBUG(0, "ASSERT: '%s' failed.", #_x_); )
#define Z_ASSERT_IF(_if_, _x_)                  Z_MOP(if (_if_) Z_ASSERT(_x_); )

#define Z_TRACE(_format_, _args_...)            Z_DEBUG(3, _format_, ##_args_)
#define Z_TRACE_IF(_if_, _format_, _args_...)   Z_DEBUG_IF(_if_, 3, _format_, ##_args_)

#define Z_TRACE_ARRAY(_f_, _e_, _i_, _n_, _a_, _args_...) \
   Z_MOP(Z_DEBUG_INTRO(3); \
          Z_DEBUG_CORE(3, _f_, ##_args_); \
          for (_i_ = 0; _i_ < _n_; ++_i_) Z_DEBUG_CORE(3, _e_, _a_[_i_]); \
          Z_DEBUG_OUTRO(3); )
#define Z_TRACE_ARRAY_IF(_if_, _f_, _e_, _i_, _n_, _a_, _args_...) Z_MOP(if (_if_) Z_TRACE_ARRAY(_f_, _e_, _i_, _n_, _a_, ##_args_); )

#endif /* Z_PACK_DEBUG */


#ifdef Z_PACK_TIME

#ifdef Z_PACK_MPI
typedef double z_time_t;
# define z_time_save(_t_)            (_t_ = MPI_Wtime())
# define z_time_diff_s(_f_, _t_)     ((_t_) - (_f_))
# define z_time_get_s()              (MPI_Wtime())
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
# include <stdlib.h>
# define Z_RAND_MIN    0
# define Z_RAND_MAX    RAND_MAX
# define z_srand(_s_)  srand(_s_)
# define z_rand()      rand()
#endif

#ifndef z_rand_minmax
# define z_rand_minmax(_min_, _max_)  (_min_ + ((double) (_max_ - _min_) * (z_rand() - Z_RAND_MIN) / (Z_RAND_MAX - Z_RAND_MIN)))
#endif

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


#endif /* __Z_PACK_H__ */
