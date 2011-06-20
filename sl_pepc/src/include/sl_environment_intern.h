/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_environment_intern.h
 *  
 */


#ifndef __SL_ENVIRONMENT_INTERN_H__
#define __SL_ENVIRONMENT_INTERN_H__


#ifndef ENV_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define ENV_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define ENV_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


#define z_alloc_hook(_n_, _s_, _p_, _fi_, _li_, _fu_)  Z_MOP( \
  Z_TRACE_IF_(_fi_, _li_, _fu_, ENV_TRACE_IF, "z_alloc: %" slint_fmt " * %" slint_fmt " -> %p", (_n_), (_s_), (_p_)); \
  rti_minc_alloc(); \
  rti_malloc((_n_) * (_s_));)

#define z_free_hook(_p_)  Z_MOP( \
  Z_TRACE_IF(ENV_TRACE_IF, "z_free: %p", (_p_)); \
  rti_minc_free();)

#define z_alloca_hook(_n_, _s_, _p_, _fi_, _li_, _fu_) \
  Z_TRACE_IF_(_fi_, _li_, _fu_, ENV_TRACE_IF, "z_alloca: %" slint_fmt " * %" slint_fmt " -> %p", (_n_), (_s_), (_p_))

#define z_freea_hook(_p_) \
  Z_TRACE_IF(ENV_TRACE_IF, "z_freea: %p", (_p_))


#endif /* __SL_ENVIRONMENT_INTERN_H__ */
