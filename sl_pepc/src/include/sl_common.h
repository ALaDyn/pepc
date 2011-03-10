/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_common.h
 *  timestamp: 2011-03-06 21:59:31 +0100
 *  
 */


#ifndef __SL_COMMON_H__
#define __SL_COMMON_H__


#ifdef HAVE_CONFIG_H
# include <config.h>
#else
# ifdef __bgp__
#  define HAVE_SPI_KERNEL_INTERFACE_H
#  define HAVE_COMMON_BGP_PERSONALITY_H
#  define HAVE_COMMON_BGP_PERSONALITY_INLINES_H
#  define HAVE__BGP_PERSONALITY_T
# endif
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>


#ifdef SLDEBUG
# ifndef SLDEBUG_OUTPUT
#  define SLDEBUG_OUTPUT  SLDEBUG
# endif
# ifndef SLDEBUG_ALLOC
#  define SLDEBUG_ALLOC   SLDEBUG
# endif
#endif

#ifdef SLDEBUG_OUTPUT_NOT
# undef SLDEBUG_OUTPUT
#endif

#ifdef SLDEBUG_ALLOC_NOT
# undef SLDEBUG_ALLOC
#endif


#include "sl_rename.h"

#include "sl_config.h"
#include "sl_config_intern.h"

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif

#include "sl_tune.h"
#include "sl_tune_intern.h"

#include "sl_deprecated.h"

#include "sl_environment.h"
#include "sl_environment_intern.h"

#include "sl_rti.h"
#include "sl_rti_intern.h"

#include "sl_elements.h"

#include "sl_pelem.h"

#include "sl_types.h"

#include "z_pack.h"

#include "sl_adds.h"

#include "sl_globals.h"

#define SL_PROTO(_f_)  _f_
#include "sl_protos.h"
#undef SL_PROTO
#undef __SL_PROTOS_H__
#define SL_PROTO(_f_)  _f_##_di
#include "sl_protos.h"
#undef SL_PROTO

#ifdef SL_USE_MPI
# define SL_PROTO(_f_)  _f_
#  include "sl_protos_mpi.h"
# undef SL_PROTO
# undef __SL_MPI_PROTOS_H__
# define SL_PROTO(_f_)  _f_##_di
#  include "sl_protos_mpi.h"
# undef SL_PROTO
#endif


#ifdef SL_USE_MPI
extern int sl_mpi_rank;
#else
 #ifdef sl_mpi_rank
 # undef sl_mpi_rank
 #endif
 #define sl_mpi_rank  -2
#endif


#endif /* __SL_COMMON_H__ */
