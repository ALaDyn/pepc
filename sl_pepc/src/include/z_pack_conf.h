/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/z_pack_conf.h
 *  timestamp: 2011-03-06 21:59:31 +0100
 *  
 */


#ifndef __Z_PACK_CONF_H__
#define __Z_PACK_CONF_H__


#include "sl_config.h"
#include "sl_config_intern.h"


#ifdef SL_PREFIX
# define Z_PREFIX  SL_PREFIX
#endif

typedef sl_int_type_c z_int_t;
#define z_int_fmt  sl_int_type_fmt

#ifdef SL_USE_MPI
# define Z_PACK_MPI
# define Z_PACK_MPI_RANK  sl_mpi_rank
#else
# define Z_PACK_MPI_RANK  -2
#endif


#define Z_PACK_NUMERIC

#define Z_PACK_ALLOC

#ifdef SLDEBUG_ALLOC
# define Z_ALLOC_DEBUG  SLDEBUG_ALLOC
#endif

#define Z_PACK_DEBUG

#ifdef SLDEBUG_OUTPUT
# define Z_DEBUG_LEVEL  SLDEBUG_OUTPUT
#endif

#define Z_PACK_TIME

#define Z_PACK_RANDOM

#define Z_PACK_DIGEST

#define Z_PACK_CRC32

typedef unsigned int z_crc32_t;
#define z_crc32_fmt  "X"


#endif /* __Z_PACK_CONF_H__ */
