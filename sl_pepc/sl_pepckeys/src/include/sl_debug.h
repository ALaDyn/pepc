/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_debug.h
 *  timestamp: 2009-11-02 16:43:48 +0100
 *  
 */


#ifndef __SL_DEBUG_H__
#define __SL_DEBUG_H__


#include <stdarg.h>
#include <stdio.h>


#ifdef DEBUG_STDERR
# define SL_DEBUG_FSTREAM  stderr
#else
# define SL_DEBUG_FSTREAM  stdout
#endif


#ifdef SL_USE_MPI
 extern int pepckeys_sl_mpi_rank;
 #define MPI_STR    "%d:"
 #define MPI_PARAM  pepckeys_sl_mpi_rank
#else
 #define MPI_STR    "%s"
 #define MPI_PARAM  ""
#endif


#define SL_NOTICE(format, args...) printf(MPI_STR format "\n", MPI_PARAM, ##args)

#define SL_ERROR(format, args...)  printf(MPI_STR "%s:%i:%s: " format "\n", MPI_PARAM, __FILE__, __LINE__, __func__, ##args)

#ifdef DEBUG
# define SL_DEBUG(level, format, args... )  \
    if (level <= DEBUG) { \
      printf(MPI_STR "%s:%i:%s " format "\n", MPI_PARAM, __FILE__, __LINE__, __func__, ##args); \
      fflush(stdout); \
    } else do {} while (0)
#else
# define SL_DEBUG(x...)  do {} while(0)
#endif

#define SL_ASSERT(x)  { if (x) {} else { SL_DEBUG(0, "assertion '%s' failed.", #x); } }

#define SL_TRACE(format, args...)  SL_DEBUG(3, format, ##args)


#endif /* __SL_DEBUG_H__ */
