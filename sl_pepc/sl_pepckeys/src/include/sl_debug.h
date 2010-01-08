/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_debug.h
 *  timestamp: 2009-12-16 18:02:58 +0100
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
 #define MPI_STR    "%d: "
 #define MPI_PARAM  pepckeys_sl_mpi_rank
#else
 #define MPI_STR    "%s"
 #define MPI_PARAM  ""
#endif


#define SL_NOTICE(_format_, _args_...)           printf(MPI_STR _format_ "\n", MPI_PARAM, ##_args_)
#define SL_NOTICE_IF(_if_, _format_, _args_...)  do { if (_if_) { SL_NOTICE(_format_, ##_args_); } } while (0)

#define SL_ERROR(_format_, _args_...)            printf(MPI_STR "%s:%i:%s: " _format_ "\n", MPI_PARAM, __FILE__, __LINE__, __func__, ##_args_)

#ifdef SLDEBUG
# define SL_DEBUG(_level_, _format_, _args_... )  \
    if (_level_ <= SLDEBUG) { \
      printf(MPI_STR "%s:%i:%s " _format_ "\n", MPI_PARAM, __FILE__, __LINE__, __func__, ##_args_); \
      fflush(stdout); \
    } else do {} while (0)
#else
# define SL_DEBUG(_x_...)                        do {} while(0)
#endif

#define SL_ASSERT(_x_)                           do { if (_x_) {} else SL_DEBUG(0, "ASSERT: '%s' failed.", #_x_); } while(0)
#define SL_ASSERT_IF(_if,_, _x_)                 do { if (_if_) { SL_ASSERT(_x_); } } while (0)

#define SL_TRACE(_format_, _args_...)            SL_DEBUG(3, _format_, ##_args_)
#define SL_TRACE_IF(_if_, _format_, _args_...)   do { if (_if_) { SL_TRACE(_format_, ##_args_); } } while (0)


#endif /* __SL_DEBUG_H__ */
