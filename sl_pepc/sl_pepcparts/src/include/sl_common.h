/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_common.h
 *  timestamp: 2010-01-05 09:03:36 +0100
 *  
 */


#ifndef __SL_COMMON_H__
#define __SL_COMMON_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>


#define xmax(a,b)           (((a)>(b))?(a):(b))
#define xmin(a,b)           (((a)<(b))?(a):(b))
#define xmax3(a,b,c)        xmax(a,xmax(b,c))
#define xmin3(a,b,c)        xmin(a,xmin(b,c))
#define xminmax(_a_, _b_, _c_)  (((_b_)<(_a_))?(_a_):(((_b_)>(_c_))?(_c_):(_b_)))

#define powof2_typed(a, t)  (((t) 1) << (a))
#define powof2(a)           powof2_typed(a, slint_t)
#define iround(x)           ((int) ((x) + 0.5))
#define xround(x, t)        ((t) ((x) + 0.5))
#define idivc(x, y)         ((x) / (y))
#define idivf(x, y)         (((x) + (y) - 1) / (y))


#define SL_MOP(_mop)        do { _mop } while (0)
#define SL_NOP()            SL_MOP()


#include "sl_config.h"
#include "sl_config_intern.h"

#ifdef SL_USE_MPI
 #include <mpi.h>
#endif

#include "sl_tune.h"
#include "sl_tune_intern.h"

#include "sl_environment.h"
#include "sl_environment_intern.h"

#include "sl_rti.h"
#include "sl_rti_intern.h"

#include "sl_elements.h"

#include "sl_pelem.h"

#include "sl_types.h"

#include "sl_adds.h"

#include "sl_debug.h"

#include "sl_protos.h"
#include "sl_protos_di.h"

#ifdef SL_USE_MPI
 #include "sl_protos_mpi.h"
 #include "sl_protos_mpi_di.h"
#endif


#ifdef SL_USE_MPI
 extern MPI_Datatype pepcparts_key_mpi_datatype;
 extern MPI_Datatype pepcparts_key_pure_mpi_datatype;
 extern MPI_Datatype pepcparts_index_mpi_datatype;
 extern MPI_Datatype pepcparts_data_mpi_datatype[data_nmax + 1];
#endif


#endif /* __SL_COMMON_H__ */
