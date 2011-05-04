
#ifndef __FORTRAN2C_TYPES_H__
#define __FORTRAN2C_TYPES_H__


/*#include "config.h"*/

#if defined (HAVE_INTTYPES_H) && (MPI_VERSION >= 2) && (MPI_SUBVERSION >= 2)

# include <inttypes.h>

# define FINT_TYPE_C       int32_t
# define FINT_TYPE_MPI     MPI_INT32_T
# define FINT_TYPE_FMT     PRId32

# define FINT8_TYPE_C      int64_t
# define FINT8_TYPE_MPI    MPI_INT64_T
# define FINT8_TYPE_FMT    PRId64

#else

# define FINT_TYPE_C       int
# define FINT_TYPE_MPI     MPI_INT
# define FINT_TYPE_FMT     "d"

# define FINT8_TYPE_C      long long
# define FINT8_TYPE_MPI    MPI_LONG_LONG
# define FINT8_TYPE_FMT    "lld"

#endif

#define FREAL8_TYPE_C      double
#define FREAL8_TYPE_MPI    MPI_DOUBLE
#define FREAL8_TYPE_FMT    "f"


#endif /* __FORTRAN2C_TYPES_H__ */
