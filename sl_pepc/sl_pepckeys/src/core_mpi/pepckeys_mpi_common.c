/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_common.c
 *  timestamp: 2009-11-19 16:59:33 +0100
 *  
 */


#include "sl_common.h"


declare_ts_temp_mpi

/* sl_var pepckeys_key_mpi_datatype pepckeys_key_pure_mpi_datatype pepckeys_index_mpi_datatype pepckeys_data_mpi_datatype */
MPI_Datatype pepckeys_key_mpi_datatype = MPI_DATATYPE_NULL;
MPI_Datatype pepckeys_key_pure_mpi_datatype = MPI_DATATYPE_NULL;
MPI_Datatype pepckeys_index_mpi_datatype = MPI_DATATYPE_NULL;
MPI_Datatype pepckeys_data_mpi_datatype[data_nmax + 1] =
{
#define xelem_call_data      MPI_DATATYPE_NULL,
#define xelem_call_data_not  MPI_DATATYPE_NULL,
#include "sl_xelem_call.h"
  MPI_DATATYPE_NULL
};

/* sl_var pepckeys_sl_mpi_rank */
int pepckeys_sl_mpi_rank = -1;


slint_t pepckeys_mpi_datatypes_init() /* pepckeys_sl_proto, sl_func pepckeys_mpi_datatypes_init */
{
#define xelem_call \
  if (xelem_size_mpi > 1) \
  { \
    MPI_Type_contiguous(xelem_size_mpi, xelem_type_mpi, &xelem_mpi_datatype); \
    MPI_Type_commit(&xelem_mpi_datatype); \
\
  } else xelem_mpi_datatype = xelem_type_mpi;
#include "sl_xelem_call.h"

  if (key_pure_size_mpi > 1)
  {
    MPI_Type_contiguous(key_pure_size_mpi, key_pure_type_mpi, &pepckeys_key_pure_mpi_datatype);
    MPI_Type_commit(&pepckeys_key_pure_mpi_datatype);

  } else pepckeys_key_pure_mpi_datatype = key_pure_type_mpi;

  return 0;
}


slint_t pepckeys_mpi_datatypes_release() /* pepckeys_sl_proto, sl_func pepckeys_mpi_datatypes_release */
{
#define xelem_call \
  if (xelem_mpi_datatype != xelem_type_mpi) MPI_Type_free(&xelem_mpi_datatype);
#include "sl_xelem_call.h"

  if (pepckeys_key_pure_mpi_datatype != key_pure_type_mpi) MPI_Type_free(&pepckeys_key_pure_mpi_datatype);

  return 0;
}
