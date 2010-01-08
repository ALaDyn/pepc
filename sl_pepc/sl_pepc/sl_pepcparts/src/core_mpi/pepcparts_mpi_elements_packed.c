/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_elements_packed.c
 *  timestamp: 2009-12-03 11:25:06 +0100
 *  
 */


#include "sl_common.h"


slint_t pepcparts_mpi_elements_packed_datatype_create(MPI_Datatype *pdt, slint_t structured) /* pepcparts_sl_proto, sl_func pepcparts_mpi_elements_packed_datatype_create */
{
  packed_element_t pe;
  
  int blengths[elem_n];
  MPI_Aint displs[elem_n];
  MPI_Datatype types[elem_n];

  int i = 0;
  
  MPI_Aint base;

  if (structured)
  {
    MPI_Get_address(&pe, &base);

#define xelem_index_not
#define xelem_call \
  blengths[i] = xelem_size_mpi; \
  MPI_Get_address(&pe.xelem_name_packed, &displs[i]); displs[i] -= base; \
  types[i] = xelem_type_mpi; \
  i++;
#include "sl_xelem_call.h"

    MPI_Type_create_struct(1 + data_n, blengths, displs, types, pdt);

  } else MPI_Type_contiguous(sizeof(packed_element_t), MPI_BYTE, pdt);

  MPI_Type_commit(pdt);

  return 0;
}


slint_t pepcparts_mpi_elements_packed_datatype_destroy(MPI_Datatype *pdt) /* pepcparts_sl_proto, sl_func pepcparts_mpi_elements_packed_datatype_destroy */
{
  MPI_Type_free(pdt);

  return 0;
}
