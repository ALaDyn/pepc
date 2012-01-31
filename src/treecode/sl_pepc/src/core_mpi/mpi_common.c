/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_common.c
 *  
 */


#include "sl_common.h"

#ifdef HAVE_SPI_KERNEL_INTERFACE_H
# include <spi/kernel_interface.h>
#endif
#ifdef HAVE_COMMON_BGP_PERSONALITY_H
# include <common/bgp_personality.h>
#endif
#ifdef HAVE_COMMON_BGP_PERSONALITY_INLINES_H
# include <common/bgp_personality_inlines.h>
#endif


MPI_Datatype int_mpi_datatype = MPI_DATATYPE_NULL;     /* sl_global MPIONLY_BEGIN, sl_var int_mpi_datatype */
MPI_Datatype key_mpi_datatype = MPI_DATATYPE_NULL;     /* sl_global, sl_var key_mpi_datatype */
MPI_Datatype pkey_mpi_datatype = MPI_DATATYPE_NULL;    /* sl_global, sl_var pkey_mpi_datatype */
MPI_Datatype pwkey_mpi_datatype = MPI_DATATYPE_NULL;   /* sl_global, sl_var pwkey_mpi_datatype */
MPI_Datatype index_mpi_datatype = MPI_DATATYPE_NULL;   /* sl_global, sl_var index_mpi_datatype */
MPI_Datatype weight_mpi_datatype = MPI_DATATYPE_NULL;  /* sl_global, sl_var weight_mpi_datatype */
MPI_Datatype data_mpi_datatype[data_nmax + 1] =        /* sl_global MPIONLY_END, sl_var data_mpi_datatype */
{
#define xelem_call_data      MPI_DATATYPE_NULL,
#define xelem_call_data_not  MPI_DATATYPE_NULL,
#include "sl_xelem_call.h"
  MPI_DATATYPE_NULL
};

int sl_mpi_rank = -2;  /* sl_global MPIONLY, sl_var sl_mpi_rank */


slint_t mpi_datatypes_init() /* sl_proto, sl_func mpi_datatypes_init */
{
  slpwkey_t pwk;
  int pwk_blenghts[2] = { key_pure_size_mpi, 1 };
  MPI_Datatype pwk_types[2] = { key_pure_type_mpi, sl_weight_type_mpi };
  MPI_Aint pwk_displs[2], base;

  /* intern integer */
  if (sl_int_size_mpi > 1)
  {
    MPI_Type_contiguous(sl_int_size_mpi, sl_int_type_mpi, &int_mpi_datatype);
    MPI_Type_commit(&int_mpi_datatype);

  } else int_mpi_datatype = sl_int_type_mpi;

#define xelem_call \
  if (xelem_size_mpi > 1) \
  { \
    MPI_Type_contiguous(xelem_size_mpi, xelem_type_mpi, &xelem_mpi_datatype); \
    MPI_Type_commit(&xelem_mpi_datatype); \
\
  } else xelem_mpi_datatype = xelem_type_mpi;
#include "sl_xelem_call.h"

  /* pure key */
  if (key_pure_size_mpi > 1)
  {
    MPI_Type_contiguous(key_pure_size_mpi, key_pure_type_mpi, &pkey_mpi_datatype);
    MPI_Type_commit(&pkey_mpi_datatype);

  } else pkey_mpi_datatype = key_pure_type_mpi;

  /* pure weighted key */
  MPI_Get_address(&pwk, &base);
  MPI_Get_address(&pwk.pkey, &pwk_displs[0]); pwk_displs[0] -= base;
  MPI_Get_address(&pwk.weight, &pwk_displs[1]); pwk_displs[1] -= base;

  MPI_Type_create_struct(2, pwk_blenghts, pwk_displs, pwk_types, &pwkey_mpi_datatype);
  MPI_Type_commit(&pwkey_mpi_datatype);

  /* weight (intern) */
  if (sl_weight_size_mpi > 1)
  {
    MPI_Type_contiguous(sl_weight_size_mpi, sl_weight_type_mpi, &weight_mpi_datatype);
    MPI_Type_commit(&weight_mpi_datatype);

  } else weight_mpi_datatype = sl_weight_type_mpi;

  return 0;
}


slint_t mpi_datatypes_release() /* sl_proto, sl_func mpi_datatypes_release */
{
  if (int_mpi_datatype != sl_int_type_mpi) MPI_Type_free(&int_mpi_datatype);

#define xelem_call \
  if (xelem_mpi_datatype != xelem_type_mpi) MPI_Type_free(&xelem_mpi_datatype);
#include "sl_xelem_call.h"

  if (pkey_mpi_datatype != key_pure_type_mpi) MPI_Type_free(&pkey_mpi_datatype);

  MPI_Type_free(&pwkey_mpi_datatype);

  if (weight_mpi_datatype != sl_weight_type_mpi) MPI_Type_free(&weight_mpi_datatype);

  return 0;
}


slint_t mpi_get_grid_properties(slint_t ndims, slint_t *dims, slint_t *pos, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_get_grid_properties */
{
  const slint_t max_ndims = 4;

  slint_t i, _dims[max_ndims], _pos[max_ndims], _rdims[max_ndims], _rpos[max_ndims];

  slint_t remap[max_ndims];

#ifdef HAVE__BGP_PERSONALITY_T

  _BGP_Personality_t personality;

  Kernel_GetPersonality(&personality, sizeof(personality));

  _dims[0] = personality.Network_Config.Xnodes;
  _dims[1] = personality.Network_Config.Ynodes;
  _dims[2] = personality.Network_Config.Znodes;

  _pos[0] = personality.Network_Config.Xcoord;
  _pos[1] = personality.Network_Config.Ycoord;
  _pos[2] = personality.Network_Config.Zcoord;

  _pos[3] = Kernel_PhysicalProcessorID();

  switch (personality.Kernel_Config.ProcessConfig)
  {
    case _BGP_PERS_PROCESSCONFIG_SMP:
      _dims[3] = 1;
      break;
    case _BGP_PERS_PROCESSCONFIG_VNM:
      _dims[3] = 4;
      break;
    case _BGP_PERS_PROCESSCONFIG_2x2:
      _dims[3] = 2;
      break;
    default:
      _dims[3] = 1;
      break;
  }

#else

/*#define STATIC*/

# ifndef STATIC

  int mpi_dims[4] = { 0, 0, 0, 0 };

  MPI_Dims_create(size, ndims, mpi_dims);

  for (i = 0; i < ndims; ++i) _dims[i] = mpi_dims[ndims - i - 1];
  for (i = ndims; i < 4; ++i) _dims[i] = 1;

# else
  switch (size)
  {
    case 4:   _dims[0] = 2;    _dims[1] = 2; _dims[2] = 1; _dims[3] = 1; break;
    case 8:   _dims[0] = 2;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 1; break;
/*    case 16:  _dims[0] = 2;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 2; break;*/
    case 16:  _dims[0] = 4;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 1; break;
    case 32:  _dims[0] = 4;    _dims[1] = 4; _dims[2] = 2; _dims[3] = 1; break;
/*    case 32:  _dims[0] = 8;    _dims[1] = 4; _dims[2] = 1; _dims[3] = 1; break;*/
    case 64:  _dims[0] = 4;    _dims[1] = 4; _dims[2] = 4; _dims[3] = 1; break;
/*    case 64:  _dims[0] = 8;    _dims[1] = 8; _dims[2] = 1; _dims[3] = 1; break;*/
    case 128: _dims[0] = 8;    _dims[1] = 4; _dims[2] = 4; _dims[3] = 1; break;
    case 512: _dims[0] = 8;    _dims[1] = 8; _dims[2] = 8; _dims[3] = 1; break;
    default:  _dims[0] = size; _dims[1] = 1; _dims[2] = 1; _dims[3] = 1; break;
  }
# endif

#undef STATIC

  _pos[3] = (rank / (1))                                  % _dims[3];
  _pos[2] = (rank / (1 * _dims[3]))                       % _dims[2];
  _pos[1] = (rank / (1 * _dims[3] * _dims[2]))            % _dims[1];
  _pos[0] = (rank / (1 * _dims[3] * _dims[2] * _dims[1])) % _dims[0];

#endif

  for (i = 0; i < max_ndims; ++i) remap[i] = i;

  /* remap dimensions */
  for (i = 0; i < max_ndims; ++i)
  {
    _rdims[i] = _dims[remap[i]];
    _rpos[i] = _pos[remap[i]];
  }

  /* reduce dimensions (lowest significant first) */
  for (i = max_ndims - 1; i >= ndims; --i)
  {
    _rpos[i - 1] = _rpos[i - 1] * _rdims[i] + _rpos[i];
    _rdims[i - 1] *= _rdims[i];

    _rpos[i] = 0;
    _rdims[i] = 1;
  }

  for (i = 0; i < max_ndims; ++i)
  {
    dims[i] = _rdims[i];
    pos[i] = _rpos[i];
  }

/*  printf("%d: <%" slint_fmt ",%" slint_fmt ",%" slint_fmt ",%" slint_fmt "> of <%" slint_fmt "x%" slint_fmt "x%" slint_fmt"x%" slint_fmt ">\n", rank, pos[0], pos[1], pos[2], pos[3], dims[0], dims[1], dims[2], dims[3]);*/

  return 0;
}


slint_t mpi_get_grid(slint_t ndims, slint_t *dims, slint_t *pos, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_get_grid */
{
  const slint_t max_ndims = 4;

  slint_t i, _dims[max_ndims], _pos[max_ndims], _rdims[max_ndims], _rpos[max_ndims];

  slint_t remap[max_ndims];

#ifdef _HAVE__BGP_PERSONALITY_T_

  _BGP_Personality_t personality;

  Kernel_GetPersonality(&personality, sizeof(personality));

  _dims[0] = personality.Network_Config.Xnodes;
  _dims[1] = personality.Network_Config.Ynodes;
  _dims[2] = personality.Network_Config.Znodes;

  _pos[0] = personality.Network_Config.Xcoord;
  _pos[1] = personality.Network_Config.Ycoord;
  _pos[2] = personality.Network_Config.Zcoord;

  _pos[3] = Kernel_PhysicalProcessorID();

  switch (personality.Kernel_Config.ProcessConfig)
  {
    case _BGP_PERS_PROCESSCONFIG_SMP:
      _dims[3] = 1;
      break;
    case _BGP_PERS_PROCESSCONFIG_VNM:
      _dims[3] = 4;
      break;
    case _BGP_PERS_PROCESSCONFIG_2x2:
      _dims[3] = 2;
      break;
    default:
      _dims[3] = 1;
      break;
  }

#ifdef BGP_TXYZ_MAPPING
  /* TXYZ mapping to ZYXT */
/*  remap[0] = 2; remap[1] = 1; remap[2] = 0; remap[3] = 3;*/
#else
  /* XYZT mapping to ZYXT (SMP mode, T = 1) or TZYX (non-SMP mode) */
  if (_dims[3] > 1) { remap[0] = 3; remap[1] = 2; remap[2] = 1; remap[3] = 0; }
  else { remap[0] = 2; remap[1] = 1; remap[2] = 0; remap[3] = 3; }
#endif

  for (i = 4; i < max_ndims; ++i) remap[i] = i;

#else

  for (i = 0; i < max_ndims; ++i) remap[i] = i;

/*#define STATIC*/

# ifndef STATIC

  int mpi_dims[4] = { 0, 0, 0, 0 };

  MPI_Dims_create(size, ndims, mpi_dims);

  for (i = 0; i < ndims; ++i) _dims[i] = mpi_dims[ndims - i - 1];
  for (i = ndims; i < 4; ++i) _dims[i] = 1;

# else
  switch (size)
  {
    case 4:   _dims[0] = 2;    _dims[1] = 2; _dims[2] = 1; _dims[3] = 1; break;
    case 8:   _dims[0] = 2;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 1; break;
/*    case 16:  _dims[0] = 2;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 2; break;*/
    case 16:  _dims[0] = 4;    _dims[1] = 2; _dims[2] = 2; _dims[3] = 1; break;
    case 32:  _dims[0] = 4;    _dims[1] = 4; _dims[2] = 2; _dims[3] = 1; break;
/*    case 32:  _dims[0] = 8;    _dims[1] = 4; _dims[2] = 1; _dims[3] = 1; break;*/
    case 64:  _dims[0] = 4;    _dims[1] = 4; _dims[2] = 4; _dims[3] = 1; break;
/*    case 64:  _dims[0] = 8;    _dims[1] = 8; _dims[2] = 1; _dims[3] = 1; break;*/
    case 128: _dims[0] = 8;    _dims[1] = 4; _dims[2] = 4; _dims[3] = 1; break;
    case 512: _dims[0] = 8;    _dims[1] = 8; _dims[2] = 8; _dims[3] = 1; break;
    default:  _dims[0] = size; _dims[1] = 1; _dims[2] = 1; _dims[3] = 1; break;
  }
# endif

#undef STATIC

  _pos[3] = (rank / (1))                                  % _dims[3];
  _pos[2] = (rank / (1 * _dims[3]))                       % _dims[2];
  _pos[1] = (rank / (1 * _dims[3] * _dims[2]))            % _dims[1];
  _pos[0] = (rank / (1 * _dims[3] * _dims[2] * _dims[1])) % _dims[0];

#endif

  /* remap dimensions */
  for (i = 0; i < max_ndims; ++i)
  {
    _rdims[i] = _dims[remap[i]];
    _rpos[i] = _pos[remap[i]];
  }

  /* reduce dimensions (lowest significant first) */
  for (i = max_ndims - 1; i >= ndims; --i)
  {
    _rpos[i - 1] = _rpos[i - 1] * _rdims[i] + _rpos[i];
    _rdims[i - 1] *= _rdims[i];

    _rpos[i] = 0;
    _rdims[i] = 1;
  }

  for (i = 0; i < max_ndims; ++i)
  {
    dims[i] = _rdims[i];
    pos[i] = _rpos[i];
  }

/*  printf("%d: <%" slint_fmt ",%" slint_fmt ",%" slint_fmt ",%" slint_fmt "> of <%" slint_fmt "x%" slint_fmt "x%" slint_fmt"x%" slint_fmt ">\n", rank, pos[0], pos[1], pos[2], pos[3], dims[0], dims[1], dims[2], dims[3]);*/

  return 0;
}


slint_t mpi_subgroups_create(slint_t nsubgroups, MPI_Comm *sub_comms, int *sub_sizes, int *sub_ranks, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_subgroups_create */
{
  slint_t i, nrtsize;


  nrtsize = (slint_t) pow((double) size, 1.0 / (double) nsubgroups);

  for (i = 0; i < nsubgroups; ++i)
  {
    if (i == 0) MPI_Comm_dup(comm, &sub_comms[i]);
    else MPI_Comm_split(sub_comms[i - 1], sub_ranks[i - 1] * nrtsize / sub_sizes[i - 1], sub_ranks[i - 1], &sub_comms[i]);

    MPI_Comm_size(sub_comms[i], &sub_sizes[i]);
    MPI_Comm_rank(sub_comms[i], &sub_ranks[i]);
  }

  return 0;
}


slint_t mpi_subgroups_delete(slint_t nsubgroups, MPI_Comm *sub_comms, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_subgroups_delete */
{
  slint_t i;

  for (i = 0; i < nsubgroups; ++i) MPI_Comm_free(&sub_comms[i]);

  return 0;
}


/* sl_macro MC_REDUCEBCAST_THRESHOLD */
/* sl_macro MC_REDUCEBCAST_ROOT */

#if !defined(MC_REDUCEBCAST_THRESHOLD) && defined(GLOBAL_REDUCEBCAST_THRESHOLD)
# define MC_REDUCEBCAST_THRESHOLD  GLOBAL_REDUCEBCAST_THRESHOLD
#endif

#ifndef MC_REDUCEBCAST_ROOT
# define MC_REDUCEBCAST_ROOT  0
#endif

int sl_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, int size, int rank) /* sl_proto, sl_func sl_MPI_Allreduce */
{
#ifdef MC_REDUCEBCAST_THRESHOLD
  if (size >= MC_REDUCEBCAST_THRESHOLD)
  {
    MPI_Reduce(sendbuf, recvbuf, count, datatype, op, MC_REDUCEBCAST_ROOT, comm);
    MPI_Bcast(recvbuf, count, datatype, MC_REDUCEBCAST_ROOT, comm);
    return 0;
  }
#endif
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}
