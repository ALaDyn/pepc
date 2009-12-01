/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_sample_ci.c
 *  timestamp: 2009-10-12 09:08:26 +0200
 *  
 */


#include "sl_common.h"


slint pepcparts_mpi_sample_ci_init(classification_info *ci, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_sample_ci_init */
{
  ci->nclasses = 0;
  ci->keys = NULL;
  ci->counts = NULL;
  ci->masks = NULL;

  ci->all_local_sizes = ci->local_lt_eq_counts = ci->all_local_lt_eq_counts = NULL;

  return 0;
}


slint pepcparts_mpi_sample_ci_duplicate(classification_info *ci_src, classification_info *ci_dup, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_sample_ci_duplicate */
{
  pepcparts_mpi_sample_ci_init(ci_dup, size, rank, comm);

  ci_dup->nclasses = ci_src->nclasses;

  if (ci_src->keys)
  {
    ci_dup->keys = sl_alloc(ci_dup->nclasses - 1, sizeof(slkey_pure_t));
    memcpy(ci_dup->keys, ci_src->keys, (ci_dup->nclasses - 1) * sizeof(slkey_pure_t));
  }

  if (ci_src->counts)
  {
    ci_dup->counts = sl_alloc(ci_dup->nclasses - 1, sizeof(slint));
    memcpy(ci_dup->counts, ci_src->counts, (ci_dup->nclasses - 1) * sizeof(slint));

  }

  return 0;
}


slint pepcparts_mpi_sample_ci_free(classification_info *ci, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_sample_ci_free */
{
  if (ci->keys) sl_free(ci->keys);
  if (ci->counts) sl_free(ci->counts);
  if (ci->masks) sl_free(ci->masks);

  ci->nclasses = 0;
  ci->keys = NULL;
  ci->counts = NULL;
  ci->masks = NULL;

  ci->all_local_sizes = ci->local_lt_eq_counts = ci->all_local_lt_eq_counts = NULL;

  return 0;
}


slint pepcparts_mpi_sample_ci_print(classification_info *ci, int size, int rank, MPI_Comm comm) /* pepcparts_sl_proto, sl_func pepcparts_mpi_sample_ci_print */
{
  slint i;

  printf("%d here: ci: size = %" sl_int_type_fmt "\n", rank, ci->nclasses);

  if (ci->keys)
  {
    printf("%d here: ci: keys = [", rank);
#ifdef key_pure_type_fmt
    for (i = 0; i < ci->nclasses - 1; i++) printf("  %" key_pure_type_fmt, ci->keys[i]);
#endif
    printf("]\n");
  }

  if (ci->counts)
  {
    printf("%d here: ci: counts = [", rank);
    for (i = 0; i < ci->nclasses - 1; i++) printf("  %" sl_int_type_fmt, ci->counts[i]);
    printf("]\n");
  }

  return 0;
}
