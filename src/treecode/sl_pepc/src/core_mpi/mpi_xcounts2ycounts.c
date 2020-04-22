/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_xcounts2ycounts.c
 *  
 */


/* sl_macro MX2Y_TRACE_IF */

#include "sl_common.h"


#ifndef MX2Y_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MX2Y_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MX2Y_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_xcounts2ycounts_all2all(int *xcounts, int *ycounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_xcounts2ycounts_all2all */
{
  MPI_Alltoall(xcounts, 1, MPI_INT, ycounts, 1, MPI_INT, comm);
  
  return 0;
}


slint_t mpi_xcounts2ycounts_sparse(int *xcounts, int *ycounts, slint_t ytotal, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_xcounts2ycounts_sparse */
{
#define MPI_XCOUNTS2YCOUNTS_SPARSE_TAG  1

  slint_t i;

  MPI_Request req;
  MPI_Status status;

  int ycount;


  for (i = 0; i < size; ++i) ycounts[i] = 0;

  for (i = 0; i < size; ++i)
  if (xcounts[i] != 0)
  {
    MPI_Isend(&xcounts[i], 1, MPI_INT, i, MPI_XCOUNTS2YCOUNTS_SPARSE_TAG, comm, &req);
    MPI_Request_free(&req);
  }
  
  while (ytotal > 0)
  {
    MPI_Recv(&ycount, 1, MPI_INT, MPI_ANY_SOURCE, MPI_XCOUNTS2YCOUNTS_SPARSE_TAG, comm, &status);
    
    ytotal -= ycount;
    ycounts[status.MPI_SOURCE] = ycount;
  }
  
  return 0;

#undef MPI_XCOUNTS2YCOUNTS_SPARSE_TAG
}


slint_t mpi_xcounts2ycounts_grouped(int *xcounts, slint_t nxcounts, int *ycounts, MPI_Comm group_comm, MPI_Comm master_comm, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_xcounts2ycounts_grouped */
{
  slint_t i;
  int group_size, group_rank;

  int *all_group_xcounts = NULL;
  MPI_Datatype rdt, rdtr;

  int group_sizes[nxcounts];

  int all_xcounts[size];

  int mscounts[nxcounts], msdispls[nxcounts];
  int mrcounts[nxcounts], mrdispls[nxcounts];

  slint_t dist_total, dist_from, dist_current;

  int my_dist_plan_size = 0;
  slint_t my_dist_plan[size * 2];
  slint_t dist_plan_size, dist_plan[size * 2];
  
  MPI_Request req;
  MPI_Status status;

  
  if (group_comm != MPI_COMM_NULL)
  {
    MPI_Comm_size(group_comm, &group_size);
    MPI_Comm_rank(group_comm, &group_rank);

  } else
  {
    group_size = 1;
    group_rank = 0;
  }

  if (group_rank == 0)
  {
    /* init space and receive type to gather xcounts (transposed) at group-root */
    all_group_xcounts = z_alloca(nxcounts * group_size, sizeof(int));

    MPI_Type_vector(nxcounts, 1, group_size, MPI_INT, &rdt);
    MPI_Type_create_resized(rdt, 0, sizeof(int), &rdtr);
    MPI_Type_commit(&rdtr);

  } else rdtr = MPI_DATATYPE_NULL; /* workaround for an OpenMPI bug (recvtype on non-root process is used in gather routine "intra_binomial") */

  /* gather xcounts (transposed) at group-root */
  MPI_Gather(xcounts, nxcounts, MPI_INT, all_group_xcounts, 1, rdtr, 0, group_comm);
  // int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, 
  //                void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)

  if (group_rank == 0)
  {
    MPI_Type_free(&rdt);
    MPI_Type_free(&rdtr);

/*    printf("%d: all_group_scounts = ", rank);
    for (i = 0; i < nparts * group_size; ++i) printf("%d  ", all_group_scounts[i]);
    printf("\n");*/

    /* distribute group sizes */
    MPI_Allgather(&group_size, 1, MPI_INT, group_sizes, 1, MPI_INT, master_comm);

    /* re-distribute xcounts among the masters */
    for (i = 0; i < nxcounts; ++i)
    {
      mscounts[i] = group_size;
      msdispls[i] = (i == 0)?0:(msdispls[i - 1] + mscounts[i - 1]);
      mrcounts[i] = group_sizes[i];
      mrdispls[i] = (i == 0)?0:(mrdispls[i - 1] + mrcounts[i - 1]);
    }

    // Superficially it looks lke all_group_xcounts could not be allocated for all ranks within master_comm. This of course
    // depends on the communicatiors used. Add a check to be sure...
    if ( all_group_xcounts == NULL ) printf("BARRFFFFF");
    MPI_Alltoallv(all_group_xcounts, mscounts, msdispls, MPI_INT, all_xcounts, mrcounts, mrdispls, MPI_INT, master_comm);
    // int MPI_Alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
    //                   const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)

    z_freea(all_group_xcounts);

/*    printf("%d: all_scounts = ", rank);
    for (i = 0; i < size; ++i) printf("%d  ", all_scounts[i]);
    printf("\n");*/
  }

  /* distribute xcounts to the group members */
  if (group_rank == 0)
  {
    /* determine total number of elements to distribute in this group */
    dist_total = 0;
    for (i = 0; i < size; ++i) dist_total += all_xcounts[i];

    /* determine elements for every group member */
    dist_from = 0;
    for (i = 0; i < group_size; ++i)
    {
      /* distribute equally */
      dist_current = (slint_t) ceil((double) dist_total / (group_size - i));

      dist_total -= dist_current;

/*      printf("%d: %" sl_int_type_fmt " receives %" sl_int_type_fmt "\n", rank, i, dist_current);*/

      /* determine how the current elements are distributed (make a distribution plan) */
      dist_plan_size = 0;
      while (dist_current > 0 && dist_from < size)
      {
        if (all_xcounts[dist_from] > 0)
        {
          dist_plan[dist_plan_size + 0] = z_min(dist_current, all_xcounts[dist_from]);
          dist_plan[dist_plan_size + 1] = dist_from;

          all_xcounts[dist_from] -= dist_plan[dist_plan_size + 0];
          dist_current -= dist_plan[dist_plan_size + 0];
          
          dist_plan_size += 2;
          
        } else
        {
          ++dist_from;
          continue;
        }
      }

/*      printf("%d: plan for %" sl_int_type_fmt ": ", rank, i);
      for (j = 0; j < dist_plan_size; j += 2) printf("(%" sl_int_type_fmt ",%" sl_int_type_fmt ")  ", dist_plan[j], dist_plan[j + 1]);
      printf("\n");*/

      /* send (or copy) a distribution plan */
      if (i == group_rank)
      {
        memcpy(my_dist_plan, dist_plan, sizeof(slint_t) * dist_plan_size);
        my_dist_plan_size = dist_plan_size;

      } else
      {
        MPI_Isend(dist_plan, dist_plan_size, int_mpi_datatype, i, 0, group_comm, &req);
        MPI_Request_free(&req);
      }

/*      printf("%d: my plan: ", rank);
      for (j = 0; j < my_dist_plan_size; j += 2) printf("(%" sl_int_type_fmt ",%" sl_int_type_fmt ")  ", my_dist_plan[j], my_dist_plan[j + 1]);
      printf("\n");*/
    }

  } else
  {
    /* receive the distribution plan */
    MPI_Recv(my_dist_plan, size * 2, int_mpi_datatype, 0, 0, group_comm, &status);
    MPI_Get_count(&status, int_mpi_datatype, &my_dist_plan_size);
  }

  /* make ycounts from the distribution plan */
  for (i = 0; i < size; ++i) ycounts[i] = 0;

  for (i = 0; i < my_dist_plan_size; i += 2) ycounts[my_dist_plan[i + 1]] = my_dist_plan[i + 0];

  return 0;
}


slint_t mpi_subxdispls2ycounts(slint_t nsubs, int *sub_xdispls, slint_t *sub_sources, slint_t *sub_sizes, MPI_Comm sub_comm, int sub_size, int *ycounts, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_subxdispls2ycounts */
{
  slint_t i, j, k;
  slint_t all_nsubs[sub_size];

  int your_sub_plans[sub_size * nsubs * 2], my_sub_plan[2 * 2 * size];
  int sub_plan_scounts[sub_size], sub_plan_sdispls[sub_size], sub_plan_rcounts[sub_size], sub_plan_rdispls[sub_size];


  Z_TRACE_IF(MX2Y_TRACE_IF, "nsubs = %" slint_fmt "", nsubs);
  Z_TRACE_ARRAY_IF(MX2Y_TRACE_IF, i, nsubs * sub_size, "%d  ", sub_xdispls[i], "sub_xdispls = ");
  Z_TRACE_ARRAY_IF(MX2Y_TRACE_IF, i, nsubs, "%" slint_fmt "  ", sub_sources[i], "sub_sources = ");
  Z_TRACE_ARRAY_IF(MX2Y_TRACE_IF, i, nsubs, "%" slint_fmt "  ", sub_sizes[i], "sub_sizes = ");

  MPI_Allgather(&nsubs, 1, int_mpi_datatype, all_nsubs, 1, int_mpi_datatype, sub_comm);
  
  for (i = 0; i < sub_size; ++i)
  {
    for (j = 0; j < nsubs; ++j)
    {
      your_sub_plans[(i * nsubs + j) * 2 + 0] = sub_sources[j];
      your_sub_plans[(i * nsubs + j) * 2 + 1] = ((i + 1 < sub_size)?sub_xdispls[(i + 1) * nsubs + j]:sub_sizes[j]) - sub_xdispls[i * nsubs + j];
    }
    sub_plan_scounts[i] = nsubs * 2;
    sub_plan_sdispls[i] = i * nsubs * 2;
    
    sub_plan_rcounts[i] = all_nsubs[i] * 2;
    sub_plan_rdispls[i] = (i == 0)?0:(sub_plan_rdispls[i - 1] + sub_plan_rcounts[i - 1]);
  }
  
  MPI_Alltoallv(your_sub_plans, sub_plan_scounts, sub_plan_sdispls, MPI_INT, my_sub_plan, sub_plan_rcounts, sub_plan_rdispls, MPI_INT, sub_comm);

  for (i = 0; i < size; ++i) ycounts[i] = 0;
  
  k = 0;
  for (i = 0; i < sub_size; i++)
  for (j = 0; j < all_nsubs[i]; ++j)
  {
    ycounts[my_sub_plan[k + 0]] += my_sub_plan[k + 1];
    k += 2;
  }
  
  return 0;
}
