/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/mpi_select_sample.c
 *  timestamp: 2011-02-10 21:20:50 +0100
 *  
 */


/* sl_macro MSS_ROOT */
/* sl_macro MSS_TRACE_IF */


#include "sl_common.h"


/*#define sl_pivot_equal(_n_, _i_, _r_)  ((slint_t) (((((double) (_i_) + 1) * (_n_)) / ((_r_) + 1)) + 0.5))*/
#define sl_pivot_equal(_n_, _i_, _r_)  (((_i_) * (_n_)) / ((_r_) + 1))

/*#define MSS_ROOT*/

#ifdef MSS_ROOT
int mss_root = -1;  /* sl_global, sl_var mss_root */
#endif

#ifndef MSS_TRACE_IF
# ifdef GLOBAL_TRACE_IF
#  define MSS_TRACE_IF  GLOBAL_TRACE_IF
# else
#  define MSS_TRACE_IF  (sl_mpi_rank == -1)
# endif
#endif


slint_t mpi_select_sample_regular(elements_t *s, slint_t nparts, partcond_t *pconds, slint_t nsamples, splitter_t *sp, int size, int rank, MPI_Comm comm) /* sl_proto, sl_func mpi_select_sample_regular */
{
  slint_t i, j;
#ifdef elem_weight
  slweight_t w, wi, wold;
#endif

  const slint_t nslocal = nsamples;
#ifdef MSS_ROOT
  const slint_t nsglobal = nslocal * size;
#endif

  const slint_t nsplitter = nparts - 1;
  slint_t splitter_skip = 0;

#ifdef elem_weight
  slpwkey_t lwskeys[nslocal];
# ifdef MSS_ROOT
  slpwkey_t gwskeys[nsglobal];
# endif
#endif
  slpkey_t lskeys[nslocal];
#ifdef MSS_ROOT
  slpkey_t gskeys[nsglobal];
#endif

  slpkey_t spkeys[nsplitter + 1];
#ifdef elem_weight
  int nlspkeys;
  int rcounts[size], rdispls[size];
  slpkey_t lspkeys[nsplitter];
#endif
  slpkey_t lspkey;

  slint_t lgcounts[2];
#ifdef elem_weight
  slweight_t lgweights[2];
#endif

  partcond_intern_t pci[nparts];

  elements_t gs, e;

#ifdef elem_weight
  slint_t doweights;
#else
# define doweights  0
#endif


  sp->displs[0] = 0;

  if (nparts < 2) return 0;

#ifdef elem_weight
  doweights = ((pconds->pcm & (SLPC_WEIGHTS_MM|SLPC_WEIGHTS_LH)) != 0);
#endif

#ifdef elem_weight
  if (doweights) mpi_elements_get_counts_and_weights(s, 1, lgcounts, lgweights, -1, size, rank, comm);
  else
#endif
    mpi_elements_get_counts(s, &lgcounts[0], &lgcounts[1], -1, size, rank, comm);

  init_partconds_intern(nparts, pci, pconds, nparts, lgcounts[1], elem_weight_ifelse(doweights?lgweights[1]:0, 0));

  Z_TRACE_IF(MSS_TRACE_IF, "counts: %" slint_fmt " / %" slint_fmt, lgcounts[0], lgcounts[1]);
#ifdef elem_weight
  if (doweights)
    Z_TRACE_IF(MSS_TRACE_IF, "weights: %" slweight_fmt " / %" slweight_fmt "", lgweights[0], lgweights[1]);
#endif

  Z_TRACE_IF(MSS_TRACE_IF, "local samples: %" slint_fmt, nslocal);

  /* select local samples */
#ifdef elem_weight
  if (doweights)
  {
    j = 0;
    w = wold = 0;

    for (i = 0; i < nslocal; ++i)
    {
      wi = (i + 1) * lgweights[0] / (nslocal + 1);

      while (w < wi && j < lgcounts[0])
      {
        w += elem_weight(s, j);
        ++j;
      }
    
      if (j < lgcounts[0]) lwskeys[i].pkey = *key_get_pure(elem_key_at(s, j));
      else lwskeys[i].pkey = *key_get_pure(elem_key_at(s, j - 1)) + 1;

      lwskeys[i].weight = w - wold;
    
      wold = w;

      Z_TRACE_IF(MSS_TRACE_IF, "local sample %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, j, lwskeys[i].pkey, lwskeys[i].weight);
    }

  } else
#endif
    for (i = 0; i < nslocal; ++i)
  {
    j = sl_pivot_equal(s->size, i + 1, nslocal);
  
    lskeys[i] = *key_get_pure(elem_key_at(s, j));

    Z_TRACE_IF(MSS_TRACE_IF, "local sample %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt, i, j, lskeys[i]);
  }

#ifdef MSS_ROOT
  /* with root-process (p^2 scaling!!!) */
  if (mss_root >= 0)
  {
    /* gather local samples at root */
#ifdef elem_weight
    if (doweights) MPI_Gather(lwskeys, nslocal, pwkey_mpi_datatype, gwskeys, nslocal, pwkey_mpi_datatype, mss_root, comm);
    else
#endif
      MPI_Gather(lskeys, nslocal, pkey_mpi_datatype, gskeys, nslocal, pkey_mpi_datatype, mss_root, comm);

    if (rank == mss_root)
    {
      Z_TRACE_IF(MSS_TRACE_IF, "global samples: %" slint_fmt, nsglobal);

      /* prepare global samples */
      elements_alloc(&gs, nsglobal, SLCM_ALL);

      gs.size = nsglobal;

      Z_TRACE_IF(MSS_TRACE_IF, "unsorted global samples");

#ifdef elem_weight
      if (doweights)
      {
        for (i = 0; i < nsglobal; ++i)
        {
          key_set_pure(elem_key_at(&gs, i), gwskeys[i].pkey);
          elem_weight(&gs, i) = gwskeys[i].weight;

          Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
        }

      } else
#endif
        for (i = 0; i < nsglobal; ++i)
      {
        key_set_pure(elem_key_at(&gs, i), gskeys[i]);

        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
      }

      /* sort global samples */
      sort_radix(&gs, NULL, -1, -1, -1);

      Z_TRACE_IF(MSS_TRACE_IF, "sorted global samples");

#ifdef elem_weight
      if (doweights)
      {
        for (i = 0; i < nsglobal; ++i)
        {
          Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
        }

      } else
#endif
        for (i = 0; i < nsglobal; ++i)
      {
        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
      }

      Z_TRACE_IF(MSS_TRACE_IF, "splitters: %" slint_fmt, nsplitter);

      /* select splitters from global samples */
#ifdef elem_weight
      if (doweights)
      {
        j = 0;
        w = 0;

        for (i = 0; i < nsplitter; ++i)
        {
          wi = (i + 1) * lgweights[1] / (nsplitter + 1) - (lgweights[1] / (nsplitter + 1) / 2);

          while (w < wi && j < nsglobal)
          {
            w += elem_weight(&gs, j);
            ++j;
          }

          if (j > 0 && (wi - w + elem_weight(&gs, j - 1)) <= (w - wi))
          {
            w -= elem_weight(&gs, j - 1);
            --j;
          }
          spkeys[i] = *key_get_pure(elem_key_at(&gs, j - 1));

          Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt " (target: %" slweight_fmt ")", i, j - 1, spkeys[i], w, wi);
        }
    
      } else
#endif
        for (i = 0; i < nsplitter; ++i)
      {
        j = sl_pivot_equal(gs.size, i + 1, nsplitter);
    
        spkeys[i] = *key_get_pure(elem_key_at(&gs, j));

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt, i, j, spkeys[i]);
      }

      elements_free(&gs);
    }

    /* broadcast splitters */
    MPI_Bcast(&spkeys, nsplitter, pkey_mpi_datatype, mss_root, comm);

  } else
#endif
  {
    /* without root-process */
    elements_t xs;

    elements_alloc(&gs, nslocal, SLCM_ALL);
    elements_alloc(&xs, nslocal / 2 + 1, SLCM_ALL);

    gs.size = nslocal;

    Z_TRACE_IF(MSS_TRACE_IF, "unsorted global samples");

#ifdef elem_weight
    if (doweights)
    {
      for (i = 0; i < gs.size; ++i)
      {
        key_set_pure(elem_key_at(&gs, i), lwskeys[i].pkey);
        elem_weight(&gs, i) = lwskeys[i].weight;

        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
      }

    } else
#endif
      for (i = 0; i < gs.size; ++i)
    {
      key_set_pure(elem_key_at(&gs, i), lskeys[i]);

      Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
    }

    /* parallel sort samples */
    mpi_mergek(&gs, sn_batcher, NULL, merge2_basic_auto_01_x, &xs, size, rank, comm);

    Z_TRACE_IF(MSS_TRACE_IF, "sorted global samples");

#ifdef elem_weight
    if (doweights)
    {
      for (i = 0; i < gs.size; ++i)
      {
        Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt, i, *key_get_pure(elem_key_at(&gs, i)), elem_weight(&gs, i));
      }

    } else
#endif
      for (i = 0; i < gs.size; ++i)
    {
      Z_TRACE_IF(MSS_TRACE_IF, "global sample %" slint_fmt ": key: %" key_pure_type_fmt, i, *key_get_pure(elem_key_at(&gs, i)));
    }

#ifdef elem_weight
    if (doweights)
    {
      nlspkeys = 0;
    
      slweight_t my_sum, my_first, next_first;
      slweight_t prev_slast, next_sfirst;
      MPI_Status status;

      my_sum = 0;
      for (i = 0; i < gs.size; ++i) my_sum += elem_weight(&gs, i);
    
      prev_slast = 0;
      MPI_Exscan(&my_sum, &prev_slast, 1, weight_mpi_datatype, MPI_SUM, comm);
    
      my_first = elem_weight(&gs, 0);
      next_first = 0;
      MPI_Sendrecv(&my_first, 1, weight_mpi_datatype, (rank > 0)?(rank - 1):MPI_PROC_NULL, 0, &next_first, 1, weight_mpi_datatype, (rank + 1 < size)?(rank + 1):MPI_PROC_NULL, 0, comm, &status);

      next_sfirst = prev_slast + my_sum + next_first;

      Z_TRACE_IF(MSS_TRACE_IF, "weights: my_sum: %" slweight_fmt ", my_first: %" slweight_fmt ", next_first: %" slweight_fmt ", prev_slast: %" slweight_fmt ", next_sfirst: %" slweight_fmt, my_sum, my_first, next_first, prev_slast, next_sfirst);

      j = 0;
      w = prev_slast;

      for (i = 0; i < nsplitter; ++i)
      {
        wi = (i + 1) * lgweights[1] / (nsplitter + 1) - (lgweights[1] / (nsplitter + 1) / 2);

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": wi: %" slweight_fmt, i, wi);

        if (wi > next_sfirst || wi < prev_slast)
        {
          Z_TRACE_IF(MSS_TRACE_IF, " -> continue");
          continue;
        }

        while (w < wi && j < gs.size)
        {
          w += elem_weight(&gs, j);
          ++j;
        }

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": A: j = %" slint_fmt ", w = %" slweight_fmt ", wi: %" slweight_fmt, i, j, w, wi);

        /* step back? */
        if (j > 0 && (w >= wi) && (wi - (w - elem_weight(&gs, j - 1))) <= (w - wi))
        {
          w -= elem_weight(&gs, j - 1);
          --j;
        }

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": B: j = %" slint_fmt ", w = %" slweight_fmt ", wi: %" slweight_fmt ", next_sfirst: %" slweight_fmt, i, j, w, wi, next_sfirst);

        if (w < wi)
        {
          if (next_sfirst < wi || (next_sfirst >= wi && (next_sfirst - wi) <= (wi - w)))
          {
            w = next_sfirst;
            ++j;
          }
        }

        Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt ": C: j = %" slint_fmt ", w = %" slweight_fmt ", wi: %" slweight_fmt, i, j, w, wi);

        if (j > 0 && j - 1 < gs.size)
        {
          lspkeys[nlspkeys] = *key_get_pure(elem_key_at(&gs, j - 1));

          Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": key: %" key_pure_type_fmt ", weight: %" slweight_fmt " (target: %" slweight_fmt ")", i, j - 1, lspkeys[nlspkeys], w, wi);

          ++nlspkeys;

        } else Z_TRACE_IF(MSS_TRACE_IF, "splitter %" slint_fmt " @ %" slint_fmt ": skip", i, j - 1);
      }
      
      Z_TRACE_IF(MSS_TRACE_IF, "local splitters: %d", nlspkeys);

      MPI_Allgather(&nlspkeys, 1, MPI_INT, rcounts, 1, MPI_INT, comm);

      rdispls[0] = 0;
      for (i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + rcounts[i - 1];

      MPI_Allgatherv(lspkeys, nlspkeys, pkey_mpi_datatype, spkeys, rcounts, rdispls, pkey_mpi_datatype, comm);

    } else
#endif
    {
      lspkey = *key_get_pure(elem_key_at(&gs, 0));

      Z_TRACE_IF(MSS_TRACE_IF, "local splitter: %" key_pure_type_fmt, lspkey);
  
      MPI_Allgather(&lspkey, 1, pkey_mpi_datatype, spkeys, 1, pkey_mpi_datatype, comm);
  
      splitter_skip = 1;
    }

    elements_free(&gs);
    elements_free(&xs);
  }

  /* determine splitting positions from splitters with binary search */
  for (i = 0; i < nsplitter; ++i)
  {
    elem_assign_at(s, sp->displs[i], &e);
    e.size = s->size - sp->displs[i];

    sp->displs[i + 1] = sp->displs[i] + sl_search_binary_lt(&e, &spkeys[splitter_skip + i]);

    Z_TRACE_IF(MSS_TRACE_IF, "displs %" slint_fmt ": %d (< %" key_pure_type_fmt ")", i + 1, sp->displs[i + 1], spkeys[splitter_skip + i]);

/*    w = 0;
    for (j = 0; j < sp->displs[i + 1] - sp->displs[i]; ++j) w += elem_weight(&e, j);

    Z_TRACE_IF(MSS_TRACE_IF, " -> count: %d, weight: %" slweight_fmt, sp->displs[i + 1] - sp->displs[i], w);*/
  }
  
  return 0;
}
