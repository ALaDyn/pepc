/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/sortnet.c
 *  timestamp: 2011-02-08 21:32:53 +0100
 *  
 */

/* a sorting_network function returns the counterpart of processor 'rank' (0..size-1) in the compare/merge-exchange operation in the 'stage' (0..)
   of the described sorting network consisting of 'size' number of participants, 'snp' may be pointing to additional preferences of the sorting-network
   - returning  < 0: there is no such 'stage' in the network
   - returning >= 0: as the rank of the counterpart, == rank or >= size means that the processor of the counterpart doesn't exist */


#include "sl_common.h"


/* simple hypercube, low -> high */
slint sn_hypercube_lh(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_hypercube_lh */
{
  slint stages = ilog2f(size);

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return rank ^ z_powof2(stage);
}

/* simple hypercube, high -> low */
slint sn_hypercube_hl(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_hypercube_hl */
{
  slint stages = ilog2f(size);

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return rank ^ z_powof2(stages - 1 - stage);
}


/* Knuth, pg. 241, 5.3.4 ex. 37 / Akl, pg. 41, 3.2 */
slint sn_odd_even_trans(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_odd_even_trans */
{
  slint stages = size;

  if (stages <= 2) --stages;

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return z_max(0, z_min(size - 1, ((stage % 2 == rank % 2)?rank + 1:rank - 1)));
}


/* single odd-stage */
slint sn_odd(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_odd */
{
  slint stages = 1;

  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return z_max(0, z_min(size - 1, ((0 == rank % 2)?rank + 1:rank - 1)));
}


/* single even-stage */
slint sn_even(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_even */
{
  slint stages = 1;

  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  return z_max(0, z_min(size - 1, ((1 == rank % 2)?rank + 1:rank - 1)));
}


/* by Batcher / Knuth, pg. 225, 5.3.4 / Akl, pg. 23, 2.3 */
slint sn_batcher(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_batcher */
{
  slint stages, _p, p, q, r, d;

  _p = ilog2f(size);
  stages = (_p * (_p + 1)) / 2;

  /* if there are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  if (up != NULL) *up = 0;

  for (_p = p = z_powof2(_p - 1); p > 0; p /= 2)
  {
    q = _p;
    r = 0;
    d = p;

    do
    {
      if (stage <= 0)
      {
        if ((rank & p) == r) return rank + d;
        if ((rank - d >= 0) && (((rank - d) & p) == r)) return rank - d;
        return size;

      } else --stage;

      d = q - p;
      q /= 2;
      r = p;
    } while (d != 0); /* the original condition was: q != p before 'd = q - p' */
  }

  return -1;
}


/* by Batcher / Knuth, pg. 232, 5.3.4 / Akl, pg. 29, 2.4 */
slint sn_bitonic(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_bitonic */
{
  slint stages, p, i = 0;

  p = ilog2f(size);
  stages = (p * (p + 1)) / 2;

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  for (i = 0; (i < p) && (stage > i); stage -= ++i);

  if (up != NULL) *up = (z_powof2(i + 1) & rank);

  return rank ^ z_powof2(i - stage);
}


/* connected */
slint sn_connected(slint size, slint rank, slint stage, void *snp, slint *up) /* sl_proto, sl_func sn_connected */
{
  sn_func *snfs = (sn_func *) snp;
  slint i = 0, stages = 0, cp = -1;

  while (snfs[i] != NULL) stages += (snfs[i++])(size, rank, -1, NULL, NULL);

  /* if there is are no stages */
  if (stages <= 0) return -1;
  /* if the rank is out of range, return 'finshed' */
  if (rank >= size) return -1;
  /* if 'stage < 0' return the number of stages */
  if (stage < 0) return stages;
  /* if the stage is to large, return 'finshed' */
  if (stage >= stages) return -1;

  i = 0;
  while ((snfs[i] != NULL) && (cp <= -1) && (stage >= 0))
  {
    cp = (snfs[i])(size, rank, stage, NULL, up);
    stage -= (snfs[i])(size, rank, -1, NULL, NULL);
    ++i;
  }

  return cp;
}
