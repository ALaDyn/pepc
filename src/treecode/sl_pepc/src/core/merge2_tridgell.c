/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/merge2_tridgell.c
 *  
 */

/* slightly optimized (performance and complexity of implementation?) merge2 from [Tridgell] */


#include "sl_common.h"


#define the_merge_x_0_1  merge2_basic_straight_x_0_1


typedef struct _block
{
  slint rank, size, quan;

  elements_t s;

} block;


slint merge2_compo_tridgell(elements_t *s0, elements_t *s1, elements_t *sx) /* sl_proto, sl_func merge2_compo_tridgell */
{
  slint blocksize = (slint) sqrt(s0->size + s1->size);

  slint s0blocks = (s0->size + blocksize - 1) / blocksize;
  slint s1blocks = (s1->size + blocksize - 1) / blocksize;

  slint s0mod = s0->size % blocksize;
  slint s1offset = (blocksize - s0mod) % blocksize;

  slint s0rank = 0;
  slint s1rank = s0blocks;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  block *blocks = NULL;
#else
  block blocks[s0blocks + s1blocks + 3];
#endif

  elements_t src0, src0b, src1, src1b, dst, _sx;

  slint nblocks, cblock, the_empty_blocks[3], empty_blocks;

  slint i, n, t, s;


/*  printf("jetzt gehts los!\n");
  printf("blocksize = %d\n", blocksize);
  printf("s0size = %d, s0blocks = %d, s0mod = %d\n", s0->size, s0blocks, s0mod);
  printf("s1size = %d, s1blocks = %d, s1offset = %d\n", s1->size, s1blocks, s1offset);*/

  /* if one list is empty, there is nothing to do */
  if (s0->size == 0 || s1->size == 0) return 0;

  if (sx == NULL)
  {
    sx = &_sx;
    elements_alloc(sx, 3 * blocksize, SLCM_ALL);
  }

  /* if the supplied eXtra-space is less than 3*blocksize */
  if (sx->size < 3 * blocksize) return -1;

#ifdef NO_VARIABLE_LENGTH_ARRAYS
  blocks = z_alloc(s0blocks + s1blocks + 3, sizeof(block));
#endif

  elem_assign(s0, &src0);
  elem_assign(s1, &src1);

/*  printf("src0.size = %d - src1.size = %d\n", src0.size, src1.size);*/

  /* skip the first 's1offset' elements in s1 when building free blocks */
  elem_assign_at(s0, 0, &src0b);
  elem_assign_at(s1, s1offset, &src1b);

  /* first free block can be undersized */
  blocks[0].rank = -3;
  elem_assign_at(sx, 0 * blocksize, &blocks[0].s);
  blocks[0].s.size = blocksize;
  blocks[0].size = 0;

  blocks[1].rank = -2;
  elem_assign_at(sx, 1 * blocksize, &blocks[1].s);
  blocks[1].s.size = blocksize;
  blocks[1].size = 0;

  blocks[2].rank = -1;
  elem_assign_at(sx, 2 * blocksize, &blocks[2].s);
  blocks[2].s.size = blocksize;
  blocks[2].size = 0;

  cblock = -1;
  nblocks = 3;

  /* merge form s0 and s1 to the destination blocks */
  while (src0.size > 0 && src1.size > 0)
  {
    /* take new free block */
    cblock++;

/*    printf("src0.size = %d - src1.size = %d\n", src0.size, src1.size);
    printf("filling %dth block @ %p\n", cblock, blocks[cblock].s.keys);*/

    /* merge to the new free block */
    blocks[cblock].size = the_merge_x_0_1(&src0, &src1, &blocks[cblock].s);

/*    printf("one merge done, |src0| = %d, |src1| = %d\n", src0.size, src1.size);
    printf("free: %d & %d\n", src0.keys - src0b.keys, src1.keys - src1b.keys);*/

    /* append a new free and complete block from s0 */
    if (src0.keys - src0b.keys >= blocksize)
    {
/*      printf("adding free #%d from s0 @ %p\n", nblocks, src0b.keys);*/

      elem_assign(&src0b, &blocks[nblocks].s);
      blocks[nblocks].s.size = blocksize;
      blocks[nblocks].rank = s0rank++;
      blocks[nblocks].size = 0;

      elem_add(&src0b, blocksize);

/*      printf("new from s0 may start @ %p\n", src0b.keys);*/

      nblocks++;
    }

    /* append a new free and complete block from s1 */
    if (src1.keys - src1b.keys >= blocksize)
    {
/*      printf("adding free #%d from s1 @ %p\n", nblocks, src1b.keys);*/

      elem_assign(&src1b, &blocks[nblocks].s);
      blocks[nblocks].s.size = blocksize;
      blocks[nblocks].rank = s1rank++;
      blocks[nblocks].size = 0;

      elem_add(&src1b, blocksize);

      nblocks++;
    }
  }

/*  for (i = 0; i < nblocks; ++i)
  {
    printf("block %d @ %p: rank = %d, content: %d of %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size);
    printf("block %d: content = %d\n", i, blocks[i].size);
    elements_print_keys(&blocks[i].s);
  }*/

/*  printf("one is empty, |src0| = %d, |src1| = %d\n", src0.size, src1.size);*/

/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);
  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* prepare the last half-filled block */
  elem_assign_at(&blocks[cblock].s, blocks[cblock].size, &dst);
  dst.size = blocks[cblock].s.size - blocks[cblock].size;

/*  printf("cblock %d (%d,%d) needs %d\n", cblock, blocks[cblock].s.size, blocks[cblock].size, dst.size);*/

  /* copy all the remaining elements from s0 to the destination blocks */
  while (src0.size > 0)
  {
    n = z_min(src0.size, dst.size);

/*    printf("filling %dth block @ %p with %d elements\n", cblock, blocks[cblock].s.keys, n);*/

    /* fill up with elements from s0 */
    elem_ncopy(&src0, &dst, n);
    elem_add(&src0, n);
    src0.size -= n;
    blocks[cblock].size += n;

/*    printf("copy done, |src0| = %d, |src1| = %d\n", src0.size, src1.size);
    printf("free: %d & %d\n", src0.keys - src0b.keys, src1.keys - src1b.keys);*/

    /* append a new free and full-sized block from s0 */
    if (src0.keys - src0b.keys >= blocksize)
    {
/*      printf("adding free #%d from s0 @ %p\n", nblocks, src0b.keys);*/

      elem_assign(&src0b, &blocks[nblocks].s);
      blocks[nblocks].s.size = blocksize;
      blocks[nblocks].rank = s0rank++;
      blocks[nblocks].size = 0;

      elem_add(&src0b, blocksize);

/*      printf("setting src0b to %p\n", src0b.keys);*/

      nblocks++;
    }

    /* take a new free block */
    cblock++;

    elem_assign(&blocks[cblock].s, &dst);
  }

  /* add an undersized free block from the end of s1, if possible */
  if (src1.keys - src1b.keys > 0)
  {
/*    printf("adding free #%d from s1 @ %p size %d\n", nblocks, src1b.keys, src1.keys - src1b.keys);*/

    elem_assign(&src1b, &blocks[nblocks].s);
    blocks[nblocks].s.size = src1.keys - src1b.keys;
    blocks[nblocks].rank = s1rank++;
    blocks[nblocks].size = 0;

    nblocks++;
  }

  empty_blocks = 0;
/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);*/

/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);
  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* if there is a (broken) block at the end of s0 and the beginning of s1 */
  if (s1offset > 0)
  {
/*    printf("handle broken block #%d of size %d, break into %d:%d, from %p to %p\n", s0blocks - 1, blocks[s0blocks - 1].size, s0mod, z_min(s1offset, src1.keys - s1->keys), blocks[s0blocks - 1].s.keys, src0b.keys);*/

    /* bring this bock in the right position */
    elem_ncopy(&blocks[s0blocks - 1].s, &src0b, s0mod);
    elem_assign_at(&blocks[s0blocks - 1].s, s0mod, &dst);
    elem_ncopy(&dst, s1, z_min(s1offset, src1.keys - s1->keys));

    /* mark the source of this block as empty */
    blocks[s0blocks - 1].size = 0;
  }

/*  printf("elements left: s0 -> %d, s1 -> %d\n", src0.size, src1.size);
  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* notice the possibly empty last 3 blocks */
  if (nblocks > 0)
  if (blocks[nblocks - 1].size == 0) the_empty_blocks[empty_blocks++] = nblocks - 1;
  if (nblocks > 1)
  if (blocks[nblocks - 2].size == 0) the_empty_blocks[empty_blocks++] = nblocks - 2;
  if (nblocks > 2)
  if (blocks[nblocks - 3].size == 0) the_empty_blocks[empty_blocks++] = nblocks - 3;

  /* notice the empty block (created by handling a broken block at the end of s0 and the beginning of s1),
     BUT only when he wasn't one of the former 3 empty block (and therefore already mentioned) */
  if (blocks[s0blocks - 1].size == 0 && (s0blocks - 1 < nblocks - 3)) the_empty_blocks[empty_blocks++] = s0blocks - 1;

/*  printf("empty_blocks at the end: %d -> %d, %d, %d\n", empty_blocks, the_empty_blocks[0], the_empty_blocks[1], the_empty_blocks[2]);

  for (i = 0; i < nblocks; ++i)
  {
    printf("block %d @ %p: rank = %d, content: %d of %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size);
    printf("block %d: content = %d\n", i, blocks[i].size);
    elements_print_keys(&blocks[i].s);
  }*/

  /* bring all blocks to their right positions */

  /* 1st, start cycling at empty blocks, this (max 3) cycles ending in the extraspace (sx) */
  for (i = 0; i < empty_blocks; ++i)
  {
    t = the_empty_blocks[i];
    s = blocks[t].rank;

/*    printf("empty cycle %d starting at %d, size = %d\n", i, the_empty_blocks[i], blocks[the_empty_blocks[i]].size);*/

    while (s >= 0)
    {
/*      printf("copying %d elements from block #%d to #%d\n", blocks[t].s.size, s, t);*/

      elem_ncopy(&blocks[s].s, &blocks[t].s, blocks[t].s.size);
      blocks[t].size = -1; /* mark the target block as 'done' */
      blocks[s].size = 0; /* mark the source block as 'empty' */

/*      printf("source: %d (%d), target: %d (%d), size = %d, first = %d, last = %d\n", s, blocks[s].rank, t, blocks[t].rank, blocks[t].s.size, blocks[t].s.keys[0], blocks[t].s.keys[blocks[t].s.size - 1]);*/

      t = s;
      s = blocks[t].rank;
    }
  }

/*  for (i = 0; i < nblocks; ++i) printf("block %d @ %p: rank = %d, content: %d of %d, first = %d, last = %d\n", i, blocks[i].s.keys, blocks[i].rank, blocks[i].size, blocks[i].s.size, blocks[i].s.keys[0], blocks[i].s.keys[blocks[i].size - 1]);*/

  /* now: the 3 blocks residing in the extraspace (sx) are at there correct positions -> extraspace (sx) is free for use */

  /* 2nd, start cycling at full blocks */
  for (i = nblocks - 1; i > 2; --i)
  if (blocks[i].size > 0) /* skip already well-placed blocks (marked as 'done' -> size == -1) */
  {
/*    printf("cycle starting at %d, size = %d\n", i, blocks[i].size);*/

    /* copy the content of the current block to an empty location, if necessary */
    elem_ncopy(&blocks[i].s, sx, blocks[i].size);

    t = i;
    s = blocks[t].rank;

    /* cycle until we stuck at the current block ('i' is where we started) */
    while (s != i)
    {
/*      printf("copying %d elements from block #%d to #%d\n", blocks[t].s.size, s, t);*/

      elem_ncopy(&blocks[s].s, &blocks[t].s, blocks[t].s.size);
      blocks[t].size = -1; /* mark the target block as 'done' */
      blocks[s].size = 0; /* mark the source block as 'empty' */

/*      printf("source: %d (%d), target: %d (%d), size = %d, first = %d, last = %d\n", s, blocks[s].rank, t, blocks[t].rank, blocks[t].s.size, blocks[t].s.keys[0], blocks[t].s.keys[blocks[t].s.size - 1]);*/

      t = s;
      s = blocks[t].rank;
    }

    /* finished this cycle by copying from the empty location */
    elem_ncopy(sx, &blocks[t].s, blocks[t].s.size);
    blocks[t].size = -1;
  }
  
#ifdef NO_VARIABLE_LENGTH_ARRAYS
  z_free(blocks);
#endif

  if (sx == &_sx) elements_free(sx);

  return 0;
}


#undef the_merge_x_0_1
