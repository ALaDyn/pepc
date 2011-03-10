/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core/z_pack.c
 *  timestamp: 2011-02-26 16:57:01 +0100
 *  
 */


#include <stdio.h>

#include "z_pack.h"


#ifdef Z_PACK_ALLOC

#endif /* Z_PACK_ALLOC */


#ifdef Z_PACK_DEBUG

/* z_var z_notice_fstream z_error_fstream z_debug_fstream */
FILE *z_notice_fstream = NULL;
FILE *z_error_fstream = NULL;
FILE *z_debug_fstream = NULL;

#endif /* Z_PACK_DEBUG */


#ifdef Z_PACK_TIME

#ifndef Z_PACK_MPI
double z_time_wtime() /* z_proto, z_func z_time_wtime */
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double) tv.tv_sec + (tv.tv_usec / 1000000.0);
}
#endif

#endif /* Z_PACK_TIME */


#ifdef Z_PACK_RANDOM

#if defined(HAVE_ESSL_H)
# include <essl.h>
#elif defined(HAVE_T4C_H)
# include <trng4c.h>
#endif


#if defined(HAVE_ESSL_H)
double z_nrandom_seed_essl = 1.0; /* z_var z_nrandom_seed_essl */
#endif


void z_nrandom_seed(unsigned long s) /* z_proto, z_func z_nrandom_seed */
{
#if defined(HAVE_ESSL_H)
  z_nrandom_seed_essl = s + 1;
#elif defined(HAVE_T4C_H)
  t4c_dnrand_seed(s);
#endif
}


double z_nrandom() /* z_proto, z_func z_nrandom */
{
  double r[2] = { 0, 0 };

#if defined(HAVE_ESSL_H)
/*  printf("dnrand\n");*/
  dnrand(&z_nrandom_seed_essl, 2, r, NULL, 0);
#elif defined(HAVE_T4C_H)
/*  printf("t4c_dnrand\n");*/
  t4c_dnrand(r);
#endif

/*  printf("z_nrandom: %f\n", r[0]);*/

  return r[0];
}


#if defined(HAVE_ESSL_H)
double z_urandom_seed_essl = 1.0;  /* z_var z_urandom_seed_essl */
#endif


void z_urandom_seed(unsigned long s) /* z_proto, z_func z_urandom_seed */
{
#if defined(HAVE_ESSL_H)
  z_urandom_seed_essl = s + 1;
#elif defined(HAVE_T4C_H)
  t4c_durand_seed(s);
#endif
}


double z_urandom() /* z_proto, z_func z_urandom */
{
  double r[1] = { 0 };

#if defined(HAVE_ESSL_H)
/*  printf("durand\n");*/
  durand(&z_urandom_seed_essl, 1, r);
#elif defined(HAVE_T4C_H)
/*  printf("t4c_durand\n");*/
  t4c_durand(r);
#endif

/*  printf("z_urandom: %f\n", r[0]);*/

  return r[0];
}

#endif /* Z_PACK_RANDOM */


#ifdef Z_PACK_DIGEST

#ifdef HAVE_GCRYPT_H
# include <gcrypt.h>
#endif


z_int_t z_digest_sum_buffer(const void *buffer, z_int_t length, void *sum) /* z_proto, z_func z_digest_sum_buffer */
{
#if defined(HAVE_GCRYPT_H)

  if (!buffer || !sum) return gcry_md_get_algo_dlen(GCRY_MD_CRC32);

  gcry_md_hash_buffer(GCRY_MD_CRC32, sum, buffer, length);

#elif defined (Z_PACK_CRC32)

  if (!buffer || !sum) return sizeof(z_crc32_t);

  *((z_crc32_t *) sum) = z_crc32_buffer(buffer, length);
  
#endif

/*  printf("z: %.8X\n", *((z_crc32_t *) sum));*/

  return 0;
}


#ifdef HAVE_GCRYPT_H
int z_digest_hash_gcrypt_algo = GCRY_MD_MD5; /* z_var, z_digest_hash_gcrypt_algo */
#endif


void z_digest_hash_open(void **hdl) /* z_proto, z_func z_digest_hash_open */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl;
  
  gcry_hdl = z_alloc(1, sizeof(gcry_md_hd_t));

  gcry_md_open(gcry_hdl, z_digest_hash_gcrypt_algo, 0);

  *hdl = gcry_hdl;
#endif
}


void z_digest_hash_close(void *hdl) /* z_proto, z_func z_digest_hash_close */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl = hdl;

  gcry_md_close(*gcry_hdl);

  z_free(gcry_hdl);
#endif
}


void z_digest_hash_write(void *hdl, const void *buffer, z_int_t length) /* z_proto, z_func z_digest_hash_write */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl = hdl;

  gcry_md_write(*gcry_hdl, buffer, length);
#endif
}


z_int_t z_digest_hash_read(void *hdl, void *hash) /* z_proto, z_func z_digest_hash_read */
{
#ifdef HAVE_GCRYPT_H
  gcry_md_hd_t *gcry_hdl = hdl;

  z_int_t dlen = gcry_md_get_algo_dlen(z_digest_hash_gcrypt_algo);

  if (!hdl || !hash) return dlen;

  memcpy(hash, gcry_md_read(*gcry_hdl, 0), dlen);
#endif

  return 0;
}

#endif /* Z_PACK_DIGEST */
