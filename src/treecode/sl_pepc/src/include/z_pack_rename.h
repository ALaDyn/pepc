/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/z_pack_rename.h
 *  
 */


#ifndef __Z_PACK_RENAME_H__
#define __Z_PACK_RENAME_H__


#define Z_CONCAT(_a_, _b_)           Z_CONCAT_(_a_, _b_)
#define Z_CONCAT_(_a_, _b_)          _a_##_b_

#define Z_CONCONCAT(_a_, _b_, _c_)   Z_CONCONCAT_(_a_, _b_, _c_)
#define Z_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#ifdef Z_PREFIX
# define Z_VAR(_v_)   Z_CONCAT(Z_PREFIX, _v_)
# define Z_FUNC(_f_)  Z_CONCAT(Z_PREFIX, _f_)
#else
# define Z_VAR(_v_)   _v_
# define Z_FUNC(_f_)  _f_
#endif


/* src/core/z_crc32.c */
#define z_crc32_table_size  Z_VAR(z_crc32_table_size)
#define z_crc32_table  Z_VAR(z_crc32_table)
#define z_crc32_make_table  Z_FUNC(z_crc32_make_table)
#define z_crc32_print_table  Z_FUNC(z_crc32_print_table)
#define z_crc32_buffer_update  Z_FUNC(z_crc32_buffer_update)
#define z_crc32_buffer  Z_FUNC(z_crc32_buffer)

/* src/core/z_pack.c */
#define z_notice_fstream  Z_VAR(z_notice_fstream)
#define z_error_fstream  Z_VAR(z_error_fstream)
#define z_debug_fstream  Z_VAR(z_debug_fstream)
#define z_time_wtime  Z_FUNC(z_time_wtime)
#define z_srandom_seed  Z_VAR(z_srandom_seed)
#define z_srandom  Z_FUNC(z_srandom)
#define z_random  Z_FUNC(z_random)
#define z_srandom64  Z_FUNC(z_srandom64)
#define z_random64  Z_FUNC(z_random64)
#define z_random64_minmax  Z_FUNC(z_random64_minmax)
#define z_random64u  Z_FUNC(z_random64u)
#define z_random64u_minmax  Z_FUNC(z_random64u_minmax)
#define z_nrandom_seed_essl  Z_VAR(z_nrandom_seed_essl)
#define z_nrandom_seed  Z_FUNC(z_nrandom_seed)
#define z_nrandom  Z_FUNC(z_nrandom)
#define z_urandom_seed_essl  Z_VAR(z_urandom_seed_essl)
#define z_urandom_seed  Z_FUNC(z_urandom_seed)
#define z_urandom  Z_FUNC(z_urandom)
#define z_digest_sum_buffer  Z_FUNC(z_digest_sum_buffer)
#define z_digest_hash_open  Z_FUNC(z_digest_hash_open)
#define z_digest_hash_close  Z_FUNC(z_digest_hash_close)
#define z_digest_hash_write  Z_FUNC(z_digest_hash_write)
#define z_digest_hash_read  Z_FUNC(z_digest_hash_read)
#define z_gmp_mpz_set_ull  Z_FUNC(z_gmp_mpz_set_ull)
#define z_gmp_mpz_set_sll  Z_FUNC(z_gmp_mpz_set_sll)
#define z_gmp_mpz_get_ull  Z_FUNC(z_gmp_mpz_get_ull)
#define z_gmp_mpz_get_sll  Z_FUNC(z_gmp_mpz_get_sll)


#endif /* __Z_PACK_RENAME_H__ */
