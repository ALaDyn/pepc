/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_environment_intern.h
 *  timestamp: 2009-11-04 18:14:36 +0100
 *  
 */


#ifndef __SL_ENVIRONMENT_INTERN_H__
#define __SL_ENVIRONMENT_INTERN_H__


#ifndef sl_alloc
 #ifdef DEBUG
  #define sl_alloc(n, s)           calloc((n), (s))
 #else
  #define sl_alloc(n, s)           malloc((n) * (s))
 #endif
#endif

#ifndef sl_free
 #define sl_free(p)                free(p)
#endif


#ifndef sl_rand

 #undef sl_srand
 #undef SL_RAND_MIN
 #undef SL_RAND_MAX

# include <stdlib.h>

# define SL_RAND_MIN         0
# define SL_RAND_MAX         RAND_MAX
# define sl_rand()           rand()
# define sl_srand(s)         srand(s)

#endif

#ifndef sl_rand01
 #define sl_rand01()         ((double) sl_rand() / (double) SL_RAND_MAX)
#endif

#ifndef sl_rand11
 #define sl_rand11()         ((sl_rand01() * 2.0) - 1.0)
#endif

#ifndef sl_rand_minmax
 #define sl_rand_minmax(_min, _max)  (_min + ((double) (_max - _min) * (sl_rand() - SL_RAND_MIN) / (SL_RAND_MAX - SL_RAND_MIN)))
#endif


#ifndef sl_ts_type

 #undef sl_ts_save
 #undef sl_ts_diff_ms
 #undef sl_ts_get_ms
 #undef sl_ts_diff_s
 #undef sl_ts_get_s

 #ifdef SL_USE_MPI

  #define sl_ts_type                 double
  #define sl_ts_save(t)              (t = MPI_Wtime())
  #define sl_ts_diff_ms(from, to)    (((to) - (from)) * 1000)
  #define sl_ts_get_ms()             (MPI_Wtime() * 1000)
  #define sl_ts_diff_s(from, to)     (((to) - (from)))
  #define sl_ts_get_s()              (MPI_Wtime())

  #define declare_ts_temp_mpi

 #else

  #include <sys/time.h>

  #define sl_ts_type                 struct timeval
  #define sl_ts_save(t)              (gettimeofday(&(t), NULL))
  #define sl_ts_diff_ms(from, to)    ((double) (((to).tv_sec - (from).tv_sec) * 1000.0 + ((to).tv_usec - (from).tv_usec) / 1000.0))
  #define sl_ts_get_ms()             (sl_ts_save(pepckeys_ts_temp), (double) (pepckeys_ts_temp.tv_sec * 1000.0 + pepckeys_ts_temp.tv_usec / 1000.0))
  #define sl_ts_diff_s(from, to)     ((double) (((to).tv_sec - (from).tv_sec) + ((to).tv_usec - (from).tv_usec) / 1000000.0))
  #define sl_ts_get_s()              (sl_ts_save(pepckeys_ts_temp), (double) (pepckeys_ts_temp.tv_sec + pepckeys_ts_temp.tv_usec / 1000000.0))

  #define declare_ts_temp            sl_ts_type pepckeys_ts_temp;
  extern sl_ts_type pepckeys_ts_temp;  /* sl_var pepckeys_ts_temp */

 #endif

#endif


#ifdef SL_DATA_IGNORE
 #undef declare_ts_temp
 #undef declare_ts_temp_mpi
#endif

#ifndef declare_ts_temp
 #define declare_ts_temp
#endif

#ifndef declare_ts_temp_mpi
 #define declare_ts_temp_mpi
#endif


#endif /* __SL_ENVIRONMENT_INTERN_H__ */
