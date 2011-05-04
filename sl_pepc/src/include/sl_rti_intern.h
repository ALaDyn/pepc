/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_rti_intern.h
 *  timestamp: 2011-03-06 21:59:31 +0100
 *  
 */


#ifndef __SL_RTI_INTERN_H__
#define __SL_RTI_INTERN_H__


#ifndef SL_USE_RTI

 #undef SL_USE_RTI_CMC  /* compare-move-counter */
 #undef SL_USE_RTI_TIM  /* timing */
 #undef SL_USE_RTI_MEM  /* memory */

#endif

#ifdef SL_USE_RTI_CMC

 /* regular commands */
 #define rti_cadd_cmp(n)          (rti_env.cmc.cmp += n)
 #define rti_cadd_movek(n)        (rti_env.cmc.movek += n)
 #define rti_cadd_moved(n)        (rti_env.cmc.moved += n)
 #define rti_cclear_cmp()         (rti_env.cmc.cmp = 0)
 #define rti_cclear_movek()       (rti_env.cmc.movek = 0)
 #define rti_cclear_moved()       (rti_env.cmc.moved = 0)
 #define rti_cclear_all()         (rti_env.cmc.cmp = rti_env.cmc.movek = rti_env.cmc.moved = 0)
 #define rti_ccmp()               my_rti_ccmp(rti_env)
 #define rti_cmovek()             my_rti_cmovek(rti_env)
 #define rti_cmoved()             my_rti_cmoved(rti_env)

 /* chained commands */
 #define cc_rti_cadd_cmp(n)       rti_cadd_cmp(n),
 #define cc_rti_cadd_movek(n)     rti_cadd_movek(n),
 #define cc_rti_cadd_moved(n)     rti_cadd_moved(n),

#else /* SL_USE_RTI_CMC */

 /* regular commands */
 #define rti_cadd_cmp(n)
 #define rti_cadd_movek(n)
 #define rti_cadd_moved(n)
 #define rti_cclear_cmp()
 #define rti_cclear_movek()
 #define rti_cclear_moved()
 #define rti_cclear_all()
 #define rti_ccmp()               0
 #define rti_cmovek()             0
 #define rti_cmoved()             0

 /* chained commands */
 #define cc_rti_cadd_cmp(n)
 #define cc_rti_cadd_movek(n)
 #define cc_rti_cadd_moved(n)

#endif /* SL_USE_RTI_CMC */


#ifdef SL_USE_RTI_TIM

 #define rti_tstart(t)            (rti_env.tim[t].start = z_time_get_s(), ++rti_env.tim[t].num)
 #define rti_tstop(t)             (rti_env.tim[t].stop = z_time_get_s(), rti_env.tim[t].cumu += (rti_env.tim[t].last = rti_env.tim[t].stop - rti_env.tim[t].start))
 #define rti_tclear(t)            (rti_env.tim[t].last = 0)
 #define rti_treset(t)            (rti_env.tim[t].last = rti_env.tim[t].cumu = 0, rti_env.tim[t].num = 0)
 #define rti_tlast(t)             my_rti_tlast(rti_env, t)
 #define rti_tcumu(t)             my_rti_tcumu(rti_env, t)
 #define rti_tnum(t)              my_rti_tnum(rti_env, t)

#else

 #define rti_tstart(t)            Z_NOP()
 #define rti_tstop(t)             Z_NOP()
 #define rti_tclear(t)            Z_NOP()
 #define rti_treset(t)            Z_NOP()
 #define rti_tlast(t)             0
 #define rti_tcumu(t)             0
 #define rti_tnum(t)              0

#endif


#ifdef SL_USE_RTI_MEM

 #define rti_minc_alloc()         ++rti_env.mem.nalloc
 #define rti_minc_free()          ++rti_env.mem.nfree
 #define rti_malloc(_s_)          (rti_env.mem.max = z_max(_s_, rti_env.mem.max), rti_env.mem.cur += _s_, rti_env.mem.cur_max = z_max(rti_env.mem.cur, rti_env.mem.cur_max))
 #define rti_mfree(_s_)           (rti_env.mem.cur -= _s_)

 #define cc_rti_minc_alloc()      rti_minc_alloc(),
 #define cc_rti_minc_free()       rti_minc_free(),
 #define cc_rti_malloc(_s_)       rti_malloc(_s_),
 #define cc_rti_mfree(_s_)        rti_mfree(_s_),

#else

 #define rti_minc_alloc()         Z_NOP()
 #define rti_minc_free()          Z_NOP()
 #define rti_malloc(_s_)          Z_NOP()
 #define rti_mfree(_s_)           Z_NOP()

 #define cc_rti_minc_alloc()
 #define cc_rti_minc_free()
 #define cc_rti_malloc(_s_)
 #define cc_rti_mfree(_s_)

#endif


#ifdef SL_USE_RTI

 #define rti_reset()              my_rti_reset(rti_env)

# ifdef DECLARE_RTI
rti rti_env;  /* sl_global, sl_var rti_env */
# endif

#else

 #define rti_reset()              Z_NOP()

#endif


#ifdef SL_DATA_IGNORE
 #undef declare_rti_env
#endif

#ifndef declare_rti_env
 #define declare_rti_env
#endif


#endif /* __SL_RTI_INTERN_H__ */
