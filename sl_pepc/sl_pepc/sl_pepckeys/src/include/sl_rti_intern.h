/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/include/sl_rti_intern.h
 *  timestamp: 2009-11-04 16:09:17 +0100
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
 #define rti_cadd_cmp(n)          (pepckeys_rti_env.cmc.cmp += n)
 #define rti_cadd_movek(n)        (pepckeys_rti_env.cmc.movek += n)
 #define rti_cadd_moved(n)        (pepckeys_rti_env.cmc.moved += n)
 #define rti_cclear_cmp()         (pepckeys_rti_env.cmc.cmp = 0)
 #define rti_cclear_movek()       (pepckeys_rti_env.cmc.movek = 0)
 #define rti_cclear_moved()       (pepckeys_rti_env.cmc.moved = 0)
 #define rti_cclear_all()         (pepckeys_rti_env.cmc.cmp = pepckeys_rti_env.cmc.movek = pepckeys_rti_env.cmc.moved = 0)
 #define rti_ccmp()               my_rti_ccmp(pepckeys_rti_env)
 #define rti_cmovek()             my_rti_cmovek(pepckeys_rti_env)
 #define rti_cmoved()             my_rti_cmoved(pepckeys_rti_env)

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

 #define rti_tstart(t)            (pepckeys_rti_env.tim[t].start = sl_ts_get_s(), 0)
 #define rti_tstop(t)             (pepckeys_rti_env.tim[t].stop = sl_ts_get_s(), pepckeys_rti_env.tim[t].cumu += (pepckeys_rti_env.tim[t].last = pepckeys_rti_env.tim[t].stop - pepckeys_rti_env.tim[t].start))
 #define rti_tclear(t)            (pepckeys_rti_env.tim[t].last = 0)
 #define rti_treset(t)            (pepckeys_rti_env.tim[t].last = pepckeys_rti_env.tim[t].cumu = 0)
 #define rti_tlast(t)             my_rti_tlast(pepckeys_rti_env, t)
 #define rti_tcumu(t)             my_rti_tcumu(pepckeys_rti_env, t)

#else

 #define rti_tstart(t)            SL_NOP()
 #define rti_tstop(t)             SL_NOP()
 #define rti_tclear(t)            SL_NOP()
 #define rti_treset(t)            SL_NOP()
 #define rti_tlast(t)             0
 #define rti_tcumu(t)             0

#endif


#ifdef SL_USE_RTI_MEM

 #define rti_minc_alloc()         ++pepckeys_rti_env.mem.nalloc
 #define rti_minc_free()          ++pepckeys_rti_env.mem.nfree
 #define rti_malloc(s)            (pepckeys_rti_env.mem.max = xmax(pepckeys_rti_env.mem.cur += s, pepckeys_rti_env.mem.max))
 #define rti_mfree(s)             (pepckeys_rti_env.mem.cur -= s)

#else

 #define rti_minc_alloc()         0
 #define rti_minc_free()          0
 #define rti_malloc(s)            0
 #define rti_mfree(s)             0

#endif


#ifdef SL_USE_RTI

 #define rti_reset()              my_rti_reset(pepckeys_rti_env)
 #define declare_rti_env          rti pepckeys_rti_env;

 extern rti pepckeys_rti_env; /* sl_var pepckeys_rti_env */

#else

 #define rti_reset()              0
 #define declare_rti_env

#endif


#ifdef SL_DATA_IGNORE
 #undef declare_rti_env
#endif

#ifndef declare_rti_env
 #define declare_rti_env
#endif


#endif /* __SL_RTI_INTERN_H__ */
