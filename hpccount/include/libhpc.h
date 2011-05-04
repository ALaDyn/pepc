/* IBM_PROLOG_BEGIN_TAG                                                   */
/* This is an automatically generated prolog.                             */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/* Licensed Materials - Property of IBM                                   */
/*                                                                        */
/* Restricted Materials of IBM                                            */
/*                                                                        */
/* (C) COPYRIGHT International Business Machines Corp. 2008               */
/* All Rights Reserved                                                    */
/*                                                                        */
/* US Government Users Restricted Rights - Use, duplication or            */
/* disclosure restricted by GSA ADP Schedule Contract with IBM Corp.      */
/*                                                                        */
/* IBM_PROLOG_END_TAG                                                     */
/* "@(#) 1.1 src/ppe/hpct/source/include/libhpc.h, ppe.hpct, ppe_rhya, rhya0837c 08/09/16 10:32:37" */
/****************************************************************************
 *
 *               include file for libhpc
 *
 *               C. Pospiech ** ACTC
 *
 *  $Id: libhpc.h,v 1.3 2009/02/01 15:20:14 cp Exp $
 *
 ****************************************************************************/

#ifndef HPM_LIBHPC_H
#define HPM_LIBHPC_H

#include <string.h>


#define HPM_NO_PARENT            -1
#define HPM_AUTO_PARENT          -2
#define HPM_ONLY_EXCLUSIVE       -3

#define HPMTHREAD_COUNTING      (1)
#define HPMPROCESS_COUNTING     (2)

extern int hpm_error_count;

/*
 * The following definitions are taken over from libhpm version 2
 * Well, mostly taken over...
 */

#define hpmInit(id, progName) \
        _hpmInit_(id, progName, strlen(progName))
#define hpmStart(id, label) _hpm_start_(id, HPM_AUTO_PARENT, \
                                         __LINE__, __FILE__, label, \
                                        strlen(__FILE__), strlen(label), \
                                        HPMPROCESS_COUNTING)
#define hpmStartx(id, par_id, label) _hpm_start_(id, par_id, \
                                         __LINE__, __FILE__, label, \
                                        strlen(__FILE__), strlen(label), \
                                        HPMPROCESS_COUNTING)
#define hpmStop(id) _hpm_stop_(id, __LINE__, HPMPROCESS_COUNTING)
#define hpmTstart(id, label) _hpm_start_(id, HPM_AUTO_PARENT, \
                                          __LINE__, __FILE__, label, \
                                          strlen(__FILE__), strlen(label), \
                                          HPMTHREAD_COUNTING)
#define hpmTstartx(id, par_id, label) _hpm_start_(id, par_id, \
                                          __LINE__, __FILE__, label, \
                                          strlen(__FILE__), strlen(label), \
                                          HPMTHREAD_COUNTING)
#define hpmTstop(id) _hpm_stop_(id, __LINE__, HPMTHREAD_COUNTING)
#define hpm_start(id) _hpm_start_(id, HPM_AUTO_PARENT, \
                                  __LINE__, __FILE__, \
                                  "&nbsp;", strlen(__FILE__), 6, \
                                   HPMPROCESS_COUNTING)
#define hpm_startx(id, par_id) _hpm_start_(id, par_id, \
                                  __LINE__, __FILE__, \
                                  "&nbsp;", strlen(__FILE__), 6, \
                                  HPMPROCESS_COUNTING)
#define hpm_stop(id) _hpm_stop_(id, __LINE__, HPMPROCESS_COUNTING)
#define hpm_tstart(id) _hpm_start_(id, HPM_AUTO_PARENT, \
                                    __LINE__, __FILE__, "&nbsp;", \
                                    strlen(__FILE__), 6, \
                                    HPMTHREAD_COUNTING)
#define hpm_tstartx(id, par_id) _hpm_start_(id, par_id, \
                                    __LINE__, __FILE__, "&nbsp;", \
                                    strlen(__FILE__), 6, \
                                    HPMTHREAD_COUNTING)
#define hpm_tstop(id) _hpm_stop_(id, __LINE__, HPMTHREAD_COUNTING)

#define hpmTerminate(id) _hpm_terminate_(id)
#define hpm_terminate(id) _hpm_terminate_(id)

#ifdef __cplusplus
extern "C" {
#endif

void _hpmInit_( int my_ID, const char* progName, int sz_progName );
void _hpm_start_( int inst_ID, int parent_ID, int line, const char* file,
                  const char* label, int sz_file, int sz_label,
		  int count_context);
void _hpm_stop_( int inst_ID, int line, int count_context);
void _hpm_terminate_( int my_ID );

#ifdef __cplusplus
}
#endif

#endif
