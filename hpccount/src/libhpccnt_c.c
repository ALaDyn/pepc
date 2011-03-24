/*==========================================================*/
/*                                                          */
/*  Interception Library for BGP_hpccount                   */
/*                                                          */
/*  Profiling Wrappers for MPI functions, Fortran and C.    */
/*  Timing data is reported when MPI_Finalize is called.    */
/*                                                          */
/*  Initial version: Sept 30, 2010                          */
/*  C. Pospiech ** ACTC  ** IBM Germany                     */
/*  © IBM 2010                                              */
/*                                                          */            
/*==========================================================*/

#include <mpi.h>
#include <libhpc.h>
#define hpm_errchk if (hpm_error_count) MPI_Abort(MPI_COMM_WORLD, 4)

/*----------------------------------------------------------*/
/*    wrapper for C: MPI_Init                               */
/*----------------------------------------------------------*/
int MPI_Init(int * argc, char *** argv) {
  int i, id, rc, bin, taskid, ntasks;
  char * ptr;

  rc = PMPI_Init(argc, argv);
  rc = PMPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  rc = PMPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  /*-----------------------------------------*/
  /* check env variables for tracing options */
  /*-----------------------------------------*/
  /* ptr = getenv("TRACE_ALL_EVENTS");
     if (ptr == NULL) trace_events = 0; */

  /*-----------------------------------------*/
  /* initialize and start the HPM counters   */
  /*-----------------------------------------*/

  hpmInit( taskid, "hpccount_auto_instrumenter" );
  hpm_errchk;
  hpmStart( 1, "All" );
  hpm_errchk;

  return rc;
}


/*----------------------------------------------------------*/
/*    wrapper for C: MPI_Finalize                           */
/*----------------------------------------------------------*/
int MPI_Finalize(void) {

  int rc;

  /*-----------------------------------------*/
  /* stop and terminate the HPM counters     */
  /*-----------------------------------------*/

  hpmStop( 1 );
  hpm_errchk;
  hpmTerminate( 0 );
  hpm_errchk;

  rc = PMPI_Finalize();
  return rc;
}
