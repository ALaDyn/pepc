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
/*    Fortran profiling interface with one underscore       */
/*----------------------------------------------------------*/
#ifdef _ONE_UNDERSCORE
#define pmpi_init              pmpi_init_
#define pmpi_init_thread       pmpi_init_thread_
#define pmpi_finalize          pmpi_finalize_
#endif

/*----------------------------------------------------------*/
/*    Fortran profiling interface with two underscores      */
/*----------------------------------------------------------*/
#ifdef _TWO_UNDERSCORE
#define pmpi_init              pmpi_init__
#define pmpi_init_thread       pmpi_init_thread__
#define pmpi_finalize          pmpi_finalize__
#endif

/*----------------------------------------------------------*/
/*    fortran profiling entry points                        */
/*----------------------------------------------------------*/
void pmpi_init(int *);
void pmpi_init_thread(int *, int *, int *);
void pmpi_finalize(int *);

/*----------------------------------------------------------*/
/*    Fortran-to-C interface with one underscore            */
/*----------------------------------------------------------*/
#ifdef _ONE_UNDERSCORE
#define mpi_init              mpi_init_
#define mpi_init_thread       mpi_init_thread_
#define mpi_finalize          mpi_finalize_
#endif

/*----------------------------------------------------------*/
/*    Fortran-to-C interface with two underscores           */
/*----------------------------------------------------------*/
#ifdef _TWO_UNDERSCORE
#define mpi_init              mpi_init__
#define mpi_init_thread       mpi_init_thread__
#define mpi_finalize          mpi_finalize__
#endif

/*----------------------------------------------------------*/
/*    Fortran-to-C interface with no underscores            */
/*----------------------------------------------------------*/

/*----------------------------------------------------------*/
/*    fortran entry points                                  */
/*----------------------------------------------------------*/
void mpi_init(int *);
void mpi_init_thread(int *, int *, int *);
void mpi_finalize(int *);

/*----------------------------------------------------------*/
/*    wrapper for Fortran: mpi_init                         */
/*----------------------------------------------------------*/
void mpi_init(int * info) {
  int i, id, rc, bin, argc, taskid, ntasks;
  int comm_world;
  char * ptr, ** argv;

  pmpi_init(&rc);

  *info = rc;

  comm_world = (int) MPI_COMM_WORLD;
  pmpi_comm_rank(&comm_world, &taskid, &rc);
  pmpi_comm_size(&comm_world, &ntasks, &rc);

  /*-----------------------------------------*/
  /* check env variables for tracing options */
  /*-----------------------------------------*/
  /*ptr = getenv("TRACE_ALL_EVENTS");
    if (ptr == NULL) trace_events = 0;*/

  /*-----------------------------------------*/
  /* initialize and start the HPM counters   */
  /*-----------------------------------------*/

  hpmInit( taskid, "hpccount_auto_instrumenter" );
  hpm_errchk;
  hpmStart( 1, "All" );
  hpm_errchk;

  return;
}

/*----------------------------------------------------------*/
/*    wrapper for Fortran: mpi_init                         */
/*----------------------------------------------------------*/
void mpi_init_thread(int * req, int * pro, int * ierr) {
  int i, id, rc, bin, argc, taskid, ntasks;
  int comm_world;
  char * ptr, ** argv;

  pmpi_init_thread(req, pro, &rc);

  *ierr = rc;

  comm_world = (int) MPI_COMM_WORLD;
  pmpi_comm_rank(&comm_world, &taskid, &rc);
  pmpi_comm_size(&comm_world, &ntasks, &rc);

  /*-----------------------------------------*/
  /* check env variables for tracing options */
  /*-----------------------------------------*/
  /*ptr = getenv("TRACE_ALL_EVENTS");
    if (ptr == NULL) trace_events = 0;*/

  /*-----------------------------------------*/
  /* initialize and start the HPM counters   */
  /*-----------------------------------------*/

  hpmInit( taskid, "hpccount_auto_instrumenter" );
  hpm_errchk;
  hpmStart( 1, "All" );
  hpm_errchk;

  return;
}


/*----------------------------------------------------------*/
/*    wrapper for Fortran: mpi_finalize                     */
/*----------------------------------------------------------*/
void mpi_finalize(int * info) {

  int rc;

  /*-----------------------------------------*/
  /* stop and terminate the HPM counters     */
  /*-----------------------------------------*/

  hpmStop( 1 );
  hpm_errchk;
  hpmTerminate( 0 );
  hpm_errchk;

  rc = PMPI_Finalize();

  *info = rc;
  return;
}

