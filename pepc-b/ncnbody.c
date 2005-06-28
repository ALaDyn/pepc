#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define ROWS 200
#define COLS 18

/* rank (number of dimensions) for each variable */
#  define RANK_data 3
#  define RANK_mintime 1


   /* variable ids */
static int data_id=-1;
static int mintime_id=-1;
static int timeid_start=0;
  /* dimension lengths */
static  size_t cols_len = -1;
static  size_t rows_len = -1;
static  size_t timeid_len = NC_UNLIMITED;


void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    exit(1);
  }
}


/*---------------------------------------------------------------------------------------*/
int _ncnbody_open(int rows, int cols, int *pncid) {
  int  ncid;			/* netCDF id */

  /* dimension ids */
  int cols_dim;
  int rows_dim;
  int timeid_dim;

  /* variable shapes */
  int data_dims[RANK_data];
  int mintime_dims[RANK_mintime];

  cols_len=cols;
  rows_len=rows;
  printf("WF ncnbody_open: rows=%d, cols=%d\n",cols_len,rows_len);

  /* enter define mode */
  int stat = nc_create("nbody.nc", NC_CLOBBER, &ncid);
  check_err(stat,__LINE__,__FILE__);

  /* define dimensions */
  stat = nc_def_dim(ncid, "cols", cols_len, &cols_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(ncid, "rows", rows_len, &rows_dim);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim(ncid, "timeid", timeid_len, &timeid_dim);
  check_err(stat,__LINE__,__FILE__);

  /* define variables */

  data_dims[0] = timeid_dim;
  data_dims[1] = rows_dim;
  data_dims[2] = cols_dim;
  stat = nc_def_var(ncid, "data", NC_FLOAT, RANK_data, data_dims, &data_id);
  check_err(stat,__LINE__,__FILE__);

  mintime_dims[0] = timeid_dim;
  stat = nc_def_var(ncid, "mintime", NC_FLOAT, RANK_mintime, mintime_dims, &mintime_id);
  check_err(stat,__LINE__,__FILE__);

  /* leave define mode */
  stat = nc_enddef (ncid);
  check_err(stat,__LINE__,__FILE__);
  printf("WF ncnbody_open: rows=%d, cols=%d ready ncid=%d\n",cols_len,rows_len,ncid);

  timeid_start=0;   
  *pncid=ncid;
  return(1); 
}

void *ncnbody_open_(int *rows, int *cols, int *ncid, int *ret) {
  *ret=_ncnbody_open(*rows,*cols,ncid);
  return;
}

void *ncnbody_open__(int *rows, int *cols, int *ncid, int *ret) {
  *ret=_ncnbody_open(*rows,*cols,ncid);
  return;
}

void *ncnbody_open(int *rows, int *cols, int *ncid, int *ret) {
  *ret=_ncnbody_open(*rows,*cols,ncid);
  return;
}

/*---------------------------------------------------------------------------------------*/
int _ncnbody_put(int ncid,float *vbuffer,int lstored,int vbufcols) {
  int stat;
  static size_t mintime_start[RANK_mintime];
  static size_t mintime_count[RANK_mintime];
  static size_t data_start[RANK_data];
  static size_t data_count[RANK_data];
  static float mintime[1];
  size_t timeid_len;

/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d\n",ncid,lstored,vbufcols); */
  
  /* store mintime */
  timeid_len = 1;		
  mintime_start[0] = timeid_start;
  mintime_count[0] = timeid_len;
  mintime[0]=vbuffer[0];	/* first value of first record */
  /* printf("parm=%d %d\n",ncid, mintime_id);  */
  stat = nc_put_vara_float(ncid, mintime_id, mintime_start, mintime_count, mintime);
  /* printf("stat=%d\n",stat); */
  check_err(stat,__LINE__,__FILE__);

  
  timeid_len = 1;			/* number of records of data data */
  data_start[0] = timeid_start;
  data_start[1] = 0;
  data_start[2] = 0;
  data_count[0] = timeid_len;
  data_count[1] = rows_len;
  data_count[2] = cols_len;
  stat = nc_put_vara_float(ncid, data_id, data_start, data_count, vbuffer);
  check_err(stat,__LINE__,__FILE__);

/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */

  timeid_start++;
  
  return(stat);
}

void* ncnbody_put_(int *ncid,float *vbuffer,int *lstored,int *vbufcols, int *rc) {
  *rc=_ncnbody_put(*ncid,vbuffer, *lstored, *vbufcols);
  return;
}

void* ncnbody_put__(int *ncid,float *vbuffer,int *lstored,int *vbufcols, int *rc) {
  *rc=_ncnbody_put(*ncid,vbuffer, *lstored, *vbufcols);
  return;
}

void* ncnbody_put(int *ncid,float *vbuffer,int *lstored,int *vbufcols, int *rc) {
  *rc=_ncnbody_put(*ncid,vbuffer, *lstored, *vbufcols);
  return;
}

/*---------------------------------------------------------------------------------------*/
int _ncnbody_close(int ncid) {
  int stat;
/*   printf("WF _ncnbody_close: ncid=%d\n",ncid); */
  stat = nc_close(ncid);     check_err(stat,__LINE__,__FILE__);
/*   printf("WF _ncnbody_close: ncid=%d ready stat=%d\n",ncid,stat); */
  return(stat);
}

void* ncnbody_close_(int *ncid, int *ret) {
  *ret=_ncnbody_close(*ncid);
  return;
}

void* ncnbody_close__(int *ncid, int *ret) {
  *ret=_ncnbody_close(*ncid);
  return;
}

void* ncnbody_close(int *ncid, int *ret) {
  *ret=_ncnbody_close(*ncid);
  return;
}
