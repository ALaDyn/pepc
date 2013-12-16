/*
* This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
* 
* Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
*                         Forschungszentrum Juelich GmbH,
*                         Germany
* 
* PEPC is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* PEPC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <time.h>

/* rank (number of dimensions) for each variable */
#  define RANK_data 3   /* time, col, row */
#  define RANK_time 2  /* time, (exec-time, sim-time, index_fielddesc, index_field1, index_field2..., index_particles) */ 
#  define RANK_field 4  /* time, x, y, z */
#  define RANK_fielddesc 3 /* time, field, desc */
#  define RANK_selfields 3
#  define RANK_vecField 5 /*time, x, y, z, components */
#  define RANK_vecFielddesc 2 /* time, desc */
#  define RANK_selVecField 2

   /* variable ids */
static int data_id=-1;
/* static int mintime_id=-1; */
static int time_id;
static int timeid_start=0;
static int field_id1 = -1;
static int field_id2 = -1;
static int field_id3 = -1;
static int field_id4 = -1;

static int field_desc_id = -1;

static int sel_fields_id = -1;

static int vecField_id = -1;
static int vecField_desc_id = -1;
static int sel_vecField_id = -1;

static int timeid_start_particles=0;
static int timeid_start_fielddesc = 0;
static int timeid_start_field1 = 0;
static int timeid_start_field2 = 0;
static int timeid_start_field3 = 0;
static int timeid_start_field4 = 0;
static int timeid_start_vecField = 0; 
static int timeid_start_vecFielddesc = 0;

  /* dimension lengths */
static  size_t cols_len = -1;
static  size_t rows_len = -1;
static  size_t timeid_len = NC_UNLIMITED;
static  size_t x_len = -1;
static  size_t y_len = -1;
static  size_t z_len = -1;

static size_t component_len = -1; 

static time_t begin;


void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    exit(1);
  }
}


/*---------------------------------------------------------------------------------------*/
int _ncnbody_open(int rows, int cols, int x, int y, int z, int *pncid) {
  
  int  ncid;			/* netCDF id */

  /* dimension ids */
  int cols_dim;
  int rows_dim;
  int timeid_dim;
  int time_and_indices_dim;

  /* field dims */
  int xdim;
  int ydim;
  int zdim;
  int compdim;

  int fielddesc_dim;
  int fieldnr_dim;
  int selfield_dim;
  
  int vecFielddesc_dim;
  int selvecField_dim;

  /* variable shapes */
  int data_dims[RANK_data];
/*   int mintime_dims[RANK_mintime]; */
  int time_dims[RANK_time];
  /*field*/
  int field_dims[RANK_field];
  int field_desc_dims[RANK_fielddesc];
  int sel_field_dims[RANK_selfields];

  /*vecField*/
  int vecField_dims[RANK_vecField];
  int vecFielddesc_dims[RANK_vecFielddesc];
  int selVecField_dims[RANK_selVecField];


  begin = time(NULL); 
  cols_len=cols;
  rows_len=rows;
  
  x_len = x;
  y_len = y;
  z_len = z;
  component_len = 3; /* for vecFields */

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
  stat =  nc_def_dim(ncid, "time_and_indices", 10, &time_and_indices_dim);

  /* field, vecField */
  stat = nc_def_dim( ncid, "xdim", x_len, &xdim );
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim( ncid, "ydim", y_len, &ydim );
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim( ncid, "zdim", z_len, &zdim );
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim( ncid, "compdim", component_len, &compdim );


  /* field desc */
  stat = nc_def_dim( ncid, "fielddesc_dim", 6, &fielddesc_dim );
  check_err(stat,__LINE__,__FILE__);
  stat = nc_def_dim( ncid, "fieldnr_dim", 4, &fieldnr_dim );
  check_err(stat,__LINE__,__FILE__);

  /* field selection */
  stat = nc_def_dim( ncid, "selfield_dim", 1, &selfield_dim );
  check_err(stat,__LINE__,__FILE__);

  /* vecField desc */
  stat = nc_def_dim( ncid, "vecFielddesc_dims", 6, &vecFielddesc_dim );
  check_err(stat,__LINE__,__FILE__);
  
  /* vecField selection */
  stat = nc_def_dim( ncid, "selVecField_dims", 1, &selvecField_dim );
  check_err(stat,__LINE__,__FILE__);
 
 
  /* define variables */

  data_dims[0] = timeid_dim;
  data_dims[1] = rows_dim;
  data_dims[2] = cols_dim;
  stat = nc_def_var(ncid, "data", NC_FLOAT, RANK_data, data_dims, &data_id);
  check_err(stat,__LINE__,__FILE__);

  field_dims[0] = timeid_dim;
  field_dims[1] = xdim;
  field_dims[2] = ydim;
  field_dims[3] = zdim;
 
  field_desc_dims[0] = timeid_dim;
  field_desc_dims[1] = fieldnr_dim;
  field_desc_dims[2] = fielddesc_dim;

  sel_field_dims[0] = timeid_dim;
  sel_field_dims[1] = fieldnr_dim;
  sel_field_dims[2] = selfield_dim;

  vecField_dims[0] = timeid_dim; 
  vecField_dims[1] = xdim; 
  vecField_dims[2] = ydim;
  vecField_dims[3] = zdim;
  vecField_dims[4] = compdim;

  vecFielddesc_dims[0] = timeid_dim;
  vecFielddesc_dims[1] = vecFielddesc_dim;
  
  selVecField_dims[0] = timeid_dim;
  selVecField_dims[1] = selvecField_dim;
 

  stat = nc_def_var( ncid, "field1", NC_FLOAT, RANK_field, field_dims, &field_id1 );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "field2", NC_FLOAT, RANK_field, field_dims, &field_id2 );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "field3", NC_FLOAT, RANK_field, field_dims, &field_id3 );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "field4", NC_FLOAT, RANK_field, field_dims, &field_id4 );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "fielddesc", NC_FLOAT, RANK_fielddesc, field_desc_dims, &field_desc_id );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "sel_fields", NC_INT, RANK_selfields, sel_field_dims, &sel_fields_id );
  check_err(stat,__LINE__,__FILE__);

  
  stat = nc_def_var( ncid, "vecField1", NC_FLOAT, RANK_vecField, vecField_dims, &vecField_id );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "vecFielddesc", NC_FLOAT, RANK_vecFielddesc, vecFielddesc_dims, &vecField_desc_id );
  check_err(stat,__LINE__,__FILE__);

  stat = nc_def_var( ncid, "sel_vecFields", NC_INT, RANK_selVecField, selVecField_dims, &sel_vecField_id );
  check_err(stat,__LINE__,__FILE__);
  
/*   mintime_dims[0] = timeid_dim; */
  time_dims[0] = timeid_dim;
  time_dims[1] = time_and_indices_dim;
  
  stat = nc_def_var(ncid, "time_and_indices", NC_FLOAT, RANK_time, time_dims, &time_id);
  check_err(stat,__LINE__,__FILE__);

  /* leave define mode */
  stat = nc_enddef (ncid);
  check_err(stat,__LINE__,__FILE__);
  printf("WF ncnbody_open: rows=%d, cols=%d ready ncid=%d\n",cols_len,rows_len,ncid);

  timeid_start=0;  
  timeid_start_particles = 0;
  timeid_start_fielddesc = 0;
  static int timeid_start_field1 = 0;
  static int timeid_start_field2 = 0;
  static int timeid_start_field3 = 0;
  static int timeid_start_field4 = 0;
  
  static int timeid_start_vecField = 0;
  timeid_start_vecFielddesc = 0;

  *pncid=ncid;
  return(1); 
}

void *ncnbody_open_(int *rows, int *cols, int *x, int *y, int *z, int *ncid, int *ret) {
  *ret=_ncnbody_open(*rows,*cols,*x, *y, *z, ncid);
  return;
}

void *ncnbody_open__(int *rows, int *cols, int *x, int *y, int *z, int *ncid, int *ret) {
  *ret=_ncnbody_open(*rows,*cols,*x, *y, *z,ncid);
  return;
}

void *ncnbody_open(int *rows, int *cols, int *x, int *y, int *z, int *ncid, int *ret) {
  *ret=_ncnbody_open(*rows,*cols,*x, *y, *z,ncid);
  return;
}

/*---------------------------------------------------------------------------------------*/
int _ncnbody_put(int ncid,float *vbuffer,int lstored,int vbufcols) {
  int stat;
  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t data_start[RANK_data];
  static size_t data_count[RANK_data];
  static float time[10];
  size_t time_len;
/*   time_t now = time(NULL); */
  
/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d\n",ncid,lstored,vbufcols); */
  
  /* store time */
  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  time[0] = vbuffer[0];	/* simtime, first value of first record */   
  time[1] = 0.1 ;/* (float)(now - begin );  exectime */
  time[2] = timeid_start_fielddesc;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);

  
  time_len = 1;			/* number of records of data data */
  data_start[0] = timeid_start_particles;
  data_start[1] = 0;
  data_start[2] = 0;
  data_count[0] = time_len;
  data_count[1] = rows_len;
  data_count[2] = cols_len;
  stat = nc_put_vara_float(ncid, data_id, data_start, data_count, vbuffer);
  check_err(stat,__LINE__,__FILE__);

/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  
  timeid_start++;
  timeid_start_particles++;

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
int _ncnbody_putfield(int ncid, float simtime, int fieldnr, int xdim, int ydim, int zdim, float *vbuffer) {
  int stat;

  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t field_start[RANK_field];
  static size_t field_count[RANK_field];
 

  static float time[10];
  size_t time_len;
/*   time_t now = time(NULL); */
  
/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d\n",ncid,lstored,vbufcols); */
  
  /* store time */
  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  time[0] = vbuffer[0];	/* simtime, first value of first record */   
  time[1] = 0.1 ;/* (float)(now - begin );  exectime */
  time[2] = timeid_start_fielddesc;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);

  field_start[1] = 0;
  field_start[2] = 0;
  field_start[3] = 0;
  
  field_count[0] = 1;
  field_count[1] = xdim;
  field_count[2] = ydim;
  field_count[3] = zdim;
  
  switch( fieldnr ){
  case 1:
    field_start[0] = timeid_start_field1; 
    stat = nc_put_vara_float(ncid, field_id1, field_start, field_count, vbuffer);
    check_err(stat,__LINE__,__FILE__);
    timeid_start_field1++;
  break;
  case 2: 
    field_start[0] = timeid_start_field2; 
    stat = nc_put_vara_float(ncid, field_id2, field_start, field_count, vbuffer);
    check_err(stat,__LINE__,__FILE__);
    timeid_start_field2++;
    break;
  case 3: 
    field_start[0] = timeid_start_field3; 
    stat = nc_put_vara_float(ncid, field_id3, field_start, field_count, vbuffer);
    check_err(stat,__LINE__,__FILE__);
    timeid_start_field3++;
    break;
  case 4: 
    field_start[0] = timeid_start_field4; 
    stat = nc_put_vara_float(ncid, field_id4, field_start, field_count, vbuffer);
    check_err(stat,__LINE__,__FILE__);
    timeid_start_field4++;
    break;

  }
  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  
  return(stat);
}

void* ncnbody_putfield_(int *ncid, float *simtime, int *fieldnr, int *xdim, int *ydim, int *zdim, float *vbuffer, int *rc) {
  *rc=_ncnbody_putfield(*ncid, *simtime, *fieldnr, *xdim, *ydim, *zdim, vbuffer);
  return;
}

void* ncnbody_putfield__(int *ncid, float *simtime, int *fieldnr, int *xdim, int *ydim, int *zdim,float *vbuffer, int *rc) {
  *rc=_ncnbody_putfield(*ncid, *simtime, *fieldnr, *xdim, *ydim, *zdim, vbuffer);
  return;
}

void* ncnbody_putfield(int *ncid, float *simtime, int *fieldnr, int *xdim, int *ydim, int *zdim,float *vbuffer, int *rc) {
  *rc=_ncnbody_putfield(*ncid, *simtime, *fieldnr, *xdim, *ydim, *zdim, vbuffer);
  return;
}
/*---------------------------------------------------------------------------------------*/
int _ncnbody_putvecfield(int ncid, float simtime, int xdim, int ydim, int zdim, float *vbuffer) {
  int stat;

  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t field_start[RANK_vecField];
  static size_t field_count[RANK_vecField];
  int compdim = 3;

  static float time[10];
  size_t time_len;
/*   time_t now = time(NULL); */
  
/*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d\n",ncid,lstored,vbufcols); */
  
  /* store time */
  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  time[0] = vbuffer[0];	/* simtime, first value of first record */   
  time[1] = 0.1 ;/* (float)(now - begin );  exectime */
  time[2] = timeid_start_fielddesc;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);

  field_start[0] = timeid_start_vecField;
  field_start[1] = 0;
  field_start[2] = 0;
  field_start[3] = 0;
  field_start[4] = 0; 
  
  field_count[0] = 1;
  field_count[1] = xdim;
  field_count[2] = ydim;
  field_count[3] = zdim;
  field_count[4] = compdim;
  
  stat = nc_put_vara_float(ncid, vecField_id, field_start, field_count, vbuffer);
  check_err(stat,__LINE__,__FILE__);
  timeid_start_vecField++;
  

  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  
  return(stat);
}

void* ncnbody_putvecfield_(int *ncid, float *simtime, int *xdim, int *ydim, int *zdim, float *vbuffer, int *rc) {
  *rc=_ncnbody_putvecfield(*ncid, *simtime, *xdim, *ydim, *zdim, vbuffer);
  return;
}

void* ncnbody_putvecfield__(int *ncid, float *simtime, int *xdim, int *ydim, int *zdim,float *vbuffer, int *rc) {
  *rc=_ncnbody_putvecfield(*ncid, *simtime, *xdim, *ydim, *zdim, vbuffer);
  return;
}

void* ncnbody_putvecfield(int *ncid, float *simtime, int *xdim, int *ydim, int *zdim,float *vbuffer, int *rc) {
  *rc=_ncnbody_putvecfield(*ncid, *simtime, *xdim, *ydim, *zdim, vbuffer);
  return;
}/*---------------------------------------------------------------------------------------*/
int _ncnbody_putvecfielddesc(int ncid, float simtime, float *vbuffer ) {
  int stat;
  int i, j;
  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t vecFielddesc_start[RANK_vecFielddesc];
  static size_t vecFielddesc_count[RANK_vecFielddesc] = {1,6};
  
  static float time[10];
  size_t time_len;
  /*   time_t now = time(NULL); */
  
  vecFielddesc_start[0] = timeid_start_vecFielddesc-1; 
  vecFielddesc_start[1] = 0;
  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  
  time[0] = simtime;
  time[1] = 0.1;   /*(float) ( now - begin );  exectime */
  time[2] = timeid_start_fielddesc-1;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);
   printf("WF _ncnbody_putvecfielddesc: ncid=%d\n",ncid);
   for( j = 0; j < vecFielddesc_count[1]; j++ ){
     printf( "%f ", vbuffer[j] );
   }
   printf( "\n" );
   
   
  stat = nc_put_vara_float(ncid, vecField_desc_id, vecFielddesc_start, vecFielddesc_count, vbuffer);
  check_err(stat,__LINE__,__FILE__);

  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  
  printf("WF _ncnbody_putvecfielddesc: ncid=%d ready, stat = %d\n",ncid, stat); 
  return(stat);
}

void* ncnbody_putvecfielddesc_(int *ncid, float *simtime, float *vbuffer, int *rc) {
  *rc=_ncnbody_putvecfielddesc(*ncid, *simtime, vbuffer);
  return;
}

void* ncnbody_putvecfielddesc__(int *ncid, float *simtime, float *vbuffer, int *rc) {
  *rc=_ncnbody_putvecfielddesc(*ncid, *simtime, vbuffer );
  return;
}

void* ncnbody_putvecfielddesc(int *ncid, float *simtime, float *vbuffer, int *rc) {
  *rc=_ncnbody_putvecfielddesc(*ncid, *simtime, vbuffer);
  return;
}
/*---------------------------------------------------------------------------------------*/
int _ncnbody_putselvecfield(int ncid, float simtime, int selVecField ) {
  int stat;

  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t selVecField_start[2];
  static size_t selVecField_count[2] = {1, 1};

  static float time[10];
  size_t time_len;
/*   time_t now = time(NULL); */

  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  
  time[0] = simtime;
  time[1] = 0.1;   /*(float) ( now - begin );  exectime */
  time[2] = timeid_start_fielddesc;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);
  
 
  selVecField_start[0] = timeid_start_fielddesc; 
  selVecField_start[1] = 0;

  printf("WF _ncnbody_putselvecfield: ncid=%d, selVecfield=%d\n",ncid,selVecField); 

  stat = nc_put_vara_int(ncid, sel_vecField_id, selVecField_start, selVecField_count, &selVecField);
  check_err(stat,__LINE__,__FILE__);
  printf("WF _ncnbody_putselvecfield: ncid=%d, selVecfield=%d, ready stat = %d\n",ncid,selVecField, stat); 
  nc_sync(ncid);
  
  timeid_start_vecFielddesc++;
  timeid_start++;
  printf("WF _ncnbody_putselvecfield ready\n" );
  return(stat);
}

void* ncnbody_putselvecfield_(int *ncid, float *simtime, int *selfield, int *rc) {
  *rc=_ncnbody_putselvecfield(*ncid, *simtime, *selfield );
  return;
}
 
void* ncnbody_putselvecfield__(int *ncid, float *simtime, int *selfield, int *rc) {
  *rc=_ncnbody_putselvecfield(*ncid, *simtime, *selfield );
  return;
}

void* ncnbody_putselvecfield(int *ncid, float *simtime, int *selfield, int *rc) {
  *rc=_ncnbody_putselvecfield(*ncid, *simtime, *selfield );
  return;
}

/*---------------------------------------------------------------------------------------*/
int _ncnbody_putfielddesc(int ncid, float simtime, float *vbuffer ) {
  int stat;
  int i, j;
  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t fielddesc_start[RANK_fielddesc];
  static size_t fielddesc_count[RANK_fielddesc] = {1, 4, 6};

  
  static float time[10];
  size_t time_len;
  /*   time_t now = time(NULL); */
  
  fielddesc_start[0] = timeid_start_fielddesc-1; 
  fielddesc_start[1] = 0;
  fielddesc_start[2] = 0;
  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  
  time[0] = simtime;
  time[1] = 0.1;   /*(float) ( now - begin );  exectime */
  time[2] = timeid_start_fielddesc-1;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);
   printf("WF _ncnbody_putfielddesc: ncid=%d\n",ncid);
   for( i = 0; i < fielddesc_count[1]; i++ ){
     printf( "field%d : ", i );
     for( j = 0; j < fielddesc_count[2]; j++ ){
       printf( "%f ", vbuffer[6*i+j] );
     }
     printf( "\n" );
   }
   
  stat = nc_put_vara_float(ncid, field_desc_id, fielddesc_start, fielddesc_count, vbuffer);
  check_err(stat,__LINE__,__FILE__);

  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  
  printf("WF _ncnbody_putfielddesc: ncid=%d ready, stat = %d\n",ncid, stat); 
  return(stat);
}

void* ncnbody_putfielddesc_(int *ncid, float *simtime, float *vbuffer, int *rc) {
  *rc=_ncnbody_putfielddesc(*ncid, *simtime, vbuffer);
  return;
}

void* ncnbody_putfielddesc__(int *ncid, float *simtime, float *vbuffer, int *rc) {
  *rc=_ncnbody_putfielddesc(*ncid, *simtime, vbuffer );
  return;
}

void* ncnbody_putfielddesc(int *ncid, float *simtime, float *vbuffer, int *rc) {
  *rc=_ncnbody_putfielddesc(*ncid, *simtime, vbuffer);
  return;
}
/*---------------------------------------------------------------------------------------*/
int _ncnbody_putselfield(int ncid, float simtime, int selfield1, int selfield2, int selfield3, int selfield4 ) {
  int stat;
  int selfield[4];

  static size_t time_start[RANK_time];
  static size_t time_count[RANK_time];
  static size_t selfield_start[3];
  static size_t selfield_count[3] = {1, 4, 1};

  static float time[10];
  size_t time_len;
/*   time_t now = time(NULL); */

  time_len = 10;		
  time_start[0] = timeid_start;
  time_start[1] = 0;
  time_count[0] = 1;
  time_count[1] = time_len;
  
  time[0] = simtime;
  time[1] = 0.1;   /*(float) ( now - begin );  exectime */
  time[2] = timeid_start_fielddesc;
  time[3] = timeid_start_field1;
  time[4] = timeid_start_field2;
  time[5] = timeid_start_field3;
  time[6] = timeid_start_field4;
  time[7] = timeid_start_particles;
  time[8] = timeid_start_vecField;
  time[9] = timeid_start_vecFielddesc;

  stat = nc_put_vara_float(ncid, time_id, time_start, time_count, time);
  check_err(stat,__LINE__,__FILE__);
  
  
  selfield[0] = selfield1;
  selfield[1] = selfield2;
  selfield[2] = selfield3;
  selfield[3] = selfield4;
  
  selfield_start[0] = timeid_start_fielddesc; 
  selfield_start[1] = 0;
  selfield_start[2] = 0;
  

  printf("WF _ncnbody_putselfield: ncid=%d, selfields=%d, %d, %d, %d\n",ncid,selfield1, selfield2, selfield3, selfield4); 

  stat = nc_put_vara_int(ncid, sel_fields_id, selfield_start, selfield_count, selfield);
  check_err(stat,__LINE__,__FILE__);
  printf("WF _ncnbody_putselfield: ncid=%d, selfields=%d, %d, %d, %d ready stat = %d\n",ncid,selfield1, selfield2, selfield3, selfield4, stat); 

 /*     printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  nc_sync(ncid);
  /*   printf("WF _ncnbody_put: ncid=%d, lstored=%d vbufcols=%d sync ready stat=%d\n",ncid,lstored,vbufcols,stat); */
  
  timeid_start_fielddesc++;
  timeid_start++;
  printf("WF _ncnbody_putselfield ready\n" );
  return(stat);
}

void* ncnbody_putselfield_(int *ncid, float *simtime, int *selfield1, int *selfield2, int *selfield3, int *selfield4, int *rc) {
  *rc=_ncnbody_putselfield(*ncid, *simtime, *selfield1, *selfield2, *selfield3, *selfield4 );
  return;
}
 
void* ncnbody_putselfield__(int *ncid, float *simtime, int *selfield1, int *selfield2, int *selfield3, int *selfield4, int *rc) {
  *rc=_ncnbody_putselfield(*ncid, *simtime, *selfield1, *selfield2, *selfield3, *selfield4 );
  return;
}

void* ncnbody_putselfield(int *ncid, float *simtime, int *selfield1, int *selfield2, int *selfield3, int *selfield4, int *rc) {
  *rc=_ncnbody_putselfield(*ncid, *simtime, *selfield1, *selfield2, *selfield3, *selfield4 );
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
