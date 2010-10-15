#include <stdio.h>
#include <assert.h>

#include <math.h>

#include "sion.h"

typedef sion_int64 integer;
typedef double     real;

const int nx = 128;
const int ny = 128;
const int nz = 128;

real xm = 180.0;
real xp = 220.0;
real ym = 40.0;
real yp = 200;
real zm = 40.0;
real zp = 200;

const int   info_header_size=7*4+24*4;
const int   read_buffer_size=8192;
      real *read_buffer;

void vtk_output_fields(integer* npp, int ntasks, int sid, FILE* vtk_file)
{
  int rcnt, fcnt;
  double *ne, *ni, *n, *jx, *jy, *jz;

  printf(" ... converting integrated fields ...\n");

  ne = (double*)(malloc(sizeof(double)*nx*ny*nz));
  assert(ne != NULL);
  ni = (double*)(malloc(sizeof(double)*nx*ny*nz));
  assert(ni != NULL);
  n  = (double*)(malloc(sizeof(double)*nx*ny*nz));
  assert(n != NULL);

  jx = (double*)(malloc(sizeof(double)*nx*ny*nz));
  assert(jx != NULL);
  jy = (double*)(malloc(sizeof(double)*nx*ny*nz));
  assert(jy != NULL);
  jz = (double*)(malloc(sizeof(double)*nx*ny*nz));
  assert(jz != NULL);

  for(fcnt=0; fcnt<nx*ny*nz; fcnt++)
    ne[fcnt] = ni[fcnt] = n[fcnt] = jx[fcnt] = jy[fcnt] = jz[fcnt] = 0;


  fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
  fprintf(vtk_file, "PEPC field vtk output\n");
  fprintf(vtk_file, "ASCII\n");
  fprintf(vtk_file, "DATASET STRUCTURED_POINTS\n"); 
  fprintf(vtk_file, "DIMENSIONS %d %d %d\n", nx+1, ny+1, nz+1); 
  fprintf(vtk_file, "ORIGIN 180.0 40.0 40.0\n");
  fprintf(vtk_file, "SPACING %f %f %f\n", (220.0-180.0)/nx, (200.0-40.0)/ny, (200.0-40.0)/nz);
  fprintf(vtk_file, "CELL_DATA %d\n", nx*ny*nz);  

  assert(NULL!=read_buffer);

  for(rcnt=0; rcnt<ntasks; rcnt++)
    {
      integer left = npp[rcnt];
      integer read_round = 0;
      integer lcnt;
      while(left>0)
	{
	  int read;
	  int bcnt;
	  for(bcnt=0; bcnt<8; bcnt++)
	    {
	      sion_int64 file_offset = info_header_size + (bcnt*npp[rcnt] + read_round*read_buffer_size)*sizeof(double);
	      int buffer_offset = bcnt*read_buffer_size;
	      
	      sion_seek(sid,rcnt,SION_ABSOLUTE_POS,file_offset);
	      read=sion_fread(read_buffer+buffer_offset,sizeof(double),read_buffer_size,sid);
	      assert(read > 0);
	    }
	  
	  read_round++;

	  left-=read;

	  if(left < 0) read+=left;
	  
	  for(lcnt=0; lcnt<read; lcnt++)
	    {
	      double x = read_buffer[lcnt];
	      double y = read_buffer[lcnt+read_buffer_size];
	      double z = read_buffer[lcnt+2*read_buffer_size];
	      
	      double ux = read_buffer[lcnt+3*read_buffer_size];
	      double uy = read_buffer[lcnt+4*read_buffer_size];
	      double uz = read_buffer[lcnt+5*read_buffer_size];
	      
	      double q = read_buffer[lcnt+6*read_buffer_size];
	      double m = read_buffer[lcnt+7*read_buffer_size];
	      
      
	      double dx = (xp - xm)/nx;
	      double dy = (yp - ym)/ny;
	      double dz = (zp - zm)/nz;
	      
	      int ix = floor((x - xm)/dx);
	      int iy = floor((y - ym)/dy);
	      int iz = floor((z - zm)/dz);
	      
	      long int ind = iz*nx*ny + iy*nx + ix;
	      
	      n[ind]++;
	      
	      if(q > 0) ni[ind]++;
	      else ne[ind]++;
	      
	      jx[ind] += ux/m*q;
	      jy[ind] += uy/m*q;
	      jz[ind] += uz/m*q;
	    }
	  
	}
    }
  
  fprintf(vtk_file, "SCALARS n float\n");
  fprintf(vtk_file, "LOOKUP_TABLE default\n");
  
  for(fcnt=0; fcnt<nx*ny*nz; fcnt++)
    fprintf(vtk_file, "%e\n", n[fcnt]);

  fprintf(vtk_file, "SCALARS ni float\n");
  fprintf(vtk_file, "LOOKUP_TABLE default\n");
  
  for(fcnt=0; fcnt<nx*ny*nz; fcnt++)
    fprintf(vtk_file, "%e\n", ni[fcnt]);

  fprintf(vtk_file, "SCALARS ne float\n");
  fprintf(vtk_file, "LOOKUP_TABLE default\n");
  
  for(fcnt=0; fcnt<nx*ny*nz; fcnt++)
    fprintf(vtk_file, "%e\n", ne[fcnt]);

  fprintf(vtk_file, "VECTORS j float\n");
  
  for(fcnt=0; fcnt<nx*ny*nz; fcnt++)
    fprintf(vtk_file, "%e %e %e\n", jx[fcnt], jy[fcnt], jz[fcnt]);
}

void vtk_output_particle_positions(integer* npp, int ntasks, integer total_particles, int sid, FILE* vtk_file)
{
  integer lcnt, rcnt;

  printf(" ... converting particle positions ...\n");

  assert(NULL!=read_buffer);

  fprintf(vtk_file, "# vtk DataFile Version 3.0\n");
  fprintf(vtk_file, "PEPC particle vtk output\n");
  fprintf(vtk_file, "ASCII\n");
  fprintf(vtk_file, "DATASET POLYDATA\n");
  fprintf(vtk_file, "POINTS %lld float\n", total_particles);

  for(rcnt=0; rcnt<ntasks; rcnt++)
    {
      integer left = npp[rcnt];
      integer read_round = 0;

      while(left>0)
	{
	  int read;
	  int bcnt;
	  for(bcnt=0; bcnt<3; bcnt++)
	    {
	      sion_int64 file_offset = info_header_size + (bcnt*npp[rcnt] + read_round*read_buffer_size)*sizeof(double);
	      int buffer_offset = bcnt*read_buffer_size;
		  
	      sion_seek(sid,rcnt,SION_ABSOLUTE_POS,file_offset);
	      read=sion_fread(read_buffer+buffer_offset,sizeof(double),read_buffer_size,sid);
	      assert(read > 0);
	    }
	      
	  read_round++;

	  left-=read;

	  if(left < 0) read+=left;
	  
	  for(lcnt=0; lcnt<read; lcnt++)
	    {
	      double x = read_buffer[lcnt];
	      double y = read_buffer[lcnt+read_buffer_size];
	      double z = read_buffer[lcnt+2*read_buffer_size];
	      
	      fprintf(vtk_file, "%e %e %e\n", x, y, z);
	    }
	  
	} /* left */
    } /* rcnt */

  fprintf(vtk_file,"VERTICES %d %d\n", total_particles, total_particles*2);
  
  for(lcnt=0; lcnt<total_particles; lcnt++)
    fprintf(vtk_file,"1 %d\n", lcnt);

}

void vtk_output_particle_vector(integer* npp, integer ntasks, integer total_particles, integer sid, FILE* vtk_file, char* name, int f1, int f2, int f3)
{
  int vector = (f2 < 0) ? 0 : 1;
  int rcnt;

  printf(" ... converting particle properties (%s) ... \n", name);

  assert(NULL!=read_buffer);
  assert(f1>=0);
  assert(f2*f3>0);

  if(1==vector)
    {
      fprintf(vtk_file,"VECTORS %s float\n", name);
    }
  else
    {
      fprintf(vtk_file, "SCALARS %s float\n", name);
      fprintf(vtk_file, "LOOKUP_TABLE default\n");
    }

  for(rcnt=0; rcnt<ntasks; rcnt++)
    {
      integer left = npp[rcnt];
      integer read_round = 0;
      integer lcnt;

      while(left>0)
	{
	  int read;
 	  sion_int64 sion_offset = info_header_size +(f1*npp[rcnt]+read_round*read_buffer_size)*sizeof(double);
	  int read_offset = 0;
	  
	  sion_seek(sid,rcnt,SION_ABSOLUTE_POS,sion_offset);
	  read=sion_fread(read_buffer,sizeof(double),read_buffer_size,sid);
	  assert(read > 0);

	  if(1==vector)
	    {
	      sion_offset = info_header_size +(f2*npp[rcnt]+read_round*read_buffer_size)*sizeof(double);
	      read_offset += read_buffer_size;
	      
	      sion_seek(sid,rcnt,SION_ABSOLUTE_POS,sion_offset);
	      read=sion_fread(read_buffer + read_offset,sizeof(double),read_buffer_size,sid);
	      assert(read > 0);

	      sion_offset = info_header_size +(f3*npp[rcnt]+read_round*read_buffer_size)*sizeof(double);
	      read_offset += read_buffer_size;

	      sion_seek(sid,rcnt,SION_ABSOLUTE_POS,sion_offset);
	      read=sion_fread(read_buffer + read_offset,sizeof(double),read_buffer_size,sid);
	      assert(read > 0);
	    }
	  
	  read_round++;
	      
	  left-=read;

	  if(left < 0) read+=left;
	  
	  if(0==vector)
	    {
	      for(lcnt=0; lcnt<read; lcnt++)
		fprintf(vtk_file, "%e\n", read_buffer[lcnt]);
	    }
	  else
	    {
	      for(lcnt=0; lcnt<read; lcnt++)
		fprintf(vtk_file, "%e %e %e\n", read_buffer[lcnt], read_buffer[lcnt+read_buffer_size], read_buffer[lcnt+2*read_buffer_size]);
	    }
	  
	}
      
    }
}

void vtk_output_particle_scalar(integer* npp, integer ntasks, integer total_particles, integer sid, FILE* vtk_file, char* name, int f1)
{
  vtk_output_particle_vector(npp, ntasks, total_particles, sid, vtk_file, name, f1, -1, -1);
}


int main(int argc, char **argv)
{
  integer acnt;

  assert(argc > 1);
  for(acnt=1; acnt<argc; acnt++)
    printf("file number %d -> %s to be processed\n", acnt, argv[acnt]);
  
  read_buffer = (double*)(malloc(sizeof(double)*8*read_buffer_size));
  assert(NULL != read_buffer);

  for(acnt=1; acnt<argc; acnt++)
    {
      integer sid;
      integer *npp;
      FILE *sion_file, *vtk_particle_file, *vtk_field_file;
      char vtk_particle_file_name[255], vtk_field_file_name[255], sion_file_name[255];
      int rcnt, ntasks, nfiles, blocks;
      sion_int64 *chunksizes=NULL;
      int        *globalranks=NULL;
      sion_int64 *sion_chunksizes;
      sion_int64 *sion_globalranks;
      sion_int64 *sion_blockcount;
      sion_int64 *sion_blocksizes;
      sion_int32 fsblksize;
      long wrote,left;
      sion_int64  globalskip, bsumread,bread;
      sion_int64  start_of_varheader;

      char info_buffer[info_header_size];

      sprintf(sion_file_name, "%s", argv[acnt]);
      sprintf(vtk_particle_file_name, "%s.particle.vtk", sion_file_name);
      sprintf(vtk_field_file_name, "%s.field.vtk", sion_file_name);
      printf("target          sion file: %s\n", argv[acnt]);
      printf("target particle vtk  file: %s\n", vtk_particle_file_name);
      printf("target     fieldvtk  file: %s\n", vtk_field_file_name);

      sid = sion_open(sion_file_name, "rb", &ntasks, &nfiles, &chunksizes, &fsblksize, &globalranks, &sion_file);

      printf(" ** sion file info - ntasks: %d\n", ntasks);
      printf(" ** sion file info - nfiles: %d\n", nfiles);
        
      npp = (sion_int64*)(malloc(sizeof(sion_int64)*ntasks));
      
      integer total_particles=0;
      /* get general information header */
      for(rcnt=0; rcnt<ntasks; rcnt++)
	{
	  sion_seek(sid,rcnt,SION_ABSOLUTE_POS,0);
	  bread=sion_fread(info_buffer,1,info_header_size,sid);
	  
	  npp[rcnt] = *((int*)(info_buffer + 4));
	  total_particles += npp[rcnt];
	  if(rcnt%35==0) printf(" ** sion file info - npp on rank %d: %lld\n", rcnt, npp[rcnt]);
	}

      /* /\* get volume information *\/ */
      /* { */
      /* 	sion_seek(sid,rcnt,SION_ABSOLUTE_POS,0); */
      /* 	bread=sion_fread(info_buffer,1,info_header_size,sid); */
	
      /* 	printf(" ** sion file info - visualization x_min: %e, x_max: %e\n", xm, xp); */
      /* } */

      printf("total particles : %lld\n", total_particles);

      vtk_field_file = fopen(vtk_field_file_name, "w");
      assert(NULL!=vtk_field_file);
      vtk_output_fields(npp, ntasks, sid, vtk_field_file);
      fclose(vtk_field_file);      

      vtk_particle_file = fopen(vtk_particle_file_name, "w");
      assert(NULL!=vtk_particle_file);
      vtk_output_particle_positions(npp, ntasks, total_particles, sid, vtk_particle_file);
      vtk_output_particle_vector(npp, ntasks, total_particles, sid, vtk_particle_file, "u", 3, 4, 5);
      vtk_output_particle_scalar(npp, ntasks, total_particles, sid, vtk_particle_file, "q", 6);
      vtk_output_particle_scalar(npp, ntasks, total_particles, sid, vtk_particle_file, "m", 7);
      fclose(vtk_particle_file);
      
      sion_close(sid);

    } /* acnt */

  free(read_buffer);

  return 0;
}


