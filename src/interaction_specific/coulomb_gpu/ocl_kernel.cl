// !!!! GLOBAL DEFINE !!!!
#define MAX_IACT_PARTNERS 2 * 256 * 2 * 2

// define off-sets in our 1-D array
#define DELTA1 (MAX_IACT_PARTNERS * (1-1))
#define DELTA2 (MAX_IACT_PARTNERS * (2-1))
#define DELTA3 (MAX_IACT_PARTNERS * (3-1))
#define CHARGE (MAX_IACT_PARTNERS * (4-1))
#define DIP1   (MAX_IACT_PARTNERS * (5-1))
#define DIP2   (MAX_IACT_PARTNERS * (6-1))
#define DIP3   (MAX_IACT_PARTNERS * (7-1))
#define QUAD1  (MAX_IACT_PARTNERS * (8-1))
#define QUAD2  (MAX_IACT_PARTNERS * (9-1))
#define QUAD3  (MAX_IACT_PARTNERS * (10-1))
#define XYQUAD (MAX_IACT_PARTNERS * (11-1))
#define YZQUAD (MAX_IACT_PARTNERS * (12-1))
#define ZXQUAD (MAX_IACT_PARTNERS * (13-1))

#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "Double precision floating point not supported by OpenCL implementation."
#endif

__kernel void ocl_gpu_kernel(int queued, double eps2, __global double* partner, __global double* pot, __global double* e_1, __global double* e_2, __global double* e_3, int dummy, __global double* debug_data)
{
   // index on GPU device to id array entry
   const int idx = get_global_id(0);

   // local temp variables
   double dist2;
   double dx, dy, dz;
   double r, rd, rd2, rd3, rd5, rd7, dx2, dy2, dz2, dx3, dy3, dz3;
   double fd1, fd2, fd3, fd4, fd5, fd6;

   if (idx < queued)
   {
      dx = partner[DELTA1+idx];
      dy = partner[DELTA2+idx];
      dz = partner[DELTA3+idx];
      dx2 = dx*dx;
      dy2 = dy*dy;
      dz2 = dz*dz;
      dx3 = dx*dx2;
      dy3 = dy*dy2;
      dz3 = dz*dz2;

      dist2 =         dx2;
      dist2 = dist2 + dy2;
      dist2 = dist2 + dz2;

      // r  = sqrt(dist2+eps2);
      //rd = 1.0 / r;
      rd = rsqrt(dist2+eps2);
      rd2 = rd *rd;
      rd3 = rd *rd2;
      rd5 = rd3*rd2;
      rd7 = rd5*rd2;

      fd1 = 3.0*dx2*rd5 - rd3;
      fd2 = 3.0*dy2*rd5 - rd3;
      fd3 = 3.0*dz2*rd5 - rd3;
      fd4 = 3.0*dx*dy*rd5;
      fd5 = 3.0*dy*dz*rd5;
      fd6 = 3.0*dx*dz*rd5;

      pot[idx] = partner[CHARGE+idx]*rd
	 + (dx*partner[DIP1+idx] + dy*partner[DIP2+idx] + dz*partner[DIP3+idx])*rd3
	 +  0.5*(fd1*partner[QUAD1+idx]  + fd2*partner[QUAD2+idx]  + fd3*partner[QUAD3+idx])
	 +       fd4*partner[XYQUAD+idx] + fd5*partner[YZQUAD+idx] + fd6*partner[ZXQUAD+idx];

      e_1[idx] = partner[CHARGE+idx]*dx*rd3
	 + fd1*partner[DIP1+idx] + fd4*partner[DIP2+idx] + fd6*partner[DIP3+idx]
	 + 3.0   * (
	       0.5 * (
		  + ( 5.0*dx3   *rd7 - 3.0*dx*rd5 )*partner[QUAD1+idx]
		  + ( 5.0*dx*dy2*rd7 -     dx*rd5 )*partner[QUAD2+idx]
		  + ( 5.0*dx*dz2*rd7 -     dx*rd5 )*partner[QUAD3+idx]
		  )
	       + ( 5.0*dy*dx2  *rd7 - dy*rd5 )*partner[XYQUAD+idx]
	       + ( 5.0*dz*dx2  *rd7 - dz*rd5 )*partner[ZXQUAD+idx]
	       + ( 5.0*dx*dy*dz*rd7          )*partner[YZQUAD+idx]
	       );

      e_2[idx] = partner[CHARGE+idx]*dy*rd3
	 + fd2*partner[DIP2+idx] + fd4*partner[DIP1+idx] + fd5*partner[DIP3+idx]
	 + 3.0 * (
	       0.5 * (
		  + ( 5.0*dy3*rd7    - 3.0*dy*rd5 )*partner[QUAD2+idx]
		  + ( 5.0*dy*dx2*rd7 -     dy*rd5 )*partner[QUAD1+idx]
		  + ( 5.0*dy*dz2*rd7 -     dy*rd5 )*partner[QUAD3+idx]
		  )
	       + ( 5.0*dx*dy2  *rd7 - dx*rd5 )*partner[XYQUAD+idx]
	       + ( 5.0*dz*dy2  *rd7 - dz*rd5 )*partner[YZQUAD+idx]
	       + ( 5.0*dx*dy*dz*rd7          )*partner[ZXQUAD+idx]
	       );

      e_3[idx] = partner[CHARGE+idx]*dz*rd3
	 + fd3*partner[DIP3+idx] + fd5*partner[DIP2+idx] + fd6*partner[DIP1+idx]
	 + 3.0 * (
	       0.5 * (
		  + ( 5.0*dz3   *rd7 - 3.0*dz*rd5 )*partner[QUAD3+idx]
		  + ( 5.0*dz*dy2*rd7 -     dz*rd5 )*partner[QUAD2+idx]
		  + ( 5.0*dz*dx2*rd7 -     dz*rd5 )*partner[QUAD1+idx]
		  )
	       + ( 5.0*dx*dz2  *rd7 - dx*rd5 )*partner[ZXQUAD+idx]
	       + ( 5.0*dy*dz2  *rd7 - dy*rd5 )*partner[YZQUAD+idx]
	       + ( 5.0*dx*dy*dz*rd7          )*partner[XYQUAD+idx]
	       );
   }
}
