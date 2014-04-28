// !!!! GLOBAL DEFINE !!!!
#define MAX_IACT_PARTNERS 1 * 256 * 1 * 2


typedef struct {                                           // type :: mpdelta
    __global double delta1[MAX_IACT_PARTNERS];             //    real*8 :: delta1(1:MAX_IACT_PARTNERS)
    __global double delta2[MAX_IACT_PARTNERS];             //    real*8 :: delta2(1:MAX_IACT_PARTNERS)
    __global double delta3[MAX_IACT_PARTNERS];             //    real*8 :: delta3(1:MAX_IACT_PARTNERS)
    __global double charge[MAX_IACT_PARTNERS];             //    real*8 :: charge(1:MAX_IACT_PARTNERS)
    __global double dip1[MAX_IACT_PARTNERS];               //    real*8 :: dip1(1:MAX_IACT_PARTNERS)
    __global double dip2[MAX_IACT_PARTNERS];               //    real*8 :: dip2(1:MAX_IACT_PARTNERS)
    __global double dip3[MAX_IACT_PARTNERS];               //    real*8 :: dip3(1:MAX_IACT_PARTNERS)
    __global double quad1[MAX_IACT_PARTNERS];              //    real*8 :: quad1(1:MAX_IACT_PARTNERS)
    __global double quad2[MAX_IACT_PARTNERS];              //    real*8 :: quad2(1:MAX_IACT_PARTNERS)
    __global double quad3[MAX_IACT_PARTNERS];              //    real*8 :: quad3(1:MAX_IACT_PARTNERS)
    __global double xyquad[MAX_IACT_PARTNERS];             //    real*8 :: xyquad(1:MAX_IACT_PARTNERS)
    __global double yzquad[MAX_IACT_PARTNERS];             //    real*8 :: yzquad(1:MAX_IACT_PARTNERS)
    __global double zxquad[MAX_IACT_PARTNERS];             //    real*8 :: zxquad(1:MAX_IACT_PARTNERS)
} mpdelta;                                                 // end type mpdelta
__kernel void force(int queued, double eps2, __global gpu_type* gpu, __global double* pot, __global double* e_1, __global double* e_2, __global double* e_3)
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
      dist2 =         gpu.delta1[idx] * gpu.delta1[idx]
      dist2 = dist2 + gpu.delta2[idx] * gpu.delta2[idx]
      dist2 = dist2 + gpu.delta3[idx] * gpu.delta3[idx]

      dx = gpu.delta1[idx];
      dy = gpu.delta2[idx];
      dz = gpu.delta3[idx];

      r  = sqrt(dist2+eps2);
      rd = 1.0/r;
      rd2 = rd *rd;
      rd3 = rd *rd2;
      rd5 = rd3*rd2;
      rd7 = rd5*rd2;

      dx2 = dx*dx;
      dy2 = dy*dy;
      dz2 = dz*dz;
      dx3 = dx*dx2;
      dy3 = dy*dy2;
      dz3 = dz*dz2;

      fd1 = 3.0*dx2*rd5 - rd3;
      fd2 = 3.0*dy2*rd5 - rd3;
      fd3 = 3.0*dz2*rd5 - rd3;
      fd4 = 3.0*dx*dy*rd5;
      fd5 = 3.0*dy*dz*rd5;
      fd6 = 3.0*dx*dz*rd5;

      // idx, gpu_id translates to idx+(gpu_i-1)*MAX_IACT_PARTNERS !!!! +/- 1 depending on start and couting...
      // however, we index the input data with gpu_id, so forget it here...

      pot[idx] = gpu.charge[idx]*rd
	 + (dx*gpu.dip1[idx] + dy*gpu.dip2[idx] + dz*gpu.dip3[idx])*rd3
	 +  0.5*(fd1*gpu.quad1[idx]  + fd2*gpu.quad2[idx]  + fd3*gpu.quad3[idx])
	 +       fd4*gpu.xyquad[idx] + fd5*gpu.yzquad[idx] + fd6*gpu.zxquad[idx];

      e_1[idx] = gpu.charge[idx]*dx*rd3
	 + fd1*gpu.dip1[idx] + fd4*gpu.dip2[idx] + fd6*gpu.dip3[idx]
	 + 3.0   * (
	       0.5 * (
		  + ( 5.0*dx3   *rd7 - 3.0*dx*rd5 )*gpu.quad1[idx]
		  + ( 5.0*dx*dy2*rd7 -     dx*rd5 )*gpu.quad2[idx]
		  + ( 5.0*dx*dz2*rd7 -     dx*rd5 )*gpu.quad3[idx]
		  )
	       + ( 5.0*dy*dx2  *rd7 - dy*rd5 )*gpu.xyquad[idx]
	       + ( 5.0*dz*dx2  *rd7 - dz*rd5 )*gpu.zxquad[idx]
	       + ( 5.0*dx*dy*dz*rd7          )*gpu.yzquad[idx]
	       );

      e_2[idx] = gpu.charge[idx]*dy*rd3
	 + fd2*gpu.dip2[idx] + fd4*gpu.dip1[idx] + fd5*gpu.dip3[idx]
	 + 3.0 * (
	       0.5 * (
		  + ( 5.0*dy3*rd7    - 3.0*dy*rd5 )*gpu.quad2[idx]
		  + ( 5.0*dy*dx2*rd7 -     dy*rd5 )*gpu.quad1[idx]
		  + ( 5.0*dy*dz2*rd7 -     dy*rd5 )*gpu.quad3[idx]
		  )
	       + ( 5.0*dx*dy2  *rd7 - dx*rd5 )*gpu.xyquad[idx]
	       + ( 5.0*dz*dy2  *rd7 - dz*rd5 )*gpu.yzquad[idx]
	       + ( 5.0*dx*dy*dz*rd7          )*gpu.zxquad[idx]
	       );

      e_3[idx] = gpu.charge[idx]*dz*rd3
	 + fd3*gpu.dip3[idx] + fd5*gpu.dip2[idx] + fd6*gpu.dip1[idx]
	 + 3.0 * (
	       0.5 * (
		  + ( 5.0*dz3   *rd7 - 3.0*dz*rd5 )*gpu.quad3[idx]
		  + ( 5.0*dz*dy2*rd7 -     dz*rd5 )*gpu.quad2[idx]
		  + ( 5.0*dz*dx2*rd7 -     dz*rd5 )*gpu.quad1[idx]
		  )
	       + ( 5.0*dx*dz2  *rd7 - dx*rd5 )*gpu.zxquad[idx]
	       + ( 5.0*dy*dz2  *rd7 - dy*rd5 )*gpu.yzquad[idx]
	       + ( 5.0*dx*dy*dz*rd7          )*gpu.xyquad[idx]
	       );
   }
}
