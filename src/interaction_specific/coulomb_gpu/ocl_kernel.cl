// !!!! GLOBAL DEFINE !!!!
#define MAX_IACT_PARTNERS 16 * 256 * 2 * 2

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

__kernel void ocl_gpu_kernel(int queued, double eps2, __global double* partner, __global double* results, int dummy)
{
   // index on GPU device to id array entry
   const int idx = get_global_id(0);
   const int local_index = get_local_id(0);
   const int local_size = get_local_size(0);

   // local buffer for global Data
   __local double pot_local[128];
   __local double e_1_local[128];
   __local double e_2_local[128];
   __local double e_3_local[128];

   if (idx >= queued)
   {
      pot_local[local_index] = 0.0;
      e_1_local[local_index] = 0.0;
      e_2_local[local_index] = 0.0;
      e_3_local[local_index] = 0.0;
   }

   // compute individal interactions (work item)
   if (idx < queued)
   {
      // private temp variables
      const double dx = partner[DELTA1+idx];
      const double dy = partner[DELTA2+idx];
      const double dz = partner[DELTA3+idx];
      const double dx2 = dx*dx;
      const double dy2 = dy*dy;
      const double dz2 = dz*dz;
      const double dx3 = dx*dx2;
      const double dy3 = dy*dy2;
      const double dz3 = dz*dz2;

      const double dist2 = dx2 + dy2 + dz2;

      // more private temp variables
      const double rd = rsqrt(dist2+eps2);
      const double rd2 = rd *rd;
      const double rd3 = rd *rd2;
      const double rd5 = rd3*rd2;
      const double rd7 = rd5*rd2;

      const double fd1 = 3.0*dx2*rd5 - rd3;
      const double fd2 = 3.0*dy2*rd5 - rd3;
      const double fd3 = 3.0*dz2*rd5 - rd3;
      const double fd4 = 3.0*dx*dy*rd5;
      const double fd5 = 3.0*dy*dz*rd5;
      const double fd6 = 3.0*dx*dz*rd5;

      const double m2rd5=2.0*rd5;
      const double m5rd7=5.0*rd7;
      const double pre1=m5rd7*dx*dy*dz;
      const double pre2x=m5rd7*dx2 - rd5;
      const double pre2y=m5rd7*dy2 - rd5;
      const double pre2z=m5rd7*dz2 - rd5;
      const double preQ1=pre2x*partner[QUAD1+idx];
      const double preQ2=pre2y*partner[QUAD2+idx];
      const double preQ3=pre2z*partner[QUAD3+idx];

      pot_local[local_index] = partner[CHARGE+idx]*rd
         + (dx*partner[DIP1+idx] + dy*partner[DIP2+idx] + dz*partner[DIP3+idx])*rd3
         +  0.5*(fd1*partner[QUAD1+idx]  + fd2*partner[QUAD2+idx]  + fd3*partner[QUAD3+idx])
         +       fd4*partner[XYQUAD+idx] + fd5*partner[YZQUAD+idx] + fd6*partner[ZXQUAD+idx];

      e_1_local[local_index] = partner[CHARGE+idx]*dx*rd3
         + fd1*partner[DIP1+idx] + fd4*partner[DIP2+idx] + fd6*partner[DIP3+idx]
         + 3.0   * (
               0.5 * dx * (
                  + ( pre2x - m2rd5 )*partner[QUAD1+idx]
                  + preQ2
                  + preQ3
                  )
               + dy * pre2x * partner[XYQUAD+idx]
               + dz * pre2x * partner[ZXQUAD+idx]
               + pre1 * partner[YZQUAD+idx]
               );

      e_2_local[local_index] = partner[CHARGE+idx]*dy*rd3
         + fd2*partner[DIP2+idx] + fd4*partner[DIP1+idx] + fd5*partner[DIP3+idx]
         + 3.0 * (
               0.5 * dy * (
                  + ( pre2y - m2rd5 )*partner[QUAD2+idx]
                  + preQ1
                  + preQ3
                  )
               + dx * pre2y * partner[XYQUAD+idx]
               + dz * pre2y * partner[YZQUAD+idx]
               + pre1 * partner[ZXQUAD+idx]
               );

      e_3_local[local_index] = partner[CHARGE+idx]*dz*rd3
         + fd3*partner[DIP3+idx] + fd5*partner[DIP2+idx] + fd6*partner[DIP1+idx]
         + 3.0 * (
               0.5 * dz * (
                  + ( pre2z - m2rd5 )*partner[QUAD3+idx]
                  + preQ2
                  + preQ1
                  )
               + dx * pre2z * partner[ZXQUAD+idx]
               + dy * pre2z * partner[YZQUAD+idx]
               + pre1 * partner[XYQUAD+idx]
               );
   }

   barrier(CLK_LOCAL_MEM_FENCE);

   // now reduce local data to global results
   for(int offset = local_size / 2;
         offset > 0;
         offset >>= 1) {
      if (local_index < offset) {
         pot_local[local_index] += pot_local[local_index + offset];
         e_1_local[local_index] += e_1_local[local_index + offset];
         e_2_local[local_index] += e_2_local[local_index + offset];
         e_3_local[local_index] += e_3_local[local_index + offset];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
   }
   if (local_index == 0) {
      results[(((queued-1)/128 + 1 ) * (1-1))+get_group_id(0)] = pot_local[0];
      results[(((queued-1)/128 + 1 ) * (2-1))+get_group_id(0)] = e_1_local[0];
      results[(((queued-1)/128 + 1 ) * (3-1))+get_group_id(0)] = e_2_local[0];
      results[(((queued-1)/128 + 1 ) * (4-1))+get_group_id(0)] = e_3_local[0];
   }

}
