#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "./Random123/threefry.h"
#include "./Random123/u01.h"

int gen_norm_double_rng(unsigned int*ctr_seed, unsigned int*key_seed, double*output) {
threefry4x32_ctr_t ctr, ran_num, ran_num2;
threefry4x32_key_t key;
ctr.v[0] = ctr_seed[0];
ctr.v[1] = ctr_seed[1];
ctr.v[2] = ctr_seed[2];
ctr.v[3] = ctr_seed[3];
key.v[0] = key_seed[0];
key.v[1] = key_seed[1];
key.v[2] = key_seed[2];
key.v[3] = key_seed[3];
//This function will generate 4 random integers
ran_num = threefry4x32(ctr, key);
//These functions will convert the above int to double in[0, 1]
output[0] = u01_closed_closed_32_53(ran_num.v[0]);
output[1] = u01_closed_closed_32_53(ran_num.v[1]);
output[2] = u01_closed_closed_32_53(ran_num.v[2]);
output[3] = u01_closed_closed_32_53(ran_num.v[3]);

ran_num2 = threefry4x32(ctr, ran_num);

output[4] = u01_closed_closed_32_53(ran_num2.v[0]);
output[5] = u01_closed_closed_32_53(ran_num2.v[1]);
output[6] = u01_closed_closed_32_53(ran_num2.v[2]);
output[7] = u01_closed_closed_32_53(ran_num2.v[3]); 
return 0;
}
