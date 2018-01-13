#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "./Random123/threefry.h"
#include "./Random123/u01.h"

int gen_norm_double_rng(unsigned int*ctr_seed, unsigned int*key_seed, double*output) {
  threefry4x32_ctr_t ctr, ran_num;
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

  return 0; 
}

//int main() {
// #pragma omp parallel num_threads(4) default(none)
//{
//int i;
//int*ctr_s, *key_s;
//double*rand_num;
//
//ctr_s = (int*) malloc(4*sizeof(int));
//key_s = (int*) malloc(4*sizeof(int));
//rand_num = (double*) malloc(4*sizeof(double));
//int id = omp_get_thread_num();
//
//for(i=0; i < 4; i + +) {
//ctr_s[i] = 123 - i*7 + id;
//key_s[i] = i*7 + id;
//rand_num[i] = 0.0;
//}
//
//gen_norm_double_rng(ctr_s, key_s, rand_num);
//
//for(i=0; i < 4; i + +) {
//printf("id: %d\tGenerated numbers: %lf\n", id, rand_num[i]);
//}
//}
//}
