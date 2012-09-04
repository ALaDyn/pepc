#!/bin/bash

echo "include definitions in this script by calling 'source $0'"

module unload mkl intel GCC parastation
module load intel/12.1.1 parastation/mpi2-intel12-mt-5.0.27-1 fftw/3.3_intel12.1.1 mkl


export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
