#!/bin/bash

module load fftw3/3.3.2_g lapack/3.3.0_g

mpixlf90_r -g -O5 -c progress.f90 -o progress.o
mpixlf90_r -g -O5 -c module_math.f90 -o module_math.o
mpixlf90_r -g -O5 -I${FFTW_INCLUDE} spatially_resolved_spherical_fields.f90 progress.o module_math.o -lfftw3 -lxlf90 -lxlopt -lxlsmp -lpthread -lm -lc -L${FFTW_LIB} -o spatially_resolved_spherical_fields.juqueen
mpixlf90_r -g -O5 -I${FFTW_INCLUDE} spatially_resolved_spherical_fields_reduced.f90 progress.o module_math.o -lfftw3 -lxlf90 -lxlopt -lxlsmp -lpthread -lm -lc -L${FFTW_LIB} -o spatially_resolved_spherical_fields_reduced.juqueen
mpixlf90_r -g -O5 -I${FFTW_INCLUDE} spatially_resolved_ccf.f90 progress.o module_math.o -llapack -lesslbg -lfftw3 -lxlf90 -lxlopt -lxlsmp -lpthread -lm -lc -L${FFTW_LIB} -L${LAPACK_LIB} -L/opt/ibmmath/essl/5.1/lib64 -o spatially_resolved_ccf.juqueen

echo "now see juqueen.sh and juqueen_template.job for usage on JUQUEEN"
