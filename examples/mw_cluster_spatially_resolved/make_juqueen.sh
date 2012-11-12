#!/bin/bash

module load fftw3/3.3.2_g

mpixlf90_r -g -O5 -c progress.f90 -o progress.o
mpixlf90_r -g -O5 -c module_math.f90 -o module_math.o
mpixlf90_r -g -O5 -I${FFTW_INCLUDE} spatially_resolved_spherical_fields.f90 progress.o module_math.o -lfftw3 -lxlf90 -lxlopt -lxlsmp -lpthread -lm -lc -L${FFTW_LIB} -o spatially_resolved_spherical_fields.juqueen

echo "now call llsubmit xterm1.job"
echo "and therein:"
echo "for dir in `ls ../../spatially_resolved_spherical_fields.juqueen /work/jzam04/jzam0415/cluster.spatially_resolved/300k/rundir-?????`; do"
echo "  runjob --np 1 : ../../spatially_resolved_spherical_fields.juqueen /work/jzam04/jzam0415/cluster.spatially_resolved/300k/rundir-?????"
echo "done"
