#!/bin/bash

# use absolute path here !
SOURCEDIRS=/lustre/jwork/jzam04/jzam0415/icosahedron_spatially_resolved/1000.fields_resolution/rundir-?????
# this may also be relative
TARGETDIR=./1000.fields_resolution

mkdir -p ${TARGETDIR}
cd ${TARGETDIR}

source modules.sh

./spatially_resolved_ccf ${SOURCEDIRS}/spatially_resolved.dat
./spatially_resolved_ccf_jackknife ${SOURCEDIRS}/spatially_resolved_momentum_eigenvalues.dat


./spatially_resolved_spherical_fields ${SOURCEDIRS}
./spatially_resolved_spherical_fields_jackknife /lustre/jwork/jzam04/jzam0415/icosahedron_spatially_resolved/1000.fields_resolution/rundir-?????


./plot_spherical_spectrum.py
./plot_spherical_fourier_coeffs.py

cd -
