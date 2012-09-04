#!/bin/bash

# use absolute path here !
SOURCEDIRS=/lustre/jwork/jzam04/jzam0415/icosahedron_spatially_resolved/1000.fields_resolution_lowdens/rundir-?????
# this may also be relative
TARGETDIR=./1000.fields_resolution_lowdens

STARTDIR=`pwd`

source modules.sh

mkdir -p ${TARGETDIR}
cd ${TARGETDIR}

echo "========================= spatially_resolved_ccf =========================="
${STARTDIR}/spatially_resolved_ccf ${SOURCEDIRS}
echo "========================= spatially_resolved_ccf_jackknife =========================="
${STARTDIR}/spatially_resolved_ccf_jackknife ${SOURCEDIRS}


echo "========================= spatially_resolved_spherical_fields =========================="
${STARTDIR}/spatially_resolved_spherical_fields ${SOURCEDIRS}
echo "========================= spatially_resolved_spherical_fields_jackknife =========================="
${STARTDIR}/spatially_resolved_spherical_fields_jackknife ${SOURCEDIRS}

echo "========================= spatially_resolved_spherical_spectrum_jackknife =========================="
${STARTDIR}/spatially_resolved_spherical_spectrum_jackknife ${SOURCEDIRS}

echo "========================= plot_spherical_spectrum =========================="
${STARTDIR}/plot_spherical_spectrum.py
echo "========================= plot_spherical_fourier_coeffs =========================="
${STARTDIR}/plot_spherical_fourier_coeffs.py

echo "========================= creating latex document =========================="
ln -sf ${STARTDIR}/allplots.tex ./allplots.tex
pdflatex allplots.tex

cd ${STARTDIR}
