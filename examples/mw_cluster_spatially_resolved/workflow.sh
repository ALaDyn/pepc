#!/bin/bash

NUMPARTS=3k

# use absolute path here !
SOURCEDIRS=/lustre/jwork/jzam04/jzam0415/PNP14/${NUMPARTS}/rundir-?????
# this may also be relative
TARGETDIR=./PNP14/${NUMPARTS}.fields_resolution_lowdens

STARTDIR=`pwd`

source modules.sh

mkdir -p ${TARGETDIR}
cd ${TARGETDIR}

echo "========================= clusterdata_jackknife =========================="
${STARTDIR}/clusterdata_jackknife ${SOURCEDIRS}
echo "========================= simparams =========================="
${STARTDIR}/simparams.py


echo "========================= momentum_electron_acf =========================="
${STARTDIR}/momentum_electron_acf ${SOURCEDIRS}
echo "========================= momentum_electron_acf_jackknife =========================="
${STARTDIR}/momentum_electron_acf_jackknife ${SOURCEDIRS}


echo "========================= spatially_resolved_spherical_fields =========================="
for dir in ${SOURCEDIRS}; do
  ${STARTDIR}/spatially_resolved_spherical_fields ${SOURCEDIRS}
done
echo "========================= spatially_resolved_spherical_fields_jackknife =========================="
${STARTDIR}/spatially_resolved_spherical_fields_jackknife ${SOURCEDIRS}


echo "========================= spatially_resolved_spherical_spectrum_jackknife =========================="
${STARTDIR}/spatially_resolved_spherical_spectrum_jackknife ${SOURCEDIRS}

echo "========================= plot_spherical_spectrum =========================="
${STARTDIR}/plot_spherical_spectrum.py
echo "========================= plot_spherical_fourier_coeffs =========================="
${STARTDIR}/plot_spherical_fourier_coeffs.py

exit

echo "========================= spatially_resolved_ccf =========================="
for dir in ${SOURCEDIRS}; do
  ${STARTDIR}/spatially_resolved_ccf ${SOURCEDIRS}
done
echo "========================= spatially_resolved_ccf_jackknife =========================="
${STARTDIR}/spatially_resolved_ccf_jackknife ${SOURCEDIRS}


echo "========================= creating latex document =========================="
ln -sf ${STARTDIR}/allplots.tex ./allplots.tex
pdflatex allplots.tex

cd ${STARTDIR}
