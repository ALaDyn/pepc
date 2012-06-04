rm -f energy.dat energies.eps

echo "Starting eqm  .."
mpiexec -l -np 2 ../../bin/pepc-b sphere.h
echo "... done" 

gle energies.gle
display energies.eps
