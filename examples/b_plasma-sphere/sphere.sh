rm -f energy.dat energies.eps

echo "Starting eqm  .."
mpirun -l -np 2 ../../pepcb sphere.h
echo "... done" 

gle energies.gle
display energies.eps