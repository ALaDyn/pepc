echo "Starting eqm  .."
cp eqm.h run.h
mpirun -np 1 ../bin/pepcb 
echo "... done" 
