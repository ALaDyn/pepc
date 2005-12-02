echo "Starting tree debug ..."
cp treetest.h run.h
mpirun -np 2 ../bin/pepcb 
echo "... done"
