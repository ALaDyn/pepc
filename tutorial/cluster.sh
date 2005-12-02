echo "Starting cluster run .."
cp cluster.h run.h
llrun -p4 ../bin/pepcb
#mpirun -np 1 ../bin/pepcb 
echo "... done" 
