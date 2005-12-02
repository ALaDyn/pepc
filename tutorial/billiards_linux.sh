echo "Starting billiard mode ..."
cp billiards.h run.h
mpirun -np 5 ../bin/pepcb 
echo "... done"
