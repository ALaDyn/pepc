echo "Laser test  ..."
cp laser.h run.h
mpirun -np 1 ../bin/pepcb 
#llrun -p 4 ../bin/pepcb 
echo "... done"
