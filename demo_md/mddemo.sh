echo "Starting MD demo  .."
cp eqm.h run.h
#llrun -p 4 ../bin/pepcevis 
llrun -p 4 ../bin/pepce_p690
#mpirun -np 1 ../bin/pepce 
echo "... done" 
