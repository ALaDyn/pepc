echo "Starting cluster run .."
cp cluster.h run.h
llrun -p4 ../pepc-b/pepcb
#mpirun -np 1 ../src/pepc 
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
