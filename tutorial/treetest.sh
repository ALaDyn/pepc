echo "Starting tree debug ..."
cp treetest.h run.h
mpirun -np 2 ../src/pepc 
cp energy.dat log/energy.$TEND
cp run.out log/run.$TEND
echo "... done"
