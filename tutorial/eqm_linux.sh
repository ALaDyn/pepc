echo "Starting eqm  .."
cp eqm.h run.h
mpirun -np 1 ../src/pepc 
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
