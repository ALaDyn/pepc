echo "Starting billiard mode ..."
cp billiards.h run.h
mpirun -np 5 ../src/pepc 
cp energy.dat log/energy.$TEND
cp run.out log/run.$TEND
echo "... done"
