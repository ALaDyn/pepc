echo "Starting PEPC run: clamped eqm ..."
cp eqm.h run.h
llrun -p4 ../src/pepc 
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
