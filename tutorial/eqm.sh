echo "Starting clamped eqm ..."
cp eqm.h run.h
llrun -p4 ../src/pepc 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat log/energy.$TEND
cp run.out log/run.$TEND
echo "... done"
