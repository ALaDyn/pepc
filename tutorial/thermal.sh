echo "Starting thermal eqm ..."
cp thermal.h run.h
llrun -p8 ../src/pepc 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat log/energy.$TEND
cp run.out log/run.$TEND
echo "... done"
