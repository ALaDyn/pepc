echo "Starting billiard mode ..."
cp billiards.h run.h
llrun -p1 ../src/pepc 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat log/energy.$TEND
cp run.out log/run.$TEND
echo "... done"
