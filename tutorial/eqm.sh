echo "Starting eqm  .."
cp eqm.h run.h
llrun -p4 ../src/pepc 
 eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
