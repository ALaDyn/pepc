echo "Starting 2nd job .."
cp run2.h run.h
llrun -p4 $HOME/tree/pepc/src/pepc2 
 eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
