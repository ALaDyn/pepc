echo "Starting clamped eqm ..."
cp clamp.h run.h
llrun -p4 $HOME/tree/pepc/src/pepc2 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat log/energy.$TEND
cp run.out log/run.$TEND
echo "... done"
