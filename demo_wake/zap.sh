echo "Firing laser .."
cp shoot.h run.h
llrun -p8 $HOME/wake/pepc/src/pepc 
 eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
