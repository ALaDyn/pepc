echo "Starting first job ..."
cp run1.h run.h
llrun -p4 $HOME/tree/pepc/src/pepc 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done"
echo "Submitting second job ..."
./block2.sh 
