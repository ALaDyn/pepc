echo "Configuring ions ..."
cp eqm.h run.h
cat eqm.h | grep "ni ="
llrun -p8 $HOME/wake/pepc/src/pepc 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done"
echo "Configuring electrons ..."
cp clamp.h run.h
cat clamp.h | grep "Te_keV ="
llrun -p8 $HOME/wake/pepc/src/pepc 
eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done"
