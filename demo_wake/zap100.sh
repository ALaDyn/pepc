echo "Setting up 100k plasma .."
cp parts_config.100k pe000/parts_dump.000000
cp parts_info.100k pe000/parts_info.000000
cp parts_info.100k parts_info.in
echo "Firing laser .."
cp shoot_100k.h run.h
llrun -p8  $HOME/wake/pepc/src/pepc 
 eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
