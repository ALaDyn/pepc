echo "Starting dust  .."
cp dust.h run.h
#cp pe000/parts_info.000400 parts_info.in
llrun -p8 ../src/pepc 
 eval `cat runstamp |  awk '{ t = $1;  printf  "TEND=%s",t}'`
cp energy.dat energy.$TEND
cp run.out run.$TEND
echo "... done" 
