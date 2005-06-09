echo "Starting KOJAK analysis ..."
cp kojak.h run.h
llrun -p4 ../src/pepc 
echo "... done"
