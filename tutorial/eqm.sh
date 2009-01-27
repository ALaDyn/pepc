echo "Starting PEPC-B: plasma sphere"
cp eqm.h run.h
#llrun  -p4 ../bin/pepcvis 
llrun -p4 ../bin/pepcb_power6_ascii 
echo "... done" 
