echo "Starting PEPC-B: plasma sphere"
cp eqm.h run.h
#llrun  -p4 ../bin/pepcvis 
llrun -p4 -T ../bin/pepcb_p690 
echo "... done" 
