echo "Starting PEPC-B: ion config "
cp ions.h run.h
llrun  -p8 ../bin/pepcb 
echo "... done" 
