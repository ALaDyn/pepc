echo "Starting PEPC-B: ion config "
cp ions.h run.h
llrun  -p6 ../bin/pepcb 
echo "... done" 
