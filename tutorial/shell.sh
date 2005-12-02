echo "Starting PEPC-B: shell target ..."
cp shell.h run.h
llrun  -p 8 ../bin/pepcb 
echo "... done" 
