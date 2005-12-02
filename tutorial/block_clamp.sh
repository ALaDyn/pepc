echo "Starting clamped eqm ..."
cp clamp.h run.h
llrun -p8 ../bin/pepcb 
echo "... done"
