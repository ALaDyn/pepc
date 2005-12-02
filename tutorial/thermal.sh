echo "Starting thermal eqm ..."
cp thermal.h run.h
llrun -p4 ../bin/pepcb 
echo "... done"
