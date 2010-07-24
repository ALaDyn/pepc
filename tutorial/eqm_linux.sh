echo "Starting eqm  .."
cp eqm.h run.h
#cp dumps/info_p0000.000040 parts_info.in
mpirun -np 1 ../bin/pepcb.ubu 
echo "... done" 
