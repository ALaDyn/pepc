echo "Starting PEPC-B: plasma sphere"
cp eqm.h run.h
#llrun  -p8 ../bin/pepcb_power6_ascii 
#llrun -np 8 ../bin/pepcb_bgp_sion 
mpiexec -n 8 ../bin/pepcb_nicole_ascii
echo "... done" 
