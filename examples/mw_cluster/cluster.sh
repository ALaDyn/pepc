echo "Starting simulation on interaction of metallic cluster with an external laser field"
mpirun -np 2 ../../bin/pepc-mw ./cluster.h
echo "..done."
echo "Visualize files particles_timeseries.pvd/-.visit using ParaView or VisIt"
