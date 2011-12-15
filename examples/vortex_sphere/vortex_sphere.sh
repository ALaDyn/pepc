echo "Starting spherical 2D vortex sheet simulation"
mpirun -np 4 ../../bin/pepc-v ./vsphere.h
echo "..done. - visualize files particles_timeseries.paraview/-.visit using ParaView or Visit"
