echo "Starting spherical 2D vortex sheet simulation"
mpirun -np 6 ../../bin/pepc-pfasst-v ./vsphere_pfasst.h
echo "..done. - visualize files particles_timeseries.paraview/-.visit using ParaView or Visit"