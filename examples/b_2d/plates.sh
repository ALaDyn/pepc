echo "Starting plates  .."
mpirun -np 1 ../../pepcb plates.h
echo "... done" 
./plot_plates.py&
./plot_particles2d.py&
