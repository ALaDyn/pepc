echo "Starting coulomb explosion"
mpirun -np 4 ../../bin/pepce ./run.h
echo "..done."

gnuplot fields_movie.gnuplot
