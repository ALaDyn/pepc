echo "Starting weibel  .."
srun --ntasks=2 --ntasks-per-node=2 --cpus-per-task=8 ../../bin/pepc-2dd ./params_weibel
echo "... done" 
